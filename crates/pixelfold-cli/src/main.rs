use std::io::{IsTerminal, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use pixelfold_core::InteractionKind;

mod analyse;
mod graph;
mod load;
mod render;
mod report;
mod resolve;

use graph::GraphFormat;
use report::Format;

/// Terminal protein structure viewer and interaction engine.
#[derive(Parser)]
#[command(
    name = "pixelfold",
    version,
    about = "Terminal protein structure viewer and interaction engine",
    subcommand_negates_reqs = true
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    /// Path to a .pdb/.cif file, or a 4-character PDB id (for example 4HHB)
    structure: Option<String>,

    /// Open the search interface to fetch structures from RCSB
    #[arg(short, long)]
    fetch: bool,

    /// Skip surface calculation for faster loading of large structures
    #[arg(long)]
    no_surface: bool,

    #[command(flatten)]
    common: Common,
}

/// Options every command that loads a structure shares.
#[derive(Clone, clap::Args)]
struct Common {
    /// Directory for cached downloads (default: the system cache dir under pixelfold/)
    #[arg(long, value_name = "DIR", global = true)]
    cache_dir: Option<PathBuf>,

    /// How to resolve atoms modelled in alternate locations
    #[arg(long, value_enum, default_value_t = AltlocArg::Occupancy, global = true)]
    altloc: AltlocArg,
}

/// Options every headless analysis shares.
#[derive(Clone, clap::Args)]
struct Analysis {
    /// Path to a .pdb/.cif file, or a 4-character PDB id
    structure: String,

    /// Restrict the report to what this selection matches, for example
    /// 'chain A and resi 10-40' or 'byres within 5 of resn STI'
    #[arg(short, long, value_name = "QUERY")]
    select: Option<String>,

    /// How to write the result
    #[arg(long, value_enum, default_value_t = Format::Table)]
    format: Format,

    #[command(flatten)]
    common: Common,
}

#[derive(Subcommand)]
enum Command {
    /// Open the interactive viewer (the default when no command is given)
    View {
        /// Path to a .pdb/.cif file, or a 4-character PDB id
        structure: String,

        /// Skip surface calculation for faster loading of large structures
        #[arg(long)]
        no_surface: bool,

        #[command(flatten)]
        common: Common,
    },
    /// Report non-covalent interactions
    Interactions {
        #[command(flatten)]
        analysis: Analysis,

        /// Report only these interaction types (repeatable)
        #[arg(short = 't', long = "type", value_enum, value_name = "TYPE")]
        types: Vec<KindArg>,

        /// Shorthand for --select 'resn CODE', the interactions of one component
        #[arg(short, long, value_name = "CODE", conflicts_with = "select")]
        ligand: Option<String>,
    },
    /// Build the residue interaction network
    Rin {
        /// Path to a .pdb/.cif file, or a 4-character PDB id
        structure: String,

        /// Restrict the network to interactions touching this selection
        #[arg(short, long, value_name = "QUERY")]
        select: Option<String>,

        /// Include only these interaction types (repeatable)
        #[arg(short = 't', long = "type", value_enum, value_name = "TYPE")]
        types: Vec<KindArg>,

        /// How to write the network
        #[arg(long, value_enum, default_value_t = GraphFormat::Json)]
        format: GraphFormat,

        /// Print a structural analysis report (components, hub residues by
        /// betweenness, cut residues) instead of exporting the graph
        #[arg(long)]
        analyze: bool,

        /// Print the shortest chain of interactions between two residues, each
        /// given as CHAIN/RESI (e.g. --path A/10 B/20)
        #[arg(long, num_args = 2, value_names = ["FROM", "TO"])]
        path: Option<Vec<String>>,

        /// Emit a MolViewSpec (.mvsj) scene for Mol*: the structure with the
        /// network's residues highlighted (scoped by --select)
        #[arg(long)]
        mvsj: bool,

        /// Write to this file instead of standard output
        #[arg(short, long, value_name = "FILE")]
        output: Option<PathBuf>,

        #[command(flatten)]
        common: Common,
    },
    /// Report secondary structure per residue
    Ss {
        #[command(flatten)]
        analysis: Analysis,
    },
    /// Report solvent-accessible surface area per residue
    Sasa {
        #[command(flatten)]
        analysis: Analysis,
    },
    /// Render a structure to a PNG file, or to the terminal when no file is given
    Render {
        /// Path to a .pdb/.cif file, or a 4-character PDB id
        structure: String,

        /// Write a PNG here; without it, render to the terminal (or pipe a PNG)
        #[arg(short, long, value_name = "FILE")]
        output: Option<PathBuf>,

        /// Image width in pixels (defaults to 1200 for a file, or the terminal
        /// width for terminal output)
        #[arg(long)]
        width: Option<u32>,

        /// Image height in pixels (defaults to 900 for a file, or the terminal
        /// height for terminal output)
        #[arg(long)]
        height: Option<u32>,

        /// How to colour atoms
        #[arg(long, value_enum, default_value_t = ColorArg::Element)]
        color: ColorArg,

        /// Render the C-alpha backbone only
        #[arg(long)]
        backbone: bool,

        /// Restrict rendering to what this selection matches
        #[arg(short, long, value_name = "QUERY")]
        select: Option<String>,

        /// Scale every atom radius (1.0 is full space-filling)
        #[arg(long, default_value_t = 1.0)]
        radius_scale: f32,

        /// Disable the depth cue (distance fog)
        #[arg(long)]
        no_depth_cue: bool,

        /// Slice the structure to a depth band: two fractions FAR NEAR in 0..1
        /// (0 is the farthest atom, 1 the nearest)
        #[arg(long, num_args = 2, value_names = ["FAR", "NEAR"])]
        slab: Option<Vec<f32>>,

        /// Terminal graphics protocol for terminal output (ignored with -o)
        #[arg(long, value_enum, default_value_t = ProtocolArg::Auto)]
        protocol: ProtocolArg,

        #[command(flatten)]
        common: Common,
    },
    /// Search RCSB and download structures
    Fetch {
        /// Search terms; omit to open the search interface empty
        query: Option<String>,

        #[command(flatten)]
        common: Common,
    },
}

/// The interaction types `--type` accepts, named as the engine labels them.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum KindArg {
    HydrogenBond,
    SaltBridge,
    Hydrophobic,
    PiStacking,
    PiCation,
    HalogenBond,
    WaterBridge,
    MetalCoordination,
    Disulfide,
}

impl From<KindArg> for InteractionKind {
    fn from(arg: KindArg) -> Self {
        match arg {
            KindArg::HydrogenBond => InteractionKind::HydrogenBond,
            KindArg::SaltBridge => InteractionKind::SaltBridge,
            KindArg::Hydrophobic => InteractionKind::Hydrophobic,
            KindArg::PiStacking => InteractionKind::PiStacking,
            KindArg::PiCation => InteractionKind::PiCation,
            KindArg::HalogenBond => InteractionKind::HalogenBond,
            KindArg::WaterBridge => InteractionKind::WaterBridge,
            KindArg::MetalCoordination => InteractionKind::MetalCoordination,
            KindArg::Disulfide => InteractionKind::Disulfide,
        }
    }
}

/// How `render` colours atoms.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum ColorArg {
    /// CPK/Jmol colour by element
    Element,
    /// Blue to red gradient over the B-factor range
    Bfactor,
    /// Secondary-structure palette
    Ss,
}

impl From<ColorArg> for pixelfold_render::Coloring {
    fn from(arg: ColorArg) -> Self {
        match arg {
            ColorArg::Element => pixelfold_render::Coloring::Element,
            ColorArg::Bfactor => pixelfold_render::Coloring::Bfactor,
            ColorArg::Ss => pixelfold_render::Coloring::SecondaryStructure,
        }
    }
}

/// The terminal graphics protocol for `render`'s terminal output.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum ProtocolArg {
    /// Detect the terminal from the environment
    Auto,
    /// Kitty graphics protocol (also Ghostty, WezTerm)
    Kitty,
    /// iTerm2 inline images
    Iterm2,
    /// Unicode half-block fallback (any truecolor terminal)
    HalfBlock,
}

impl From<ProtocolArg> for render::TerminalTarget {
    fn from(arg: ProtocolArg) -> Self {
        match arg {
            ProtocolArg::Auto => render::TerminalTarget::Auto,
            ProtocolArg::Kitty => render::TerminalTarget::Kitty,
            ProtocolArg::Iterm2 => render::TerminalTarget::Iterm2,
            ProtocolArg::HalfBlock => render::TerminalTarget::HalfBlock,
        }
    }
}

/// Alternate-location handling policy for the CLI.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum AltlocArg {
    /// Keep the highest-occupancy conformer per atom (default)
    Occupancy,
    /// Keep only altloc A (and atoms with no altloc)
    A,
    /// Keep only altloc B
    B,
    /// Keep every conformer
    All,
}

impl From<AltlocArg> for pixelfold_core::AltlocPolicy {
    fn from(arg: AltlocArg) -> Self {
        match arg {
            AltlocArg::Occupancy => pixelfold_core::AltlocPolicy::Occupancy,
            AltlocArg::A => pixelfold_core::AltlocPolicy::A,
            AltlocArg::B => pixelfold_core::AltlocPolicy::B,
            AltlocArg::All => pixelfold_core::AltlocPolicy::All,
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Some(command) => run(command),
        None if cli.fetch => search(cli.structure, &cli.common),
        None => match cli.structure {
            Some(structure) => view(&structure, cli.no_surface, &cli.common),
            None => {
                Cli::command().print_help()?;
                println!();
                Ok(())
            }
        },
    }
}

fn run(command: Command) -> Result<()> {
    match command {
        Command::View {
            structure,
            no_surface,
            common,
        } => view(&structure, no_surface, &common),

        Command::Fetch { query, common } => search(query, &common),

        Command::Interactions {
            analysis,
            types,
            ligand,
        } => {
            let query = ligand
                .map(|code| format!("resn {code}"))
                .or(analysis.select);
            let protein = loaded(&analysis.structure, &analysis.common)?;
            let chosen = load::chosen(&protein, query.as_deref())?;
            let kinds: Vec<InteractionKind> = types.into_iter().map(Into::into).collect();

            emit(
                &analyse::interactions(&protein, &chosen, &kinds),
                analysis.format,
            )
        }

        Command::Rin {
            structure,
            select,
            types,
            format,
            analyze,
            path,
            mvsj,
            output,
            common,
        } => {
            let protein = loaded(&structure, &common)?;
            let chosen = load::chosen(&protein, select.as_deref())?;
            let kinds: Vec<InteractionKind> = types.into_iter().map(Into::into).collect();

            let interactions = analyse::selected(&protein, &chosen, &kinds);
            let network = pixelfold_core::rin::build(&protein, &interactions);

            if let Some(endpoints) = path {
                let node_of = |id: &str| {
                    network
                        .nodes
                        .iter()
                        .position(|node| node.id == id)
                        .with_context(|| {
                            format!("residue {id} is not in the network (does it interact?)")
                        })
                };
                let from = node_of(&endpoints[0])?;
                let to = node_of(&endpoints[1])?;
                let route = pixelfold_core::rin::shortest_path(&network, from, to);

                emit_to(output.as_deref(), |out| {
                    graph::write_path(
                        out,
                        &network,
                        &endpoints[0],
                        &endpoints[1],
                        route.as_deref(),
                    )
                })
            } else if mvsj {
                let (url, parse_format) = mvsj_source(&structure)?;
                let title = format!("PixelFold — {structure} interaction network");
                emit_to(output.as_deref(), |out| {
                    graph::write_mvsj(out, &network, &url, parse_format, &title)
                })
            } else if analyze {
                let report = pixelfold_core::rin::analyze(&network);
                emit_to(output.as_deref(), |out| {
                    graph::write_analysis(out, &network, &report)
                })
            } else if format == GraphFormat::Tsv
                && output.is_none()
                && std::io::stdout().is_terminal()
            {
                emit_to(None, |out| graph::write_tsv_aligned(out, &network))
            } else {
                emit_to(output.as_deref(), |out| graph::write(out, &network, format))
            }
        }

        Command::Ss { analysis } => {
            let protein = loaded(&analysis.structure, &analysis.common)?;
            let chosen = load::chosen(&protein, analysis.select.as_deref())?;

            emit(&analyse::secondary(&protein, &chosen), analysis.format)
        }

        Command::Sasa { analysis } => {
            let protein = loaded(&analysis.structure, &analysis.common)?;
            let chosen = load::chosen(&protein, analysis.select.as_deref())?;

            emit(&analyse::sasa(&protein, &chosen), analysis.format)
        }

        Command::Render {
            structure,
            output,
            width,
            height,
            color,
            backbone,
            select,
            radius_scale,
            no_depth_cue,
            slab,
            protocol,
            common,
        } => {
            let protein = loaded(&structure, &common)?;

            let selected = match select {
                Some(query) => Some(
                    load::chosen(&protein, Some(&query))?
                        .iter()
                        .collect::<Vec<usize>>(),
                ),
                None => None,
            };

            let slab = match slab.as_deref() {
                None => None,
                Some([far, near]) if (0.0..=1.0).contains(far) && far < near && *near <= 1.0 => {
                    Some((*far, *near))
                }
                Some(_) => {
                    anyhow::bail!("--slab takes two fractions FAR NEAR with 0 <= FAR < NEAR <= 1")
                }
            };

            let params = render::RenderParams {
                width,
                height,
                coloring: color.into(),
                backbone,
                radius_scale,
                depth_cue: !no_depth_cue,
                slab,
            };

            render::render(
                &protein,
                selected.as_deref(),
                &params,
                output.as_deref(),
                protocol.into(),
            )
        }
    }
}

fn search(query: Option<String>, common: &Common) -> Result<()> {
    let dir = cache_dir(common)?;
    std::fs::create_dir_all(&dir)?;

    pixelfold_tui::search(query, dir)
}

fn view(structure: &str, no_surface: bool, common: &Common) -> Result<()> {
    let dir = cache_dir(common)?;
    let path = load::structure_path(structure, &dir)?;

    pixelfold_tui::view(&path, no_surface, common.altloc.into())
}

fn loaded(structure: &str, common: &Common) -> Result<pixelfold_core::Protein> {
    let dir = cache_dir(common)?;

    load::analysis_structure(structure, &dir, common.altloc.into())
}

/// The URL and parse format a MolViewSpec scene should load the structure from.
///
/// A 4-character PDB id resolves to its canonical RCSB mmCIF URL, which any Mol*
/// can fetch. A local path becomes a `file://` URL, which a desktop or
/// locally-served Mol* can read.
fn mvsj_source(structure: &str) -> Result<(String, &'static str)> {
    let is_pdb_id = structure.len() == 4 && structure.chars().all(|c| c.is_ascii_alphanumeric());
    if is_pdb_id {
        return Ok((
            format!(
                "https://files.rcsb.org/download/{}.cif",
                structure.to_uppercase()
            ),
            "mmcif",
        ));
    }

    let path = Path::new(structure)
        .canonicalize()
        .with_context(|| format!("cannot resolve {structure} for the scene URL"))?;
    let format = match path.extension().and_then(|ext| ext.to_str()) {
        Some("pdb" | "ent") => "pdb",
        _ => "mmcif",
    };

    Ok((format!("file://{}", path.display()), format))
}

/// Write records to standard output.
///
/// Closed pipes are valid, e.g. `pixelfold interactions ... | head`.
fn emit<T>(records: &[T], format: Format) -> Result<()>
where
    T: report::Row + serde::Serialize,
{
    emit_to(None, |out| report::write(out, records, format))
}

/// Run `render` against a file when one is named, or standard output otherwise.
///
/// A closed pipe on standard output is valid (`... | head`); a write to a file
/// that fails is a real error and is reported.
fn emit_to(
    output: Option<&Path>,
    render: impl FnOnce(&mut dyn std::io::Write) -> Result<()>,
) -> Result<()> {
    match output {
        Some(path) => {
            let mut file = std::io::BufWriter::new(
                std::fs::File::create(path)
                    .with_context(|| format!("could not write {}", path.display()))?,
            );
            render(&mut file)?;
            file.flush()?;

            Ok(())
        }
        None => {
            let mut out = std::io::stdout().lock();
            match render(&mut out) {
                Err(error) if is_broken_pipe(&error) => Ok(()),
                other => other,
            }
        }
    }
}

/// Whether an error is the reader at the other end of the pipe.
///
/// Two spellings reach here: a write through `writeln!` fails as an
/// `io::Error`, while a write through serde fails as a `serde_json::Error`.
fn is_broken_pipe(error: &anyhow::Error) -> bool {
    let pipe = |kind| kind == std::io::ErrorKind::BrokenPipe;

    error
        .downcast_ref::<std::io::Error>()
        .is_some_and(|io| pipe(io.kind()))
        || error
            .downcast_ref::<serde_json::Error>()
            .and_then(serde_json::Error::io_error_kind)
            .is_some_and(pipe)
}

/// Where downloads are cached.
///
/// The directory is not created here.
fn cache_dir(common: &Common) -> Result<PathBuf> {
    match &common.cache_dir {
        Some(dir) => Ok(dir.clone()),
        None => Ok(dirs::cache_dir()
            .ok_or_else(|| {
                anyhow::anyhow!("could not determine a cache directory; pass --cache-dir")
            })?
            .join("pixelfold")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_command_line_is_well_formed() {
        Cli::command().debug_assert();
    }

    #[test]
    fn a_bare_structure_still_opens_the_viewer() {
        let cli = Cli::try_parse_from(["pixelfold", "4HHB"]).expect("parses");

        assert!(cli.command.is_none());
        assert_eq!(cli.structure.as_deref(), Some("4HHB"));
    }

    #[test]
    fn subcommands_take_their_own_structure() {
        let cli =
            Cli::try_parse_from(["pixelfold", "sasa", "1UBQ", "--format", "tsv"]).expect("parses");

        let Some(Command::Sasa { analysis }) = cli.command else {
            panic!("expected the sasa command");
        };
        assert_eq!(analysis.structure, "1UBQ");
        assert_eq!(analysis.format, Format::Tsv);
    }

    #[test]
    fn a_ligand_shorthand_and_a_selection_are_mutually_exclusive() {
        assert!(
            Cli::try_parse_from([
                "pixelfold",
                "interactions",
                "1IEP",
                "--ligand",
                "STI",
                "--select",
                "chain A",
            ])
            .is_err()
        );
    }

    #[test]
    fn interaction_types_repeat() {
        let cli = Cli::try_parse_from([
            "pixelfold",
            "interactions",
            "1IEP",
            "--type",
            "hydrogen-bond",
            "--type",
            "pi-stacking",
        ])
        .expect("parses");

        let Some(Command::Interactions { types, .. }) = cli.command else {
            panic!("expected the interactions command");
        };
        let kinds: Vec<InteractionKind> = types.into_iter().map(Into::into).collect();
        assert_eq!(
            kinds,
            vec![InteractionKind::HydrogenBond, InteractionKind::PiStacking]
        );
    }

    #[test]
    fn rin_takes_a_graph_format_and_an_output_file() {
        let cli = Cli::try_parse_from([
            "pixelfold",
            "rin",
            "4HHB",
            "--format",
            "graphml",
            "-o",
            "net.graphml",
        ])
        .expect("parses");

        let Some(Command::Rin {
            structure,
            format,
            output,
            ..
        }) = cli.command
        else {
            panic!("expected the rin command");
        };
        assert_eq!(structure, "4HHB");
        assert_eq!(format, GraphFormat::Graphml);
        assert_eq!(output.as_deref(), Some(Path::new("net.graphml")));
    }
}
