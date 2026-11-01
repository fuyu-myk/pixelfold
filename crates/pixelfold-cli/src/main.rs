use std::io::{IsTerminal, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use pixelfold_core::parser::{AssemblyChoice, Load};
use pixelfold_core::rin::{self, Ensemble};
use pixelfold_core::{AtomSet, InteractionKind, Network, Protein};

mod analyse;
mod graph;
mod load;
mod poses;
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

    /// Terminal graphics protocol for the structure pane. Detection is
    /// automatic; force half-block if the pane stays blank
    #[arg(long, value_enum, default_value_t = ProtocolArg::Auto)]
    protocol: ProtocolArg,

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

    /// Which model of a multi-model file (an NMR ensemble, say) to load: its
    /// number, or `all` for ss, sasa and interactions to report on every model
    /// and rin to build one network over them (default: the first)
    #[arg(long, value_name = "N|all", value_parser = model_choice, global = true)]
    model: Option<ModelChoice>,

    /// Build a biological assembly instead of the deposited coordinates: the
    /// primary one with a bare --assembly, another with --assembly=ID. Copies of
    /// a chain are named A, A-2, A-3, ...
    #[arg(
        long,
        value_name = "ID",
        num_args = 0..=1,
        require_equals = true,
        global = true
    )]
    assembly: Option<Option<String>>,
}

impl Common {
    /// What to build: model `model` and the assembly `--assembly` names.
    fn load(&self, model: Option<usize>) -> Load {
        Load {
            model,
            assembly: self.assembly.as_ref().map(|id| match id {
                Some(id) => AssemblyChoice::Id(id.clone()),
                None => AssemblyChoice::Primary,
            }),
            altloc: self.altloc.into(),
        }
    }

    /// The model a command that works on one model at a time loads.
    fn one_model(&self) -> Result<Option<usize>> {
        match self.model {
            Some(ModelChoice::All) => {
                anyhow::bail!("this command works on one model at a time; pass --model N")
            }
            Some(ModelChoice::One(number)) => Ok(Some(number)),
            None => Ok(None),
        }
    }
}

/// Which model of a multi-model file `--model` asks for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ModelChoice {
    One(usize),
    All,
}

fn fraction(value: &str) -> Result<f32, String> {
    value
        .parse()
        .ok()
        .filter(|fraction| (0.0..=1.0).contains(fraction))
        .ok_or_else(|| format!("expected a fraction between 0 and 1, not `{value}`"))
}

fn model_choice(value: &str) -> Result<ModelChoice, String> {
    if value.eq_ignore_ascii_case("all") {
        return Ok(ModelChoice::All);
    }

    value
        .parse()
        .map(ModelChoice::One)
        .map_err(|_| format!("expected a model number or `all`, not `{value}`"))
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

        /// Terminal graphics protocol for the structure pane. Detection is
        /// automatic; force half-block if the pane stays blank
        #[arg(long, value_enum, default_value_t = ProtocolArg::Auto)]
        protocol: ProtocolArg,

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
    /// Build the residue interaction network.
    /// With --model all, the network of the whole ensemble: every edge any model
    /// has, weighted by the fraction of models it is present in
    Rin {
        /// Path to a .pdb/.cif file, or a 4-character PDB id
        structure: String,

        /// Restrict the network to interactions touching this selection
        #[arg(short, long, value_name = "QUERY")]
        select: Option<String>,

        /// Include only these interaction types (repeatable)
        #[arg(short = 't', long = "type", value_enum, value_name = "TYPE")]
        types: Vec<KindArg>,

        /// With --model all, keep only the edges present in at least this
        /// fraction of the models, from 0 to 1 (default: every edge). The
        /// fraction compared is exact, an edge's models over the network's, not
        /// the frequency as rounded for display
        #[arg(long, value_name = "FRACTION", value_parser = fraction)]
        min_frequency: Option<f32>,

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
    /// Rank ligand poses by the reference contacts they keep, with their
    /// pocket-aligned RMSD, their clashes and the reference stereocentres they
    /// invert.
    /// Every model of each pose file is a pose (or only --model N); the
    /// reference is its first model. Residues are matched with the reference's
    /// by aligning chain sequences, so a pose may number them as it likes. A
    /// pose holding several copies of the ligand is scored by the copy most
    /// like the reference's, named in the copy column
    Poses {
        /// The pose files: .pdb/.cif paths, 4-character PDB ids, or .sdf files
        /// of docked ligands, each record a pose
        #[arg(required = true, value_name = "POSE")]
        poses: Vec<String>,

        /// The structure whose ligand contacts the poses are compared with
        #[arg(short, long, value_name = "STRUCTURE")]
        reference: String,

        /// Which atoms are the ligand, in the reference and in every pose
        #[arg(short, long, value_name = "QUERY", default_value = "ligand")]
        select: String,

        /// Shorthand for --select 'resn CODE'
        #[arg(short, long, value_name = "CODE", conflicts_with = "select")]
        ligand: Option<String>,

        /// Which atoms are the ligand in the reference, when they are not the
        /// ones --select or --ligand picks there, such as one copy of a ligand
        /// the reference holds several of
        #[arg(long, value_name = "QUERY")]
        reference_select: Option<String>,

        /// The protein the .sdf poses were docked into, its own copies of the
        /// reference ligand left out [default: the reference without its
        /// ligand]
        #[arg(long, value_name = "STRUCTURE")]
        receptor: Option<String>,

        /// How to write the result
        #[arg(long, value_enum, default_value_t = Format::Table)]
        format: Format,

        /// Open the viewer on the ranked poses instead: a table to sort and
        /// step through, each pose drawn on the reference's pocket
        #[arg(long, conflicts_with = "format")]
        view: bool,

        /// Terminal graphics protocol for the viewer's structure pane
        #[arg(long, value_enum, requires = "view")]
        protocol: Option<ProtocolArg>,

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
    PiHydrogenBond,
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
            KindArg::PiHydrogenBond => InteractionKind::PiHydrogenBond,
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

impl From<ProtocolArg> for pixelfold_tui::ViewProtocol {
    fn from(arg: ProtocolArg) -> Self {
        match arg {
            ProtocolArg::Auto => pixelfold_tui::ViewProtocol::Auto,
            ProtocolArg::Kitty => pixelfold_tui::ViewProtocol::Kitty,
            ProtocolArg::Iterm2 => pixelfold_tui::ViewProtocol::Iterm2,
            ProtocolArg::HalfBlock => pixelfold_tui::ViewProtocol::HalfBlock,
        }
    }
}

/// Alternate-location handling policy for the CLI.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum AltlocArg {
    /// Keep one conformer per residue, the most occupied (default)
    Occupancy,
    /// Keep only altloc A (and atoms with no altloc)
    A,
    /// Keep only altloc B (and atoms with no altloc)
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
            Some(structure) => view(&structure, cli.protocol, &cli.common),
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
            protocol,
            common,
        } => view(&structure, protocol, &common),

        Command::Fetch { query, common } => search(query, &common),

        Command::Interactions {
            analysis,
            types,
            ligand,
        } => {
            let query = ligand
                .map(|code| format!("resn {code}"))
                .or_else(|| analysis.select.clone());
            let kinds: Vec<InteractionKind> = types.into_iter().map(Into::into).collect();

            emit_analysis(&analysis, query.as_deref(), |protein, chosen| {
                analyse::interactions(protein, chosen, &kinds)
            })
        }

        Command::Rin {
            structure,
            select,
            types,
            min_frequency,
            format,
            analyze,
            path,
            mvsj,
            output,
            common,
        } => {
            if mvsj && common.assembly.is_some() {
                anyhow::bail!("--mvsj cannot show a biological assembly yet; drop --assembly");
            }
            let kinds: Vec<InteractionKind> = types.into_iter().map(Into::into).collect();
            let network = network_of(
                &structure,
                select.as_deref(),
                &kinds,
                min_frequency,
                &common,
            )?;

            if let Some(endpoints) = path {
                let node_of = |id: &str| {
                    network
                        .nodes
                        .iter()
                        .position(|node| node.id == id)
                        .with_context(|| {
                            let often = min_frequency.map_or(String::new(), |fraction| {
                                format!(" in at least {fraction} of the models")
                            });
                            format!("residue {id} is not in the network (does it interact{often}?)")
                        })
                };
                let from = node_of(&endpoints[0])?;
                let to = node_of(&endpoints[1])?;
                let route = rin::shortest_path(&network, from, to);

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
                let report = rin::analyze(&network);
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

        Command::Poses {
            poses,
            reference,
            select,
            ligand,
            reference_select,
            receptor,
            format,
            view,
            protocol,
            common,
        } => {
            let query = ligand.map_or(select, |code| format!("resn {code}"));
            let reference_query = reference_select.unwrap_or_else(|| query.clone());
            let shown = view.then(|| protocol.unwrap_or(ProtocolArg::Auto));

            compare_poses(
                &reference,
                &reference_query,
                receptor.as_deref(),
                &poses,
                &query,
                shown.map_or(Output::Format(format), Output::View),
                &common,
            )
        }

        Command::Ss { analysis } => {
            emit_analysis(&analysis, analysis.select.as_deref(), analyse::secondary)
        }

        Command::Sasa { analysis } => {
            emit_analysis(&analysis, analysis.select.as_deref(), analyse::sasa)
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

fn view(structure: &str, protocol: ProtocolArg, common: &Common) -> Result<()> {
    let model = common.one_model()?;
    let dir = cache_dir(common)?;
    let path = load::structure_path(structure, &dir)?;

    pixelfold_tui::view(&path, common.load(model), protocol.into())
}

fn loaded(structure: &str, common: &Common) -> Result<Protein> {
    let model = common.one_model()?;
    let dir = cache_dir(common)?;

    load::analysis_structure(structure, &dir, &common.load(model))
}

/// Write `analysis` for the model `--model` names, or for every model, with a
/// `model` column, when it names `all`. Under `all` the file is parsed once and
/// each model built in turn, so an ensemble is never held whole, and a
/// selection is an error only when it matches nothing in any model.
fn emit_analysis<R, F>(analysis: &Analysis, query: Option<&str>, report: F) -> Result<()>
where
    R: report::Row + serde::Serialize,
    F: Fn(&Protein, &AtomSet) -> Vec<R>,
{
    if analysis.common.model != Some(ModelChoice::All) {
        let protein = loaded(&analysis.structure, &analysis.common)?;
        let chosen = load::chosen(&protein, query)?;
        return emit(&report(&protein, &chosen), analysis.format);
    }

    let dir = cache_dir(&analysis.common)?;
    let parsed = load::parsed(&analysis.structure, &dir)?;
    let mut records = Vec::new();
    let mut matched = false;
    for model in parsed.model_numbers() {
        let protein = parsed.protein(&analysis.common.load(Some(model)))?;
        let chosen = load::selection(&protein, query)?;
        matched |= !chosen.is_empty();
        records.extend(
            report(&protein, &chosen)
                .into_iter()
                .map(|record| report::InModel { model, record }),
        );
    }
    if !matched {
        return Err(matched_no_model(query));
    }

    emit(&records, analysis.format)
}

/// The residue interaction network of the `kinds` of interaction touching what
/// `query` selects: in the model `--model` names, or, when it names `all`, in
/// every model folded into one, keeping the edges present in at least
/// `min_frequency` of them.
fn network_of(
    structure: &str,
    query: Option<&str>,
    kinds: &[InteractionKind],
    min_frequency: Option<f32>,
    common: &Common,
) -> Result<Network> {
    if common.model != Some(ModelChoice::All) {
        if min_frequency.is_some() {
            anyhow::bail!("--min-frequency weighs edges across models; pass --model all");
        }
        let protein = loaded(structure, common)?;
        let chosen = load::chosen(&protein, query)?;
        return Ok(rin::build(
            &protein,
            &analyse::selected(&protein, &chosen, kinds),
        ));
    }

    let dir = cache_dir(common)?;
    let parsed = load::parsed(structure, &dir)?;
    let mut ensemble = Ensemble::default();
    let mut matched = false;
    for model in parsed.model_numbers() {
        let protein = parsed
            .protein(&common.load(Some(model)))
            .with_context(|| format!("model {model}"))?;
        let chosen = load::selection(&protein, query)?;
        matched |= !chosen.is_empty();
        ensemble.add(&protein, &analyse::selected(&protein, &chosen, kinds));
    }
    if !matched {
        return Err(matched_no_model(query));
    }

    Ok(ensemble.network(min_frequency.unwrap_or(0.0)))
}

fn matched_no_model(query: Option<&str>) -> anyhow::Error {
    anyhow::anyhow!(
        "the selection `{}` matched no atoms in any model",
        query.unwrap_or("all")
    )
}

/// Where ranked poses go: written out, or listed in the viewer.
enum Output {
    Format(Format),
    View(ProtocolArg),
}

/// Score the `query` ligand of every pose against the `reference_query` ligand
/// of the reference, best first. Each pose file is parsed once and its models
/// built in turn; each record of an SD file is a pose, its ligand laid into
/// `receptor`, or into the reference without its ligand.
fn compare_poses(
    reference: &str,
    reference_query: &str,
    receptor: Option<&str>,
    poses: &[String],
    query: &str,
    output: Output,
    common: &Common,
) -> Result<()> {
    let docking = poses
        .iter()
        .any(|pose| pixelfold_core::sdf::is_sdf(Path::new(pose)));
    anyhow::ensure!(
        docking || receptor.is_none(),
        "--receptor is the protein .sdf poses were docked into, and no pose is an .sdf file"
    );
    let dir = cache_dir(common)?;
    let mut warnings = poses::Warnings::default();
    let structure = load::analysis_structure(reference, &dir, &common.load(None))
        .with_context(|| format!("the reference {reference}"))?;
    let scoring = poses::Reference::new(&structure, reference_query, &mut warnings)
        .with_context(|| format!("the reference {reference}"))?;
    let receptor = match (docking, receptor) {
        (false, _) => None,
        (true, Some(given)) => Some(
            load::analysis_structure(given, &dir, &common.load(None))
                .map(|protein| scoring.receptor(protein, &mut warnings))
                .with_context(|| format!("the receptor {given}"))?,
        ),
        (true, None) => Some(poses::reference_receptor(&structure, reference_query)?),
    };
    let chosen = match common.model {
        Some(ModelChoice::One(number)) => Some(number),
        _ => None,
    };

    let mut records = Vec::new();
    for pose in poses {
        if let Some(receptor) = receptor
            .as_ref()
            .filter(|_| pixelfold_core::sdf::is_sdf(Path::new(pose)))
        {
            records.extend(poses::docked(
                &scoring,
                receptor,
                pose,
                chosen,
                &mut warnings,
            )?);
            continue;
        }
        let parsed = load::parsed(pose, &dir)?;
        let models = match chosen {
            Some(number) => vec![number],
            None => parsed.model_numbers(),
        };
        for model in models {
            let record = parsed
                .protein(&common.load(Some(model)))
                .and_then(|protein| scoring.assess(protein, query, pose, model, &mut warnings))
                .with_context(|| format!("{pose}, model {model}"))?;
            records.push(record);
        }
    }

    let ranked = poses::ranked(records);
    match output {
        Output::Format(format) => emit(&ranked, format),
        Output::View(protocol) => {
            anyhow::ensure!(!ranked.is_empty(), "no pose file holds a model to view");
            pixelfold_tui::view_poses(
                poses::viewed(ranked, |pose| load::structure_path(pose, &dir))?,
                receptor,
                common.load(None),
                protocol.into(),
            )
        }
    }
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

    #[test]
    fn poses_take_a_reference_and_any_number_of_pose_files() {
        let cli = Cli::try_parse_from([
            "pixelfold",
            "poses",
            "--reference",
            "ref.cif",
            "a.cif",
            "b.cif",
        ])
        .expect("parses");

        let Some(Command::Poses {
            poses,
            reference,
            select,
            ligand,
            reference_select,
            ..
        }) = cli.command
        else {
            panic!("expected the poses command");
        };
        assert_eq!(poses, ["a.cif", "b.cif"]);
        assert_eq!(reference, "ref.cif");
        assert_eq!(select, "ligand");
        assert_eq!(ligand, None);
        assert_eq!(reference_select, None);

        let apart = [
            "pixelfold",
            "poses",
            "-r",
            "ref.cif",
            "--ligand",
            "STI",
            "--reference-select",
            "chain A and resn STI",
            "a.cif",
        ];
        let Some(Command::Poses {
            reference_select, ..
        }) = Cli::try_parse_from(apart).expect("parses").command
        else {
            panic!("expected the poses command");
        };
        assert_eq!(reference_select.as_deref(), Some("chain A and resn STI"));

        let no_poses = ["pixelfold", "poses", "--reference", "ref.cif"];
        assert!(Cli::try_parse_from(no_poses).is_err());
        let both = [
            "pixelfold",
            "poses",
            "-r",
            "ref.cif",
            "--ligand",
            "STI",
            "--select",
            "chain B",
            "a.cif",
        ];
        assert!(Cli::try_parse_from(both).is_err());
    }
}
