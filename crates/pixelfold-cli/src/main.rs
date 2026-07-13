use std::path::PathBuf;

use anyhow::Result;
use clap::{CommandFactory, Parser};

mod resolve;

/// Terminal 3D protein structure viewer.
#[derive(Parser)]
#[command(
    name = "pixelfold",
    version,
    about = "Terminal 3D protein structure viewer"
)]
struct Cli {
    /// Path to a .pdb/.cif file, or a 4-character PDB id (for example 4HHB)
    structure: Option<String>,

    /// Open the search interface to fetch structures from RCSB
    #[arg(short, long)]
    fetch: bool,

    /// Skip surface calculation for faster loading of large structures
    #[arg(long)]
    no_surface: bool,

    /// Directory for cached downloads (default: the system cache dir under pixelfold/)
    #[arg(long, value_name = "DIR")]
    cache_dir: Option<PathBuf>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let cache_dir = match cli.cache_dir {
        Some(dir) => dir,
        None => default_cache_dir()?,
    };
    std::fs::create_dir_all(&cache_dir)?;

    if cli.fetch {
        return pixelfold_tui::search(cli.structure, cache_dir);
    }

    let Some(input) = cli.structure else {
        Cli::command().print_help()?;
        println!();
        return Ok(());
    };

    let path = resolve::resolve(&input, &cache_dir)?;
    pixelfold_tui::view(&path, cli.no_surface)
}

/// The default cache directory: the system cache dir with a `pixelfold` subdirectory.
fn default_cache_dir() -> Result<PathBuf> {
    let base = dirs::cache_dir().ok_or_else(|| {
        anyhow::anyhow!("could not determine a cache directory; pass --cache-dir")
    })?;

    Ok(base.join("pixelfold"))
}
