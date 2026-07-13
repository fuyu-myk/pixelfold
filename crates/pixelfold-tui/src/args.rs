use std::env;

use anyhow::Result;

use crate::App;
use crate::app::run_app;

pub fn run() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    // Parse command-line arguments
    let mut protein_name: Option<String> = None;
    let mut skip_surface = false;

    for arg in args.iter().skip(1) {
        if arg == "--no-surface" {
            skip_surface = true;
        } else if arg == "--help" || arg == "-h" {
            println!("PixelFold - Terminal-based 3D protein structure viewer");
            println!();
            println!("--- Searching and fetching structures ---");
            println!("Usage: pixelfold [OPTIONS] <optional_protein_name>");
            println!();
            println!("Options (compulsory):");
            println!(
                "  --fetch, -f     Opens a search interface to fetch protein structures by name"
            );
            println!();
            println!("--- Visualization ---");
            println!("Usage: pixelfold <path/to/protein.pdb|protein.cif> [OPTIONS]");
            println!("       pixelfold <protein> [OPTIONS]");
            println!();
            println!("Options:");
            println!(
                "  --no-surface    Skip surface calculation for faster loading (large proteins)"
            );
            println!("  --help, -h      Show this help message");
            println!();
            println!("Controls:");
            println!("  WASDZX    Rotate structure");
            println!("  +/-       Zoom in/out");
            println!("  Arrows    Pan view (or adjust settings when mode active)");
            println!("  1         Show all atoms");
            println!("  2         Show backbone atoms");
            println!("  V         Toggle surface visualization");
            println!("    ↑↓      Adjust surface density when surface visible");
            println!("  C         Toggle backbone connections");
            println!("  B         Toggle B-factor coloring");
            println!("  H         Toggle hydrogen bond display");
            println!("    ↑↓      Adjust H-bond energy threshold when H-bonds visible");
            println!("  N         Toggle H-bond network analysis overlay");
            println!("  F         Frame the protein in view");
            println!("  I         Inspect mode");
            println!("    ↑↓      Cycle through nearby atoms");
            println!("  R         Toggle residue highlighting (in inspect mode)");
            println!("  Q         Quit");
            return Ok(());
        } else if !arg.starts_with('-') {
            protein_name = Some(arg.clone());
        }
    }

    let mut app = App::new();
    let mut terminal = ratatui::init();

    crossterm::execute!(std::io::stdout(), crossterm::event::EnableMouseCapture)?;

    let result = if args.iter().any(|arg| arg == "--fetch" || arg == "-f") {
        crate::search::fetch_structures(&mut terminal, protein_name)
    } else {
        if let Some(name) = protein_name {
            // Initial terminal size
            let size = terminal.size()?;
            let width = size.width as f32 * 2.0; // Braille canvas width
            let height = size.height as f32 * 4.0; // Braille canvas height

            let project_dir = env!("CARGO_MANIFEST_DIR");

            let path = if name.starts_with("data") {
                format!("{}/{}", project_dir, name)
            } else {
                if name.ends_with(".cif") {
                    format!("{}/data/{}", project_dir, name)
                } else if name.ends_with(".pdb") {
                    format!("{}/data/{}", project_dir, name)
                } else {
                    format!("{}/data/{}.cif", project_dir, name)
                }
            };

            match app.load_protein(&path, width, height, skip_surface) {
                Ok(_) => {}
                Err(e) => eprintln!("Failed to load protein: {}", e),
            }
        }

        run_app(&mut terminal, &mut app)
    };

    crossterm::execute!(std::io::stdout(), crossterm::event::DisableMouseCapture)?;

    ratatui::restore();
    result
}
