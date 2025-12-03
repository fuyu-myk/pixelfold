# PixelFold

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A terminal-based 3D protein structure viewer using Braille Unicode characters for high-resolution visualization.

## Features

- **Braille Rendering**: Uses Unicode Braille characters (2×4 pixel resolution per character cell) for detailed protein structure visualization
- **PDB & mmCIF Support**: Parses both PDB and mmCIF/PDBx protein structure formats using `pdbtbx`
- **Interactive 3D Controls**: Rotate, zoom, and pan the protein structure in real-time
- **Performance Modes**: Multiple loading options for different protein sizes

## Installation

Clone this repository and run the following command:

```bash
cargo install --path .
```

The binary will be available at `target/release/pixelfold`.

## Usage

```bash
# Load and visualize a protein structure
pixelfold path/to/protein.pdb

# Or with mmCIF format
pixelfold path/to/protein.cif
```

### Controls

- **WASD**: Rotate the structure (W/S = pitch, A/D = yaw)
- **Z/X**: Roll rotation
- **+/-**: Zoom in/out
- **Arrow Keys**: Pan the view
- **1/2**: Toggle display modes (1: All atoms (default); 2: Alpha carbon backbone)
- **F**: Auto-frame (reset and fit protein to view)
- **I**: Inspect-mode (interactive clicking enabled to view more information)
- **R**: Toggle highlight amino acid (in inspect mode only)
- **C**: Toggle backbone connections
- **B**: Toggle b-factor coloring
- **Q**: Quit

## Implementation Details

### Architecture

The project consists of several modules:

- **`parser`**: Loads PDB/mmCIF files using `pdbtbx` and converts to internal data structures
  - `load_protein()`: Load all atoms
  - `load_protein_backbone()`: Load only backbone atoms (CA, C, N, O)
  - `load_protein_ca_only()`: Load only alpha carbons for large proteins

- **`renderer`**: Handles 3D-to-2D projection and camera controls
  - Orthographic projection with rotation/zoom/pan
  - Depth sorting for proper atom occlusion
  - Camera utilities (auto-framing, bounds calculation)

- **`main`**: TUI application using `ratatui` with Braille canvas rendering

### Why Braille?

Braille Unicode characters provide **8× higher resolution** compared to ASCII:

- Each terminal character cell can display 2×4 sub-pixels (8 dots)
- Excellent for visualizing protein density and structure detail
- Native support in `ratatui::widgets::canvas::Canvas` with `Marker::Braille`

### Dependencies

- **`ratatui`**: Terminal UI framework
- **`crossterm`**: Terminal manipulation
- **`pdbtbx`**: PDB and mmCIF parsing
- **`glam`**: 3D math library (vectors, matrices)
- **`anyhow`**: Error handling

### Current limitations

**Color bleeding**: Due to how braille is rendered and colored, the colors of atoms in close proximity to each other may bleed into surrounding ones

- This is especially prevalent in inspect mode, where the atom is marked
- One workaround is to zoom in as much as possible to clearly distinguish which atom is marked

## Example Workflow

1. Download a protein structure from [RCSB PDB](https://www.rcsb.org/)
2. Run PixelFold with the file: `pixelfold 1crn.pdb`
3. Use WASD to rotate and explore the structure
4. Press F to auto-frame the view
5. Use +/- to zoom in on specific regions

## Future TODO

- [x] Bond rendering between atoms (alpha carbons)
- [] Multiple selection and filtering modes
- [] Save/load camera positions
- [] Animation and structure comparison
