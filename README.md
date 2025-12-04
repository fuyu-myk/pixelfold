# PixelFold

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A terminal-based 3D protein structure viewer using Braille Unicode characters for high-resolution visualization.

## Features

- **Braille Rendering**: Uses Unicode Braille characters (2×4 pixel resolution per character cell) for detailed protein structure visualization
- **PDB & mmCIF Support**: Parses both PDB and mmCIF/PDBx protein structure formats using `pdbtbx`
- **Interactive 3D Controls**: Rotate, zoom, and pan the protein structure in real-time
- **Surface Visualization**: Solvent-accessible surface with Kyte-Doolittle hydrophobicity coloring using the Shrake-Rupley algorithm
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
- **Arrow Keys**: Pan the view (or cycle through candidate atoms in inspect mode, or adjust surface density when surface is visible)
- **1/2**: Toggle display modes (1: All atoms (default); 2: Alpha carbon backbone)
- **F**: Auto-frame (reset and fit protein to view)
- **I**: Inspect-mode (interactive clicking enabled to view more information)
  - **Up/Down arrows** (upon clicking on an atom): Cycle between closest 5 atoms
- **R**: Toggle highlight amino acid (in inspect mode only)
- **C**: Toggle backbone connections
- **B**: Toggle b-factor coloring
- **V**: Toggle solvent-accessible surface visualization
  - **Note**: When surface is visible, atom rendering is disabled to show clear hydrophobicity patterns
  - **Up/Down arrows** (when surface visible): Adjust point density (100-500 points/atom)
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

- **`surface`**: Solvent-accessible surface calculation using Shrake-Rupley algorithm
  - Powered by `rust-sasa` library for optimized performance
  - Fibonacci spiral sphere point generation for uniform sampling
  - Van der Waals radii lookup for common atoms
  - Solvent accessibility testing with 1.4Å water probe radius
  - Kyte-Doolittle hydrophobicity scale for surface coloring
  - Parallel computation for large proteins

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

**Surface computation performance**: For very large proteins (>5000 atoms), surface calculation may take several seconds on load

- Surface computation uses the optimized `rust-sasa` library with parallel computation
- Use `--no-surface` flag to skip surface calculation and compute on-demand with 'V' key
- Surface points are computed once and cached for the session

## Example Workflow

1. Download a protein structure from [RCSB PDB](https://www.rcsb.org/)
2. Run PixelFold with the file: `pixelfold 1crn.pdb`
3. Use WASD to rotate and explore the structure
4. Press F to auto-frame the view
5. Press V to toggle surface visualization
6. Use up/down arrows to adjust surface point density
7. Use +/- to zoom in on specific regions

## Future TODO

- [x] Bond rendering between atoms (alpha carbons)
- [x] Solvent-accessible surface visualization with hydrophobicity coloring
- [] Multiple selection and filtering modes
- [] Save/load camera positions
- [] Animation and structure comparison
- [] Electrostatic potential surface coloring (via APBS integration)

## Citations

```bibtex
@article{Kabsch1983,
  author = {Kabsch, Wolfgang and Sander, Christian},
  title = {Dictionary of protein secondary structure: Pattern recognition of hydrogen-bonded and geometrical features},
  journal = {Biopolymers},
  volume = {22},
  number = {12},
  pages = {2577--2637},
  year = {1983},
  doi = {10.1002/bip.360221211}
}

@article{Kyte1982,
  author = {Kyte, Jack and Doolittle, Russell F.},
  title = {A simple method for displaying the hydropathic character of a protein},
  journal = {Journal of Molecular Biology},
  volume = {157},
  number = {1},
  pages = {105--132},
  year = {1982},
  doi = {10.1016/0022-2836(82)90515-0}
}

@article{Bondi1964,
  author = {Bondi, A.},
  title = {van der Waals Volumes and Radii},
  journal = {The Journal of Physical Chemistry},
  volume = {68},
  number = {3},
  pages = {441--451},
  year = {1964},
  doi = {10.1021/j100785a001}
}

@article{Shrake1973,
  author = {Shrake, A. and Rupley, J. A.},
  title = {Environment and exposure to solvent of protein atoms. Lysozyme and insulin},
  journal = {Journal of Molecular Biology},
  volume = {79},
  number = {2},
  pages = {351--371},
  year = {1973},
  doi = {10.1016/0022-2836(73)90011-9}
}
```
