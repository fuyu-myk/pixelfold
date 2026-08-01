//! Struct-of-arrays scene the rasteriser consumes.
//!
//! Projection touches these arrays every frame, so they hold only the numbers
//! the rasteriser needs and never a `String`. `atom_idx` carries each sphere's
//! index in the source protein, which the ID buffer writes through for picking.

use glam::Vec3;

use pixelfold_core::{Protein, calculate_bfactor_range, get_calpha_indices, sasa::get_vdw_radius};

use crate::draw::{bfactor_to_color, element_to_color};

/// How to colour atoms.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Coloring {
    /// CPK/Jmol colour by element.
    #[default]
    Element,
    /// Blue to red gradient over the B-factor percentile range.
    Bfactor,
    /// Helix/sheet/turn/coil palette.
    SecondaryStructure,
}

/// Atoms flattened to parallel arrays for rasterisation.
pub struct Scene {
    pub positions: Vec<Vec3>,
    /// Van der Waals radii in Angstroms; the rasteriser scales to pixels.
    pub radii: Vec<f32>,
    pub colors: Vec<[u8; 3]>,
    pub atom_idx: Vec<u32>,
}

impl Scene {
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Build a scene from a protein.
    ///
    /// `indices` restricts the scene to those atoms (a selection, or the
    /// C-alpha trace); `None` renders every atom. B-factor colouring normalises
    /// over the whole protein so the palette is stable regardless of the subset.
    pub fn from_protein(protein: &Protein, coloring: Coloring, indices: Option<&[usize]>) -> Self {
        let (b_min, b_max) = match coloring {
            Coloring::Bfactor => calculate_bfactor_range(protein),
            _ => (0.0, 1.0),
        };

        let color_of = |idx: usize| -> [u8; 3] {
            let atom = &protein.atoms[idx];
            let (r, g, b) = match coloring {
                Coloring::Element => element_to_color(atom.element.as_str()),
                Coloring::Bfactor => bfactor_to_color(atom.b_factor, b_min, b_max),
                Coloring::SecondaryStructure => atom.secondary_structure.color_rgb(),
            };
            [r, g, b]
        };

        let build = |indices: &[usize]| {
            let mut positions = Vec::with_capacity(indices.len());
            let mut radii = Vec::with_capacity(indices.len());
            let mut colors = Vec::with_capacity(indices.len());
            let mut atom_idx = Vec::with_capacity(indices.len());
            for &idx in indices {
                let atom = &protein.atoms[idx];
                positions.push(atom.position);
                radii.push(get_vdw_radius(atom.element.as_str()));
                colors.push(color_of(idx));
                atom_idx.push(idx as u32);
            }

            Scene {
                positions,
                radii,
                colors,
                atom_idx,
            }
        };

        match indices {
            Some(indices) => build(indices),
            None => build(&(0..protein.atoms.len()).collect::<Vec<_>>()),
        }
    }

    /// C-alpha-only scene, for the backbone display mode.
    pub fn backbone(protein: &Protein, coloring: Coloring) -> Self {
        Self::from_protein(protein, coloring, Some(&get_calpha_indices(protein)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pixelfold_core::{Atom, FixedStr, SecondaryStructure};

    fn atom(element: &str, ss: SecondaryStructure, b: f32, pos: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new("CA"),
            element: FixedStr::new(element),
            residue_name: FixedStr::new("ALA"),
            residue_seq: 1,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: pos,
            b_factor: b,
            secondary_structure: ss,
        }
    }

    fn protein(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    #[test]
    fn element_coloring_uses_cpk_and_vdw_radius() {
        let p = protein(vec![
            atom("N", SecondaryStructure::Helix, 10.0, Vec3::ZERO),
            atom("O", SecondaryStructure::Coil, 20.0, Vec3::X),
        ]);
        let scene = Scene::from_protein(&p, Coloring::Element, None);

        assert_eq!(scene.len(), 2);
        assert_eq!(scene.colors[0], [48, 80, 248]); // nitrogen blue
        assert_eq!(scene.colors[1], [255, 13, 13]); // oxygen red
        assert_eq!(scene.radii[0], 1.55); // nitrogen vdW
        assert_eq!(scene.radii[1], 1.52); // oxygen vdW
        assert_eq!(scene.atom_idx, vec![0, 1]);
    }

    #[test]
    fn secondary_structure_coloring() {
        let p = protein(vec![atom("C", SecondaryStructure::Sheet, 0.0, Vec3::ZERO)]);
        let scene = Scene::from_protein(&p, Coloring::SecondaryStructure, None);
        assert_eq!(scene.colors[0], [100, 100, 255]); // sheet blue
    }

    #[test]
    fn indices_restrict_the_scene_but_keep_source_atom_ids() {
        let p = protein(vec![
            atom("C", SecondaryStructure::Coil, 0.0, Vec3::ZERO),
            atom("N", SecondaryStructure::Coil, 0.0, Vec3::X),
            atom("O", SecondaryStructure::Coil, 0.0, Vec3::Y),
        ]);
        let scene = Scene::from_protein(&p, Coloring::Element, Some(&[2]));

        assert_eq!(scene.len(), 1);
        assert_eq!(scene.atom_idx, vec![2], "ID buffer must map back to atom 2");
        assert_eq!(scene.colors[0], [255, 13, 13]); // the oxygen
    }

    #[test]
    fn bfactor_coloring_normalises_over_the_whole_protein() {
        // Colours must be identical whether or not the subset is passed, because
        // the palette is normalised over the full protein, not the subset.
        let p = protein(vec![
            atom("C", SecondaryStructure::Coil, 0.0, Vec3::ZERO),
            atom("C", SecondaryStructure::Coil, 50.0, Vec3::X),
            atom("C", SecondaryStructure::Coil, 100.0, Vec3::Y),
        ]);
        let full = Scene::from_protein(&p, Coloring::Bfactor, None);
        let subset = Scene::from_protein(&p, Coloring::Bfactor, Some(&[1]));
        assert_eq!(subset.colors[0], full.colors[1]);
    }
}
