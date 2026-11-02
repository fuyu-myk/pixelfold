//! How near a ligand sits to the rest of its structure, as PoseBusters
//! measures a docked pose (Buttenschoen, Morris and Deane, Chem. Sci.
//! 15:3130, 2024): its clashes with the protein and organic cofactors, with
//! metal and inorganic ions, and with waters, each a pair of heavy atoms
//! closer than 0.75 of their radii summed, and whether any protein lies within
//! 5 A of it.
//!
//! Ions are measured on covalent radii, since a metal's coordination bonds sit
//! well inside the van der Waals sum, and everything else on van der Waals
//! radii, both as RDKit gives them. PoseBusters takes for ions 17 elements
//! (chlorine and bromine among them, in any residue) and eight ion residues
//! (sulfate, phosphate, iron-sulfur clusters and the like); here every metal
//! the interaction detector knows is an ion too, so a mercury or platinum bond
//! is no clash.

use glam::Vec3;

use crate::fixed_str::FixedStr;
use crate::geometry::vdw_radius;
use crate::interactions::classify::{is_metal, is_water};
use crate::pose::{Attachments, CLASH_RATIO, CLASH_REACH, attached_pairs, is_hydrogen};
use crate::select::AtomSet;
use crate::spatial::Grid;
use crate::structure::{Atom, Protein};

/// The farthest the protein may lie from a ligand that binds it, PoseBusters'
/// maximum distance.
pub const MAX_PROTEIN_DISTANCE: f32 = 5.0;

/// The pairs that clash by each of PoseBusters' distance checks, and how near
/// the protein comes.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct Distances {
    /// Pairs with protein or organic cofactor atoms closer than 0.75 of their
    /// van der Waals radii summed.
    pub clashes: usize,
    /// Pairs with metal or inorganic ion atoms closer than 0.75 of their
    /// covalent radii summed.
    pub metal_clashes: usize,
    /// Pairs with water oxygens closer than 0.75 of their van der Waals radii
    /// summed.
    pub water_clashes: usize,
    /// The shortest distance from a ligand heavy atom to a polymer one, in
    /// angstroms; `None` when the structure holds no polymer.
    pub protein_distance: Option<f32>,
}

impl Distances {
    /// Whether no polymer atom lies within [`MAX_PROTEIN_DISTANCE`] of the
    /// ligand.
    pub fn too_far(&self) -> bool {
        self.protein_distance
            .is_none_or(|distance| distance > MAX_PROTEIN_DISTANCE)
    }
}

/// Residues PoseBusters takes for inorganic ions whatever their elements.
const ION_RESIDUES: &[&str] = &["FES", "MOS", "PO3", "PO4", "PPK", "SO3", "SO4", "VO4"];

/// Elements PoseBusters takes for inorganic ions in any residue, besides
/// metals.
const ION_ELEMENTS: &[&str] = &["CL", "BR"];

/// Metals the interaction detector's list leaves out: beryllium and
/// molybdenum, which PoseBusters counts, and vanadium, which it counts only
/// as vanadate.
const MORE_METALS: &[&str] = &["BE", "MO", "V"];

/// What a heavy atom beside the ligand is to PoseBusters.
#[derive(Clone, Copy)]
enum Partner {
    /// The protein, or an organic cofactor.
    Organic,
    Ion,
    Water,
}

/// How near the `ligand` atoms of `protein` sit to the rest of it, the pairs
/// `attached` names left out of the clashes. Hydrogens never count, and a
/// metal of the ligand's own clashes with nothing.
pub fn distances(protein: &Protein, ligand: &AtomSet, attached: &Attachments) -> Distances {
    let heavy: Vec<usize> = ligand
        .iter()
        .filter(|&atom| !is_hydrogen(&protein.atoms[atom]))
        .collect();
    let excluded = attached_pairs(protein, ligand, attached);
    let positions: Vec<Vec3> = protein.atoms.iter().map(|atom| atom.position).collect();
    let grid = Grid::build(&positions, CLASH_REACH);

    let mut found = Distances::default();
    for &atom in heavy
        .iter()
        .filter(|&&atom| !is_any_metal(&protein.atoms[atom]))
    {
        let at = &protein.atoms[atom];
        grid.for_each_within(positions[atom], CLASH_REACH, |other| {
            let partner = &protein.atoms[other];
            if ligand.contains(other) || is_hydrogen(partner) || excluded.contains(&(atom, other)) {
                return;
            }

            let (radii, count) = match partner_kind(partner) {
                Partner::Organic => (vdw(at) + vdw(partner), &mut found.clashes),
                Partner::Ion => (covalent(at) + covalent(partner), &mut found.metal_clashes),
                Partner::Water => (vdw(at) + vdw(partner), &mut found.water_clashes),
            };
            if positions[atom].distance(positions[other]) < CLASH_RATIO * radii {
                *count += 1;
            }
        });
    }

    let polymer: Vec<Vec3> = protein
        .atoms
        .iter()
        .enumerate()
        .filter(|&(index, atom)| {
            atom.is_polymer() && !ligand.contains(index) && !is_hydrogen(atom) && !is_ion(atom)
        })
        .map(|(_, atom)| atom.position)
        .collect();
    found.protein_distance = heavy
        .iter()
        .flat_map(|&atom| {
            let from = positions[atom];
            polymer.iter().map(move |at| at.distance(from))
        })
        .min_by(f32::total_cmp);

    found
}

fn partner_kind(atom: &Atom) -> Partner {
    if is_water(atom.residue_name.as_str()) {
        Partner::Water
    } else if is_ion(atom) {
        Partner::Ion
    } else {
        Partner::Organic
    }
}

fn is_ion(atom: &Atom) -> bool {
    is_any_metal(atom)
        || ION_RESIDUES.contains(&atom.residue_name.as_str())
        || ION_ELEMENTS.contains(&element(atom).as_str())
}

fn is_any_metal(atom: &Atom) -> bool {
    is_metal(atom) || MORE_METALS.contains(&element(atom).as_str())
}

fn element(atom: &Atom) -> FixedStr<2> {
    FixedStr::uppercased(atom.element.as_str())
}

fn vdw(atom: &Atom) -> f32 {
    vdw_radius(element(atom))
}

/// An atom's covalent radius as RDKit gives it (Cordero et al., Dalton Trans.
/// 2008, 2832), in angstroms.
fn covalent(atom: &Atom) -> f32 {
    match element(atom).as_str() {
        "B" => 0.84,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "SI" => 1.11,
        "P" => 1.07,
        "S" => 1.05,
        "CL" => 1.02,
        "BR" | "SE" | "GE" => 1.2,
        "TE" => 1.38,
        "I" => 1.39,
        "AS" => 1.19,
        "LI" => 1.28,
        "BE" => 0.96,
        "NA" => 1.66,
        "MG" => 1.41,
        "AL" => 1.21,
        "K" => 2.03,
        "CA" => 1.76,
        "V" => 1.52,
        "CR" | "MN" | "PD" | "SB" => 1.39,
        "FE" | "CU" | "HG" => 1.32,
        "CO" => 1.26,
        "NI" => 1.24,
        "ZN" | "GA" => 1.22,
        "RB" => 2.2,
        "SR" => 1.95,
        "MO" => 1.54,
        "RU" => 1.46,
        "RH" | "IN" => 1.42,
        "IR" => 1.41,
        "AG" | "TL" => 1.45,
        "CD" | "OS" => 1.44,
        "CS" => 2.44,
        "BA" => 2.15,
        "LA" => 2.07,
        "CE" => 2.04,
        "PR" => 2.03,
        "SM" | "EU" => 1.98,
        "GD" => 1.96,
        "TB" => 1.94,
        "YB" | "LU" => 1.87,
        "W" => 1.62,
        "PT" | "AU" => 1.36,
        "PB" => 1.46,
        _ => 1.5,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::select::select;
    use crate::structure::{Molecule, SecondaryStructure};

    fn atom(residue: &str, name: &str, element: &str, position: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue),
            residue_seq: if residue == "LIG" { 900 } else { 10 },
            chain_id: FixedStr::new("A"),
            is_hetatm: residue != "LEU",
            molecule: Molecule::Unknown,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    /// A ligand nitrogen at the origin and one atom `distance` from it.
    fn measured(residue: &str, name: &str, element: &str, distance: f32) -> Distances {
        let structure = Protein {
            atoms: vec![
                atom("LIG", "N1", "N", Vec3::ZERO),
                atom(residue, name, element, Vec3::X * distance),
            ],
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        };
        let ligand = select("resn LIG", &structure).expect("a valid query");

        distances(&structure, &ligand, &Attachments::default())
    }

    #[test]
    fn a_metal_clashes_inside_three_quarters_of_the_covalent_radii_summed() {
        // N and Zn: 0.75 x (0.71 + 1.22) = 1.45 A, well short of a 2.0 A
        // coordination bond.
        for (distance, expected) in [(1.4, 1), (1.5, 0), (2.0, 0)] {
            let found = measured("ZN", "ZN", "ZN", distance);

            assert_eq!(found.metal_clashes, expected, "at {distance} A");
            assert_eq!(found.clashes, 0, "at {distance} A");
        }
    }

    #[test]
    fn an_ion_residue_is_measured_as_an_ion_whatever_its_elements() {
        // A sulfate oxygen at 2.0 A is inside N and O's van der Waals limit
        // (2.36 A) and outside their covalent one (1.03 A).
        let found = measured("SO4", "O1", "O", 2.0);
        assert_eq!((found.clashes, found.metal_clashes), (0, 0));

        assert_eq!(measured("SO4", "O1", "O", 1.0).metal_clashes, 1);
    }

    #[test]
    fn a_chlorine_or_bromine_is_measured_as_an_ion_whatever_its_residue() {
        // N and Cl: 0.75 x (1.6 + 1.8) = 2.55 A on van der Waals radii, 0.75 x
        // (0.71 + 1.02) = 1.30 A on covalent ones.
        let found = measured("XYZ", "CL1", "CL", 2.5);
        assert_eq!((found.clashes, found.metal_clashes), (0, 0));

        assert_eq!(measured("XYZ", "BR1", "BR", 1.2).metal_clashes, 1);
    }

    #[test]
    fn a_water_clashes_on_van_der_waals_radii_apart_from_the_protein() {
        // N and O: 0.75 x (1.6 + 1.55) = 2.36 A.
        for (distance, expected) in [(2.3, 1), (2.4, 0)] {
            let found = measured("HOH", "O", "O", distance);

            assert_eq!(found.water_clashes, expected, "at {distance} A");
            assert_eq!(found.clashes, 0, "at {distance} A");
        }
    }

    #[test]
    fn the_protein_and_organic_cofactors_clash_alike() {
        // N and C: 0.75 x (1.6 + 1.7) = 2.475 A.
        for residue in ["LEU", "HEM"] {
            for (distance, expected) in [(2.4, 1), (2.5, 0)] {
                let found = measured(residue, "C1", "C", distance);

                assert_eq!(found.clashes, expected, "{residue} at {distance} A");
            }
        }
    }

    #[test]
    fn hydrogens_never_clash() {
        let found = measured("LEU", "HD11", "H", 0.5);

        assert_eq!(found, Distances::default());
    }

    #[test]
    fn the_nearest_polymer_atom_says_whether_the_ligand_is_too_far() {
        for (distance, too_far) in [(4.9, false), (5.1, true)] {
            let found = measured("LEU", "CD1", "C", distance);

            let nearest = found.protein_distance.expect("a polymer atom");
            assert!((nearest - distance).abs() < 1e-5, "{nearest}");
            assert_eq!(found.too_far(), too_far, "at {distance} A");
        }
        // A water, an ion or a cofactor is no protein.
        for (residue, name, element) in [("HOH", "O", "O"), ("ZN", "ZN", "ZN"), ("HEM", "C1", "C")]
        {
            let found = measured(residue, name, element, 3.0);

            assert_eq!(found.protein_distance, None, "{residue}");
            assert!(found.too_far(), "{residue}");
        }
    }
}
