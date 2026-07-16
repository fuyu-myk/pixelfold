//! Per-type detectors. Each one narrows the atom set to the roles it needs,
//! grids those positions, and only tests candidates inside the cutoff.

use glam::Vec3;

use crate::fixed_str::FixedStr;
use crate::spatial::Grid;
use crate::structure::{Atom, Protein};

use super::classify::{Charge, charged_group, is_cysteine_sulfur, is_metal, is_metal_donor};
use super::params::{DISULFIDE_DIST_MAX, METAL_DIST_MAX, MIN_DIST, SALTBRIDGE_DIST_MAX};
use super::{Interaction, InteractionKind};

type ResidueKey = (FixedStr<4>, u32, Option<char>);

fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

/// Indices of the atoms satisfying `role`, with their positions.
fn collect(protein: &Protein, role: impl Fn(&Atom) -> bool) -> (Vec<usize>, Vec<Vec3>) {
    let indices: Vec<usize> = protein
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| role(atom))
        .map(|(index, _)| index)
        .collect();
    let positions = indices.iter().map(|&i| protein.atoms[i].position).collect();

    (indices, positions)
}

/// Cysteine SG-SG disulfide bridges. Detection is distance-only; the CB-SG-SG-CB
/// dihedral clusters near 90 degrees in real bridges but is not a gate.
pub fn disulfides(protein: &Protein) -> Vec<Interaction> {
    let (sulfurs, positions) = collect(protein, is_cysteine_sulfur);
    if sulfurs.len() < 2 {
        return Vec::new();
    }

    let grid = Grid::build(&positions, DISULFIDE_DIST_MAX);

    let mut out = Vec::new();
    for (a, &atom_a) in sulfurs.iter().enumerate() {
        let center = positions[a];
        grid.for_each_within(center, DISULFIDE_DIST_MAX, |b| {
            if b <= a {
                return; // report each pair once, never a sulfur against itself
            }

            let distance = (center - positions[b]).length();
            if distance <= MIN_DIST {
                return;
            }

            out.push(Interaction {
                kind: InteractionKind::Disulfide,
                atoms_a: vec![atom_a],
                atoms_b: vec![sulfurs[b]],
                distance,
                angle: None,
            });
        });
    }

    out
}

/// Metal ions coordinated by an O, N, or S donor.
pub fn metal_coordination(protein: &Protein) -> Vec<Interaction> {
    let (metals, metal_positions) = collect(protein, is_metal);
    if metals.is_empty() {
        return Vec::new();
    }

    let (donors, donor_positions) = collect(protein, is_metal_donor);
    if donors.is_empty() {
        return Vec::new();
    }

    let grid = Grid::build(&donor_positions, METAL_DIST_MAX);

    let mut out = Vec::new();
    for (m, &metal) in metals.iter().enumerate() {
        let center = metal_positions[m];
        grid.for_each_within(center, METAL_DIST_MAX, |d| {
            let distance = (center - donor_positions[d]).length();
            if !(MIN_DIST < distance && distance < METAL_DIST_MAX) {
                return;
            }

            out.push(Interaction {
                kind: InteractionKind::MetalCoordination,
                atoms_a: vec![metal],
                atoms_b: vec![donors[d]],
                distance,
                angle: None,
            });
        });
    }

    out
}

/// A residue's charged group: the centroid of its charged atoms.
struct ChargedGroup {
    sign: Charge,
    centroid: Vec3,
    atoms: Vec<usize>,
}

/// One charged group per charged residue, built from the atoms actually present
/// (a residue missing its carboxylate oxygens contributes nothing).
fn charged_groups(protein: &Protein) -> Vec<ChargedGroup> {
    let mut groups = Vec::new();
    let mut index = 0;

    while index < protein.atoms.len() {
        let start = index;
        let key = residue_key(&protein.atoms[start]);
        while index < protein.atoms.len() && residue_key(&protein.atoms[index]) == key {
            index += 1;
        }

        let Some((sign, names)) = charged_group(protein.atoms[start].residue_name.as_str()) else {
            continue;
        };
        let atoms: Vec<usize> = (start..index)
            .filter(|&i| names.contains(&protein.atoms[i].name.as_str()))
            .collect();
        if atoms.is_empty() {
            continue;
        }

        let sum = atoms
            .iter()
            .fold(Vec3::ZERO, |acc, &i| acc + protein.atoms[i].position);
        groups.push(ChargedGroup {
            sign,
            centroid: sum / atoms.len() as f32,
            atoms,
        });
    }

    groups
}

/// Salt bridges between oppositely charged groups, measured centroid to
/// centroid.
pub fn salt_bridges(protein: &Protein) -> Vec<Interaction> {
    let groups = charged_groups(protein);
    let (positive, negative): (Vec<&ChargedGroup>, Vec<&ChargedGroup>) = groups
        .iter()
        .partition(|group| group.sign == Charge::Positive);
    if positive.is_empty() || negative.is_empty() {
        return Vec::new();
    }

    let negative_centroids: Vec<Vec3> = negative.iter().map(|group| group.centroid).collect();
    let grid = Grid::build(&negative_centroids, SALTBRIDGE_DIST_MAX);

    let mut out = Vec::new();
    for cation in &positive {
        grid.for_each_within(cation.centroid, SALTBRIDGE_DIST_MAX, |n| {
            let anion = negative[n];
            let distance = (cation.centroid - anion.centroid).length();
            if !(MIN_DIST < distance && distance < SALTBRIDGE_DIST_MAX) {
                return;
            }

            out.push(Interaction {
                kind: InteractionKind::SaltBridge,
                atoms_a: cation.atoms.clone(),
                atoms_b: anion.atoms.clone(),
                distance,
                angle: None,
            });
        });
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::SecondaryStructure;

    fn atom(residue: &str, seq: u32, name: &str, element: &str, position: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn protein(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
        }
    }

    #[test]
    fn disulfide_found_at_bond_length_and_rejected_when_too_far() {
        // Two Cys gamma sulfurs at the nominal 2.05 A S-S bond length.
        let bonded = protein(vec![
            atom("CYS", 1, "SG", "S", Vec3::ZERO),
            atom("CYS", 2, "SG", "S", Vec3::new(2.05, 0.0, 0.0)),
        ]);
        let found = disulfides(&bonded);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::Disulfide);
        assert!((found[0].distance - 2.05).abs() < 1e-4);

        // Not a bridge beyond the 2.5 A cutoff.
        let apart = protein(vec![
            atom("CYS", 1, "SG", "S", Vec3::ZERO),
            atom("CYS", 2, "SG", "S", Vec3::new(3.0, 0.0, 0.0)),
        ]);
        assert!(disulfides(&apart).is_empty());
    }

    #[test]
    fn disulfide_reports_each_pair_once() {
        let structure = protein(vec![
            atom("CYS", 1, "SG", "S", Vec3::ZERO),
            atom("CYS", 2, "SG", "S", Vec3::new(2.05, 0.0, 0.0)),
        ]);
        assert_eq!(disulfides(&structure).len(), 1); // not two, mirrored
    }

    #[test]
    fn metal_coordinated_by_donor_but_not_by_carbon() {
        // Zinc with an Asp carboxylate oxygen in range, plus an apolar carbon
        // equally close that must not count.
        let structure = protein(vec![
            atom("ZN", 1, "ZN", "Zn", Vec3::ZERO),
            atom("ASP", 2, "OD1", "O", Vec3::new(2.1, 0.0, 0.0)),
            atom("ALA", 3, "CB", "C", Vec3::new(0.0, 2.1, 0.0)),
        ]);
        let found = metal_coordination(&structure);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::MetalCoordination);
        assert_eq!(found[0].atoms_b, vec![1]); // the Asp oxygen

        // Not coordinated beyond 3.0 A.
        let far = protein(vec![
            atom("ZN", 1, "ZN", "Zn", Vec3::ZERO),
            atom("ASP", 2, "OD1", "O", Vec3::new(4.0, 0.0, 0.0)),
        ]);
        assert!(metal_coordination(&far).is_empty());
    }

    #[test]
    fn calcium_ion_coordinates_but_c_alpha_never_does() {
        // Both are named "CA"; only the element distinguishes them.
        let ion = protein(vec![
            atom("CA", 1, "CA", "Ca", Vec3::ZERO),
            atom("GLU", 2, "OE1", "O", Vec3::new(2.4, 0.0, 0.0)),
        ]);
        assert_eq!(metal_coordination(&ion).len(), 1);

        let backbone = protein(vec![
            atom("ALA", 1, "CA", "C", Vec3::ZERO),
            atom("GLU", 2, "OE1", "O", Vec3::new(2.4, 0.0, 0.0)),
        ]);
        assert!(metal_coordination(&backbone).is_empty());
    }

    #[test]
    fn salt_bridge_between_carboxylate_and_guanidinium_centroids() {
        // Asp carboxylate centroid at x=0, Arg guanidinium centroid at x=4.
        let structure = protein(vec![
            atom("ASP", 1, "OD1", "O", Vec3::new(0.0, 1.0, 0.0)),
            atom("ASP", 1, "OD2", "O", Vec3::new(0.0, -1.0, 0.0)),
            atom("ARG", 2, "NE", "N", Vec3::new(4.0, 1.0, 0.0)),
            atom("ARG", 2, "NH1", "N", Vec3::new(4.0, 0.0, 0.0)),
            atom("ARG", 2, "NH2", "N", Vec3::new(4.0, -1.0, 0.0)),
        ]);
        let found = salt_bridges(&structure);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::SaltBridge);
        assert!((found[0].distance - 4.0).abs() < 1e-4); // centroid separation
        assert_eq!(found[0].atoms_a.len(), 3); // whole guanidinium
        assert_eq!(found[0].atoms_b.len(), 2); // whole carboxylate
    }

    #[test]
    fn salt_bridge_rejected_beyond_cutoff_and_between_like_charges() {
        let far = protein(vec![
            atom("ASP", 1, "OD1", "O", Vec3::ZERO),
            atom("ASP", 1, "OD2", "O", Vec3::ZERO),
            atom("LYS", 2, "NZ", "N", Vec3::new(7.0, 0.0, 0.0)),
        ]);
        assert!(salt_bridges(&far).is_empty());

        // Two anions close together must not pair.
        let like = protein(vec![
            atom("ASP", 1, "OD1", "O", Vec3::ZERO),
            atom("GLU", 2, "OE1", "O", Vec3::new(3.0, 0.0, 0.0)),
        ]);
        assert!(salt_bridges(&like).is_empty());
    }
}
