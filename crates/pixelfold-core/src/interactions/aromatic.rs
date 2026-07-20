//! Aromatic ring perception.
//!
//! Pi interactions are between ring faces, so each aromatic ring is reduced to a
//! centre and a plane normal. The four aromatic residues have known rings, so
//! their members are tabulated rather than perceived: phenylalanine and tyrosine
//! share one six-membered ring, histidine has one five-membered imidazole, and
//! tryptophan has two, a five-membered pyrrole fused to a six-membered benzene.
//! PLIP reaches the same set by taking the smallest set of smallest rings and
//! keeping those OpenBabel calls aromatic or that belong to one of these four
//! residue names.
//!
//! The normal is fitted from every ring atom by Newell's method rather than from
//! three probe atoms as PLIP does. For a planar ring the two agree; Newell's is
//! steadier when thermal motion leaves a real ring slightly puckered, which is
//! the case that decides a stacking angle near its threshold. Every downstream
//! use of the normal is insensitive to its sign, so the winding direction does
//! not matter.
//!
//! Ligand rings wait on the chemical component dictionary; until then only these
//! residues carry rings.
//!
//! Source: PLIP `Mol.find_rings` in `plip/structure/preparation.py`.

use glam::Vec3;

use crate::fixed_str::FixedStr;
use crate::structure::{Atom, Protein};

type ResidueKey = (FixedStr<4>, u32, Option<char>);

fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

/// An aromatic ring reduced to the geometry the pi detectors need.
pub struct Ring {
    /// The ring's atom indices.
    pub atoms: Vec<usize>,
    /// The residue the ring belongs to.
    pub residue: ResidueKey,
    pub center: Vec3,
    /// Unit normal to the ring plane.
    pub normal: Vec3,
}

impl Ring {
    /// The lateral offset of `point`: how far its foot in the ring plane lies
    /// from the ring centre. Zero when `point` is directly over the centre.
    pub fn planar_offset(&self, point: Vec3) -> f32 {
        let from_center = point - self.center;
        let in_plane = from_center - self.normal * from_center.dot(self.normal);

        in_plane.length()
    }
}

/// The acute angle between two ring planes, in degrees (0 to 90).
///
/// A ring normal has no preferred sign, so the raw angle between two normals and
/// its supplement describe the same pair of planes; the smaller is the answer.
pub fn plane_angle(a: Vec3, b: Vec3) -> f32 {
    let (Some(a), Some(b)) = (a.try_normalize(), b.try_normalize()) else {
        return 0.0;
    };
    let angle = a.dot(b).clamp(-1.0, 1.0).acos().to_degrees();

    angle.min(180.0 - angle)
}

/// The ring atom names of a residue, each ring in traversal order so a polygon
/// normal is well defined. Phenylalanine and tyrosine share the six-membered
/// ring; the hydroxyl of tyrosine is not part of it.
fn ring_members(residue_name: &str) -> &'static [&'static [&'static str]] {
    match residue_name {
        "PHE" | "TYR" => &[&["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
        "HIS" => &[&["CG", "ND1", "CE1", "NE2", "CD2"]],
        "TRP" => &[
            &["CG", "CD1", "NE1", "CE2", "CD2"],
            &["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"],
        ],
        _ => &[],
    }
}

/// Every aromatic ring in the structure.
pub fn rings(protein: &Protein) -> Vec<Ring> {
    let mut out = Vec::new();
    let mut index = 0;

    while index < protein.atoms.len() {
        let start = index;
        let key = residue_key(&protein.atoms[start]);
        while index < protein.atoms.len() && residue_key(&protein.atoms[index]) == key {
            index += 1;
        }
        let residue = start..index;

        for names in ring_members(protein.atoms[start].residue_name.as_str()) {
            if let Some(ring) = build_ring(protein, residue.clone(), key, names) {
                out.push(ring);
            }
        }
    }

    out
}

/// Assemble one ring, or nothing if any of its atoms is missing.
fn build_ring(
    protein: &Protein,
    residue: std::ops::Range<usize>,
    key: ResidueKey,
    names: &[&str],
) -> Option<Ring> {
    let atoms: Vec<usize> = names
        .iter()
        .map(|&name| residue.clone().find(|&i| protein.atoms[i].name == name))
        .collect::<Option<_>>()?;

    let points: Vec<Vec3> = atoms.iter().map(|&i| protein.atoms[i].position).collect();
    let center = points.iter().copied().sum::<Vec3>() / points.len() as f32;

    Some(Ring {
        atoms,
        residue: key,
        center,
        normal: newell_normal(&points)?,
    })
}

/// The plane normal of a polygon by Newell's method, which sums signed edge
/// contributions over every vertex and so tolerates a ring that is not perfectly
/// flat. `None` for a degenerate ring with no area.
fn newell_normal(points: &[Vec3]) -> Option<Vec3> {
    let mut normal = Vec3::ZERO;
    for i in 0..points.len() {
        let current = points[i];
        let next = points[(i + 1) % points.len()];
        normal.x += (current.y - next.y) * (current.z + next.z);
        normal.y += (current.z - next.z) * (current.x + next.x);
        normal.z += (current.x - next.x) * (current.y + next.y);
    }

    normal.try_normalize()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::SecondaryStructure;

    fn atom(residue: &str, seq: u32, name: &str, position: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(&name[..1]),
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

    /// A flat regular hexagon in the z = h plane, named as a phenylalanine ring.
    fn phenylalanine(seq: u32, height: f32) -> Vec<Atom> {
        let names = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"];
        (0..6)
            .map(|k| {
                let theta = std::f32::consts::PI / 3.0 * k as f32;
                atom(
                    "PHE",
                    seq,
                    names[k],
                    Vec3::new(1.4 * theta.cos(), 1.4 * theta.sin(), height),
                )
            })
            .collect()
    }

    #[test]
    fn a_flat_hexagon_has_its_centre_and_normal() {
        let ring = &rings(&protein(phenylalanine(1, 0.0)))[0];

        assert_eq!(ring.atoms.len(), 6);
        assert!(ring.center.length() < 1e-5, "centre at the origin");
        // The ring lies in the xy plane, so its normal is the z axis.
        assert!(plane_angle(ring.normal, Vec3::Z) < 1e-3);
    }

    #[test]
    fn tryptophan_perceives_both_of_its_rings() {
        let mut atoms = vec![
            atom("TRP", 1, "CG", Vec3::new(0.0, 0.0, 0.0)),
            atom("TRP", 1, "CD1", Vec3::new(1.2, 0.5, 0.0)),
            atom("TRP", 1, "NE1", Vec3::new(2.0, -0.5, 0.0)),
            atom("TRP", 1, "CE2", Vec3::new(1.3, -1.6, 0.0)),
            atom("TRP", 1, "CD2", Vec3::new(0.0, -1.4, 0.0)),
            atom("TRP", 1, "CZ2", Vec3::new(1.6, -2.9, 0.0)),
            atom("TRP", 1, "CH2", Vec3::new(0.6, -3.9, 0.0)),
            atom("TRP", 1, "CZ3", Vec3::new(-0.7, -3.7, 0.0)),
            atom("TRP", 1, "CE3", Vec3::new(-1.0, -2.4, 0.0)),
        ];
        atoms.reverse(); // ring assembly must not depend on atom order

        let found = rings(&protein(atoms));
        assert_eq!(found.len(), 2);
        assert_eq!(found[0].atoms.len(), 5); // pyrrole
        assert_eq!(found[1].atoms.len(), 6); // benzene
    }

    #[test]
    fn a_ring_missing_an_atom_is_not_perceived() {
        let mut atoms = phenylalanine(1, 0.0);
        atoms.pop(); // drop CD2
        assert!(rings(&protein(atoms)).is_empty());
    }

    #[test]
    fn a_non_aromatic_residue_has_no_ring() {
        let alanine = protein(vec![
            atom("ALA", 1, "CA", Vec3::ZERO),
            atom("ALA", 1, "CB", Vec3::new(1.5, 0.0, 0.0)),
        ]);
        assert!(rings(&alanine).is_empty());
    }

    #[test]
    fn parallel_and_perpendicular_planes_read_as_zero_and_ninety() {
        assert!(plane_angle(Vec3::Z, Vec3::Z) < 1e-4);
        // A normal and its reverse describe the same plane: still parallel.
        assert!(plane_angle(Vec3::Z, -Vec3::Z) < 1e-4);
        assert!((plane_angle(Vec3::Z, Vec3::X) - 90.0).abs() < 1e-3);
    }

    #[test]
    fn planar_offset_is_the_in_plane_displacement() {
        let ring = &rings(&protein(phenylalanine(1, 0.0)))[0];

        // A point straight above the centre has no lateral offset.
        assert!(ring.planar_offset(Vec3::new(0.0, 0.0, 5.0)) < 1e-4);
        // Shifted sideways by 2, whatever its height.
        assert!((ring.planar_offset(Vec3::new(2.0, 0.0, 5.0)) - 2.0).abs() < 1e-4);
    }
}
