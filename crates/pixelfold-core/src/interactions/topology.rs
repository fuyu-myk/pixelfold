//! Covalent topology of the standard amino acids.
//!
//! A crystal structure carries coordinates, not bonds, and pixelfold has no
//! chemical component dictionary yet, so the heavy-atom connectivity of the 20
//! standard residues is tabulated here.
//!
//! Tabulating bonds rather than roles keeps the roles derivable. PLIP defines a
//! hydrophobic atom as a carbon whose neighbours are all carbon or hydrogen
//! (`Mol.hydrophobic_atoms`), a question it answers through OpenBabel's
//! perceived connectivity. Asking the same question of this table gives the same
//! answer from data that can be checked bond by bond against a textbook, rather
//! than from a hand-copied list of atom names whose derivation is invisible.
//!
//! Hydrogens are absent from the coordinates and so from this table.

/// Bonds present in every standard residue's backbone. `OXT` is the C-terminal
/// oxygen, absent from interior residues; a bond naming an atom the residue does
/// not have simply never matches.
///
/// The peptide bond C(i)-N(i+1) is deliberately omitted: it is the one bond that
/// crosses a residue boundary, and is unneeded. Both of its atoms already carry
/// an intra-residue non-carbon neighbour.
const BACKBONE: &[(&str, &str)] = &[("N", "CA"), ("CA", "C"), ("C", "O"), ("C", "OXT")];

/// Side-chain heavy-atom bonds, from CB outward. Glycine has no side chain, and
/// proline's CD closes the ring back onto the backbone nitrogen.
fn side_chain_bonds(residue_name: &str) -> Option<&'static [(&'static str, &'static str)]> {
    let bonds: &'static [(&'static str, &'static str)] = match residue_name {
        "ALA" => &[("CA", "CB")],
        "ARG" => &[
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("CD", "NE"),
            ("NE", "CZ"),
            ("CZ", "NH1"),
            ("CZ", "NH2"),
        ],
        "ASN" => &[("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")],
        "ASP" => &[("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")],
        "CYS" => &[("CA", "CB"), ("CB", "SG")],
        "GLN" => &[
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("CD", "OE1"),
            ("CD", "NE2"),
        ],
        "GLU" => &[
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("CD", "OE1"),
            ("CD", "OE2"),
        ],
        "GLY" => &[],
        // Imidazole: CG-ND1-CE1-NE2-CD2 closing back to CG.
        "HIS" => &[
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "ND1"),
            ("CG", "CD2"),
            ("ND1", "CE1"),
            ("CE1", "NE2"),
            ("NE2", "CD2"),
        ],
        "ILE" => &[("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1")],
        "LEU" => &[("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")],
        "LYS" => &[
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("CD", "CE"),
            ("CE", "NZ"),
        ],
        "MET" => &[("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
        // Benzene: CG-CD1-CE1-CZ-CE2-CD2 closing back to CG.
        "PHE" => &[
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD1"),
            ("CG", "CD2"),
            ("CD1", "CE1"),
            ("CD2", "CE2"),
            ("CE1", "CZ"),
            ("CE2", "CZ"),
        ],
        "PRO" => &[("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")],
        "SER" => &[("CA", "CB"), ("CB", "OG")],
        "THR" => &[("CA", "CB"), ("CB", "OG1"), ("CB", "CG2")],
        // Indole: a five-membered CG-CD1-NE1-CE2-CD2 ring fused to a six-membered
        // CD2-CE2-CZ2-CH2-CZ3-CE3 ring across the CD2-CE2 bond.
        "TRP" => &[
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD1"),
            ("CG", "CD2"),
            ("CD1", "NE1"),
            ("NE1", "CE2"),
            ("CD2", "CE2"),
            ("CD2", "CE3"),
            ("CE2", "CZ2"),
            ("CE3", "CZ3"),
            ("CZ2", "CH2"),
            ("CZ3", "CH2"),
        ],
        "TYR" => &[
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD1"),
            ("CG", "CD2"),
            ("CD1", "CE1"),
            ("CD2", "CE2"),
            ("CE1", "CZ"),
            ("CE2", "CZ"),
            ("CZ", "OH"),
        ],
        "VAL" => &[("CA", "CB"), ("CB", "CG1"), ("CB", "CG2")],
        _ => return None,
    };

    Some(bonds)
}

/// The heavy atoms bonded to `atom_name` within its own residue.
///
/// A residue outside the standard 20 has no known connectivity and so no
/// neighbours: its backbone is not assumed, because a ligand atom named `O` is
/// not a carbonyl oxygen just because a residue's would be.
pub fn heavy_neighbours<'a>(
    residue_name: &str,
    atom_name: &'a str,
) -> impl Iterator<Item = &'static str> + 'a {
    let side_chain = side_chain_bonds(residue_name);
    let backbone: &'static [(&'static str, &'static str)] = match side_chain {
        Some(_) => BACKBONE,
        None => &[],
    };

    backbone
        .iter()
        .chain(side_chain.unwrap_or(&[]).iter())
        .filter_map(move |&(a, b)| match () {
            _ if a == atom_name => Some(b),
            _ if b == atom_name => Some(a),
            _ => None,
        })
}

/// The element of a standard-residue atom, read from the first character of its
/// name. Every atom of the 20 standard residues is C, N, O, or S, and none of
/// those collide with the two-letter element symbols a ligand could carry.
fn element_of(atom_name: &str) -> Option<u8> {
    atom_name.bytes().next()
}

/// True when `atom_name` has at least one bonded heavy neighbour and every one
/// of them is a carbon.
///
/// This is the bond-level half of PLIP's hydrophobic test; the caller supplies
/// the element check. An atom with no neighbours in the table is not a known
/// atom of the residue, so it answers false rather than vacuously true.
pub fn bonded_only_to_carbon(residue_name: &str, atom_name: &str) -> bool {
    let mut neighbours = heavy_neighbours(residue_name, atom_name).peekable();

    neighbours.peek().is_some() && neighbours.all(|n| element_of(n) == Some(b'C'))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every atom named by a bond must be reachable from the atom on its other
    /// end: the table is a symmetric relation, not a directed one.
    #[test]
    fn every_bond_is_symmetric() {
        for residue in STANDARD_RESIDUES {
            for &(a, b) in side_chain_bonds(residue).unwrap() {
                assert!(
                    heavy_neighbours(residue, a).any(|n| n == b),
                    "{residue}: {a} does not list {b}"
                );
                assert!(
                    heavy_neighbours(residue, b).any(|n| n == a),
                    "{residue}: {b} does not list {a}"
                );
            }
        }
    }

    const STANDARD_RESIDUES: &[&str] = &[
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ];

    #[test]
    fn all_twenty_standard_residues_are_known_and_nothing_else_is() {
        assert_eq!(STANDARD_RESIDUES.len(), 20);
        for residue in STANDARD_RESIDUES {
            assert!(side_chain_bonds(residue).is_some(), "{residue} missing");
        }

        for other in ["HOH", "ZN", "STI", "MSE"] {
            assert!(side_chain_bonds(other).is_none(), "{other} is not standard");
            assert_eq!(heavy_neighbours(other, "CA").count(), 0);
        }

        // A water oxygen must not inherit the backbone carbonyl's carbon.
        assert!(!bonded_only_to_carbon("HOH", "O"));
    }

    #[test]
    fn rings_close() {
        // Each ring atom carries two ring bonds, so a walk returns to its start.
        assert!(heavy_neighbours("PRO", "N").any(|n| n == "CD"));
        assert!(heavy_neighbours("PRO", "CD").any(|n| n == "N"));

        // Phenylalanine's CZ is para to CG, reached either way round the ring.
        assert!(heavy_neighbours("PHE", "CZ").any(|n| n == "CE1"));
        assert!(heavy_neighbours("PHE", "CZ").any(|n| n == "CE2"));

        // Tryptophan's two rings are fused across the CD2-CE2 bond, so both
        // atoms carry three neighbours.
        assert_eq!(heavy_neighbours("TRP", "CD2").count(), 3);
        assert_eq!(heavy_neighbours("TRP", "CE2").count(), 3);
    }

    #[test]
    fn backbone_is_shared_by_every_residue() {
        for residue in STANDARD_RESIDUES {
            assert!(heavy_neighbours(residue, "N").any(|n| n == "CA"));
            assert!(heavy_neighbours(residue, "CA").any(|n| n == "C"));
            assert!(heavy_neighbours(residue, "C").any(|n| n == "O"));
        }

        // Glycine's CA is the backbone and nothing else.
        assert_eq!(heavy_neighbours("GLY", "CA").count(), 2);
        assert!(heavy_neighbours("ALA", "CA").any(|n| n == "CB"));
    }

    #[test]
    fn no_backbone_carbon_is_apolar() {
        // CA always carries the amide nitrogen and C always carries the carbonyl
        // oxygen, so hydrophobic contacts are a side-chain phenomenon.
        for residue in STANDARD_RESIDUES {
            assert!(!bonded_only_to_carbon(residue, "CA"), "{residue} CA");
            assert!(!bonded_only_to_carbon(residue, "C"), "{residue} C");
        }
    }

    #[test]
    fn a_carbon_bonded_to_sulfur_is_not_apolar() {
        // Methionine's CG and CE flank SD, and cysteine's CB carries SG, so PLIP
        // excludes them despite the residues being hydrophobic overall.
        assert!(!bonded_only_to_carbon("MET", "CG"));
        assert!(!bonded_only_to_carbon("MET", "CE"));
        assert!(!bonded_only_to_carbon("CYS", "CB"));

        // Methionine's CB is two bonds from the sulfur and so survives.
        assert!(bonded_only_to_carbon("MET", "CB"));
    }

    #[test]
    fn a_carbon_bonded_to_oxygen_or_nitrogen_is_not_apolar() {
        assert!(!bonded_only_to_carbon("SER", "CB")); // OG
        assert!(!bonded_only_to_carbon("THR", "CB")); // OG1
        assert!(!bonded_only_to_carbon("TYR", "CZ")); // OH
        assert!(!bonded_only_to_carbon("ASP", "CG")); // carboxylate
        assert!(!bonded_only_to_carbon("ARG", "CZ")); // guanidinium
        assert!(!bonded_only_to_carbon("LYS", "CE")); // NZ
        assert!(!bonded_only_to_carbon("TRP", "CD1")); // NE1
        assert!(!bonded_only_to_carbon("HIS", "CG")); // ND1
    }

    #[test]
    fn unknown_atoms_and_residues_are_not_apolar() {
        assert!(!bonded_only_to_carbon("ALA", "CZ9")); // not an atom of alanine
        assert!(!bonded_only_to_carbon("STI", "C1")); // ligand: unknown chemistry
    }

    /// The apolar carbons of every standard residue, derived by hand from the
    /// bond table and pinned here so a change to the table has to be deliberate.
    #[test]
    fn apolar_carbons_match_the_derivation() {
        let expected: &[(&str, &[&str])] = &[
            ("ALA", &["CB"]),
            ("ARG", &["CB", "CG"]),
            ("ASN", &["CB"]),
            ("ASP", &["CB"]),
            ("CYS", &[]),
            ("GLN", &["CB", "CG"]),
            ("GLU", &["CB", "CG"]),
            ("GLY", &[]),
            ("HIS", &["CB"]),
            ("ILE", &["CB", "CG1", "CG2", "CD1"]),
            ("LEU", &["CB", "CG", "CD1", "CD2"]),
            ("LYS", &["CB", "CG", "CD"]),
            ("MET", &["CB"]),
            ("PHE", &["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
            ("PRO", &["CB", "CG"]),
            ("SER", &[]),
            ("THR", &["CG2"]),
            ("TRP", &["CB", "CG", "CD2", "CE3", "CZ2", "CZ3", "CH2"]),
            ("TYR", &["CB", "CG", "CD1", "CD2", "CE1", "CE2"]),
            ("VAL", &["CB", "CG1", "CG2"]),
        ];

        for &(residue, apolar) in expected {
            let mut found: Vec<&str> = atoms_of(residue)
                .into_iter()
                .filter(|name| element_of(name) == Some(b'C'))
                .filter(|name| bonded_only_to_carbon(residue, name))
                .collect();
            found.sort_unstable();

            let mut want = apolar.to_vec();
            want.sort_unstable();

            assert_eq!(found, want, "{residue}");
        }
    }

    /// Every atom named anywhere in a residue's table.
    fn atoms_of(residue: &str) -> Vec<&'static str> {
        let mut atoms: Vec<&'static str> = BACKBONE
            .iter()
            .chain(side_chain_bonds(residue).unwrap().iter())
            .flat_map(|&(a, b)| [a, b])
            .collect();
        atoms.sort_unstable();
        atoms.dedup();

        atoms
    }
}
