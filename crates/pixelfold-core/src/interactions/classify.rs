//! Per-atom functional roles: which atoms can coordinate a metal, carry a
//! formal charge, or form a disulfide.
//!
//! PLIP derives these roles at runtime from OpenBabel typing; pixelfold uses
//! explicit residue tables instead (the approach RING takes). The memberships
//! below follow PLIP's enumeration so the two agree on what an atom *is*,
//! leaving only the geometric thresholds to differ.
//!
//! Protonation is the standard pH 7.4 model: Asp and Glu deprotonated (-1), Lys
//! and Arg protonated (+1), His admitted to the cationic set, Cys and Tyr
//! neutral.

use crate::structure::Atom;

use super::params::METAL_ELEMENTS;
use super::topology;

/// Formal charge sign of a residue's charged group.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Charge {
    Positive,
    Negative,
}

/// True when the residue name denotes a water molecule.
pub fn is_water(residue_name: &str) -> bool {
    matches!(residue_name, "HOH" | "WAT" | "H2O" | "DOD" | "TIP")
}

/// True when the atom is a coordinating metal ion, decided on the element
/// column rather than the atom name.
pub fn is_metal(atom: &Atom) -> bool {
    let element = atom.element.as_str().to_uppercase();
    METAL_ELEMENTS.contains(&element.as_str())
}

/// True when the atom can donate a lone pair to a metal ion.
///
/// Follows PLIP's enumeration: every backbone carbonyl oxygen, the side-chain
/// oxygens of Asp/Glu/Ser/Thr/Tyr, the His side-chain nitrogens, the Cys
/// sulfur, and water oxygens. Ligand donor groups arrive with the chemical
/// component dictionary.
pub fn is_metal_donor(atom: &Atom) -> bool {
    let residue = atom.residue_name.as_str();
    let name = atom.name.as_str();

    if is_water(residue) {
        return atom.element.as_str().eq_ignore_ascii_case("O");
    }

    // Backbone carbonyl, including the C-terminal oxygen.
    if name == "O" || name == "OXT" {
        return true;
    }

    matches!(
        (residue, name),
        ("ASP", "OD1" | "OD2")
            | ("GLU", "OE1" | "OE2")
            | ("SER", "OG")
            | ("THR", "OG1")
            | ("TYR", "OH")
            | ("HIS", "ND1" | "NE2")
            | ("CYS", "SG")
    )
}

/// The charged group of a residue for salt-bridge centroids: its sign and the
/// atom names whose centroid defines the charge centre.
///
/// Negative: Asp and Glu carboxylates. Positive: the Arg guanidinium, the Lys
/// ammonium, and the His imidazole. PLIP's protein path adds no explicit
/// terminus groups, so neither does this.
pub fn charged_group(residue_name: &str) -> Option<(Charge, &'static [&'static str])> {
    match residue_name {
        "ASP" => Some((Charge::Negative, &["OD1", "OD2"])),
        "GLU" => Some((Charge::Negative, &["OE1", "OE2"])),
        "ARG" => Some((Charge::Positive, &["NE", "NH1", "NH2"])),
        "LYS" => Some((Charge::Positive, &["NZ"])),
        "HIS" => Some((Charge::Positive, &["ND1", "NE2"])),
        _ => None,
    }
}

/// True when the atom is a cysteine gamma sulfur, the disulfide partner.
pub fn is_cysteine_sulfur(atom: &Atom) -> bool {
    atom.residue_name == "CYS" && atom.name == "SG"
}

/// True when the atom is an apolar carbon: a carbon bonded only to carbons and
/// hydrogens.
///
/// This is PLIP's hydrophobic atom, verbatim (`Mol.hydrophobic_atoms`: carbons
/// whose neighbours' atomic numbers are a subset of {1, 6}). PLIP reads the
/// neighbours from OpenBabel's perceived bonds; pixelfold reads them from
/// [`super::topology`], which gives the same answer for a standard residue and
/// no answer at all for a ligand until the chemical component dictionary lands.
///
/// Sulfur is excluded on both counts, which surprises chemical intuition and is
/// worth stating: methionine's SD is not an apolar atom, and its flanking CG and
/// CE are disqualified for bonding to it, so methionine contributes only CB.
/// RING instead admits sulfur to its van der Waals contacts.
pub fn is_apolar_carbon(atom: &Atom) -> bool {
    atom.element.as_str().eq_ignore_ascii_case("C")
        && topology::bonded_only_to_carbon(atom.residue_name.as_str(), atom.name.as_str())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixed_str::FixedStr;
    use crate::structure::SecondaryStructure;
    use glam::Vec3;

    fn atom(residue: &str, name: &str, element: &str) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue),
            residue_seq: 1,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    #[test]
    fn metal_is_decided_by_element_not_atom_name() {
        // The whole reason Atom carries an element column: a calcium ion and a
        // C-alpha carbon are both named "CA".
        assert!(is_metal(&atom("CA", "CA", "Ca")));
        assert!(!is_metal(&atom("ALA", "CA", "C")));
        assert!(is_metal(&atom("ZN", "ZN", "Zn")));
        assert!(!is_metal(&atom("ALA", "N", "N")));
    }

    #[test]
    fn metal_donors_cover_backbone_sidechain_and_water() {
        assert!(is_metal_donor(&atom("ALA", "O", "O"))); // backbone carbonyl
        assert!(is_metal_donor(&atom("GLY", "OXT", "O"))); // C-terminal oxygen
        assert!(is_metal_donor(&atom("ASP", "OD1", "O")));
        assert!(is_metal_donor(&atom("HIS", "NE2", "N")));
        assert!(is_metal_donor(&atom("CYS", "SG", "S")));
        assert!(is_metal_donor(&atom("HOH", "O", "O"))); // water

        // Not donors: backbone amide N, apolar carbons, a non-listed sidechain.
        assert!(!is_metal_donor(&atom("ALA", "N", "N")));
        assert!(!is_metal_donor(&atom("ALA", "CB", "C")));
        assert!(!is_metal_donor(&atom("LYS", "NZ", "N")));
    }

    #[test]
    fn charged_groups_match_the_ph_74_model() {
        assert_eq!(
            charged_group("ASP"),
            Some((Charge::Negative, &["OD1", "OD2"][..]))
        );
        assert_eq!(
            charged_group("ARG"),
            Some((Charge::Positive, &["NE", "NH1", "NH2"][..]))
        );
        assert_eq!(charged_group("LYS"), Some((Charge::Positive, &["NZ"][..])));
        assert!(charged_group("ALA").is_none());
        assert!(charged_group("SER").is_none()); // neutral hydroxyl
    }

    #[test]
    fn only_cysteine_gamma_sulfur_forms_disulfides() {
        assert!(is_cysteine_sulfur(&atom("CYS", "SG", "S")));
        assert!(!is_cysteine_sulfur(&atom("MET", "SD", "S")));
        assert!(!is_cysteine_sulfur(&atom("CYS", "CB", "C")));
    }

    #[test]
    fn apolar_carbons_are_side_chain_carbons_with_no_polar_neighbour() {
        assert!(is_apolar_carbon(&atom("LEU", "CD1", "C")));
        assert!(is_apolar_carbon(&atom("PHE", "CZ", "C")));
        assert!(is_apolar_carbon(&atom("ALA", "CB", "C")));

        // Backbone carbons always reach a nitrogen or an oxygen.
        assert!(!is_apolar_carbon(&atom("LEU", "CA", "C")));
        assert!(!is_apolar_carbon(&atom("LEU", "C", "C")));

        // Carbons one bond from a heteroatom.
        assert!(!is_apolar_carbon(&atom("SER", "CB", "C"))); // OG
        assert!(!is_apolar_carbon(&atom("TYR", "CZ", "C"))); // OH
        assert!(!is_apolar_carbon(&atom("MET", "CG", "C"))); // SD
    }

    #[test]
    fn sulfur_is_never_apolar_and_nor_is_a_ligand_carbon() {
        assert!(!is_apolar_carbon(&atom("MET", "SD", "S")));
        assert!(!is_apolar_carbon(&atom("CYS", "SG", "S")));

        // Ligand chemistry needs the chemical component dictionary.
        assert!(!is_apolar_carbon(&atom("STI", "C1", "C")));
    }
}
