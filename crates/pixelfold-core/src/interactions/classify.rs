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

use crate::structure::{Atom, Protein};

use crate::components::{BondOrder, Component};

use super::connectivity::Bonds;
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

/// True when the atom can accept a hydrogen bond.
///
/// PLIP delegates this to OpenBabel's `IsHbondAcceptor`, whose rules reduce to
/// something small once only the standard residues can reach them:
///
/// - **Every oxygen accepts.** OpenBabel rejects an oxygen only for being nitro,
///   aromatic, a sulfone, a diaryl ether, or an sp3 ester bridge, and a protein
///   residue has none of those. Backbone carbonyls survive the ester rule
///   because that test needs a single bond, and carboxylates survive it by an
///   explicit exception.
/// - **Almost no nitrogen accepts.** OpenBabel rejects a nitrogen that is either
///   four-coordinate sp3 (an ammonium) or three-coordinate sp2 (an amide or a
///   pyrrole-type ring). Every nitrogen in a protein is one or the other, except
///   an unprotonated histidine ring nitrogen, which is two-coordinate.
/// - **No sulfur accepts**, since OpenBabel admits that sulfur is only at a formal
///   charge of -1. Cysteine's SG and methionine's SD are therefore excluded.
///   HBPLUS, CHARMM and ChimeraX disagree about those two atoms, and the
///   spectroscopy says sulfur hydrogen bonds are rare rather than weak, so this
///   is a real omission rather than a settled one.
///
/// Histidine's ND1 and NE2 both appear here and both also donate: only one
/// carries a hydrogen in a given tautomer, and the coordinates do not say which.
pub fn is_hbond_acceptor(protein: &Protein, bonds: &Bonds, index: usize) -> bool {
    let atom = &protein.atoms[index];
    let residue = atom.residue_name.as_str();

    if is_water(residue) {
        // A water is not an ordinary acceptor here.
        return false;
    }

    if topology::is_standard_residue(residue) {
        if atom.element.as_str().eq_ignore_ascii_case("O") {
            return true;
        }

        return matches!((residue, atom.name.as_str()), ("HIS", "ND1" | "NE2"));
    }

    let name = atom.name.as_str();
    let Some(component) = protein.components.get(residue) else {
        return false;
    };

    // A nitrogen the coordinates bond outside its own definition is the nitrogen
    // of a linkage, which is an amide in every common case.
    if atom.element.as_str().eq_ignore_ascii_case("N")
        && is_linked(protein, bonds, component, index)
    {
        return false;
    }

    accepts(component, name)
}

/// Whether the coordinates bond an atom to something its own component
/// definition does not list.
fn is_linked(protein: &Protein, bonds: &Bonds, component: &Component, index: usize) -> bool {
    let name = protein.atoms[index].name.as_str();

    bonds.of(index).iter().any(|&n| {
        let neighbour = protein.atoms[n].name.as_str();
        !component.neighbours(name).any(|other| other == neighbour)
    })
}

/// True when an oxygen belongs to an acid group that is ionised at pH 7.4: a
/// carboxylate or a phosphate.
pub(super) fn is_acid_oxygen(component: &Component, name: &str) -> bool {
    component
        .element_of(name)
        .is_some_and(|element| element.eq_ignore_ascii_case("O"))
        && component.neighbours(name).count() == 1
        && component
            .neighbours(name)
            .any(|centre| free_oxygens(component, centre) >= 2)
}

/// Whether an atom of a component accepts a hydrogen bond, on OpenBabel's rules.
///
/// The standard residues are answered by table above, because their chemistry is
/// known and fixed. Everything else is decided here from what the file says the
/// component is: its elements, its bond orders, its aromatic flags, and how many
/// hydrogens each atom carries.
///
/// Sulfur is the one branch left out. OpenBabel admits it only at a formal
/// charge of -1, and the file's component definitions carry no charge here, so a
/// thiolate is missed. Neutral sulfur would not have accepted anyway.
fn accepts(component: &Component, name: &str) -> bool {
    let Some(element) = component.element_of(name) else {
        return false;
    };

    match element.to_ascii_uppercase().as_str() {
        "O" => accepts_at_oxygen(component, name),
        "N" => accepts_at_nitrogen(component, name),
        // Organic fluorine does not accept; a bare fluoride would.
        "F" => !component
            .neighbours(name)
            .any(|n| component.element_of(n).is_some_and(is_carbon)),
        _ => false,
    }
}

/// An oxygen accepts unless it is nitro, aromatic, a sulfone oxygen, a bridge
/// between two aromatics, or the single-bonded oxygen of an ester.
fn accepts_at_oxygen(component: &Component, name: &str) -> bool {
    if is_nitro_oxygen(component, name)
        || component.is_aromatic(name)
        || is_sulfone_oxygen(component, name)
    {
        return false;
    }

    let mut aromatic_neighbours = 0;
    for (neighbour, order) in component.bonded(name) {
        if component.is_aromatic(neighbour) {
            aromatic_neighbours += 1;
            if aromatic_neighbours == 2 {
                return false; // an ether between two aromatic rings
            }
            continue;
        }

        if order == BondOrder::Single
            && component.element_of(neighbour).is_some_and(is_carbon)
            && has_carbonyl(component, neighbour)
            && !is_carboxyl_oxygen(component, name)
        {
            return false; // the bridging oxygen of an ester
        }
    }

    true
}

/// A nitrogen accepts unless its lone pair is spoken for: a four-coordinate
/// ammonium, or a three-coordinate sp2 nitrogen such as an amide or a
/// pyrrole-type ring.
fn accepts_at_nitrogen(component: &Component, name: &str) -> bool {
    let degree = component.neighbours(name).count() + component.hydrogen_count(name);

    !((degree == 4 && hybridisation(component, name) == 3)
        || (degree == 3 && hybridisation(component, name) == 2))
}

/// The hybridisation OpenBabel's typing rules assign, decoded from its SMARTS
/// table into the graph terms available here.
pub(super) fn hybridisation(component: &Component, name: &str) -> u8 {
    let bonded: Vec<(&str, BondOrder)> = component.bonded(name).collect();

    if bonded.iter().any(|&(_, order)| order == BondOrder::Triple)
        || bonded
            .iter()
            .filter(|&&(_, order)| order == BondOrder::Double)
            .count()
            >= 2
    {
        return 1;
    }

    let unsaturated_neighbour = bonded.iter().any(|&(neighbour, order)| {
        if order != BondOrder::Single && order != BondOrder::Aromatic {
            return false;
        }
        let Some(element) = component.element_of(neighbour) else {
            return false;
        };
        let organic = matches!(element.to_ascii_uppercase().as_str(), "C" | "N" | "O");

        organic
            && component
                .bonded(neighbour)
                .any(|(_, o)| o != BondOrder::Single)
    });

    let sp2 = component.is_aromatic(name)
        || bonded.iter().any(|&(_, order)| order == BondOrder::Double)
        || unsaturated_neighbour
        || bonded
            .iter()
            .any(|&(neighbour, _)| component.is_aromatic(neighbour));

    if sp2 { 2 } else { 3 }
}

/// A terminal oxygen hanging off `name`, which is how OpenBabel counts the
/// oxygens of a nitro, sulfone, or carboxyl group.
fn free_oxygens(component: &Component, name: &str) -> usize {
    component
        .neighbours(name)
        .filter(|&n| {
            component
                .element_of(n)
                .is_some_and(|e| e.eq_ignore_ascii_case("O"))
                && component.neighbours(n).count() == 1
        })
        .count()
}

fn free_sulfurs(component: &Component, name: &str) -> usize {
    component
        .neighbours(name)
        .filter(|&n| {
            component
                .element_of(n)
                .is_some_and(|e| e.eq_ignore_ascii_case("S"))
                && component.neighbours(n).count() == 1
        })
        .count()
}

/// A terminal oxygen whose only heavy neighbour is `element`, and the name of
/// that neighbour.
fn terminal_on<'a>(component: &'a Component, name: &'a str, element: &str) -> Option<&'a str> {
    if component.neighbours(name).count() != 1 {
        return None;
    }

    component.neighbours(name).find(|&n| {
        component
            .element_of(n)
            .is_some_and(|e| e.eq_ignore_ascii_case(element))
    })
}

fn is_nitro_oxygen(component: &Component, name: &str) -> bool {
    terminal_on(component, name, "N").is_some_and(|n| free_oxygens(component, n) == 2)
}

/// A sulfone oxygen: OpenBabel lets any nitrogen on the sulfur rescue the oxygen
/// as an acceptor.
fn is_sulfone_oxygen(component: &Component, name: &str) -> bool {
    terminal_on(component, name, "S").is_some_and(|s| {
        free_oxygens(component, s) == 2
            && !component.neighbours(s).any(|n| {
                component
                    .element_of(n)
                    .is_some_and(|e| e.eq_ignore_ascii_case("N"))
            })
    })
}

fn is_carboxyl_oxygen(component: &Component, name: &str) -> bool {
    terminal_on(component, name, "C").is_some_and(|c| {
        free_oxygens(component, c) == 2
            || (free_oxygens(component, c) == 1 && free_sulfurs(component, c) == 1)
    })
}

/// Whether a carbon carries a double bond to an oxygen, which is what makes a
/// single C-O bond an ester bond.
fn has_carbonyl(component: &Component, carbon: &str) -> bool {
    component.bonded(carbon).any(|(neighbour, order)| {
        order == BondOrder::Double
            && component
                .element_of(neighbour)
                .is_some_and(|e| e.eq_ignore_ascii_case("O"))
    })
}

fn is_carbon(element: &str) -> bool {
    element.eq_ignore_ascii_case("C")
}

/// True when the atom is an apolar carbon: a carbon bonded only to carbons and
/// hydrogens.
///
/// This is PLIP's hydrophobic atom, verbatim (`Mol.hydrophobic_atoms`: carbons
/// whose neighbours' atomic numbers are a subset of {1, 6}). PLIP reads the
/// neighbours from OpenBabel's perceived bonds; these come from the bond graph,
/// which draws on the residue table for a standard residue and on the file's own
/// component definitions for a ligand. Hydrogens are absent from the graph, so
/// "every heavy neighbour is a carbon" is the same rule.
///
/// An atom with no bonds at all is not apolar.which keeps a lone ion out of the
/// hydrophobic set.
///
/// Sulfur is excluded on both counts, which surprises chemical intuition and is
/// worth stating: methionine's SD is not an apolar atom, and its flanking CG and
/// CE are disqualified for bonding to it, so methionine contributes only CB.
/// RING instead admits sulfur to its van der Waals contacts.
pub fn is_apolar_carbon(protein: &Protein, bonds: &Bonds, index: usize) -> bool {
    let atom = &protein.atoms[index];
    if !atom.element.as_str().eq_ignore_ascii_case("C") {
        return false;
    }

    let mut neighbours = bonds.of(index).iter().peekable();

    neighbours.peek().is_some()
        && neighbours.all(|&n| protein.atoms[n].element.as_str().eq_ignore_ascii_case("C"))
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

    /// A residue laid out so its side chain can be walked: the bond graph,
    /// decides whether a carbon is apolar.
    fn residue(name: &str, atoms: &[(&str, &str)]) -> Protein {
        let atoms = atoms
            .iter()
            .map(|(atom_name, element)| {
                let mut a = atom(name, atom_name, element);
                a.residue_seq = 1;
                a
            })
            .collect();

        Protein {
            atoms,
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    fn apolar(residue_name: &str, atoms: &[(&str, &str)], wanted: &str) -> bool {
        let structure = residue(residue_name, atoms);
        let bonds = super::super::connectivity::bonds(&structure);
        let index = structure
            .atoms
            .iter()
            .position(|a| a.name == wanted)
            .expect("the atom under test");

        is_apolar_carbon(&structure, &bonds, index)
    }

    const LEUCINE: &[(&str, &str)] = &[
        ("N", "N"),
        ("CA", "C"),
        ("C", "C"),
        ("O", "O"),
        ("CB", "C"),
        ("CG", "C"),
        ("CD1", "C"),
        ("CD2", "C"),
    ];
    const SERINE: &[(&str, &str)] = &[
        ("N", "N"),
        ("CA", "C"),
        ("C", "C"),
        ("O", "O"),
        ("CB", "C"),
        ("OG", "O"),
    ];
    const METHIONINE: &[(&str, &str)] = &[
        ("N", "N"),
        ("CA", "C"),
        ("C", "C"),
        ("O", "O"),
        ("CB", "C"),
        ("CG", "C"),
        ("SD", "S"),
        ("CE", "C"),
    ];

    /// One component exercising each branch of the acceptor rules, written the
    /// way a structure file states them.
    fn typed(atom_name: &str) -> bool {
        const DEFINITION: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
LIG C1   C  N
LIG O2   O  N
LIG N3   N  N
LIG H3   H  N
LIG C4   C  N
LIG N5   N  N
LIG H5A  H  N
LIG H5B  H  N
LIG H5C  H  N
LIG C6   C  N
LIG N7   N  Y
LIG C8   C  Y
LIG NO   N  N
LIG ONA  O  N
LIG ONB  O  N
LIG CE1  C  N
LIG OE2  O  N
LIG CE3  C  N
LIG OE4  O  N
LIG F1   F  N
LIG CF   C  N
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG C1  O2   doub
LIG C1  N3   sing
LIG N3  H3   sing
LIG N3  C4   sing
LIG N5  H5A  sing
LIG N5  H5B  sing
LIG N5  H5C  sing
LIG N5  C6   sing
LIG N7  C8   arom
LIG NO  ONA  doub
LIG NO  ONB  doub
LIG NO  C1   sing
LIG CE1 OE2  doub
LIG CE1 OE4  sing
LIG OE4 CE3  sing
LIG F1  CF   sing
#
";
        let dictionary = crate::components::Dictionary::parse(DEFINITION);
        let component = dictionary.get("LIG").expect("the component");

        accepts(component, atom_name)
    }

    #[test]
    fn a_ligand_oxygen_accepts_unless_its_lone_pair_is_spoken_for() {
        assert!(typed("O2"), "a carbonyl oxygen accepts");
        assert!(typed("OE2"), "an ester carbonyl oxygen accepts");

        assert!(!typed("ONA"), "a nitro oxygen does not");
        assert!(!typed("OE4"), "the bridging oxygen of an ester does not");
    }

    #[test]
    fn a_ligand_nitrogen_accepts_only_with_a_free_lone_pair() {
        // An amide nitrogen is three-coordinate counting its hydrogen, and sp2
        // through the neighbouring carbonyl.
        assert!(!typed("N3"), "an amide nitrogen does not accept");
        // An ammonium is four-coordinate.
        assert!(!typed("N5"), "an ammonium does not accept");
        // An aromatic nitrogen with no hydrogen is two-coordinate.
        assert!(typed("N7"), "a pyridine-type nitrogen accepts");
    }

    #[test]
    fn organic_fluorine_never_accepts() {
        assert!(!typed("F1"));
    }

    #[test]
    fn apolar_carbons_are_side_chain_carbons_with_no_polar_neighbour() {
        assert!(apolar("LEU", LEUCINE, "CD1"));
        assert!(apolar("LEU", LEUCINE, "CG"));
        assert!(apolar("LEU", LEUCINE, "CB"));

        // Backbone carbons always reach a nitrogen or an oxygen.
        assert!(!apolar("LEU", LEUCINE, "CA"));
        assert!(!apolar("LEU", LEUCINE, "C"));

        // Carbons one bond from a heteroatom.
        assert!(!apolar("SER", SERINE, "CB")); // reaches OG
        assert!(!apolar("MET", METHIONINE, "CG")); // reaches SD
        assert!(!apolar("MET", METHIONINE, "CE")); // reaches SD
        assert!(apolar("MET", METHIONINE, "CB"));
    }

    /// A modified residue's backbone nitrogen is an amide nitrogen and accepts
    /// nothing, but its own definition cannot show that.
    #[test]
    fn a_nitrogen_bonded_outside_its_own_component_does_not_accept() {
        const HYDROXYPROLINE: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
HYP N   N
HYP H   H
HYP CA  C
HYP CD  C
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
HYP N  H
HYP N  CA
HYP N  CD
#
";
        let mut structure = residue("HYP", &[("N", "N"), ("CA", "C"), ("CD", "C")]);
        structure.components = crate::components::Dictionary::parse(HYDROXYPROLINE);
        structure.atoms[1].position = Vec3::new(1.5, 0.0, 0.0);
        structure.atoms[2].position = Vec3::new(-0.7, 1.3, 0.0);

        // On its own, the free component reads as an acceptor.
        let alone = super::super::connectivity::bonds(&structure);
        assert!(is_hbond_acceptor(&structure, &alone, 0));

        // Bonded into a chain, the same nitrogen is a tertiary amide.
        let mut linked = structure.atoms.clone();
        let mut carbon = atom("GLY", "C", "C");
        carbon.residue_seq = 2;
        carbon.position = Vec3::new(-0.5, -1.3, 0.0);
        linked.push(carbon);
        structure.atoms = linked;

        let bonds = super::super::connectivity::bonds(&structure);
        assert!(bonds.of(0).contains(&3), "the linkage was perceived");
        assert!(!is_hbond_acceptor(&structure, &bonds, 0));
    }

    /// The definitions write the neutral acid, and at pH 7.4 it is not.
    #[test]
    fn phosphate_and_carboxylate_oxygens_are_acid_oxygens_but_a_hydroxyl_is_not() {
        const ACIDS: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
LIG P1  P
LIG OP1 O
LIG OP2 O
LIG OP3 O
LIG C1  C
LIG O1  O
LIG O2  O
LIG C2  C
LIG O3  O
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG P1 OP1 doub
LIG P1 OP2 sing
LIG P1 OP3 sing
LIG C1 O1  doub
LIG C1 O2  sing
LIG C2 O3  sing
#
";
        let dictionary = crate::components::Dictionary::parse(ACIDS);
        let component = dictionary.get("LIG").expect("the component");

        for name in ["OP1", "OP2", "OP3", "O1", "O2"] {
            assert!(is_acid_oxygen(component, name), "{name} is an acid oxygen");
        }
        assert!(!is_acid_oxygen(component, "O3"), "a lone hydroxyl is not");
        assert!(!is_acid_oxygen(component, "C1"));
    }

    #[test]
    fn sulfur_is_never_apolar_and_nor_is_an_unbonded_atom() {
        assert!(!apolar("MET", METHIONINE, "SD"));

        // A component with no definition and no residue table has no bonds (apolar).
        assert!(!apolar("STI", &[("C1", "C"), ("C2", "C")], "C1"));
    }
}
