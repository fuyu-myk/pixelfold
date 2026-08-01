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

use std::collections::HashMap;

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

/// One charged group of a component: its sign and the atom names whose centroid
/// is the charge centre.
pub(super) struct ChargedSite {
    pub sign: Charge,
    pub atoms: Vec<String>,
    /// The atom the group hangs off, when a bond from it to a neighbouring residue
    /// would neutralise the group. `None` when no such linkage neutralises it.
    pub centre: Option<String>,
}

/// The charged groups a component carries at pH 7.4, perceived from its own
/// definition the way the residue tables fix them for the standard residues.
///
/// This is the same protonation model, extended to a ligand or modified residue
/// from its own connectivity:
///
/// - **Negative** are the deprotonated oxyacids: a carboxylate (a carbon bearing
///   two terminal oxygens), a phosphate or phosphonate (a phosphorus with two or
///   more), a sulfonate or sulfate (a sulfur with three or more).
/// - **Positive** are the protonated bases: a guanidinium (a carbon with three
///   nitrogens, centred on them like an arginine), an amidinium (a carbon with
///   two nitrogens and no carbonyl, which separates it from a urea), and an
///   aliphatic amine (an sp3 nitrogen bonded only to carbon, primary through
///   quaternary).
///
/// Note: The exact set OpenBabel would protonate differs at the margins (a weak
/// base such as an aniline or an imidazole sits near the pH where the call is a
/// judgement), and those edges are where pixelfold and PLIP are expected to
/// disagree rather than fault.
///
/// A component definition describes the free molecule.
pub(super) fn component_charged_groups(component: &Component) -> Vec<ChargedSite> {
    let mut names: Vec<&str> = component.atoms().collect();
    names.sort_unstable();

    let mut sites = Vec::new();
    for name in names {
        let Some(element) = component.element_of(name) else {
            continue;
        };

        let site = match element.to_ascii_uppercase().as_str() {
            "C" => carbon_cation(component, name).or_else(|| carboxylate(component, name)),
            "P" => oxyanion(component, name, 2),
            "S" => oxyanion(component, name, 3),
            "N" => amine(component, name),
            _ => None,
        };

        if let Some(site) = site {
            sites.push(site);
        }
    }

    sites
}

/// A guanidinium or amidinium centred on `carbon`.
///
/// The carbon must lie outside any ring. A guanidinium or amidinium is a discrete
/// group on an sp2 carbon that stands apart from any ring; a ring carbon bearing
/// several nitrogens is an atom of a neutral heteroaromatic base, such as the C2
/// of a guanine or the C4 of a cytosine, and is not a cation.
fn carbon_cation(component: &Component, carbon: &str) -> Option<ChargedSite> {
    if in_ring(component, carbon) {
        return None;
    }

    let nitrogens: Vec<String> = component
        .neighbours(carbon)
        .filter(|&n| component.element_of(n).is_some_and(is_nitrogen))
        .map(String::from)
        .collect();

    // Three nitrogens is a guanidinium; two with no carbonyl is an amidinium,
    // the carbonyl being what would make it a neutral urea or carbamate instead.
    let charged = nitrogens.len() >= 3
        || (nitrogens.len() == 2
            && hybridisation(component, carbon) == 2
            && !has_carbonyl(component, carbon));

    charged.then_some(ChargedSite {
        sign: Charge::Positive,
        atoms: nitrogens,
        centre: None,
    })
}

/// A deprotonated oxyacid centred on `atom`, needing at least `minimum` terminal
/// oxygens: the terminal oxygens themselves.
///
/// A phosphate or sulfonate stays charged when it bridges two residues, as a
/// nucleic-acid phosphodiester does, so it carries no neutralising centre.
fn oxyanion(component: &Component, atom: &str, minimum: usize) -> Option<ChargedSite> {
    let oxygens: Vec<String> = terminal_oxygens(component, atom)
        .map(String::from)
        .collect();

    (oxygens.len() >= minimum).then_some(ChargedSite {
        sign: Charge::Negative,
        atoms: oxygens,
        centre: None,
    })
}

/// A carboxylate centred on `carbon`: its two terminal oxygens, or a
/// thiocarboxylate's oxygen and sulfur.
fn carboxylate(component: &Component, carbon: &str) -> Option<ChargedSite> {
    let oxygens: Vec<String> = terminal_oxygens(component, carbon)
        .map(String::from)
        .collect();
    if oxygens.len() >= 2 {
        return Some(ChargedSite {
            sign: Charge::Negative,
            atoms: oxygens,
            centre: Some(carbon.to_string()),
        });
    }

    (free_oxygens(component, carbon) == 1 && free_sulfurs(component, carbon) == 1).then(|| {
        let atoms = component
            .neighbours(carbon)
            .filter(|&n| {
                component.neighbours(n).count() == 1
                    && component
                        .element_of(n)
                        .is_some_and(|e| e.eq_ignore_ascii_case("O") || e.eq_ignore_ascii_case("S"))
            })
            .map(String::from)
            .collect();

        ChargedSite {
            sign: Charge::Negative,
            atoms,
            centre: Some(carbon.to_string()),
        }
    })
}

/// A protonated aliphatic amine at `nitrogen`: an sp3 nitrogen where every
/// heavy neighbour of which is a carbon, from a primary amine to a quaternary
/// ammonium. The nitrogen is its own charge centre.
fn amine(component: &Component, nitrogen: &str) -> Option<ChargedSite> {
    let mut neighbours = component.neighbours(nitrogen).peekable();
    let all_carbon = neighbours.peek().is_some()
        && neighbours.all(|n| component.element_of(n).is_some_and(is_carbon));

    (all_carbon && hybridisation(component, nitrogen) == 3).then(|| ChargedSite {
        sign: Charge::Positive,
        atoms: vec![nitrogen.to_string()],
        centre: Some(nitrogen.to_string()),
    })
}

/// The terminal oxygens hanging off `atom`, the ones with no other heavy bond.
fn terminal_oxygens<'a>(component: &'a Component, atom: &'a str) -> impl Iterator<Item = &'a str> {
    component.neighbours(atom).filter(|&n| {
        component.neighbours(n).count() == 1
            && component
                .element_of(n)
                .is_some_and(|e| e.eq_ignore_ascii_case("O"))
    })
}

fn is_nitrogen(element: &str) -> bool {
    element.eq_ignore_ascii_case("N")
}

/// Whether `atom` lies on a ring within its component, found by asking whether
/// any two of its neighbours stay connected once it is removed.
fn in_ring(component: &Component, atom: &str) -> bool {
    let neighbours: Vec<&str> = component.neighbours(atom).collect();
    if neighbours.len() < 2 {
        return false;
    }

    let mut reached: HashMap<&str, usize> = HashMap::new();
    for (source, &seed) in neighbours.iter().enumerate() {
        let mut stack = vec![seed];
        while let Some(current) = stack.pop() {
            match reached.get(current) {
                // Already reached from another neighbour: a ring closes here.
                Some(&other) if other != source => return true,
                Some(_) => continue,
                None => {
                    reached.insert(current, source);
                    stack.extend(component.neighbours(current).filter(|&n| n != atom));
                }
            }
        }
    }

    false
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

    /// A component carrying one of every charged group the pH 7.4 model
    /// recognises, plus the neutral groups it must leave out: a nitro, a urea,
    /// and a two-nitrogen carbon that is a ring atom of a base.
    fn charges(atom_name: &str) -> Option<(Charge, Vec<String>)> {
        const DEFINITION: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
LIG CC   C  N
LIG OC1  O  N
LIG OC2  O  N
LIG PP   P  N
LIG OP1  O  N
LIG OP2  O  N
LIG OP3  O  N
LIG SS   S  N
LIG OS1  O  N
LIG OS2  O  N
LIG OS3  O  N
LIG CG   C  N
LIG NG1  N  N
LIG NG2  N  N
LIG NG3  N  N
LIG CM   C  N
LIG NM1  N  N
LIG NM2  N  N
LIG NA   N  N
LIG CA1  C  N
LIG CU   C  N
LIG OU   O  N
LIG NU   N  N
LIG NU2  N  N
LIG CT   C  N
LIG NT   N  N
LIG OT1  O  N
LIG OT2  O  N
LIG CR1  C  Y
LIG NR1  N  Y
LIG CR2  C  Y
LIG NR2  N  Y
LIG CR3  C  Y
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG CC  OC1  doub
LIG CC  OC2  sing
LIG PP  OP1  doub
LIG PP  OP2  sing
LIG PP  OP3  sing
LIG SS  OS1  doub
LIG SS  OS2  doub
LIG SS  OS3  sing
LIG CG  NG1  doub
LIG CG  NG2  sing
LIG CG  NG3  sing
LIG CM  NM1  doub
LIG CM  NM2  sing
LIG NA  CA1  sing
LIG CU  OU   doub
LIG CU  NU   sing
LIG CU  NU2  sing
LIG CT  NT   sing
LIG NT  OT1  doub
LIG NT  OT2  doub
LIG CR1 NR1  arom
LIG NR1 CR2  arom
LIG CR2 NR2  arom
LIG NR2 CR3  arom
LIG CR3 CR1  arom
#
";
        let dictionary = crate::components::Dictionary::parse(DEFINITION);
        let component = dictionary.get("LIG").expect("the component");

        component_charged_groups(component)
            .into_iter()
            .find(|site| site.atoms.iter().any(|a| a == atom_name))
            .map(|site| {
                let mut atoms = site.atoms;
                atoms.sort();
                (site.sign, atoms)
            })
    }

    #[test]
    fn the_deprotonated_oxyacids_are_negative() {
        assert_eq!(
            charges("OC1"),
            Some((Charge::Negative, vec!["OC1".into(), "OC2".into()])),
            "a carboxylate"
        );
        assert_eq!(
            charges("OP1"),
            Some((
                Charge::Negative,
                vec!["OP1".into(), "OP2".into(), "OP3".into()]
            )),
            "a phosphate"
        );
        assert_eq!(
            charges("OS1"),
            Some((
                Charge::Negative,
                vec!["OS1".into(), "OS2".into(), "OS3".into()]
            )),
            "a sulfonate"
        );
    }

    #[test]
    fn the_protonated_bases_are_positive() {
        assert_eq!(
            charges("NG1"),
            Some((
                Charge::Positive,
                vec!["NG1".into(), "NG2".into(), "NG3".into()]
            )),
            "a guanidinium"
        );
        assert_eq!(
            charges("NM1"),
            Some((Charge::Positive, vec!["NM1".into(), "NM2".into()])),
            "an amidinium"
        );
        assert_eq!(
            charges("NA"),
            Some((Charge::Positive, vec!["NA".into()])),
            "an aliphatic amine"
        );
    }

    #[test]
    fn the_neutral_groups_carry_no_charge() {
        // A urea carbon has two nitrogens but a carbonyl, so it is not an
        // amidinium.
        assert_eq!(charges("NU"), None, "a urea nitrogen");
        // A nitro group is net neutral; its nitrogen and oxygens are left out.
        assert_eq!(charges("NT"), None, "a nitro nitrogen");
        assert_eq!(charges("OT1"), None, "a nitro oxygen");
        // A ring carbon bearing two ring nitrogens is a base, not an amidinium.
        assert_eq!(charges("NR1"), None, "an aromatic ring nitrogen");
    }

    /// A ring carbon whose exocyclic substituent is bonded first must still be
    /// found to be in the ring, so the imidazolinone carbon of a chromophore is
    /// not mistaken for an amidinium. This is the shape of the GFP chromophore's
    /// CA1-C1(-N2)(-N3) centre, with the exocyclic CA1 listed first.
    #[test]
    fn a_ring_carbon_with_an_exocyclic_first_neighbour_is_still_in_the_ring() {
        const CHROMOPHORE: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
LIG CA1 C
LIG C1  C
LIG N2  N
LIG CA2 C
LIG C2  C
LIG N3  N
LIG O2  O
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG CA1 C1  sing
LIG C1  N2  doub
LIG N2  CA2 sing
LIG CA2 C2  sing
LIG C2  N3  sing
LIG N3  C1  sing
LIG C2  O2  doub
#
";
        let dictionary = crate::components::Dictionary::parse(CHROMOPHORE);
        let component = dictionary.get("LIG").expect("the component");

        // C1 carries N2 and N3 and would read as an amidinium if the ring were
        // missed. It closes C1-N2-CA2-C2-N3-C1, so it is neutral.
        assert!(in_ring(component, "C1"));
        assert!(
            component_charged_groups(component).is_empty(),
            "the imidazolinone carbon is not a cation"
        );
    }
}
