//! Chemical component definitions: what a residue or ligand is made of.
//!
//! A ligand's chemistry is not derivable from its coordinates. It can be
//! derived in an mmCIF structure file as it already carries the definitions
//! of every component it uses, in the `_chem_comp_atom` and `_chem_comp_bond`
//! categories, complete with element and aromatic flags.
//!
//! A legacy PDB file has no such categories, and a handful of mmCIF files omit
//! them; both yield an empty dictionary, and callers fall back to the built-in
//! tables for the standard residues.
//!
//! Hydrogens are not kept as atoms, since the coordinates are heavy-atom only
//! and there would be nothing to match them against. They are counted instead:
//! whether an atom bears a hydrogen is what decides if it can donate, and the
//! count is part of the degree that decides if a nitrogen can accept.

use std::collections::{HashMap, HashSet};

use glam::Vec3;
use mmcif_cif::{Block, Column, Document};

/// How strongly two atoms are bonded, as the file writes it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondOrder {
    /// The bond's order counted in half bonds, so that an aromatic bond is
    /// three: an atom's sum is then the same whether a file writes its ring
    /// aromatic or in a Kekule structure.
    pub(crate) fn halves(self) -> u8 {
        match self {
            BondOrder::Single => 2,
            BondOrder::Double => 4,
            BondOrder::Triple => 6,
            BondOrder::Aromatic => 3,
        }
    }

    fn parse(value: &str) -> Self {
        match value.to_ascii_uppercase().as_str() {
            "DOUB" => BondOrder::Double,
            "TRIP" => BondOrder::Triple,
            "AROM" => BondOrder::Aromatic,
            _ => BondOrder::Single,
        }
    }
}

/// One heavy-atom bond.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Bond {
    first: String,
    second: String,
    order: BondOrder,
    /// Whether the definition calls the bond E or Z.
    stereo: bool,
}

/// One component's atoms and bonds, by atom name.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Component {
    /// Element symbol per atom name, as written in the mmCIF file.
    elements: HashMap<String, String>,
    /// The heavy atoms' names in the order the file lists them.
    order: Vec<String>,
    /// Atom names flagged aromatic.
    aromatic: Vec<String>,
    /// Heavy-atom bonds.
    bonds: Vec<Bond>,
    /// How many hydrogens each heavy atom carries.
    hydrogens: HashMap<String, usize>,
    /// Each heavy atom's formal charge, where the file gives charges.
    charges: HashMap<String, i8>,
    /// Each heavy atom's ideal position, where the file gives one: the Chemical
    /// Component Dictionary's model of the free molecule.
    ideal: HashMap<String, Vec3>,
    /// The heavy atoms the definition calls R or S.
    stereocentres: HashSet<String>,
}

impl Component {
    /// The element of a named atom.
    pub fn element_of(&self, atom_name: &str) -> Option<&str> {
        self.elements.get(atom_name).map(String::as_str)
    }

    /// Whether the file flags this atom as aromatic.
    pub fn is_aromatic(&self, atom_name: &str) -> bool {
        self.aromatic.iter().any(|a| a == atom_name)
    }

    /// Every atom the file flags as aromatic, in the order it lists them.
    pub fn aromatic_atoms(&self) -> impl Iterator<Item = &str> {
        self.aromatic.iter().map(String::as_str)
    }

    /// Every heavy atom the component defines, in no particular order.
    pub fn atoms(&self) -> impl Iterator<Item = &str> {
        self.elements.keys().map(String::as_str)
    }

    /// Every heavy atom the component defines, in the order the file lists
    /// them.
    pub fn atoms_in_order(&self) -> impl Iterator<Item = &str> {
        self.order.iter().map(String::as_str)
    }

    /// A heavy atom's formal charge, or `None` when the file gives no charges,
    /// as a structure file's own definitions do not.
    pub fn charge(&self, atom_name: &str) -> Option<i8> {
        self.charges.get(atom_name).copied()
    }

    /// A heavy atom's position in the component's ideal model, which the
    /// Chemical Component Dictionary gives and a structure file does not.
    pub fn ideal(&self, atom_name: &str) -> Option<Vec3> {
        self.ideal.get(atom_name).copied()
    }

    /// Whether the definition calls this atom an R or S stereocentre.
    pub fn is_stereocentre(&self, atom_name: &str) -> bool {
        self.stereocentres.contains(atom_name)
    }

    /// Whether the definition calls the bond between these atoms E or Z.
    pub fn has_stereo_bond(&self, first: &str, second: &str) -> bool {
        self.bonds.iter().any(|bond| {
            bond.stereo
                && ((bond.first == first && bond.second == second)
                    || (bond.first == second && bond.second == first))
        })
    }

    /// How many hydrogens the definition gives this heavy atom.
    pub fn hydrogen_count(&self, atom_name: &str) -> usize {
        self.hydrogens.get(atom_name).copied().unwrap_or(0)
    }

    /// The heavy atoms bonded to `atom_name`.
    pub fn neighbours<'a>(&'a self, atom_name: &'a str) -> impl Iterator<Item = &'a str> + 'a {
        self.bonded(atom_name).map(|(name, _)| name)
    }

    /// The heavy atoms bonded to `atom_name`, each with the bond's order.
    pub fn bonded<'a>(
        &'a self,
        atom_name: &'a str,
    ) -> impl Iterator<Item = (&'a str, BondOrder)> + 'a {
        self.bonds.iter().filter_map(move |bond| {
            if bond.first == atom_name {
                Some((bond.second.as_str(), bond.order))
            } else if bond.second == atom_name {
                Some((bond.first.as_str(), bond.order))
            } else {
                None
            }
        })
    }

    /// A component made of the given heavy atoms, in order, and the bonds
    /// between them, as a molecule's own bond table defines it.
    pub(crate) fn of_parts(
        atoms: Vec<DefinedAtom>,
        bonds: Vec<(String, String, BondOrder)>,
    ) -> Component {
        let mut component = Component::default();
        for atom in atoms {
            component.order.push(atom.name.clone());
            component.elements.insert(atom.name.clone(), atom.element);
            if atom.aromatic {
                component.aromatic.push(atom.name.clone());
            }
            component
                .hydrogens
                .insert(atom.name.clone(), atom.hydrogens);
            component.charges.insert(atom.name, atom.charge);
        }
        component.bonds = bonds
            .into_iter()
            .map(|(first, second, order)| Bond {
                first,
                second,
                order,
                stereo: false,
            })
            .collect();

        component
    }

    /// How many atoms the component defines.
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
}

// Every ideal coordinate is finite, refused at reading otherwise, so equality
// is reflexive.
impl Eq for Component {}

/// A heavy atom as [`Component::of_parts`] takes it.
pub(crate) struct DefinedAtom {
    pub(crate) name: String,
    pub(crate) element: String,
    pub(crate) aromatic: bool,
    pub(crate) hydrogens: usize,
    pub(crate) charge: i8,
}

/// Every component definition a structure file carries, by component id.
#[derive(Clone, Debug, Default)]
pub struct Dictionary {
    components: HashMap<String, Component>,
    /// The SMILES each component is given (`_chem_comp.pdbx_smiles`), which a
    /// co-folding model's file carries for a ligand it defines no atoms of.
    smiles: HashMap<String, String>,
}

impl Dictionary {
    /// Read the `_chem_comp_atom` and `_chem_comp_bond` categories out of the
    /// first data block of mmCIF text. Text that is not CIF, or that carries no
    /// such categories, yields an empty dictionary.
    pub fn parse(file_text: &str) -> Self {
        Document::parse(file_text)
            .ok()
            .and_then(|document| document.block(0).map(|block| Self::from_block(&block)))
            .unwrap_or_default()
    }

    /// Read the `_chem_comp_atom` and `_chem_comp_bond` categories of one data
    /// block. A block without them yields an empty dictionary.
    pub(crate) fn from_block(block: &Block<'_>) -> Self {
        let smiles = smiles_of(block);
        let mut components: HashMap<String, Component> = HashMap::new();
        // Hydrogen names per component.
        let mut hydrogen_names: HashMap<String, HashSet<String>> = HashMap::new();

        if let Some(atoms) = block.category("chem_comp_atom") {
            let (Some(ids), Some(names), Some(symbols)) = (
                atoms.column("comp_id"),
                atoms.column("atom_id"),
                atoms.column("type_symbol"),
            ) else {
                return Self {
                    smiles,
                    ..Self::default()
                };
            };
            let aromatic = atoms.column("pdbx_aromatic_flag");
            let charges = atoms.column("charge");
            let stereo = atoms.column("pdbx_stereo_config");
            let ideal = [
                atoms.column("pdbx_model_Cartn_x_ideal"),
                atoms.column("pdbx_model_Cartn_y_ideal"),
                atoms.column("pdbx_model_Cartn_z_ideal"),
            ];

            for row in 0..atoms.row_count() {
                let (Some(id), Some(name), Some(symbol)) =
                    (text(&ids, row), text(&names, row), text(&symbols, row))
                else {
                    continue;
                };
                if symbol.eq_ignore_ascii_case("H") || symbol.eq_ignore_ascii_case("D") {
                    hydrogen_names
                        .entry(id.to_string())
                        .or_default()
                        .insert(name.to_string());
                    continue; // the coordinates hold no hydrogens to match
                }

                let component = components.entry(id.to_string()).or_default();
                component
                    .elements
                    .insert(name.to_string(), symbol.to_string());
                component.order.push(name.to_string());
                if aromatic.as_ref().and_then(|flags| text(flags, row)) == Some("Y") {
                    component.aromatic.push(name.to_string());
                }
                if let Some(charge) = charges
                    .as_ref()
                    .and_then(|charges| text(charges, row)?.parse().ok())
                {
                    component.charges.insert(name.to_string(), charge);
                }
                if matches!(
                    stereo.as_ref().and_then(|flags| text(flags, row)),
                    Some("R" | "S")
                ) {
                    component.stereocentres.insert(name.to_string());
                }
                let coordinate = |axis: &Option<Column<'_>>| -> Option<f32> {
                    text(axis.as_ref()?, row)?
                        .parse()
                        .ok()
                        .filter(|value: &f32| value.is_finite())
                };
                if let [Some(x), Some(y), Some(z)] = ideal.each_ref().map(coordinate) {
                    component.ideal.insert(name.to_string(), Vec3::new(x, y, z));
                }
            }
        }

        if let Some(bonds) = block.category("chem_comp_bond") {
            let (Some(ids), Some(firsts), Some(seconds)) = (
                bonds.column("comp_id"),
                bonds.column("atom_id_1"),
                bonds.column("atom_id_2"),
            ) else {
                return Self { components, smiles };
            };
            let orders = bonds.column("value_order");
            let stereo = bonds.column("pdbx_stereo_config");

            for row in 0..bonds.row_count() {
                let (Some(id), Some(first), Some(second)) =
                    (text(&ids, row), text(&firsts, row), text(&seconds, row))
                else {
                    continue;
                };
                let hydrogens = hydrogen_names.get(id);
                let Some(component) = components.get_mut(id) else {
                    continue;
                };

                let first_heavy = component.elements.contains_key(first);
                let second_heavy = component.elements.contains_key(second);
                let is_hydrogen = |name: &str| hydrogens.is_some_and(|names| names.contains(name));

                if first_heavy && second_heavy {
                    component.bonds.push(Bond {
                        first: first.to_string(),
                        second: second.to_string(),
                        order: orders
                            .as_ref()
                            .and_then(|orders| text(orders, row))
                            .map_or(BondOrder::Single, BondOrder::parse),
                        stereo: matches!(
                            stereo.as_ref().and_then(|flags| text(flags, row)),
                            Some("E" | "Z")
                        ),
                    });
                } else if first_heavy && is_hydrogen(second) {
                    *component.hydrogens.entry(first.to_string()).or_default() += 1;
                } else if second_heavy && is_hydrogen(first) {
                    *component.hydrogens.entry(second.to_string()).or_default() += 1;
                }
            }
        }

        Self { components, smiles }
    }

    /// The SMILES the file gives a component, if any.
    pub fn smiles(&self, component_id: &str) -> Option<&str> {
        self.smiles.get(component_id).map(String::as_str)
    }

    /// The definition of a component, if the file carried one.
    pub fn get(&self, component_id: &str) -> Option<&Component> {
        self.components.get(component_id)
    }

    /// This dictionary with `component` defining `id`, in place of any
    /// definition it had.
    pub fn with(mut self, id: &str, component: Component) -> Dictionary {
        self.components.insert(id.to_string(), component);
        self
    }

    /// This dictionary, plus `fallback`'s definition of every component it lacks.
    pub fn with_fallback(mut self, fallback: &Dictionary) -> Dictionary {
        for (id, component) in &fallback.components {
            self.components
                .entry(id.clone())
                .or_insert_with(|| component.clone());
        }
        for (id, string) in &fallback.smiles {
            self.smiles
                .entry(id.clone())
                .or_insert_with(|| string.clone());
        }

        self
    }

    /// How many components are defined.
    pub fn len(&self) -> usize {
        self.components.len()
    }

    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }
}

/// Each component's SMILES from a block's `_chem_comp`.
fn smiles_of(block: &Block<'_>) -> HashMap<String, String> {
    let Some(category) = block.category("chem_comp") else {
        return HashMap::new();
    };
    let (Some(ids), Some(strings)) = (category.column("id"), category.column("pdbx_smiles")) else {
        return HashMap::new();
    };

    (0..category.row_count())
        .filter_map(|row| {
            let string = text(&strings, row).filter(|string| !matches!(*string, "?" | "."))?;
            Some((text(&ids, row)?.to_string(), string.to_string()))
        })
        .collect()
}

fn text<'a>(column: &Column<'a>, row: usize) -> Option<&'a str> {
    column.get(row)?.as_str()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An alanine and a small halogenated ligand, in the shape a structure file
    /// writes them.
    const SAMPLE: &str = "\
data_TEST
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_ordinal
ALA N   N  N 1
ALA CA  C  N 2
ALA CB  C  N 3
ALA H   H  N 4
LIG C1  C  Y 5
LIG C2  C  Y 6
LIG BR3 BR N 7
LIG H4  H  N 8
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_ordinal
ALA N  CA  sing N 1
ALA CA CB  sing N 2
ALA N  H   sing N 3
LIG C1 C2  arom Y 4
LIG C1 BR3 sing N 5
LIG C2 H4  sing N 6
#
";

    #[test]
    fn a_file_yields_its_components() {
        let dictionary = Dictionary::parse(SAMPLE);

        assert_eq!(dictionary.len(), 2);
        assert!(dictionary.get("ALA").is_some());
        assert!(dictionary.get("LIG").is_some());
        assert!(dictionary.get("ABSENT").is_none());
    }

    #[test]
    fn a_fallback_fills_only_the_components_a_dictionary_lacks() {
        // A second file defining LIG differently and adding GOL.
        let other = Dictionary::parse(
            "\
data_OTHER
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_ordinal
LIG X1 C N 1
GOL O1 O N 1
GOL C1 C N 2
#
",
        );
        let own = Dictionary::parse(SAMPLE);
        let merged = Dictionary::parse(SAMPLE).with_fallback(&other);

        assert_eq!(merged.len(), 3);
        assert_eq!(merged.get("LIG"), own.get("LIG"));
        assert_eq!(merged.get("GOL"), other.get("GOL"));
    }

    #[test]
    fn hydrogens_are_dropped_along_with_their_bonds() {
        let dictionary = Dictionary::parse(SAMPLE);
        let alanine = dictionary.get("ALA").unwrap();

        assert_eq!(alanine.len(), 3); // N, CA, CB and not H
        assert!(alanine.element_of("H").is_none());
        // The N-H bond went with it, so N reaches only CA.
        let neighbours: Vec<&str> = alanine.neighbours("N").collect();
        assert_eq!(neighbours, vec!["CA"]);
    }

    #[test]
    fn bonds_are_symmetric() {
        let dictionary = Dictionary::parse(SAMPLE);
        let alanine = dictionary.get("ALA").unwrap();

        assert!(alanine.neighbours("CA").any(|n| n == "N"));
        assert!(alanine.neighbours("N").any(|n| n == "CA"));
        assert!(alanine.neighbours("CA").any(|n| n == "CB"));
    }

    #[test]
    fn a_ligand_keeps_its_elements_and_aromatic_flags() {
        let dictionary = Dictionary::parse(SAMPLE);
        let ligand = dictionary.get("LIG").unwrap();

        assert_eq!(ligand.element_of("BR3"), Some("BR"));
        assert_eq!(ligand.element_of("C1"), Some("C"));
        assert!(ligand.is_aromatic("C1"));
        assert!(!ligand.is_aromatic("BR3"));

        // The bromine hangs off the ring carbon, which is what a halogen bond
        // needs to know.
        let attached: Vec<&str> = ligand.neighbours("BR3").collect();
        assert_eq!(attached, vec!["C1"]);
    }

    #[test]
    fn hydrogens_are_counted_against_the_atom_they_hang_from() {
        let dictionary = Dictionary::parse(SAMPLE);
        let alanine = dictionary.get("ALA").unwrap();
        let ligand = dictionary.get("LIG").unwrap();

        // The alanine nitrogen carries one hydrogen in the definition.
        assert_eq!(alanine.hydrogen_count("N"), 1);
        assert_eq!(alanine.hydrogen_count("CA"), 0);
        // The ligand's second ring carbon carries one; the bromine none.
        assert_eq!(ligand.hydrogen_count("C2"), 1);
        assert_eq!(ligand.hydrogen_count("BR3"), 0);
        // An atom the definition never mentions carries none.
        assert_eq!(alanine.hydrogen_count("ABSENT"), 0);
    }

    #[test]
    fn bond_order_is_read_from_the_file() {
        let dictionary = Dictionary::parse(SAMPLE);
        let ligand = dictionary.get("LIG").unwrap();

        let orders: Vec<(&str, BondOrder)> = ligand.bonded("C1").collect();
        assert!(orders.contains(&("C2", BondOrder::Aromatic)));
        assert!(orders.contains(&("BR3", BondOrder::Single)));
    }

    /// A component written as the Chemical Component Dictionary writes one: a
    /// charged oxygen, an R carbon, an E double bond, ideal coordinates.
    const DICTIONARY_ENTRY: &str = "\
data_TST
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
_chem_comp_atom.pdbx_ordinal
TST C3  C3  C 0  N N N 0.0 0.0 0.0  1.500  0.000 0.000 1
TST C1  C1  C 0  N N R 0.0 0.0 0.0  0.000  0.000 0.000 2
TST O4  O4  O -1 N N N 0.0 0.0 0.0 -0.700  1.200 0.000 3
TST C2  C2  C 0  N N N ? ? ?        2.200  1.100 0.000 4
TST H1  H1  H 0  N N N 0.0 0.0 0.0 -0.300 -0.500 0.900 5
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
TST C1 C3 SING N N 1
TST C3 C2 DOUB N E 2
TST C1 O4 SING N N 3
TST C1 H1 SING N N 4
#
";

    #[test]
    fn a_dictionary_entry_gives_charges_ideal_coordinates_and_stereo() {
        let dictionary = Dictionary::parse(DICTIONARY_ENTRY);
        let component = dictionary.get("TST").unwrap();

        assert_eq!(component.charge("O4"), Some(-1));
        assert_eq!(component.charge("C1"), Some(0));
        assert_eq!(component.ideal("C1"), Some(Vec3::ZERO));
        assert_eq!(component.ideal("O4"), Some(Vec3::new(-0.7, 1.2, 0.0)));
        assert!(component.is_stereocentre("C1"));
        assert!(!component.is_stereocentre("C3"));
        assert!(component.has_stereo_bond("C2", "C3"));
        assert!(!component.has_stereo_bond("C1", "C3"));
    }

    #[test]
    fn an_ideal_coordinate_that_is_not_a_finite_number_is_none() {
        let dictionary = Dictionary::parse(&DICTIONARY_ENTRY.replace("2.200  1.100", "nan  1.100"));
        let component = dictionary.get("TST").unwrap();

        assert_eq!(component.ideal("C2"), None);
        assert_eq!(component.clone(), *component);
    }

    #[test]
    fn a_co_folding_model_s_smiles_are_read_from_its_components() {
        // AlphaFold 3 writes each component's SMILES and no atoms for it.
        let dictionary = Dictionary::parse(
            "data_MODEL\nloop_\n_chem_comp.id\n_chem_comp.pdbx_smiles\n_chem_comp.type\n\
             ALA C[C@H](N)C(O)=O \"L-PEPTIDE LINKING\"\nLIG_B CCO non-polymer\nUNK ? other\n",
        );

        assert_eq!(dictionary.smiles("LIG_B"), Some("CCO"));
        assert_eq!(dictionary.smiles("ALA"), Some("C[C@H](N)C(O)=O"));
        assert_eq!(dictionary.smiles("UNK"), None);
        assert!(dictionary.get("LIG_B").is_none());
    }

    #[test]
    fn heavy_atoms_keep_the_order_the_file_lists_them_in() {
        let dictionary = Dictionary::parse(DICTIONARY_ENTRY);
        let component = dictionary.get("TST").unwrap();

        let order: Vec<&str> = component.atoms_in_order().collect();
        assert_eq!(order, ["C3", "C1", "O4", "C2"]);
    }

    #[test]
    fn a_structure_file_without_charges_or_ideal_coordinates_gives_none() {
        let dictionary = Dictionary::parse(SAMPLE);
        let ligand = dictionary.get("LIG").unwrap();

        assert_eq!(ligand.charge("C1"), None);
        assert_eq!(ligand.ideal("C1"), None);
        assert!(!ligand.is_stereocentre("C1"));
    }

    #[test]
    fn a_file_without_component_categories_is_empty() {
        let dictionary = Dictionary::parse("data_TEST\n#\nloop_\n_atom_site.id\n1\n#\n");

        assert!(dictionary.is_empty());
        assert!(dictionary.get("ALA").is_none());
    }
}
