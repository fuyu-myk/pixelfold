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
}

/// One component's atoms and bonds, by atom name.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Component {
    /// Element symbol per atom name, as written in the mmCIF file.
    elements: HashMap<String, String>,
    /// Atom names flagged aromatic.
    aromatic: Vec<String>,
    /// Heavy-atom bonds.
    bonds: Vec<Bond>,
    /// How many hydrogens each heavy atom carries.
    hydrogens: HashMap<String, usize>,
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

    /// How many atoms the component defines.
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
}

/// Every component definition a structure file carries, by component id.
#[derive(Clone, Debug, Default)]
pub struct Dictionary {
    components: HashMap<String, Component>,
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
        let mut components: HashMap<String, Component> = HashMap::new();
        // Hydrogen names per component.
        let mut hydrogen_names: HashMap<String, HashSet<String>> = HashMap::new();

        if let Some(atoms) = block.category("chem_comp_atom") {
            let (Some(ids), Some(names), Some(symbols)) = (
                atoms.column("comp_id"),
                atoms.column("atom_id"),
                atoms.column("type_symbol"),
            ) else {
                return Self::default();
            };
            let aromatic = atoms.column("pdbx_aromatic_flag");

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
                if aromatic.as_ref().and_then(|flags| text(flags, row)) == Some("Y") {
                    component.aromatic.push(name.to_string());
                }
            }
        }

        if let Some(bonds) = block.category("chem_comp_bond") {
            let (Some(ids), Some(firsts), Some(seconds)) = (
                bonds.column("comp_id"),
                bonds.column("atom_id_1"),
                bonds.column("atom_id_2"),
            ) else {
                return Self { components };
            };
            let orders = bonds.column("value_order");

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
                    });
                } else if first_heavy && is_hydrogen(second) {
                    *component.hydrogens.entry(first.to_string()).or_default() += 1;
                } else if second_heavy && is_hydrogen(first) {
                    *component.hydrogens.entry(second.to_string()).or_default() += 1;
                }
            }
        }

        Self { components }
    }

    /// The definition of a component, if the file carried one.
    pub fn get(&self, component_id: &str) -> Option<&Component> {
        self.components.get(component_id)
    }

    /// How many components are defined.
    pub fn len(&self) -> usize {
        self.components.len()
    }

    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }
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

    #[test]
    fn a_file_without_component_categories_is_empty() {
        let dictionary = Dictionary::parse("data_TEST\n#\nloop_\n_atom_site.id\n1\n#\n");

        assert!(dictionary.is_empty());
        assert!(dictionary.get("ALA").is_none());
    }
}
