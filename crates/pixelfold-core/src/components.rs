//! Chemical component definitions: what a residue or ligand is made of.
//!
//! A ligand's chemistry is not derivable from its coordinates. It can be
//! derived in an mmCIF structure file as it already carries the definitions
//! of every component it uses, in the `_chem_comp_atom` and `_chem_comp_bond`
//! categories, complete with element and aromatic flags.
//!
//! A legacy PDB file has no such categories, and a handful of mmCIF files omit
//! them; both yield an empty dictionary, and callers fall back to the built-in
//! tables for the standard residues. Hydrogens are dropped, since the coordinates
//! are heavy-atom only.

use std::collections::HashMap;

use crate::mmcif::{column_of, read_loop};

/// One component's atoms and bonds, by atom name.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Component {
    /// Element symbol per atom name, as written in the mmCIF file.
    elements: HashMap<String, String>,
    /// Atom names flagged aromatic.
    aromatic: Vec<String>,
    /// Heavy-atom bonds, as name pairs.
    bonds: Vec<(String, String)>,
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

    /// The heavy atoms bonded to `atom_name`.
    pub fn neighbours<'a>(&'a self, atom_name: &'a str) -> impl Iterator<Item = &'a str> + 'a {
        self.bonds.iter().filter_map(move |(a, b)| {
            if a == atom_name {
                Some(b.as_str())
            } else if b == atom_name {
                Some(a.as_str())
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
    /// Read the `_chem_comp_atom` and `_chem_comp_bond` categories out of an
    /// mmCIF file. A file without them yields an empty dictionary.
    pub fn parse(file_text: &str) -> Self {
        let lines: Vec<&str> = file_text.lines().collect();
        let mut components: HashMap<String, Component> = HashMap::new();

        if let Some((columns, rows)) = read_loop(&lines, "_chem_comp_atom.") {
            let (Some(id), Some(name), Some(symbol)) = (
                column_of(&columns, "comp_id"),
                column_of(&columns, "atom_id"),
                column_of(&columns, "type_symbol"),
            ) else {
                return Self::default();
            };
            let aromatic = column_of(&columns, "pdbx_aromatic_flag");

            for row in rows {
                let (Some(id), Some(name), Some(symbol)) =
                    (row.get(id), row.get(name), row.get(symbol))
                else {
                    continue;
                };
                if symbol.eq_ignore_ascii_case("H") || symbol.eq_ignore_ascii_case("D") {
                    continue; // the coordinates hold no hydrogens to match
                }

                let component = components.entry(id.clone()).or_default();
                component.elements.insert(name.clone(), symbol.clone());
                if aromatic.and_then(|a| row.get(a)).is_some_and(|f| f == "Y") {
                    component.aromatic.push(name.clone());
                }
            }
        }

        if let Some((columns, rows)) = read_loop(&lines, "_chem_comp_bond.") {
            let (Some(id), Some(first), Some(second)) = (
                column_of(&columns, "comp_id"),
                column_of(&columns, "atom_id_1"),
                column_of(&columns, "atom_id_2"),
            ) else {
                return Self { components };
            };

            for row in rows {
                let (Some(id), Some(first), Some(second)) =
                    (row.get(id), row.get(first), row.get(second))
                else {
                    continue;
                };
                let Some(component) = components.get_mut(id) else {
                    continue;
                };

                // A bond naming an atom the component did not define is a bond
                // to a hydrogen, which was dropped above.
                if component.elements.contains_key(first) && component.elements.contains_key(second)
                {
                    component.bonds.push((first.clone(), second.clone()));
                }
            }
        }

        Self::with_numeric_aliases(components)
    }

    /// Register the spelling a numeric-looking component id takes on in the
    /// coordinates, so a definition that is present cannot be unreachable.
    fn with_numeric_aliases(mut components: HashMap<String, Component>) -> Self {
        let aliases: Vec<(String, Component)> = components
            .iter()
            .filter_map(|(id, component)| {
                let alias = numeric_alias(id)?;
                (!components.contains_key(&alias)).then(|| (alias, component.clone()))
            })
            .collect();
        components.extend(aliases);

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

/// How a component id reads back after a numeric parse, or nothing when it is
/// not a number or already reads as itself.
fn numeric_alias(component_id: &str) -> Option<String> {
    let value: f64 = component_id.parse().ok()?;
    let rendered = format!("{value}");

    (rendered != component_id).then_some(rendered)
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

    /// A component whose id reads as a number is reachable under both spellings,
    /// because the coordinates and the definitions are read by different parsers
    /// that disagree about it.
    #[test]
    fn a_numeric_looking_component_id_is_reachable_from_the_coordinates() {
        let text = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
4E2 C1  C
4E2 CL2 CL
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
4E2 C1 CL2
#
";
        let dictionary = Dictionary::parse(text);

        // Under the name the file writes.
        assert!(dictionary.get("4E2").is_some());
        // And under the name a numeric parse leaves, which is what the
        // coordinates carry.
        let aliased = dictionary.get("400").expect("reachable as 400 too");
        assert_eq!(aliased.neighbours("CL2").collect::<Vec<_>>(), vec!["C1"]);
    }

    #[test]
    fn a_real_component_keeps_a_name_an_alias_would_take() {
        let text = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
4E2 C1  C
400 O1  O
#
";
        let dictionary = Dictionary::parse(text);

        // The alias must not overwrite a component genuinely named 400.
        assert_eq!(dictionary.get("400").unwrap().element_of("O1"), Some("O"));
    }

    #[test]
    fn an_ordinary_component_id_gains_no_alias() {
        let dictionary = Dictionary::parse(SAMPLE);

        // ALA and LIG read as no number, so nothing extra is registered.
        assert_eq!(dictionary.len(), 2);
    }

    #[test]
    fn a_file_without_component_categories_is_empty() {
        let dictionary = Dictionary::parse("data_TEST\n#\nloop_\n_atom_site.id\n1\n#\n");

        assert!(dictionary.is_empty());
        assert!(dictionary.get("ALA").is_none());
    }
}
