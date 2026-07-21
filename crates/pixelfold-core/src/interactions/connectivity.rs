//! Which atoms are bonded to which.
//!
//! A standard residue's bonds come from [`super::topology`], a hand-checked table.
//! Everything else comes from the definitions the structure file carries.
//!
//! An atom whose chemistry no source describes simply has no bonds.

use std::collections::HashMap;

use crate::structure::{Atom, Protein};

use super::hydrogens::MAX_PEPTIDE_BOND;
use super::topology;

/// The bonded heavy-atom neighbours of every atom, indexed by atom.
///
/// This is built once and shared.
pub struct Bonds {
    neighbours: Vec<Vec<usize>>,
}

impl Bonds {
    /// The atoms bonded to `atom`.
    pub fn of(&self, atom: usize) -> &[usize] {
        self.neighbours.get(atom).map_or(&[], Vec::as_slice)
    }

    /// The bonded neighbours of `atom` whose element is one of `elements`,
    /// matched case-insensitively.
    pub fn of_elements<'a>(
        &'a self,
        protein: &'a Protein,
        atom: usize,
        elements: &'a [&'a str],
    ) -> impl Iterator<Item = usize> + 'a {
        self.of(atom).iter().copied().filter(move |&n| {
            let element = protein.atoms[n].element.as_str();
            elements.iter().any(|e| element.eq_ignore_ascii_case(e))
        })
    }

    /// The one bonded neighbour of `atom` drawn from `elements`, or nothing when
    /// it has none or several.
    ///
    /// A single antecedent is what the halogen bond needs on both sides, as it
    /// fixes the direction the interaction is measured against. An atom with two
    /// of them, a ring nitrogen or a thioether sulfur, has no single one.
    pub fn sole_neighbour(
        &self,
        protein: &Protein,
        atom: usize,
        elements: &[&str],
    ) -> Option<usize> {
        let mut found = self.of_elements(protein, atom, elements);
        let first = found.next()?;
        found.next().is_none().then_some(first)
    }
}

/// Build the bond graph of a structure.
pub fn bonds(protein: &Protein) -> Bonds {
    let mut neighbours = vec![Vec::new(); protein.atoms.len()];
    let mut link = |a: usize, b: usize| {
        if a != b {
            neighbours[a].push(b);
            neighbours[b].push(a);
        }
    };

    // The carbonyl carbon of the residue before, for the peptide bond.
    let mut preceding: Option<(usize, crate::fixed_str::FixedStr<4>)> = None;

    let mut index = 0;
    while index < protein.atoms.len() {
        let start = index;
        let key = residue_key(&protein.atoms[start]);
        while index < protein.atoms.len() && residue_key(&protein.atoms[index]) == key {
            index += 1;
        }
        let residue = start..index;

        let by_name: HashMap<&str, usize> = residue
            .clone()
            .map(|i| (protein.atoms[i].name.as_str(), i))
            .collect();

        let residue_name = protein.atoms[start].residue_name.as_str();
        for i in residue.clone() {
            let name = protein.atoms[i].name.as_str();
            for other in bonded_names(protein, residue_name, name) {
                if let Some(&j) = by_name.get(other.as_str())
                    && i < j
                {
                    link(i, j);
                }
            }
        }

        // The peptide bond: this residue's N to the previous residue's C.
        if let (Some((carbon, chain)), Some(&nitrogen)) = (preceding, by_name.get("N"))
            && chain == protein.atoms[start].chain_id
            && (protein.atoms[carbon].position - protein.atoms[nitrogen].position).length()
                <= MAX_PEPTIDE_BOND
        {
            link(carbon, nitrogen);
        }
        if let Some(&carbon) = by_name.get("C") {
            preceding = Some((carbon, protein.atoms[start].chain_id));
        }

        index = residue.end;
    }

    Bonds { neighbours }
}

/// The names bonded to `name` within its component, from the residue table for a
/// standard residue and from the file's own definitions otherwise.
fn bonded_names(protein: &Protein, residue_name: &str, name: &str) -> Vec<String> {
    if topology::is_standard_residue(residue_name) {
        return topology::heavy_neighbours(residue_name, name)
            .map(String::from)
            .collect();
    }

    protein
        .components
        .get(residue_name)
        .map(|component| component.neighbours(name).map(String::from).collect())
        .unwrap_or_default()
}

type ResidueKey = (crate::fixed_str::FixedStr<4>, u32, Option<char>);

fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixed_str::FixedStr;
    use crate::structure::SecondaryStructure;
    use glam::Vec3;

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
            components: Default::default(),
        }
    }

    /// Two bonded alanines, laid out so the peptide bond is in range.
    fn dipeptide() -> Protein {
        protein(vec![
            atom("ALA", 1, "N", "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("ALA", 1, "CA", "C", Vec3::new(1.5, 0.0, 0.0)),
            atom("ALA", 1, "C", "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("ALA", 1, "O", "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("ALA", 2, "N", "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("ALA", 2, "CA", "C", Vec3::new(4.0, 2.8, 0.0)),
        ])
    }

    #[test]
    fn a_standard_residue_uses_the_residue_table() {
        let structure = dipeptide();
        let bonds = bonds(&structure);

        // Residue 1's CA reaches N, C and CB; CB is absent here, so N and C.
        let mut ca: Vec<usize> = bonds.of(1).to_vec();
        ca.sort_unstable();
        assert_eq!(ca, vec![0, 2]);
        assert!(bonds.of(2).contains(&3)); // C to O
    }

    #[test]
    fn the_peptide_bond_joins_consecutive_residues() {
        let structure = dipeptide();
        let bonds = bonds(&structure);

        // Residue 1's C bonds residue 2's N.
        assert!(bonds.of(2).contains(&4), "C(1) reaches N(2)");
        assert!(bonds.of(4).contains(&2), "and the bond is symmetric");
    }

    #[test]
    fn a_backbone_nitrogen_in_a_chain_has_two_antecedents() {
        // This is why a backbone amide nitrogen is not a halogen acceptor: it
        // carries the peptide bond as well as its own C-alpha.
        let structure = dipeptide();
        let bonds = bonds(&structure);

        assert!(
            bonds
                .sole_neighbour(&structure, 4, &["C", "N", "P", "S"])
                .is_none(),
            "the interior amide nitrogen has no single antecedent"
        );
        // The first residue's nitrogen has only its C-alpha, so it does.
        assert_eq!(
            bonds.sole_neighbour(&structure, 0, &["C", "N", "P", "S"]),
            Some(1)
        );
    }

    #[test]
    fn a_break_too_long_for_a_peptide_bond_is_not_joined() {
        let mut atoms = dipeptide().atoms;
        atoms[4].position = Vec3::new(30.0, 1.5, 0.0); // residue 2's N, far away
        let structure = protein(atoms);
        let bonds = bonds(&structure);

        assert!(!bonds.of(2).contains(&4));
    }

    #[test]
    fn a_carbonyl_oxygen_has_one_antecedent() {
        let structure = dipeptide();
        let bonds = bonds(&structure);

        assert_eq!(
            bonds.sole_neighbour(&structure, 3, &["C", "N", "P", "S"]),
            Some(2),
            "the carbonyl oxygen hangs off its carbon alone"
        );
    }

    #[test]
    fn an_unknown_component_has_no_bonds() {
        // No residue table and no definition in the file: nothing is invented.
        let structure = protein(vec![
            atom("LIG", 1, "C1", "C", Vec3::ZERO),
            atom("LIG", 1, "BR2", "BR", Vec3::new(1.9, 0.0, 0.0)),
        ]);
        let bonds = bonds(&structure);

        assert!(bonds.of(0).is_empty());
        assert!(bonds.of(1).is_empty());
    }
}
