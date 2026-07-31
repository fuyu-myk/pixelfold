//! The headless analyses: what each subcommand computes and the records it
//! reports.
//!
//! Each analysis narrows the structure to the chosen atoms, runs the core
//! engine, and turns the result into records [`crate::report`] knows how to
//! write.

use pixelfold_core::interactions::{Interaction, InteractionKind, detect};
use pixelfold_core::sasa::{SurfaceCalculator, relative_sasa};
use pixelfold_core::{Atom, AtomSet, Protein};
use serde::Serialize;

use crate::report::Row;

/// An insertion code as a cell: the letter, or empty where there is none.
fn code(insertion: Option<char>) -> String {
    insertion.map(String::from).unwrap_or_default()
}

/// Round to `places` decimals.
fn rounded(value: f32, places: u32) -> f32 {
    let scale = 10f32.powi(places as i32);

    (value * scale).round() / scale
}

/// One detected interaction, flattened to a row.
#[derive(Serialize)]
pub struct InteractionRecord {
    pub kind: String,
    pub chain_a: String,
    pub resi_a: u32,
    pub icode_a: Option<char>,
    pub resn_a: String,
    pub atoms_a: String,
    pub chain_b: String,
    pub resi_b: u32,
    pub icode_b: Option<char>,
    pub resn_b: String,
    pub atoms_b: String,
    pub distance: f32,
    pub angle: Option<f32>,
}

impl Row for InteractionRecord {
    fn columns() -> &'static [&'static str] {
        &[
            "kind", "chain_a", "resi_a", "icode_a", "resn_a", "atoms_a", "chain_b", "resi_b",
            "icode_b", "resn_b", "atoms_b", "distance", "angle",
        ]
    }

    fn cells(&self) -> Vec<String> {
        vec![
            self.kind.clone(),
            self.chain_a.clone(),
            self.resi_a.to_string(),
            code(self.icode_a),
            self.resn_a.clone(),
            self.atoms_a.clone(),
            self.chain_b.clone(),
            self.resi_b.to_string(),
            code(self.icode_b),
            self.resn_b.clone(),
            self.atoms_b.clone(),
            format!("{:.2}", self.distance),
            self.angle.map_or_else(String::new, |a| format!("{a:.1}")),
        ]
    }
}

/// Every detected interaction touching the chosen atoms, of the wanted kinds.
///
/// An interaction is kept when *either* side is chosen. An empty `kinds` keeps
/// every kind.
///
/// This is the shared filter behind both the flat report and the network.
pub fn selected(
    protein: &Protein,
    chosen: &AtomSet,
    kinds: &[InteractionKind],
) -> Vec<Interaction> {
    detect(protein)
        .into_iter()
        .filter(|found| kinds.is_empty() || kinds.contains(&found.kind))
        .filter(|found| touches(found, chosen))
        .collect()
}

/// The chosen interactions as flat records, in a deterministic order.
pub fn interactions(
    protein: &Protein,
    chosen: &AtomSet,
    kinds: &[InteractionKind],
) -> Vec<InteractionRecord> {
    let mut records: Vec<InteractionRecord> = selected(protein, chosen, kinds)
        .iter()
        .map(|found| record(protein, found))
        .collect();

    records.sort_by(|a, b| {
        (
            &a.kind, &a.chain_a, a.resi_a, a.icode_a, &a.chain_b, a.resi_b, a.icode_b,
        )
            .cmp(&(
                &b.kind, &b.chain_a, b.resi_a, b.icode_a, &b.chain_b, b.resi_b, b.icode_b,
            ))
    });

    records
}

fn touches(interaction: &Interaction, chosen: &AtomSet) -> bool {
    interaction
        .atoms_a
        .iter()
        .chain(&interaction.atoms_b)
        .any(|&atom| chosen.contains(atom))
}

fn record(protein: &Protein, interaction: &Interaction) -> InteractionRecord {
    let side = |atoms: &[usize]| {
        let first = &protein.atoms[atoms[0]];
        let names: Vec<&str> = atoms
            .iter()
            .map(|&i| protein.atoms[i].name.as_str())
            .collect();

        (
            first.chain_id.as_str().to_string(),
            first.residue_seq,
            first.insertion_code,
            first.residue_name.as_str().to_string(),
            names.join("+"),
        )
    };

    let (chain_a, resi_a, icode_a, resn_a, atoms_a) = side(&interaction.atoms_a);
    let (chain_b, resi_b, icode_b, resn_b, atoms_b) = side(&interaction.atoms_b);

    InteractionRecord {
        kind: interaction.kind.label().to_string(),
        chain_a,
        resi_a,
        icode_a,
        resn_a,
        atoms_a,
        chain_b,
        resi_b,
        icode_b,
        resn_b,
        atoms_b,
        distance: rounded(interaction.distance, 2),
        angle: interaction.angle.map(|angle| rounded(angle, 1)),
    }
}

/// One residue's secondary structure.
#[derive(Serialize)]
pub struct SecondaryRecord {
    pub chain: String,
    pub resi: u32,
    pub icode: Option<char>,
    pub resn: String,
    pub ss: char,
}

impl Row for SecondaryRecord {
    fn columns() -> &'static [&'static str] {
        &["chain", "resi", "icode", "resn", "ss"]
    }

    fn cells(&self) -> Vec<String> {
        vec![
            self.chain.clone(),
            self.resi.to_string(),
            code(self.icode),
            self.resn.clone(),
            self.ss.to_string(),
        ]
    }
}

/// The secondary structure of every chosen residue.
///
/// The assignment is made at load. A residue's letter comes from its C-alpha,
/// the atom the assignment is anchored on.
pub fn secondary(protein: &Protein, chosen: &AtomSet) -> Vec<SecondaryRecord> {
    residues(protein, chosen)
        .map(|(range, first)| {
            let anchor = range
                .clone()
                .find(|&i| protein.atoms[i].name == "CA")
                .unwrap_or(range.start);

            SecondaryRecord {
                chain: first.chain_id.as_str().to_string(),
                resi: first.residue_seq,
                icode: first.insertion_code,
                resn: first.residue_name.as_str().to_string(),
                ss: protein.atoms[anchor].secondary_structure.code(),
            }
        })
        .collect()
}

/// One residue's solvent-accessible surface area, absolute and relative.
#[derive(Serialize)]
pub struct SasaRecord {
    pub chain: String,
    pub resi: u32,
    pub icode: Option<char>,
    pub resn: String,
    pub sasa: f32,
    /// SASA over the residue's theoretical maximum; `None` for non-standard
    /// residues with no reference (ligands, modified residues).
    pub rsa: Option<f32>,
}

impl Row for SasaRecord {
    fn columns() -> &'static [&'static str] {
        &["chain", "resi", "icode", "resn", "sasa", "rsa"]
    }

    fn cells(&self) -> Vec<String> {
        vec![
            self.chain.clone(),
            self.resi.to_string(),
            code(self.icode),
            self.resn.clone(),
            format!("{:.2}", self.sasa),
            self.rsa.map_or(String::new(), |rsa| format!("{rsa:.3}")),
        ]
    }
}

/// Per-residue solvent-accessible surface area, in square Angstroms.
///
/// The areas are computed over the whole structure and only then narrowed to the
/// chosen residues.
///
/// The structure here is the polymer alone. Solvent and hydrogens are dropped,
/// mirroring the atom set FreeSASA is run on to produce the benchmark golden
/// files. It is therefore the set pixelfold's published SASA numbers reflect.
pub fn sasa(protein: &Protein, chosen: &AtomSet) -> Vec<SasaRecord> {
    let surface: Vec<usize> = (0..protein.atoms.len())
        .filter(|&i| bears_surface(&protein.atoms[i]))
        .collect();
    let atoms: Vec<Atom> = surface.iter().map(|&i| protein.atoms[i].clone()).collect();

    let computed = SurfaceCalculator::default().calculate_atom_sasa(&atoms);
    let mut areas = vec![0.0f32; protein.atoms.len()];
    for (slot, &atom) in surface.iter().enumerate() {
        areas[atom] = computed.get(slot).copied().unwrap_or(0.0);
    }

    let has_surface =
        |range: &std::ops::Range<usize>| range.clone().any(|i| bears_surface(&protein.atoms[i]));

    residues(protein, chosen)
        .filter(|(range, _)| has_surface(range))
        .map(|(range, first)| {
            let resn = first.residue_name.as_str().to_string();
            let absolute: f32 = range.map(|i| areas.get(i).copied().unwrap_or(0.0)).sum();
            let rsa = relative_sasa(absolute, &resn).map(|value| rounded(value, 3));

            SasaRecord {
                chain: first.chain_id.as_str().to_string(),
                resi: first.residue_seq,
                icode: first.insertion_code,
                resn,
                sasa: rounded(absolute, 2),
                rsa,
            }
        })
        .collect()
}

/// Whether an atom is part of the surface being measured: a polymer heavy atom.
///
/// This is the filter `pixelfold_validate` applies before comparing against FreeSASA.
fn bears_surface(atom: &Atom) -> bool {
    !atom.is_hetatm && !atom.element.as_str().eq_ignore_ascii_case("H")
}

/// Each chosen residue as its atom range and first atom, in file order.
fn residues<'a>(
    protein: &'a Protein,
    chosen: &'a AtomSet,
) -> impl Iterator<Item = (std::ops::Range<usize>, &'a Atom)> + 'a {
    let key = |atom: &Atom| (atom.chain_id, atom.residue_seq, atom.insertion_code);

    let mut index = 0;
    std::iter::from_fn(move || {
        while index < protein.atoms.len() {
            let start = index;
            let current = key(&protein.atoms[start]);
            while index < protein.atoms.len() && key(&protein.atoms[index]) == current {
                index += 1;
            }

            let range = start..index;
            if range.clone().any(|i| chosen.contains(i)) {
                return Some((range, &protein.atoms[start]));
            }
        }

        None
    })
}
