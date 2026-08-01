use glam::Vec3;

use crate::dssp::HBond;
use crate::fixed_str::FixedStr;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DisplayMode {
    AllAtoms,
    Backbone,
}

#[derive(Clone)]
pub struct Atom {
    pub serial: u32,
    pub name: FixedStr<4>,
    /// Element symbol (e.g. "C", "Ca", "Fe"), parsed from the structure's
    /// element column so a calcium ion is distinct from a C-alpha carbon.
    pub element: FixedStr<2>,
    pub residue_name: FixedStr<6>,
    pub residue_seq: u32,
    pub chain_id: FixedStr<4>,
    pub is_hetatm: bool,
    pub altloc: Option<char>,
    pub occupancy: f32,
    pub insertion_code: Option<char>,
    pub model_number: usize,
    pub position: Vec3,
    pub b_factor: f32,
    pub secondary_structure: SecondaryStructure,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SecondaryStructure {
    Helix,
    Sheet,
    Turn,
    Coil,
}

impl SecondaryStructure {
    /// The DSSP-style letter, which is also what the selection language's `ss`
    /// predicate matches on.
    pub fn code(&self) -> char {
        match self {
            SecondaryStructure::Helix => 'H',
            SecondaryStructure::Sheet => 'E',
            SecondaryStructure::Turn => 'T',
            SecondaryStructure::Coil => 'C',
        }
    }

    pub fn color_rgb(&self) -> (u8, u8, u8) {
        match self {
            SecondaryStructure::Helix => (255, 100, 100), // Red
            SecondaryStructure::Sheet => (100, 100, 255), // Blue
            SecondaryStructure::Turn => (255, 255, 100),  // Yellow
            SecondaryStructure::Coil => (100, 255, 100),  // Green
        }
    }

    pub fn color_256(&self) -> u8 {
        match self {
            SecondaryStructure::Helix => 196, // Red
            SecondaryStructure::Sheet => 33,  // Blue
            SecondaryStructure::Turn => 226,  // Yellow
            SecondaryStructure::Coil => 46,   // Green
        }
    }
}

#[derive(Clone)]
pub struct Protein {
    pub atoms: Vec<Atom>,
    pub title: String,
    pub hbonds: Vec<HBond>,
    /// Set when the file declares a biological assembly that differs from the
    /// deposited asymmetric unit (the coordinates shown). `None` when the
    /// deposited coordinates already are the biological unit, or when no
    /// assembly is declared.
    pub assembly: Option<BiologicalAssembly>,
    /// The chemical definitions the file carried for its own components,
    /// informing a ligand's connectivity. Empty when the file declared, as a
    /// legacy PDB file does not.
    pub components: crate::components::Dictionary,
}

/// A biological assembly whose symmetry expansion differs from the deposited
/// asymmetric unit. Present only when the primary assembly applies more than the
/// identity operator, so the deposited coordinates are a partial view of the
/// functional unit. Assembly generation itself is not yet performed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BiologicalAssembly {
    /// Number of symmetry operators the primary assembly applies. Greater than
    /// one means the deposited coordinates are only part of the biological unit.
    pub operator_count: usize,
    /// Human-readable oligomeric state stated in the file, if any (mmCIF
    /// `oligomeric_details` or PDB `REMARK 350` biological-unit line, for
    /// example "tetrameric"). `None` when the file does not state one.
    pub oligomer: Option<String>,
}

/// Calculate B-factor percentiles for normalization
pub fn calculate_bfactor_range(protein: &Protein) -> (f32, f32) {
    let mut b_factors: Vec<f32> = protein
        .atoms
        .iter()
        .map(|a| a.b_factor)
        .filter(|b| b.is_finite())
        .collect();
    if b_factors.is_empty() {
        return (0.0, 100.0);
    }
    b_factors.sort_by(f32::total_cmp);

    // Avoiding outliers
    let idx_min = (b_factors.len() as f32 * 0.02) as usize;
    let idx_max = (b_factors.len() as f32 * 0.98) as usize;

    let b_min = b_factors[idx_min.min(b_factors.len() - 1)];
    let b_max = b_factors[idx_max.min(b_factors.len() - 1)];

    (b_min, b_max)
}

/// Get indices of C-alpha atoms in the protein
pub fn get_calpha_indices(protein: &Protein) -> Vec<usize> {
    protein
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.name == "CA")
        .map(|(idx, _)| idx)
        .collect()
}

/// Get pairs of C-alpha indices that should be connected
/// Only connects sequential C-alphas in the same chain within distance threshold
pub fn get_calpha_connections(protein: &Protein, ca_indices: &[usize]) -> Vec<(usize, usize)> {
    let mut connections = Vec::new();
    let distance_threshold = 4.2; // Angstroms (typical CA-CA distance is ~3.8Å)

    // Sort C-alphas by chain first, then by residue sequence
    let mut sorted_cas: Vec<(FixedStr<4>, u32, usize)> = ca_indices
        .iter()
        .map(|&idx| {
            let atom = &protein.atoms[idx];
            (atom.chain_id, atom.residue_seq, idx)
        })
        .collect();

    sorted_cas.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    for i in 0..sorted_cas.len().saturating_sub(1) {
        let (chain1, res1, idx1) = &sorted_cas[i];
        let (chain2, res2, idx2) = &sorted_cas[i + 1];

        // Only connect if same chain
        if chain1 == chain2 {
            // Exactly sequential (res2 = res1 + 1)
            if *res2 == *res1 + 1 {
                let atom1 = &protein.atoms[*idx1];
                let atom2 = &protein.atoms[*idx2];
                let distance = (atom2.position - atom1.position).length();

                if distance >= 2.5 && distance <= distance_threshold {
                    connections.push((*idx1, *idx2));
                }
            }
        }
    }

    connections
}

/// Policy for resolving atoms modelled in multiple alternate locations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum AltlocPolicy {
    /// Keep the highest-occupancy conformer per atom; ties prefer altloc 'A'.
    #[default]
    Occupancy,
    /// Keep only altloc 'A' and atoms with no altloc.
    A,
    /// Keep only altloc 'B' and atoms with no altloc.
    B,
    /// Keep every conformer.
    All,
}

/// Apply an alternate-location policy. Atoms with no altloc are always kept.
pub fn filter_altlocs(atoms: Vec<Atom>, policy: AltlocPolicy) -> Vec<Atom> {
    match policy {
        AltlocPolicy::All => atoms,
        AltlocPolicy::A => atoms
            .into_iter()
            .filter(|a| a.altloc.is_none() || a.altloc == Some('A'))
            .collect(),
        AltlocPolicy::B => atoms
            .into_iter()
            .filter(|a| a.altloc.is_none() || a.altloc == Some('B'))
            .collect(),
        AltlocPolicy::Occupancy => resolve_occupancy(atoms),
    }
}

type ResidueKey = (FixedStr<4>, u32, Option<char>);

fn residue_of(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

fn is_backbone(name: &FixedStr<4>) -> bool {
    matches!(name.as_str(), "N" | "CA" | "C" | "O")
}

/// The label with the greatest summed occupancy, ties preferring the earlier
/// label. `labels` is never empty at the call sites.
fn pick_label(labels: &std::collections::HashMap<char, f32>) -> char {
    labels
        .iter()
        .max_by(|(la, oa), (lb, ob)| oa.total_cmp(ob).then_with(|| lb.cmp(la)))
        .map(|(&label, _)| label)
        .expect("a residue in the occupancy map has at least one label")
}

/// Resolve alternate locations to a single conformer per residue.
///
/// Choosing a winner independently per atom can leave a residue chimeric (its
/// atoms drawn from different conformers) which places a sidechain rotor about an
/// axis that is a bond in neither conformer. Instead one label wins the whole
/// residue, decided by the backbone, by greatest summed occupancy (ties prefer
/// the earlier label).
fn resolve_occupancy(atoms: Vec<Atom>) -> Vec<Atom> {
    use std::collections::HashMap;

    let mut overall: HashMap<ResidueKey, HashMap<char, f32>> = HashMap::new();
    let mut backbone: HashMap<ResidueKey, HashMap<char, f32>> = HashMap::new();
    for atom in &atoms {
        let Some(label) = atom.altloc else { continue };
        let residue = residue_of(atom);
        *overall
            .entry(residue)
            .or_default()
            .entry(label)
            .or_default() += atom.occupancy;
        if is_backbone(&atom.name) {
            *backbone
                .entry(residue)
                .or_default()
                .entry(label)
                .or_default() += atom.occupancy;
        }
    }
    let chosen: HashMap<ResidueKey, char> = overall
        .iter()
        .map(|(residue, all_labels)| {
            (
                *residue,
                pick_label(backbone.get(residue).unwrap_or(all_labels)),
            )
        })
        .collect();

    let mut best: HashMap<(ResidueKey, FixedStr<4>), usize> = HashMap::new();
    let mut kept: Vec<Atom> = Vec::new();
    for atom in atoms {
        if atom.altloc.is_none() {
            kept.push(atom);
            continue;
        }

        let want = chosen.get(&residue_of(&atom)).copied();
        let key = (residue_of(&atom), atom.name);
        match best.get(&key).copied() {
            None => {
                best.insert(key, kept.len());
                kept.push(atom);
            }
            Some(idx) if prefer(&atom, &kept[idx], want) => kept[idx] = atom,
            Some(_) => {}
        }
    }

    kept
}

/// Whether `new` should replace `existing` as a residue's kept copy of an atom.
/// The chosen conformer always wins; where neither copy is the chosen label (the
/// atom is only modelled in other conformers) the more-occupied copy wins, ties
/// preferring the earlier label.
fn prefer(new: &Atom, existing: &Atom, want: Option<char>) -> bool {
    match (new.altloc == want, existing.altloc == want) {
        (true, false) => true,
        (false, true) => false,
        _ => {
            new.occupancy > existing.occupancy
                || (new.occupancy == existing.occupancy && new.altloc < existing.altloc)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom_with_bfactor(b: f32) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new("CA"),
            element: FixedStr::new("C"),
            residue_name: FixedStr::new("ALA"),
            residue_seq: 1,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: b,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    #[test]
    fn bfactor_range_does_not_panic_on_nan() {
        // Real PDB files contain NaN coordinates and B-factors; the depth/value
        // sorts must not panic on them.
        let protein = Protein {
            atoms: vec![
                atom_with_bfactor(10.0),
                atom_with_bfactor(f32::NAN),
                atom_with_bfactor(50.0),
            ],
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        };

        let (b_min, b_max) = calculate_bfactor_range(&protein);
        assert!(b_min.is_finite() && b_max.is_finite());
    }

    fn altloc_atom(name: &str, altloc: Option<char>, occupancy: f32) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new("C"),
            residue_name: FixedStr::new("SER"),
            residue_seq: 1,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc,
            occupancy,
            insertion_code: None,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    #[test]
    fn altloc_occupancy_keeps_highest_and_breaks_ties_with_a() {
        // Same atom, two conformers: B has higher occupancy so it wins.
        let kept = filter_altlocs(
            vec![
                altloc_atom("OG", Some('A'), 0.4),
                altloc_atom("OG", Some('B'), 0.6),
            ],
            AltlocPolicy::Occupancy,
        );
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].altloc, Some('B'));

        // On an occupancy tie, altloc 'A' wins even if listed second.
        let tie = filter_altlocs(
            vec![
                altloc_atom("OG", Some('B'), 0.5),
                altloc_atom("OG", Some('A'), 0.5),
            ],
            AltlocPolicy::Occupancy,
        );
        assert_eq!(tie.len(), 1);
        assert_eq!(tie[0].altloc, Some('A'));
    }

    #[test]
    fn occupancy_keeps_one_consistent_conformer_per_residue() {
        // CB alone prefers A (0.6 > 0.4) while OG alone prefers C (0.7 > 0.3):
        // a per-atom winner would keep CB(A) with OG(C), a chimera. The residue's
        // summed occupancy is A = 0.9 vs C = 1.1, so it resolves wholly to C.
        let kept = filter_altlocs(
            vec![
                altloc_atom("CB", Some('A'), 0.6),
                altloc_atom("CB", Some('C'), 0.4),
                altloc_atom("OG", Some('A'), 0.3),
                altloc_atom("OG", Some('C'), 0.7),
            ],
            AltlocPolicy::Occupancy,
        );
        assert_eq!(kept.len(), 2);
        assert!(
            kept.iter().all(|a| a.altloc == Some('C')),
            "the whole residue takes one conformer"
        );
    }

    #[test]
    fn occupancy_lets_the_backbone_decide_over_the_sidechain() {
        // Summed over all atoms, conformer B wins (1.55 vs 1.45); but the
        // backbone CA favours A (0.55 vs 0.45). The residue must follow the
        // backbone so a disordered sidechain never shifts the main chain.
        let kept = filter_altlocs(
            vec![
                altloc_atom("CA", Some('A'), 0.55),
                altloc_atom("CA", Some('B'), 0.45),
                altloc_atom("CB", Some('A'), 0.45),
                altloc_atom("CB", Some('B'), 0.55),
                altloc_atom("CG", Some('A'), 0.45),
                altloc_atom("CG", Some('B'), 0.55),
            ],
            AltlocPolicy::Occupancy,
        );
        assert!(
            kept.iter().all(|a| a.altloc == Some('A')),
            "the backbone decides the residue conformer"
        );
    }

    #[test]
    fn occupancy_keeps_atoms_absent_from_the_chosen_conformer() {
        // The residue resolves to conformer A (0.7 > 0.3), but its guanidinium
        // tip is modelled only in B; that atom must survive rather than vanish.
        let kept = filter_altlocs(
            vec![
                altloc_atom("CB", Some('A'), 0.7),
                altloc_atom("CB", Some('B'), 0.3),
                altloc_atom("NH1", Some('B'), 0.3),
            ],
            AltlocPolicy::Occupancy,
        );
        let cb = kept.iter().find(|a| a.name.as_str() == "CB").unwrap();
        let nh1 = kept.iter().find(|a| a.name.as_str() == "NH1");
        assert_eq!(cb.altloc, Some('A'), "CB takes the chosen conformer");
        assert_eq!(
            nh1.map(|a| a.altloc),
            Some(Some('B')),
            "an atom the chosen conformer omits falls back to its own copy"
        );
    }

    #[test]
    fn altloc_a_keeps_a_and_none_drops_b() {
        let kept = filter_altlocs(
            vec![
                altloc_atom("OG", Some('A'), 0.5),
                altloc_atom("OG", Some('B'), 0.5),
                altloc_atom("N", None, 1.0),
            ],
            AltlocPolicy::A,
        );
        assert_eq!(kept.len(), 2);
        assert!(kept.iter().all(|a| a.altloc != Some('B')));
    }

    #[test]
    fn altloc_all_keeps_every_conformer() {
        let kept = filter_altlocs(
            vec![
                altloc_atom("OG", Some('A'), 0.5),
                altloc_atom("OG", Some('B'), 0.5),
            ],
            AltlocPolicy::All,
        );
        assert_eq!(kept.len(), 2);
    }
}
