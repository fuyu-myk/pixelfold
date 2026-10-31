//! Matching the residues of two structures of one molecule by sequence.
//!
//! A crystal structure numbers its residues as its authors chose, insertion
//! codes included, and leaves out those it did not resolve; co-folding output
//! numbers every residue from 1 and names its chains A, B, C. The same residue is
//! found in both by aligning the residue names of each pair of chains. Two chains
//! pair when at least half the residues of the shorter one, and no fewer than
//! [`MIN_IDENTICAL`], align with residues of the same name, so a point mutant or
//! a modified residue still pairs and an unrelated chain does not. A chain
//! shorter than that, a peptide, pairs only where it matches whole and
//! ungapped, since a few residues align somewhere in any long chain.

use std::collections::HashMap;
use std::ops::Add;

use crate::fixed_str::FixedStr;
use crate::structure::{Atom, ChainId, Protein};

/// A residue's identity: chain, number and insertion code.
pub type ResidueKey = (ChainId, i32, Option<char>);

/// A residue: its identity and its name.
pub type Residue = (ResidueKey, FixedStr<6>);

/// A polymer chain: its residues in the order the structure gives them.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Chain {
    pub id: ChainId,
    residues: Vec<Residue>,
}

impl Chain {
    pub fn contains(&self, residue: &ResidueKey) -> bool {
        self.residues.iter().any(|(key, _)| key == residue)
    }

    fn names(&self) -> Vec<FixedStr<6>> {
        self.residues.iter().map(|&(_, name)| name).collect()
    }
}

/// The fewest identical residues that pair two chains, unless the shorter one
/// is shorter still.
pub const MIN_IDENTICAL: usize = 10;

/// How far below a pose chain's best alignment, as the share of the shorter
/// chain that is identical, another reference chain may fall and still
/// compete for it. Copies of one chain in a crystal align alike, however many
/// residues each resolves; a paralog falls well below.
const IDENTITY_MARGIN: f32 = 0.1;

/// The polymer chains of `protein`, in the order they first appear.
pub fn polymer_chains(protein: &Protein) -> Vec<Chain> {
    let mut chains: Vec<Chain> = Vec::new();
    let mut index: HashMap<ChainId, usize> = HashMap::new();
    for atom in protein.atoms.iter().filter(|atom| atom.is_polymer()) {
        let at = *index.entry(atom.chain_id).or_insert_with(|| {
            chains.push(Chain {
                id: atom.chain_id,
                residues: Vec::new(),
            });
            chains.len() - 1
        });
        let residues = &mut chains[at].residues;
        let key = residue_key(atom);
        if residues.last().is_none_or(|(last, _)| *last != key) {
            residues.push((key, atom.residue_name));
        }
    }

    chains
}

/// One chain of a pose aligned with one chain of the reference.
#[derive(Clone, Debug)]
pub struct ChainPair {
    pose: usize,
    reference: usize,
    residues: HashMap<ResidueKey, Residue>,
    identical: usize,
    /// The share of the shorter chain that aligns identically.
    identity: f32,
}

impl ChainPair {
    /// Chain `pose_index` of a pose paired with chain `reference_index` of the
    /// reference through `alignment`, or `None` when the two are not the same
    /// molecule.
    fn new(
        (pose_index, pose): (usize, &Chain),
        (reference_index, reference): (usize, &Chain),
        alignment: &Alignment,
    ) -> Option<ChainPair> {
        let shorter = pose.residues.len().min(reference.residues.len());
        let same = if shorter < MIN_IDENTICAL {
            alignment.identical == shorter
                && alignment
                    .pairs
                    .windows(2)
                    .all(|step| step[1].0 == step[0].0 + 1 && step[1].1 == step[0].1 + 1)
        } else {
            alignment.identical >= MIN_IDENTICAL && alignment.identical * 2 >= shorter
        };
        if !same || shorter == 0 {
            return None;
        }

        Some(ChainPair {
            pose: pose_index,
            reference: reference_index,
            residues: alignment
                .pairs
                .iter()
                .map(|&(p, r)| (pose.residues[p].0, reference.residues[r]))
                .collect(),
            identical: alignment.identical,
            identity: alignment.identical as f32 / shorter as f32,
        })
    }

    /// The reference residue a pose residue of this chain aligns with.
    pub fn reference_residue(&self, residue: &ResidueKey) -> Option<Residue> {
        self.residues.get(residue).copied()
    }
}

/// What a chain pairing is worth, compared field by field.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
struct Worth {
    kept: usize,
    preferred: usize,
    identical: usize,
    same_name: usize,
}

impl Add for Worth {
    type Output = Worth;

    fn add(self, other: Worth) -> Worth {
        Worth {
            kept: self.kept + other.kept,
            preferred: self.preferred + other.preferred,
            identical: self.identical + other.identical,
            same_name: self.same_name + other.same_name,
        }
    }
}

/// Which residue of the reference each polymer residue of a pose stands for.
#[derive(Clone, Debug, Default)]
pub struct ResidueMap {
    /// Every polymer residue of the pose, with its counterpart if it has one.
    to_reference: HashMap<ResidueKey, Option<Residue>>,
    /// Every polymer residue of the reference, with its counterpart if it has
    /// one.
    to_pose: HashMap<ResidueKey, Option<ResidueKey>>,
    /// The chains paired, pose first.
    chains: Vec<(ChainId, ChainId)>,
    /// The reference chains some pose chain aligns with, paired or not.
    alignable: Vec<ChainId>,
}

impl ResidueMap {
    /// The reference residue the pose residue `residue` aligns with.
    pub fn reference_residue(&self, residue: &ResidueKey) -> Option<Residue> {
        self.to_reference.get(residue).copied().flatten()
    }

    /// The pose residue the reference residue `residue` aligns with.
    pub fn pose_residue(&self, residue: &ResidueKey) -> Option<ResidueKey> {
        self.to_pose.get(residue).copied().flatten()
    }

    /// True when `residue` is a polymer residue of the pose.
    pub fn is_pose_polymer(&self, residue: &ResidueKey) -> bool {
        self.to_reference.contains_key(residue)
    }

    /// True when `residue` is a polymer residue of the reference.
    pub fn is_reference_polymer(&self, residue: &ResidueKey) -> bool {
        self.to_pose.contains_key(residue)
    }

    /// The chains paired, each as its pose chain and its reference chain.
    pub fn chains(&self) -> impl Iterator<Item = &(ChainId, ChainId)> {
        self.chains.iter()
    }

    /// True when some pose chain aligns with the reference chain `chain`, even
    /// if another reference chain took it.
    pub fn is_alignable(&self, chain: &ChainId) -> bool {
        self.alignable.contains(chain)
    }
}

/// The pose's chains paired with the reference's they align with, each
/// reference chain used once.
///
/// Where a pose chain could pair with several reference chains, as the chains
/// of a homodimer can, the pairing that `kept` scores highest wins: `kept` counts
/// what a pairing keeps of the reference, such as the contacts a ligand makes.
/// Only reference chains that align about as well as the pose chain's best
/// compete, so `kept` chooses among copies of one chain and never favours a
/// paralog. Ties go to the `preferred` reference chains, then to the better
/// alignment, then to the chain of the same name.
pub fn residue_map(
    pose: &[Chain],
    reference: &[Chain],
    preferred: &[ChainId],
    kept: impl Fn(&ChainPair) -> usize,
) -> ResidueMap {
    let pairs = candidate_pairs(pose, reference);
    let worth: Vec<Worth> = pairs
        .iter()
        .map(|pair| Worth {
            kept: kept(pair),
            preferred: usize::from(preferred.contains(&reference[pair.reference].id)),
            identical: pair.identical,
            same_name: usize::from(pose[pair.pose].id == reference[pair.reference].id),
        })
        .collect();
    let find = |p: usize, r: usize| {
        pairs
            .iter()
            .position(|pair| pair.pose == p && pair.reference == r)
    };

    // Only the chains whose pairing changes what is kept need the full search;
    // the rest pair greedily, best alignment first.
    let contested: Vec<usize> = (0..pose.len())
        .filter(|&p| (0..pairs.len()).any(|i| pairs[i].pose == p && worth[i].kept > 0))
        .collect();
    let mut chosen: Vec<usize> = best_pairing(contested.len(), reference.len(), |row, r| {
        find(contested[row], r).map(|i| worth[i])
    })
    .into_iter()
    .filter_map(|(row, r)| find(contested[row], r))
    .collect();

    let mut rest: Vec<usize> = (0..pairs.len())
        .filter(|&i| !contested.contains(&pairs[i].pose))
        .collect();
    rest.sort_by(|&a, &b| worth[b].cmp(&worth[a]));
    for i in rest {
        let taken = chosen
            .iter()
            .any(|&c| pairs[c].pose == pairs[i].pose || pairs[c].reference == pairs[i].reference);
        if !taken {
            chosen.push(i);
        }
    }
    chosen.sort_by_key(|&i| pairs[i].pose);

    let mut alignable: Vec<ChainId> = pairs
        .iter()
        .map(|pair| reference[pair.reference].id)
        .collect();
    alignable.sort();
    alignable.dedup();
    let mut map = ResidueMap {
        to_reference: residue_keys(pose).map(|key| (key, None)).collect(),
        to_pose: residue_keys(reference).map(|key| (key, None)).collect(),
        chains: Vec::new(),
        alignable,
    };
    for pair in chosen.into_iter().map(|i| &pairs[i]) {
        for (&pose_residue, &reference_residue) in &pair.residues {
            map.to_reference
                .insert(pose_residue, Some(reference_residue));
            map.to_pose.insert(reference_residue.0, Some(pose_residue));
        }
        map.chains
            .push((pose[pair.pose].id, reference[pair.reference].id));
    }

    map
}

/// Every pair of a pose chain and a reference chain that are the same
/// molecule, less those aligning well below the pose chain's best. Each
/// distinct pair of sequences is aligned once, which is what keeps a homomer's
/// many pairs cheap.
fn candidate_pairs(pose: &[Chain], reference: &[Chain]) -> Vec<ChainPair> {
    let pose_names: Vec<Vec<FixedStr<6>>> = pose.iter().map(Chain::names).collect();
    let reference_names: Vec<Vec<FixedStr<6>>> = reference.iter().map(Chain::names).collect();
    let mut alignments: HashMap<(&[FixedStr<6>], &[FixedStr<6>]), Alignment> = HashMap::new();

    let mut pairs = Vec::new();
    for (p, chain) in pose.iter().enumerate() {
        for (r, other) in reference.iter().enumerate() {
            let alignment = alignments
                .entry((pose_names[p].as_slice(), reference_names[r].as_slice()))
                .or_insert_with(|| align(&pose_names[p], &reference_names[r]));
            pairs.extend(ChainPair::new((p, chain), (r, other), alignment));
        }
    }

    let best: Vec<f32> = (0..pose.len())
        .map(|p| {
            pairs
                .iter()
                .filter(|pair| pair.pose == p)
                .map(|pair| pair.identity)
                .fold(0.0, f32::max)
        })
        .collect();
    pairs.retain(|pair| pair.identity >= best[pair.pose] - IDENTITY_MARGIN);

    pairs
}

fn residue_keys(chains: &[Chain]) -> impl Iterator<Item = ResidueKey> + '_ {
    chains
        .iter()
        .flat_map(|chain| chain.residues.iter().map(|&(key, _)| key))
}

/// The pairs of rows and columns, each row and each column in at most one,
/// whose weights sum highest; `weight` is `None` where a row and a column cannot
/// pair. The search tries each row's columns heaviest first, so of equal sums
/// the first found wins, and abandons a branch that cannot beat the best so far
/// even if every row left took its heaviest column. It stops at
/// [`SEARCH_LIMIT`] complete pairings.
pub(crate) fn best_pairing<W>(
    rows: usize,
    columns: usize,
    weight: impl Fn(usize, usize) -> Option<W>,
) -> Vec<(usize, usize)>
where
    W: Copy + Ord + Default + Add<Output = W>,
{
    let candidates: Vec<Vec<(usize, W)>> = (0..rows)
        .map(|row| {
            let mut options: Vec<(usize, W)> = (0..columns)
                .filter_map(|column| Some((column, weight(row, column)?)))
                .collect();
            options.sort_by_key(|&(_, weight)| std::cmp::Reverse(weight));
            options
        })
        .collect();
    // Adding one weight to two sums keeps their order, so the heaviest weights
    // of the rows left bound what a branch can still reach.
    let mut bounds = vec![W::default(); rows + 1];
    for row in (0..rows).rev() {
        let heaviest = candidates[row]
            .first()
            .map_or(W::default(), |&(_, weight)| weight.max(W::default()));
        bounds[row] = bounds[row + 1] + heaviest;
    }
    let mut search = Search {
        candidates: &candidates,
        bounds: &bounds,
        used: vec![false; columns],
        chosen: Vec::new(),
        best: None,
        budget: SEARCH_LIMIT,
    };
    search.run(0, W::default());

    search.best.map(|(_, pairs)| pairs).unwrap_or_default()
}

/// How many complete pairings [`best_pairing`] weighs before settling for the
/// best so far.
const SEARCH_LIMIT: usize = 100_000;

/// A depth-first search over pairings, one row at a time.
struct Search<'a, W> {
    candidates: &'a [Vec<(usize, W)>],
    /// The most each row onwards can add.
    bounds: &'a [W],
    used: Vec<bool>,
    chosen: Vec<(usize, usize)>,
    best: Option<(W, Vec<(usize, usize)>)>,
    budget: usize,
}

impl<W: Copy + Ord + Add<Output = W>> Search<'_, W> {
    fn run(&mut self, row: usize, total: W) {
        let hopeless = self
            .best
            .as_ref()
            .is_some_and(|(best, _)| total + self.bounds[row] <= *best);
        if self.budget == 0 || hopeless {
            return;
        }
        let candidates = self.candidates;
        let Some(options) = candidates.get(row) else {
            self.budget -= 1;
            if self.best.as_ref().is_none_or(|(best, _)| total > *best) {
                self.best = Some((total, self.chosen.clone()));
            }
            return;
        };

        for &(column, weight) in options {
            if self.used[column] {
                continue;
            }
            self.used[column] = true;
            self.chosen.push((row, column));
            self.run(row + 1, total + weight);
            self.chosen.pop();
            self.used[column] = false;
        }
        self.run(row + 1, total);
    }
}

const MATCH: i32 = 2;
const MISMATCH: i32 = -1;
/// The cost of opening a gap, charged on top of extending it by its first
/// residue.
const GAP_OPEN: i32 = -4;
/// Nothing: a crystal may leave out a loop of any length, and a long gap
/// charged by the residue costs more than the few residues resolved after it
/// are worth, so they would be aligned with the loop instead.
const GAP_EXTEND: i32 = 0;
const UNREACHABLE: i32 = i32::MIN / 2;

/// The sequences an alignment step advances.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Step {
    Both,
    OnlyA,
    OnlyB,
}

const STEPS: [Step; 3] = [Step::Both, Step::OnlyA, Step::OnlyB];

/// The steps each state of a cell came from, two bits apiece.
fn packed(steps: [Step; 3]) -> u8 {
    steps
        .iter()
        .enumerate()
        .fold(0, |bits, (state, &step)| bits | (step as u8) << (2 * state))
}

fn unpacked(bits: u8, state: Step) -> Step {
    STEPS[usize::from((bits >> (2 * state as u8)) & 0b11)]
}

/// The best of `scores`, one per step ending there, each plus its cost: the
/// score and the step it came from, the earliest of equals.
fn best_of(scores: [i32; 3], costs: [i32; 3]) -> (i32, Step) {
    (0..3).fold((UNREACHABLE, Step::Both), |best, s| {
        let score = scores[s] + costs[s];
        if score > best.0 {
            (score, STEPS[s])
        } else {
            best
        }
    })
}

/// Two sequences aligned: the index pairs aligned, in order, how many of them
/// are identical, and the score.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Alignment {
    pairs: Vec<(usize, usize)>,
    identical: usize,
    score: i32,
}

/// The best alignment of `a` and `b` with affine gap costs (Gotoh, J Mol Biol
/// 162:705, 1982), charging nothing for a gap at either end, since one structure
/// may resolve or construct more of a chain than the other.
fn align<T: PartialEq>(a: &[T], b: &[T]) -> Alignment {
    let width = b.len() + 1;
    let open = GAP_OPEN + GAP_EXTEND;
    // Every cell of the top row and the left column starts an alignment, which
    // is what leaves a leading gap free.
    let start = [0, UNREACHABLE, UNREACHABLE];
    let mut previous = vec![start; width];
    let mut current = vec![start; width];
    let mut came_from = vec![0u8; (a.len() + 1) * width];
    let mut best = (0, 0, 0, Step::Both);

    for i in 1..=a.len() {
        for j in 1..width {
            let substitution = if a[i - 1] == b[j - 1] {
                MATCH
            } else {
                MISMATCH
            };
            let (both, from_both) = best_of(previous[j - 1], [substitution; 3]);
            let (only_a, from_a) = best_of(previous[j], [open, GAP_EXTEND, open]);
            let (only_b, from_b) = best_of(current[j - 1], [open, open, GAP_EXTEND]);
            current[j] = [both, only_a, only_b];
            came_from[i * width + j] = packed([from_both, from_a, from_b]);

            // A trailing gap is free too: the alignment may end on the last row
            // or the last column.
            if i == a.len() || j == b.len() {
                for (s, &score) in current[j].iter().enumerate() {
                    if score > best.0 {
                        best = (score, i, j, STEPS[s]);
                    }
                }
            }
        }
        std::mem::swap(&mut previous, &mut current);
    }

    let (score, mut i, mut j, mut step) = best;
    let mut pairs = Vec::new();
    while i > 0 && j > 0 {
        let from = unpacked(came_from[i * width + j], step);
        match step {
            Step::Both => {
                pairs.push((i - 1, j - 1));
                i -= 1;
                j -= 1;
            }
            Step::OnlyA => i -= 1,
            Step::OnlyB => j -= 1,
        }
        step = from;
    }
    pairs.reverse();
    let identical = pairs.iter().filter(|&&(i, j)| a[i] == b[j]).count();

    Alignment {
        pairs,
        identical,
        score,
    }
}

pub(crate) fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Molecule, SecondaryStructure};
    use glam::Vec3;

    fn letters(text: &str) -> Vec<char> {
        text.chars().collect()
    }

    fn aligned(a: &str, b: &str) -> Alignment {
        align(&letters(a), &letters(b))
    }

    #[test]
    fn identical_sequences_align_residue_for_residue() {
        let found = aligned("ACDEF", "ACDEF");

        assert_eq!(found.pairs, (0..5).map(|i| (i, i)).collect::<Vec<_>>());
        assert_eq!(found.identical, 5);
    }

    #[test]
    fn residues_one_structure_leaves_unresolved_open_one_gap() {
        // The second lacks G and H, a loop its crystal did not resolve.
        let found = aligned("ACDEFGHIKL", "ACDEFIKL");

        assert_eq!(found.identical, 8);
        assert!(found.pairs.contains(&(4, 4)));
        assert!(found.pairs.contains(&(7, 5)));
    }

    #[test]
    fn residues_beyond_either_end_cost_nothing() {
        // A full-length prediction against a crystal of its middle.
        let found = aligned("MSKACDEFGHIKLMNPQ", "ACDEFGHIK");

        assert_eq!(found.pairs.first(), Some(&(3, 0)));
        assert_eq!(found.identical, 9);
        assert_eq!(found.score, aligned("ACDEFGHIK", "ACDEFGHIK").score);
    }

    #[test]
    fn a_substitution_aligns_rather_than_opening_gaps() {
        let found = aligned("ACDEFGHIK", "ACDWFGHIK");

        assert_eq!(found.pairs.len(), 9);
        assert_eq!(found.identical, 8);
    }

    fn atom(chain: &str, seq: i32, icode: Option<char>, residue: &str, name: &str) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(&name[..1]),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new(chain),
            is_hetatm: matches!(residue, "HOH" | "ZN"),
            molecule: Molecule::Unknown,
            altloc: None,
            occupancy: 1.0,
            insertion_code: icode,
            model_number: 1,
            position: Vec3::ZERO,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    /// A structure of one alpha carbon per residue: each chain as its name, its
    /// first number and its residues, the numbers running on from the first.
    fn chains_of(chains: &[(&str, i32, &str)]) -> Protein {
        let atoms = chains
            .iter()
            .flat_map(|&(chain, first, residues)| {
                residues
                    .split_whitespace()
                    .zip(first..)
                    .map(move |(residue, seq)| atom(chain, seq, None, residue, "CA"))
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

    fn key(chain: &str, seq: i32, icode: Option<char>) -> ResidueKey {
        (FixedStr::new(chain), seq, icode)
    }

    #[test]
    fn polymer_chains_hold_each_residue_once_in_order() {
        let mut structure = chains_of(&[("H", 16, "ILE VAL GLY")]);
        // A second atom of the first residue, a residue under an insertion code,
        // a zinc ion and a water: only the polymer residues count.
        structure.atoms.insert(1, atom("H", 16, None, "ILE", "CB"));
        structure.atoms.push(atom("H", 18, Some('A'), "ASP", "CA"));
        structure.atoms.push(atom("H", 301, None, "ZN", "ZN"));
        structure.atoms.push(atom("H", 401, None, "HOH", "O"));

        let chains = polymer_chains(&structure);
        assert_eq!(chains.len(), 1);
        let residues: Vec<ResidueKey> = chains[0].residues.iter().map(|(key, _)| *key).collect();
        assert_eq!(
            residues,
            [
                key("H", 16, None),
                key("H", 17, None),
                key("H", 18, None),
                key("H", 18, Some('A'))
            ]
        );
    }

    #[test]
    fn chains_pair_by_sequence_whatever_their_names_and_numbers() {
        let reference = chains_of(&[
            ("L", 1, "THR PHE GLY SER GLY GLU ALA ASP CYS"),
            (
                "H",
                16,
                "ILE VAL GLU GLY SER ASP ALA GLU ILE GLY MET SER PRO TRP",
            ),
        ]);
        // Co-folding order and numbering: the heavy chain first, as A, and both
        // from 1, the light chain one residue longer at its start.
        let pose = chains_of(&[
            (
                "A",
                1,
                "ILE VAL GLU GLY SER ASP ALA GLU ILE GLY MET SER PRO TRP",
            ),
            ("B", 1, "MET THR PHE GLY SER GLY GLU ALA ASP CYS"),
        ]);

        let map = residue_map(
            &polymer_chains(&pose),
            &polymer_chains(&reference),
            &[],
            |_| 0,
        );

        let name = FixedStr::new;
        assert_eq!(
            map.reference_residue(&key("A", 1, None)),
            Some((key("H", 16, None), name("ILE")))
        );
        assert_eq!(
            map.reference_residue(&key("B", 2, None)),
            Some((key("L", 1, None), name("THR")))
        );
        assert_eq!(map.reference_residue(&key("B", 1, None)), None);
        assert!(map.is_pose_polymer(&key("B", 1, None)));
        assert_eq!(
            map.pose_residue(&key("H", 29, None)),
            Some(key("A", 14, None))
        );
        let chains: Vec<(ChainId, ChainId)> = map.chains().copied().collect();
        let chain = FixedStr::new;
        assert_eq!(chains, [(chain("A"), chain("H")), (chain("B"), chain("L"))]);
    }

    #[test]
    fn an_unrelated_chain_pairs_with_nothing() {
        let reference = chains_of(&[("A", 1, "ALA CYS ASP GLU PHE GLY HIS ILE")]);
        let pose = chains_of(&[("A", 1, "TRP TRP TYR TYR MET MET PRO PRO")]);

        let map = residue_map(
            &polymer_chains(&pose),
            &polymer_chains(&reference),
            &[],
            |_| 0,
        );

        assert_eq!(map.chains().count(), 0);
        assert_eq!(map.reference_residue(&key("A", 1, None)), None);
        assert!(map.is_reference_polymer(&key("A", 1, None)));
    }

    const MONOMER: &str = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU";

    #[test]
    fn of_identical_chains_the_pairing_that_keeps_the_most_wins() {
        let reference = chains_of(&[("A", 1, MONOMER), ("B", 1, MONOMER)]);
        let pose = chains_of(&[("A", 1, MONOMER), ("B", 1, MONOMER)]);
        let (pose_chains, reference_chains) = (polymer_chains(&pose), polymer_chains(&reference));

        // Only the crossed pairing keeps anything, as when a co-folding tool
        // names the chains of a homodimer the other way round.
        let crossed = |pair: &ChainPair| usize::from(pair.pose != pair.reference);
        let map = residue_map(&pose_chains, &reference_chains, &[], crossed);
        assert_eq!(
            map.pose_residue(&key("A", 3, None)),
            Some(key("B", 3, None))
        );
        assert_eq!(
            map.pose_residue(&key("B", 3, None)),
            Some(key("A", 3, None))
        );

        // With nothing to tell them apart, each keeps its own name.
        let map = residue_map(&pose_chains, &reference_chains, &[], |_| 0);
        assert_eq!(
            map.pose_residue(&key("A", 3, None)),
            Some(key("A", 3, None))
        );
        assert_eq!(
            map.pose_residue(&key("B", 3, None)),
            Some(key("B", 3, None))
        );
    }

    #[test]
    fn a_pose_with_more_copies_than_the_reference_leaves_one_unpaired() {
        let reference = chains_of(&[("A", 1, MONOMER)]);
        let pose = chains_of(&[("A", 1, MONOMER), ("B", 1, MONOMER)]);
        let kept_by_b = |pair: &ChainPair| usize::from(pair.pose == 1);

        let map = residue_map(
            &polymer_chains(&pose),
            &polymer_chains(&reference),
            &[],
            kept_by_b,
        );

        let chains: Vec<(ChainId, ChainId)> = map.chains().copied().collect();
        assert_eq!(chains, [(FixedStr::new("B"), FixedStr::new("A"))]);
        assert_eq!(map.reference_residue(&key("A", 1, None)), None);
    }

    #[test]
    fn the_best_pairing_maximises_the_sum_not_each_row() {
        // Row 0 prefers column 0, but giving it column 1 frees column 0 for row
        // 1, which can take nothing else.
        let weights = [[Some(5), Some(4)], [Some(3), None]];

        let found = best_pairing(2, 2, |row, column| weights[row][column]);

        assert_eq!(found, [(0, 1), (1, 0)]);
    }

    /// `n` letters drawn from the twenty amino-acid codes by a fixed generator.
    fn random_letters(seed: u64, n: usize) -> Vec<char> {
        let alphabet: Vec<char> = "ACDEFGHIKLMNPQRSTVWY".chars().collect();
        let mut state = seed;
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                alphabet[(state >> 33) as usize % alphabet.len()]
            })
            .collect()
    }

    /// Letters as residue names, for [`chains_of`].
    fn spaced(letters: &[char]) -> String {
        letters
            .iter()
            .map(char::to_string)
            .collect::<Vec<_>>()
            .join(" ")
    }

    #[test]
    fn a_short_segment_after_a_long_missing_loop_still_aligns() {
        // A full-length chain against a crystal that left out 150 residues, as a
        // GPCR structure leaves out its third intracellular loop, with only 15
        // resolved after it.
        let (head, missing, tail) = (
            random_letters(1, 60),
            random_letters(2, 150),
            random_letters(3, 15),
        );
        let full = [head.clone(), missing, tail.clone()].concat();
        let resolved = [head, tail].concat();

        let found = align(&full, &resolved);

        assert_eq!(found.identical, 75);
        assert!(found.pairs.contains(&(210, 60)), "{:?}", found.pairs.last());
    }

    #[test]
    fn a_short_chain_pairs_only_with_an_exact_match() {
        let peptide = letters("GKDALS");
        let around = random_letters(4, 100);
        let with = |inside: &str| [around.clone(), letters(inside), around.clone()].concat();
        let pose = chains_of(&[("P", 1, &spaced(&peptide))]);

        for (inside, pairs) in [("GKDALS", true), ("GKDWLS", false)] {
            let reference = chains_of(&[("A", 1, &spaced(&with(inside)))]);
            let map = residue_map(
                &polymer_chains(&pose),
                &polymer_chains(&reference),
                &[],
                |_| 0,
            );

            assert_eq!(map.chains().count() == 1, pairs, "{inside}");
        }
    }

    #[test]
    fn unrelated_chains_do_not_pair_whatever_their_lengths() {
        for (seed, lengths) in [(5, (100, 1000)), (6, (250, 250)), (7, (30, 600))] {
            let pose = chains_of(&[("A", 1, &spaced(&random_letters(seed, lengths.0)))]);
            let reference = chains_of(&[("A", 1, &spaced(&random_letters(seed + 10, lengths.1)))]);

            let map = residue_map(
                &polymer_chains(&pose),
                &polymer_chains(&reference),
                &[],
                |_| 0,
            );

            assert_eq!(map.chains().count(), 0, "{lengths:?}");
            assert!(!map.is_alignable(&FixedStr::new("A")));
        }
    }

    #[test]
    fn a_paralog_that_keeps_more_does_not_displace_the_true_partner() {
        let own = random_letters(8, 40);
        // Two residues in five substituted: aligned, but plainly another protein.
        let paralog: Vec<char> = own
            .iter()
            .enumerate()
            .map(|(i, &letter)| if i % 5 < 2 { 'W' } else { letter })
            .collect();
        let reference = chains_of(&[("A", 1, &spaced(&own)), ("B", 1, &spaced(&paralog))]);
        let pose = chains_of(&[("X", 1, &spaced(&own))]);
        let kept_by_paralog = |pair: &ChainPair| usize::from(pair.reference == 1);

        let map = residue_map(
            &polymer_chains(&pose),
            &polymer_chains(&reference),
            &[],
            kept_by_paralog,
        );

        let chains: Vec<(ChainId, ChainId)> = map.chains().copied().collect();
        assert_eq!(chains, [(FixedStr::new("X"), FixedStr::new("A"))]);
    }

    #[test]
    fn a_preferred_reference_chain_wins_when_nothing_else_decides() {
        // A monomer pose against a dimer whose ligand sits on chain B: the pose
        // chain stands for B, though it shares A's name.
        let reference = chains_of(&[("A", 1, MONOMER), ("B", 1, MONOMER)]);
        let pose = chains_of(&[("A", 1, MONOMER)]);
        let preferred = [FixedStr::new("B")];

        let map = residue_map(
            &polymer_chains(&pose),
            &polymer_chains(&reference),
            &preferred,
            |_| 0,
        );

        assert_eq!(
            map.reference_residue(&key("A", 1, None)),
            Some((key("B", 1, None), FixedStr::new("ALA")))
        );
        // Either reference chain could have paired; only one did.
        assert!(map.is_alignable(&FixedStr::new("A")) && map.is_alignable(&FixedStr::new("B")));
        assert_eq!(map.chains().count(), 1);
    }
}
