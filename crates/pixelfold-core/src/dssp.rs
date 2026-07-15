//! DSSP secondary-structure assignment.
//!
//! Reference: Kabsch, W. & Sander, C. (1983). Dictionary of protein secondary
//! structure. Biopolymers 22(12), 2577-2637.
//!
//! The backbone hydrogen bond drives everything downstream, so its geometry must
//! be right. The amide hydrogen is almost never present in a crystal structure;
//! it is inferred 1.0 A from the amide nitrogen, opposite the C=O bond of the
//! preceding residue. Proline has no amide hydrogen, and the first residue of a
//! chain (or one after a chain break) has no preceding carbonyl, so neither can
//! donate.

use std::collections::{HashMap, HashSet};

use glam::Vec3;

use crate::fixed_str::FixedStr;
use crate::spatial::Grid;
use crate::structure::{Atom, SecondaryStructure};

/// Below this energy (kcal/mol) a backbone N-H...O=C pair counts as bonded.
const HBOND_ENERGY_THRESHOLD: f32 = -0.5;
/// DSSP clamps the energy of near-overlapping atoms to this floor.
const MIN_ENERGY: f32 = -9.9;
/// Interatomic distances below this (A) are treated as coincident.
const MIN_DISTANCE: f32 = 0.1;
/// A C(prev)-N separation beyond this (A) is a chain break, not a peptide bond.
const MAX_PEPTIDE_BOND: f32 = 2.5;
/// Donor N atoms beyond this (A) from an acceptor O cannot reach the -0.5
/// kcal/mol threshold, so the neighbour search need not look further.
const HBOND_SEARCH_RADIUS: f32 = 6.0;

/// Backbone hydrogen bond: the C=O of `acceptor_residue` to the N-H of
/// `donor_residue`. Residue fields are indices into the ordered residue list.
#[derive(Clone, Copy, Debug)]
pub struct HBond {
    pub donor_residue: usize,
    pub acceptor_residue: usize,
    pub donor_atom_idx: usize,
    pub acceptor_atom_idx: usize,
    pub energy: f32,
}

/// One residue's backbone geometry, in chain then sequence order.
struct Residue {
    chain_id: FixedStr<4>,
    residue_seq: u32,
    insertion_code: Option<char>,
    is_proline: bool,
    n: Option<Vec3>,
    c: Option<Vec3>,
    o: Option<Vec3>,
    /// Inferred amide hydrogen; `None` when the residue cannot donate.
    h: Option<Vec3>,
    n_atom_idx: Option<usize>,
    o_atom_idx: Option<usize>,
}

/// Assign secondary structure to `atoms` in place and return the backbone
/// hydrogen bonds (for visualisation and the residue interaction network).
pub fn assign(atoms: &mut [Atom]) -> Vec<HBond> {
    let mut residues = build_residues(atoms);
    infer_amide_hydrogens(&mut residues);
    let hbonds = compute_hbonds(&residues);
    let ss = assign_from_hbonds(&hbonds, &residues);

    let mut ss_by_key: HashMap<(FixedStr<4>, u32, Option<char>), SecondaryStructure> =
        HashMap::new();
    for (res, &s) in residues.iter().zip(ss.iter()) {
        ss_by_key.insert((res.chain_id, res.residue_seq, res.insertion_code), s);
    }
    for atom in atoms.iter_mut() {
        if let Some(&s) = ss_by_key.get(&(atom.chain_id, atom.residue_seq, atom.insertion_code)) {
            atom.secondary_structure = s;
        }
    }

    hbonds
}

/// Group atoms into residues, preserving file order (chain then sequence) and
/// capturing backbone atom positions. A new residue starts whenever the
/// (chain, sequence, insertion code) key changes.
fn build_residues(atoms: &[Atom]) -> Vec<Residue> {
    let mut residues: Vec<Residue> = Vec::new();
    let mut current: Option<(FixedStr<4>, u32, Option<char>)> = None;

    for (idx, atom) in atoms.iter().enumerate() {
        let key = (atom.chain_id, atom.residue_seq, atom.insertion_code);
        if current != Some(key) {
            current = Some(key);
            residues.push(Residue {
                chain_id: atom.chain_id,
                residue_seq: atom.residue_seq,
                insertion_code: atom.insertion_code,
                is_proline: atom.residue_name == "PRO",
                n: None,
                c: None,
                o: None,
                h: None,
                n_atom_idx: None,
                o_atom_idx: None,
            });
        }
        if let Some(res) = residues.last_mut() {
            match atom.name.as_str() {
                "N" => {
                    res.n = Some(atom.position);
                    res.n_atom_idx = Some(idx);
                }
                "C" => res.c = Some(atom.position),
                "O" => {
                    res.o = Some(atom.position);
                    res.o_atom_idx = Some(idx);
                }
                _ => {}
            }
        }
    }

    residues
}

/// Place the amide hydrogen for every residue that can donate one:
///
/// `H = N + 1.0 A * normalize(C_prev - O_prev)`.
fn infer_amide_hydrogens(residues: &mut [Residue]) {
    for i in 0..residues.len() {
        if i == 0 || residues[i].is_proline {
            continue; // proline never donates; residue 0 has no predecessor
        }

        let Some(n) = residues[i].n else { continue };
        if residues[i - 1].chain_id != residues[i].chain_id {
            continue; // start of a new chain
        }

        let (Some(c_prev), Some(o_prev)) = (residues[i - 1].c, residues[i - 1].o) else {
            continue;
        };

        if (c_prev - n).length() > MAX_PEPTIDE_BOND {
            continue; // chain break: no genuine peptide bond to the predecessor
        }

        residues[i].h = Some(n + (c_prev - o_prev).normalize());
    }
}

/// The Kabsch-Sander electrostatic hydrogen bond energy (kcal/mol):
///
/// `E = 0.084 * 332 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)`.
fn hbond_energy(c: Vec3, o: Vec3, n: Vec3, h: Vec3) -> f32 {
    let r_on = (o - n).length();
    let r_ch = (c - h).length();
    let r_oh = (o - h).length();
    let r_cn = (c - n).length();

    if r_on < MIN_DISTANCE || r_ch < MIN_DISTANCE || r_oh < MIN_DISTANCE || r_cn < MIN_DISTANCE {
        return 0.0;
    }

    let energy = 0.084 * 332.0 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn);
    energy.max(MIN_ENERGY)
}

/// Every backbone hydrogen bond below the energy threshold. The acceptor
/// contributes C=O, the donor contributes N-H; sequential same-chain neighbours
/// are peptide-bonded rather than hydrogen-bonded and are skipped.
///
/// A uniform grid over the donor nitrogens turns the acceptor-donor search from
/// O(N^2) into O(N): each acceptor only tests donors within the search radius.
fn compute_hbonds(residues: &[Residue]) -> Vec<HBond> {
    let mut donor_positions: Vec<Vec3> = Vec::new();
    let mut donor_residue: Vec<usize> = Vec::new();
    for (j, res) in residues.iter().enumerate() {
        if let (Some(n), Some(_)) = (res.n, res.h) {
            donor_positions.push(n);
            donor_residue.push(j);
        }
    }
    let grid = Grid::build(&donor_positions, HBOND_SEARCH_RADIUS);

    let mut hbonds = Vec::new();
    for (i, acc) in residues.iter().enumerate() {
        let (Some(c), Some(o), Some(o_idx)) = (acc.c, acc.o, acc.o_atom_idx) else {
            continue;
        };
        let acc_chain = acc.chain_id;

        grid.for_each_within(o, HBOND_SEARCH_RADIUS, |d| {
            let j = donor_residue[d];
            if i == j || (acc_chain == residues[j].chain_id && i.abs_diff(j) == 1) {
                return;
            }

            let (Some(n), Some(h), Some(n_idx)) =
                (residues[j].n, residues[j].h, residues[j].n_atom_idx)
            else {
                return;
            };

            let energy = hbond_energy(c, o, n, h);
            if energy < HBOND_ENERGY_THRESHOLD {
                hbonds.push(HBond {
                    donor_residue: j,
                    acceptor_residue: i,
                    donor_atom_idx: n_idx,
                    acceptor_atom_idx: o_idx,
                    energy,
                });
            }
        });
    }

    hbonds
}

/// Assign secondary structure from the hydrogen bond set using the DSSP pattern
/// rules, reduced to helix / sheet / turn / coil.
///
/// `bonded(a, b)` reads "the C=O of residue a hydrogen-bonds the N-H of residue
/// b". An n-turn is `bonded(i, i + n)`; a helix needs two consecutive n-turns; a
/// bridge is the parallel or antiparallel pairing of two strands. Priority
/// follows DSSP: alpha helix, then sheet, then 3-10 / pi helix, then turn.
fn assign_from_hbonds(hbonds: &[HBond], residues: &[Residue]) -> Vec<SecondaryStructure> {
    let n = residues.len();
    let mut ss = vec![SecondaryStructure::Coil; n];
    if n == 0 {
        return ss;
    }

    let bonds: HashSet<(usize, usize)> = hbonds
        .iter()
        .map(|hb| (hb.acceptor_residue, hb.donor_residue))
        .collect();
    let bonded = |a: usize, b: usize| bonds.contains(&(a, b));
    let same_chain = |a: usize, b: usize| residues[a].chain_id == residues[b].chain_id;

    // n-turn: C=O(i) ... N-H(i + step), within one chain.
    let turn = |i: usize, step: usize| {
        let j = i + step;
        j < n && same_chain(i, j) && bonded(i, j)
    };

    // Alpha helix (highest priority): two consecutive 4-turns cover i..=i+3.
    for i in 1..n {
        if turn(i - 1, 4) && turn(i, 4) {
            let end = (i + 3).min(n - 1);
            ss[i..=end].fill(SecondaryStructure::Helix);
        }
    }

    // Beta bridges -> sheet, over coil only (alpha wins ties).
    let bridge = |i: usize, j: usize| -> bool {
        let prev = |x: usize| (x > 0 && same_chain(x - 1, x)).then(|| x - 1);
        let next = |x: usize| (x + 1 < n && same_chain(x, x + 1)).then_some(x + 1);
        let (pi, ni, pj, nj) = (prev(i), next(i), prev(j), next(j));

        let parallel = matches!((pi, ni), (Some(p), Some(q)) if bonded(p, j) && bonded(j, q))
            || matches!((pj, nj), (Some(p), Some(q)) if bonded(p, i) && bonded(i, q));
        let antiparallel = (bonded(i, j) && bonded(j, i))
            || matches!((pi, nj, pj, ni), (Some(a), Some(b), Some(c), Some(d))
                if bonded(a, b) && bonded(c, d));

        parallel || antiparallel
    };
    for i in 0..n {
        for j in (i + 1)..n {
            if same_chain(i, j) && j - i < 3 {
                continue; // too close in one chain to be a genuine bridge
            }

            if bridge(i, j) {
                if ss[i] == SecondaryStructure::Coil {
                    ss[i] = SecondaryStructure::Sheet;
                }
                if ss[j] == SecondaryStructure::Coil {
                    ss[j] = SecondaryStructure::Sheet;
                }
            }
        }
    }

    // 3-10 and pi helices, over coil only (sheet wins ties).
    for (step, span) in [(3usize, 2usize), (5, 4)] {
        for i in 1..n {
            if turn(i - 1, step) && turn(i, step) {
                let end = (i + span).min(n - 1);
                for s in ss[i..=end].iter_mut() {
                    if *s == SecondaryStructure::Coil {
                        *s = SecondaryStructure::Helix;
                    }
                }
            }
        }
    }

    // Turns: residues inside an n-turn that are still coil.
    for i in 0..n {
        for step in [3usize, 4, 5] {
            if turn(i, step) {
                let end = (i + step - 1).min(n - 1);
                for s in ss[(i + 1)..=end].iter_mut() {
                    if *s == SecondaryStructure::Coil {
                        *s = SecondaryStructure::Turn;
                    }
                }
            }
        }
    }

    ss
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residue(chain: &str, seq: u32, proline: bool) -> Residue {
        Residue {
            chain_id: FixedStr::new(chain),
            residue_seq: seq,
            insertion_code: None,
            is_proline: proline,
            n: None,
            c: None,
            o: None,
            h: None,
            n_atom_idx: None,
            o_atom_idx: None,
        }
    }

    /// `n` residues, all in chain A, for pattern tests driven by synthetic bonds.
    fn chain_a(n: usize) -> Vec<Residue> {
        (0..n).map(|i| residue("A", i as u32 + 1, false)).collect()
    }

    /// A hydrogen bond whose C=O side is `acceptor` and N-H side is `donor`.
    fn hbond(acceptor: usize, donor: usize) -> HBond {
        HBond {
            donor_residue: donor,
            acceptor_residue: acceptor,
            donor_atom_idx: 0,
            acceptor_atom_idx: 0,
            energy: -1.0,
        }
    }

    #[test]
    fn consecutive_four_turns_form_an_alpha_helix() {
        let residues = chain_a(8);
        let hbonds = [hbond(0, 4), hbond(1, 5), hbond(2, 6), hbond(3, 7)];
        let ss = assign_from_hbonds(&hbonds, &residues);
        // The overlap of the 4-turns covers residues 1..=6.
        for (offset, s) in ss[1..=6].iter().enumerate() {
            assert_eq!(*s, SecondaryStructure::Helix, "residue {}", offset + 1);
        }
        assert_ne!(ss[0], SecondaryStructure::Helix);
    }

    #[test]
    fn a_single_four_turn_is_a_turn_not_a_helix() {
        let residues = chain_a(6);
        let ss = assign_from_hbonds(&[hbond(0, 4)], &residues);
        assert!(ss.iter().all(|&s| s != SecondaryStructure::Helix));
        assert_eq!(ss[2], SecondaryStructure::Turn); // inside the lone turn
    }

    #[test]
    fn antiparallel_bridge_forms_a_sheet() {
        let residues = chain_a(9);
        // Two strands with mutual C=O/N-H pairs at (0,8) and (2,6).
        let hbonds = [hbond(2, 6), hbond(6, 2), hbond(0, 8), hbond(8, 0)];
        let ss = assign_from_hbonds(&hbonds, &residues);
        for i in [0usize, 2, 6, 8] {
            assert_eq!(ss[i], SecondaryStructure::Sheet, "residue {i}");
        }
    }

    #[test]
    fn parallel_bridge_forms_a_sheet() {
        let residues = chain_a(9);
        // Parallel: C=O(i-1)->N-H(j) and C=O(j)->N-H(i+1) for i=3, j=7.
        let ss = assign_from_hbonds(&[hbond(2, 7), hbond(7, 4)], &residues);
        assert_eq!(ss[3], SecondaryStructure::Sheet);
        assert_eq!(ss[7], SecondaryStructure::Sheet);
    }

    #[test]
    fn amide_hydrogen_is_placed_opposite_the_preceding_carbonyl() {
        let mut prev = residue("A", 1, false);
        prev.c = Some(Vec3::new(1.0, 0.0, 0.0));
        prev.o = Some(Vec3::new(2.0, 0.0, 0.0));
        let mut cur = residue("A", 2, false);
        cur.n = Some(Vec3::new(0.0, 0.0, 0.0));

        let mut residues = vec![prev, cur];
        infer_amide_hydrogens(&mut residues);

        // H = N + normalize(C_prev - O_prev) = (0,0,0) + normalize((-1,0,0)).
        let h = residues[1]
            .h
            .expect("residue with a valid predecessor donates");
        assert!((h - Vec3::new(-1.0, 0.0, 0.0)).length() < 1e-5);
    }

    #[test]
    fn proline_and_chain_start_do_not_donate() {
        let mut prev = residue("A", 1, false);
        prev.c = Some(Vec3::new(1.0, 0.0, 0.0));
        prev.o = Some(Vec3::new(2.0, 0.0, 0.0));
        let mut proline = residue("A", 2, true); // proline
        proline.n = Some(Vec3::new(0.0, 0.0, 0.0));
        let mut chain_start = residue("B", 1, false); // first of a new chain
        chain_start.n = Some(Vec3::new(0.0, 0.0, 0.0));

        let mut residues = vec![prev, proline, chain_start];
        infer_amide_hydrogens(&mut residues);

        assert!(residues[0].h.is_none()); // residue 0: no predecessor
        assert!(residues[1].h.is_none()); // proline
        assert!(residues[2].h.is_none()); // different chain than its predecessor
    }

    #[test]
    fn chain_break_prevents_donation() {
        let mut prev = residue("A", 1, false);
        prev.c = Some(Vec3::new(0.0, 0.0, 0.0));
        prev.o = Some(Vec3::new(1.0, 0.0, 0.0));
        let mut cur = residue("A", 9, false); // same chain but far away (missing loop)
        cur.n = Some(Vec3::new(30.0, 0.0, 0.0));

        let mut residues = vec![prev, cur];
        infer_amide_hydrogens(&mut residues);
        assert!(residues[1].h.is_none());
    }

    #[test]
    fn energy_is_favourable_for_a_linear_hydrogen_bond() {
        // C=O ... H-N roughly collinear, O...H ~1.8 A.
        let c = Vec3::new(-1.23, 0.0, 0.0);
        let o = Vec3::new(0.0, 0.0, 0.0);
        let n = Vec3::new(2.8, 0.0, 0.0);
        let h = Vec3::new(1.8, 0.0, 0.0);
        assert!(hbond_energy(c, o, n, h) < HBOND_ENERGY_THRESHOLD);

        // Far apart: no bond.
        let far_n = Vec3::new(100.0, 0.0, 0.0);
        let far_h = Vec3::new(99.0, 0.0, 0.0);
        assert!(hbond_energy(c, o, far_n, far_h) > HBOND_ENERGY_THRESHOLD);
    }

    #[test]
    fn compute_hbonds_finds_the_favourable_pair_via_the_grid() {
        // Acceptor at residue 0, donor at residue 3 (non-adjacent), positioned
        // for a favourable bond; residues 1 and 2 carry no backbone atoms.
        let mut acceptor = residue("A", 1, false);
        acceptor.c = Some(Vec3::new(-1.23, 0.0, 0.0));
        acceptor.o = Some(Vec3::new(0.0, 0.0, 0.0));
        acceptor.o_atom_idx = Some(10);
        let mut donor = residue("A", 4, false);
        donor.n = Some(Vec3::new(2.8, 0.0, 0.0));
        donor.h = Some(Vec3::new(1.8, 0.0, 0.0));
        donor.n_atom_idx = Some(20);

        let residues = vec![
            acceptor,
            residue("A", 2, false),
            residue("A", 3, false),
            donor,
        ];
        let hbonds = compute_hbonds(&residues);

        assert_eq!(hbonds.len(), 1);
        assert_eq!(hbonds[0].acceptor_residue, 0);
        assert_eq!(hbonds[0].donor_residue, 3);
        assert_eq!(hbonds[0].acceptor_atom_idx, 10);
        assert_eq!(hbonds[0].donor_atom_idx, 20);
    }

    #[test]
    fn build_residues_keeps_chains_separate() {
        // Same residue_seq (1) in two chains must not merge.
        let atoms = vec![
            atom("A", 1, "N"),
            atom("A", 1, "CA"),
            atom("B", 1, "N"),
            atom("B", 1, "CA"),
        ];
        let residues = build_residues(&atoms);
        assert_eq!(residues.len(), 2);
        assert_eq!(residues[0].chain_id.as_str(), "A");
        assert_eq!(residues[1].chain_id.as_str(), "B");
    }

    fn atom(chain: &str, seq: u32, name: &str) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new("C"),
            residue_name: FixedStr::new("ALA"),
            residue_seq: seq,
            chain_id: FixedStr::new(chain),
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
}
