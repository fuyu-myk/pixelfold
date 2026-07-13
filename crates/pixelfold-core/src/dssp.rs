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

use std::collections::HashMap;

use glam::Vec3;

use crate::fixed_str::FixedStr;
use crate::structure::{Atom, SecondaryStructure};

/// Below this energy (kcal/mol) a backbone N-H...O=C pair counts as bonded.
const HBOND_ENERGY_THRESHOLD: f32 = -0.5;
/// DSSP clamps the energy of near-overlapping atoms to this floor.
const MIN_ENERGY: f32 = -9.9;
/// Interatomic distances below this (A) are treated as coincident.
const MIN_DISTANCE: f32 = 0.1;
/// A C(prev)-N separation beyond this (A) is a chain break, not a peptide bond.
const MAX_PEPTIDE_BOND: f32 = 2.5;

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
fn compute_hbonds(residues: &[Residue]) -> Vec<HBond> {
    let mut hbonds = Vec::new();

    for (i, acc) in residues.iter().enumerate() {
        let (Some(c), Some(o), Some(o_idx)) = (acc.c, acc.o, acc.o_atom_idx) else {
            continue;
        };

        for (j, don) in residues.iter().enumerate() {
            if i == j {
                continue;
            }
            if acc.chain_id == don.chain_id && i.abs_diff(j) == 1 {
                continue;
            }
            let (Some(n), Some(h), Some(n_idx)) = (don.n, don.h, don.n_atom_idx) else {
                continue;
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
        }
    }

    hbonds
}

/// Interim secondary-structure assignment from the hydrogen bond set. This is a
/// chain-aware version of the original offset heuristic.
fn assign_from_hbonds(hbonds: &[HBond], residues: &[Residue]) -> Vec<SecondaryStructure> {
    let n = residues.len();
    let mut ss = vec![SecondaryStructure::Coil; n];

    let mut acceptor_to_donors: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut donor_to_acceptors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for hb in hbonds {
        acceptor_to_donors[hb.acceptor_residue].push(hb.donor_residue);
        donor_to_acceptors[hb.donor_residue].push(hb.acceptor_residue);
    }

    let same_chain = |a: usize, b: usize| residues[a].chain_id == residues[b].chain_id;

    // Helices: a residue donates to i+3/i+4/i+5 within the same chain.
    for i in 0..n {
        for &acc in &donor_to_acceptors[i] {
            if !same_chain(i, acc) {
                continue;
            }

            let offset = i.abs_diff(acc);
            if offset == 4 {
                let (lo, hi) = (i.min(acc), i.max(acc));
                ss[lo..=hi].fill(SecondaryStructure::Helix);
            } else if (offset == 3 || offset == 5) && ss[i] == SecondaryStructure::Coil {
                ss[i] = SecondaryStructure::Helix;
                ss[acc] = SecondaryStructure::Helix;
            }
        }
    }

    // Sheets: non-local hydrogen bonds (different chain, or far in sequence).
    for i in 0..n {
        if ss[i] == SecondaryStructure::Helix {
            continue;
        }

        for &partner in acceptor_to_donors[i]
            .iter()
            .chain(donor_to_acceptors[i].iter())
        {
            let nonlocal = !same_chain(i, partner) || i.abs_diff(partner) > 4;
            if nonlocal {
                ss[i] = SecondaryStructure::Sheet;
                if ss[partner] != SecondaryStructure::Helix {
                    ss[partner] = SecondaryStructure::Sheet;
                }
            }
        }
    }

    // Turns: coil residues flanked by structured residues in the same chain.
    for i in 1..n.saturating_sub(1) {
        if ss[i] != SecondaryStructure::Coil {
            continue;
        }

        let prev = if same_chain(i, i - 1) {
            ss[i - 1]
        } else {
            SecondaryStructure::Coil
        };
        let next = if same_chain(i, i + 1) {
            ss[i + 1]
        } else {
            SecondaryStructure::Coil
        };
        if prev != SecondaryStructure::Coil || next != SecondaryStructure::Coil {
            ss[i] = SecondaryStructure::Turn;
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
