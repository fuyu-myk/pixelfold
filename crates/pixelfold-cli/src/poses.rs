//! `pixelfold poses`: ligand poses ranked against a reference by the contacts
//! the ligand makes, with the clashes each suffers.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::PathBuf;

use anyhow::Result;
use pixelfold_core::components::Dictionary;
use pixelfold_core::interactions::{Interaction, detect};
use pixelfold_core::numbering::{self, Chain, ResidueKey, ResidueMap};
use pixelfold_core::pose::{self, Attachments, Contact, Contents, Fingerprint};
use pixelfold_core::rmsd::{AtomId, NoRmsd, POCKET_REACH, RmsdReference, Superposition};
use pixelfold_core::{Atom, AtomSet, ChainId, Protein};
use serde::Serialize;

use crate::analyse::rounded;
use crate::load;
use crate::report::Row;

/// One pose against the reference, flattened to a row.
#[derive(Serialize)]
pub struct PoseRecord {
    pub pose: String,
    pub model: usize,
    /// The copy of the ligand scored, named by its first residue, when the
    /// pose holds several and the one most like the reference's is kept.
    pub copy: Option<String>,
    pub contacts: usize,
    pub shared: usize,
    pub recovery: f32,
    pub similarity: f32,
    /// The ligand's heavy-atom RMSD in angstroms once the pose is superposed
    /// on the reference's pocket, or `None` when it cannot be measured.
    pub rmsd: Option<f32>,
    pub clashes: usize,
    /// The reference ligand's stereocentres the pose inverts, or `None` when
    /// its atoms could not be paired with the reference's.
    pub inverted: Option<Vec<String>>,
    pub lost: Vec<String>,
    /// The scored copy's atoms, which the viewer tints.
    #[serde(skip)]
    pub ligand: Vec<usize>,
    /// What lays the pose on the reference's pocket, for the viewer.
    #[serde(skip)]
    pub superposition: Option<Superposition>,
}

impl Row for PoseRecord {
    fn columns() -> Vec<&'static str> {
        vec![
            "pose",
            "model",
            "copy",
            "contacts",
            "shared",
            "recovery",
            "similarity",
            "rmsd",
            "clashes",
            "inverted",
            "lost",
        ]
    }

    fn cells(&self) -> Vec<String> {
        vec![
            self.pose.clone(),
            self.model.to_string(),
            self.copy.clone().unwrap_or_default(),
            self.contacts.to_string(),
            self.shared.to_string(),
            format!("{:.3}", self.recovery),
            format!("{:.3}", self.similarity),
            self.rmsd.map_or(String::new(), |rmsd| format!("{rmsd:.2}")),
            self.clashes.to_string(),
            self.inverted
                .as_ref()
                .map_or(String::new(), |inverted| inverted.join(", ")),
            self.lost.join(", "),
        ]
    }
}

/// Warnings for standard error, each printed once however many poses raise it.
#[derive(Default)]
pub struct Warnings {
    seen: HashSet<String>,
}

impl Warnings {
    fn warn(&mut self, message: String) {
        if !self.seen.contains(&message) {
            eprintln!("warning: {message}");
            self.seen.insert(message);
        }
    }
}

/// What every pose is compared with: the reference ligand's contacts, its
/// covalent attachments, the chains a pose's are aligned with, what its RMSD
/// is measured against, what one copy of it holds, and the chemistry a pose
/// may borrow.
pub struct Reference {
    components: Dictionary,
    fingerprint: Fingerprint,
    attachments: Attachments,
    chains: Vec<Chain>,
    rmsd: RmsdReference,
    /// What the reference ligand holds, when it is one covalent piece: a pose
    /// ligand in several pieces that each hold as much holds copies of it.
    copy: Option<Contents>,
}

/// One copy of a pose's ligand scored, with what the warnings about it need.
struct Scored {
    record: PoseRecord,
    map: ResidueMap,
    no_rmsd: Option<NoRmsd>,
}

/// A model of a pose file, with what scoring any copy of its ligand needs.
struct Model<'a> {
    pose: &'a str,
    number: usize,
    protein: &'a Protein,
    interactions: Vec<Interaction>,
    chains: Vec<Chain>,
}

impl Reference {
    /// The reference made of the `query` atoms of `protein`.
    pub fn new(protein: &Protein, query: &str, warnings: &mut Warnings) -> Result<Self> {
        let ligand = load::chosen(protein, Some(query))?;
        let residues = ligand_residues(protein, &ligand, |_| true);
        if residues.len() > 1 {
            let names: Vec<&str> = residues.keys().copied().collect();
            warnings.warn(format!(
                "the ligand `{query}` covers {} components ({}); pass --ligand CODE if only \
                 one is meant",
                names.len(),
                names.join(", ")
            ));
        }
        // A peptide repeats its residue names by nature; only a small molecule
        // repeated is a copy.
        let small = ligand_residues(protein, &ligand, |atom| !atom.is_polymer());
        for (name, copies) in small.iter().filter(|(_, copies)| copies.len() > 1) {
            let ids: Vec<String> = copies.iter().map(|&residue| residue_id(residue)).collect();
            warnings.warn(format!(
                "the ligand `{query}` covers {} residues of {name} ({}); if they are copies of \
                 one ligand, pick one with --reference-select",
                copies.len(),
                ids.join(", ")
            ));
        }
        warn_undefined(protein, &ligand, warnings);
        let rmsd = RmsdReference::new(protein, &ligand);
        if !rmsd.has_pocket() {
            warnings.warn(format!(
                "the ligand `{query}` has fewer than three protein backbone atoms within \
                 {POCKET_REACH} A, so no pose gets an RMSD"
            ));
        }

        let fingerprint = pose::fingerprint(protein, &detect(protein), &ligand);
        if fingerprint.is_empty() {
            anyhow::bail!("the ligand `{query}` makes no contacts to compare poses with");
        }

        Ok(Self {
            components: protein.components.clone(),
            fingerprint,
            attachments: pose::attachments(protein, &ligand),
            chains: numbering::polymer_chains(protein),
            rmsd,
            copy: (pose::pieces(protein, &ligand).len() == 1)
                .then(|| Contents::of(protein, &ligand)),
        })
    }

    /// Score model `model` of the pose file `pose` against the reference.
    ///
    /// A ligand the pose file does not define takes the reference's definition
    /// of the same component, matched by atom name. The pose's residues are
    /// matched with the reference's by aligning chain sequences. A pose holding
    /// several copies of the ligand is scored by the copy most like the
    /// reference's.
    pub fn assess(
        &self,
        protein: Protein,
        query: &str,
        pose: &str,
        model: usize,
        warnings: &mut Warnings,
    ) -> Result<PoseRecord> {
        let protein = Protein {
            components: protein.components.with_fallback(&self.components),
            ..protein
        };
        let ligand = load::chosen(&protein, Some(query))?;
        warn_undefined(&protein, &ligand, warnings);

        let copies = self.copies(&protein, &ligand, query, pose, warnings);
        let model = Model {
            pose,
            number: model,
            protein: &protein,
            interactions: detect(&protein),
            chains: numbering::polymer_chains(&protein),
        };
        let mut best = copies
            .iter()
            .map(|copy| self.score(&model, copy, copies.len() > 1))
            .min_by(|a, b| preferred(&a.record, &b.record))
            .expect("copies holds at least the whole ligand");
        best.record.superposition = self.rmsd.superposition(&protein, &best.map).ok();

        self.check_pairing(&best.map, pose, warnings)?;
        if let Some(reason) = &best.no_rmsd
            && self.rmsd.has_pocket()
        {
            warnings.warn(format!("no RMSD for {pose}: {reason}"));
        }

        Ok(best.record)
    }

    /// The copies of the reference ligand among the `ligand` atoms of the pose
    /// file `pose`: its covalent pieces when there are several and each may be
    /// a copy, else the whole, with a warning when the whole is in pieces.
    fn copies(
        &self,
        protein: &Protein,
        ligand: &AtomSet,
        query: &str,
        pose: &str,
        warnings: &mut Warnings,
    ) -> Vec<AtomSet> {
        let Some(copy) = &self.copy else {
            return vec![ligand.clone()];
        };
        let pieces = pose::pieces(protein, ligand);
        if pieces.len() > 1 {
            if pieces
                .iter()
                .all(|piece| copy.admits(&Contents::of(protein, piece)))
            {
                return pieces;
            }
            warnings.warn(format!(
                "the ligand `{query}` is {} separate pieces in {pose}, not all copies of the \
                 reference's, so they are scored as one; narrow --select to score a copy",
                pieces.len()
            ));
        }

        vec![ligand.clone()]
    }

    /// One `copy` of the ligand of `model` scored, the copy named in the
    /// record when `named`.
    fn score(&self, model: &Model, copy: &AtomSet, named: bool) -> Scored {
        let reference = &self.fingerprint;
        let found = pose::fingerprint(model.protein, &model.interactions, copy);
        let map = pose::residue_map(&found, &model.chains, reference, &self.chains);
        let fingerprint = pose::renumbered(&found, reference, &map);
        let fit = self.rmsd.fit(model.protein, copy, &map);

        let record = PoseRecord {
            pose: model.pose.to_string(),
            model: model.number,
            copy: named.then(|| copy_label(model.protein, copy)),
            contacts: fingerprint.len(),
            shared: fingerprint.shared(reference),
            recovery: rounded(fingerprint.recovery(reference), 3),
            similarity: rounded(fingerprint.similarity(reference), 3),
            rmsd: fit.as_ref().ok().map(|fit| rounded(fit.rmsd, 2)),
            clashes: pose::clashes(model.protein, copy, &self.attachments.renumbered(&map)),
            inverted: fit
                .as_ref()
                .ok()
                .map(|fit| fit.inverted.iter().map(atom_label).collect()),
            lost: fingerprint.lost(reference).map(contact_label).collect(),
            ligand: copy.iter().collect(),
            superposition: None,
        };

        Scored {
            record,
            map,
            no_rmsd: fit.err(),
        }
    }

    /// Fail when no chain of the pose aligns with any reference chain the
    /// ligand touches, and warn of each such chain left without a counterpart.
    fn check_pairing(&self, map: &ResidueMap, pose: &str, warnings: &mut Warnings) -> Result<()> {
        let touched: BTreeSet<ChainId> = self
            .fingerprint
            .contacts()
            .filter(|contact| map.is_reference_polymer(&contact.residue))
            .map(|contact| contact.residue.0)
            .collect();
        if !touched.is_empty() && !touched.iter().any(|chain| map.is_alignable(chain)) {
            let names: Vec<String> = touched.iter().map(|chain| chain.to_string()).collect();
            let chains = if names.len() == 1 { "chain" } else { "chains" };
            anyhow::bail!(
                "no chain here aligns by sequence with the reference's {chains} {}, which the \
                 ligand touches",
                names.join(", ")
            );
        }

        let paired: BTreeSet<ChainId> = map.chains().map(|&(_, reference)| reference).collect();
        for chain in touched.difference(&paired) {
            warnings.warn(format!(
                "chain {chain} of the reference has no counterpart in {pose}, so the contacts \
                 the ligand makes with it are lost"
            ));
        }

        Ok(())
    }
}

/// The residues of the ligand atoms that `include` admits, by component.
fn ligand_residues<'a>(
    protein: &'a Protein,
    ligand: &AtomSet,
    include: impl Fn(&Atom) -> bool,
) -> BTreeMap<&'a str, BTreeSet<ResidueKey>> {
    let mut residues: BTreeMap<&str, BTreeSet<ResidueKey>> = BTreeMap::new();
    for atom in ligand
        .iter()
        .map(|atom| &protein.atoms[atom])
        .filter(|atom| include(atom))
    {
        residues
            .entry(atom.residue_name.as_str())
            .or_default()
            .insert((atom.chain_id, atom.residue_seq, atom.insertion_code));
    }

    residues
}

/// A copy of the ligand named by its first residue.
fn copy_label(protein: &Protein, copy: &AtomSet) -> String {
    copy.iter()
        .next()
        .map(|atom| {
            let atom = &protein.atoms[atom];
            residue_id((atom.chain_id, atom.residue_seq, atom.insertion_code))
        })
        .unwrap_or_default()
}

/// The records best first: most similar, then most recovered, then fewest
/// clashes, ties kept in the order given.
pub fn ranked(mut records: Vec<PoseRecord>) -> Vec<PoseRecord> {
    records.sort_by(|a, b| {
        b.similarity
            .total_cmp(&a.similarity)
            .then(b.recovery.total_cmp(&a.recovery))
            .then(a.clashes.cmp(&b.clashes))
    });

    records
}

/// The ranked `records` as the viewer lists them, each pose file found by
/// `path_of`.
pub fn viewed(
    records: Vec<PoseRecord>,
    path_of: impl Fn(&str) -> Result<PathBuf>,
) -> Result<Vec<pixelfold_tui::Pose>> {
    records
        .into_iter()
        .map(|record| {
            Ok(pixelfold_tui::Pose {
                path: path_of(&record.pose)?,
                label: record.pose,
                model: record.model,
                copy: record.copy,
                similarity: record.similarity,
                recovery: record.recovery,
                rmsd: record.rmsd,
                clashes: record.clashes,
                inverted: record.inverted,
                lost: record.lost,
                ligand: record.ligand,
                superposition: record.superposition,
            })
        })
        .collect()
}

/// Which of two copies of a pose's ligand to keep: the more similar, then the
/// more recovered, then the nearer the reference ligand, then the one with
/// fewer clashes, then the first. Nearness comes before clashes because the
/// copy kept stands for the reference's, whose own crystal contacts can clash.
fn preferred(a: &PoseRecord, b: &PoseRecord) -> Ordering {
    let rmsd = |record: &PoseRecord| record.rmsd.unwrap_or(f32::INFINITY);
    b.similarity
        .total_cmp(&a.similarity)
        .then(b.recovery.total_cmp(&a.recovery))
        .then(rmsd(a).total_cmp(&rmsd(b)))
        .then(a.clashes.cmp(&b.clashes))
}

fn warn_undefined(protein: &Protein, ligand: &AtomSet, warnings: &mut Warnings) {
    for (residue, name) in pose::undefined_ligand_residues(protein, ligand) {
        warnings.warn(format!(
            "ligand {name} {} has atoms without a chemical definition, so its contacts are \
             incomplete",
            residue_id(residue)
        ));
    }
}

/// A residue as `chain/seq`, with the insertion code appended when it has one.
fn residue_id((chain, seq, icode): ResidueKey) -> String {
    match icode {
        Some(code) => format!("{chain}/{seq}{code}"),
        None => format!("{chain}/{seq}"),
    }
}

/// An atom as its residue and name, e.g. `0G6 H/1 CA`.
fn atom_label((residue, residue_name, name): &AtomId) -> String {
    format!("{residue_name} {} {name}", residue_id(*residue))
}

/// A contact as its residue and kind, e.g. `HIS H/57 hydrogen-bond`.
fn contact_label(contact: &Contact) -> String {
    format!(
        "{} {} {}",
        contact.residue_name,
        residue_id(contact.residue),
        contact.kind.label()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(pose: &str, similarity: f32, recovery: f32, clashes: usize) -> PoseRecord {
        PoseRecord {
            pose: pose.to_string(),
            model: 1,
            copy: None,
            contacts: 0,
            shared: 0,
            recovery,
            similarity,
            rmsd: None,
            clashes,
            inverted: None,
            lost: Vec::new(),
            ligand: Vec::new(),
            superposition: None,
        }
    }

    #[test]
    fn of_copies_alike_in_contacts_the_nearer_is_kept_before_the_unclashed() {
        let copy = |pose: &str, similarity: f32, rmsd: Option<f32>, clashes: usize| PoseRecord {
            rmsd,
            ..record(pose, similarity, 1.0, clashes)
        };
        let copies = [
            copy("unmeasured", 1.0, None, 0),
            copy("far", 1.0, Some(1.5), 0),
            copy("fewer-contacts", 0.9, Some(0.0), 0),
            copy("in-place", 1.0, Some(0.0), 1),
            copy("in-place-clashing-more", 1.0, Some(0.0), 2),
            copy("in-place-too", 1.0, Some(0.0), 1),
        ];

        let best = copies.into_iter().min_by(preferred).expect("copies");
        assert_eq!(best.pose, "in-place");
    }

    #[test]
    fn poses_rank_by_similarity_then_recovery_then_fewest_clashes() {
        let given = vec![
            record("tied-first", 0.5, 0.5, 2),
            record("fewer-clashes", 0.5, 0.5, 1),
            record("more-recovered", 0.5, 0.8, 9),
            record("most-similar", 0.9, 0.1, 9),
            record("tied-second", 0.5, 0.5, 2),
        ];

        let order: Vec<String> = ranked(given).into_iter().map(|r| r.pose).collect();
        assert_eq!(
            order,
            [
                "most-similar",
                "more-recovered",
                "fewer-clashes",
                "tied-first",
                "tied-second"
            ]
        );
    }
}
