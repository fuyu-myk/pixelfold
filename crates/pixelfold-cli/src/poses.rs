//! `pixelfold poses`: ligand poses ranked against a reference by the contacts
//! the ligand makes, with the clashes each suffers.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use anyhow::Result;
use pixelfold_core::components::Dictionary;
use pixelfold_core::interactions::detect;
use pixelfold_core::numbering::{self, Chain, ResidueKey, ResidueMap};
use pixelfold_core::pose::{self, Attachments, Contact, Fingerprint};
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
    pub contacts: usize,
    pub shared: usize,
    pub recovery: f32,
    pub similarity: f32,
    pub clashes: usize,
    pub lost: Vec<String>,
}

impl Row for PoseRecord {
    fn columns() -> Vec<&'static str> {
        vec![
            "pose",
            "model",
            "contacts",
            "shared",
            "recovery",
            "similarity",
            "clashes",
            "lost",
        ]
    }

    fn cells(&self) -> Vec<String> {
        vec![
            self.pose.clone(),
            self.model.to_string(),
            self.contacts.to_string(),
            self.shared.to_string(),
            format!("{:.3}", self.recovery),
            format!("{:.3}", self.similarity),
            self.clashes.to_string(),
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
/// covalent attachments, the chains a pose's are aligned with, and the
/// chemistry a pose may borrow.
pub struct Reference {
    components: Dictionary,
    fingerprint: Fingerprint,
    attachments: Attachments,
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

        let fingerprint = pose::fingerprint(protein, &detect(protein), &ligand);
        if fingerprint.is_empty() {
            anyhow::bail!("the ligand `{query}` makes no contacts to compare poses with");
        }

        Ok(Self {
            components: protein.components.clone(),
            fingerprint,
            attachments: pose::attachments(protein, &ligand),
            chains: numbering::polymer_chains(protein),
        })
    }

    /// Score model `model` of the pose file `pose` against the reference.
    ///
    /// A ligand the pose file does not define takes the reference's definition
    /// of the same component, matched by atom name. The pose's residues are
    /// matched with the reference's by aligning chain sequences.
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

        let reference = &self.fingerprint;
        let found = pose::fingerprint(&protein, &detect(&protein), &ligand);
        let chains = numbering::polymer_chains(&protein);
        let map = pose::residue_map(&found, &chains, reference, &self.chains);
        self.check_pairing(&map, pose, warnings)?;
        let fingerprint = pose::renumbered(&found, reference, &map);

        Ok(PoseRecord {
            pose: pose.to_string(),
            model,
            contacts: fingerprint.len(),
            shared: fingerprint.shared(reference),
            recovery: rounded(fingerprint.recovery(reference), 3),
            similarity: rounded(fingerprint.similarity(reference), 3),
            clashes: pose::clashes(&protein, &ligand, &self.attachments.renumbered(&map)),
            lost: fingerprint.lost(reference).map(contact_label).collect(),
        })
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
            contacts: 0,
            shared: 0,
            recovery,
            similarity,
            clashes,
            lost: Vec::new(),
        }
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
