//! `pixelfold poses` without a reference: each pose checked as PoseBusters
//! checks a docked pose, its clashes with the protein, its own geometry and
//! its stereochemistry against the Chemical Component Dictionary, and ranked
//! by the checks it fails, then by the confidence its co-folding model gives
//! it.

use std::cmp::Ordering;
use std::collections::{BTreeSet, HashMap};
use std::path::Path;

use anyhow::Result;
use pixelfold_core::checks;
use pixelfold_core::interactions::detect;
use pixelfold_core::pose::{self, Attachments, Link};
use pixelfold_core::sdf;
use pixelfold_core::smiles::Smiles;
use pixelfold_core::stated::Stated;
use pixelfold_core::{Atom, AtomSet, Component, Dictionary, Protein};
use serde::Serialize;

use crate::confidence::{self, Confidence};
use crate::load;
use crate::poses::{GeometryRecord, Warnings, atom_label, bond_label};
use crate::report::Row;

/// One pose checked, flattened to a row.
#[derive(Serialize)]
pub struct CheckRecord {
    pub pose: String,
    pub model: usize,
    /// The checks the pose fails: `clashes`, `bonds`, `angles`,
    /// `internal_clashes`, `flatness`, `flat_rings`, `stereo`.
    pub failed: Vec<&'static str>,
    /// The checks that could not be made: `geometry` without the ligand's
    /// bond table, `stereo` without a dictionary entry its atoms match.
    pub unchecked: Vec<&'static str>,
    pub contacts: usize,
    pub clashes: usize,
    pub geometry: Option<GeometryRecord>,
    /// The labelled stereocentres the pose inverts and E or Z double bonds it
    /// turns the other way, or `None` when no residue could be checked.
    pub stereo: Option<Vec<String>>,
    pub confidence: Option<Confidence>,
}

impl Row for CheckRecord {
    fn columns() -> Vec<&'static str> {
        vec![
            "pose",
            "model",
            "failed",
            "clashes",
            "geometry",
            "stereo",
            "confidence",
            "contacts",
        ]
    }

    fn cells(&self) -> Vec<String> {
        let unchecked = || "unchecked".to_string();
        vec![
            self.pose.clone(),
            self.model.to_string(),
            self.failed.join(", "),
            self.clashes.to_string(),
            self.geometry
                .as_ref()
                .map_or_else(unchecked, |geometry| geometry.failed().join(", ")),
            self.stereo
                .as_ref()
                .map_or_else(unchecked, |stereo| stereo.join(", ")),
            self.confidence
                .as_ref()
                .map_or(String::new(), |confidence| match confidence.ligand_iptm {
                    Some(iptm) => format!("{:.2} (ligand ipTM {iptm:.2})", confidence.score),
                    None => format!("{:.2}", confidence.score),
                }),
            self.contacts.to_string(),
        ]
    }
}

/// What checking every pose shares: the ligand picked, where dictionary
/// entries are cached, and the entries already read.
pub struct Checker<'a> {
    query: &'a str,
    /// The component `--ligand` names, in capitals, as the dictionary has it.
    code: Option<String>,
    /// The ligand's SMILES, as `--smiles` gives it.
    smiles: Option<Smiles>,
    cache_dir: &'a Path,
    entries: HashMap<String, Option<Component>>,
}

impl<'a> Checker<'a> {
    /// A checker of the `query` atoms of each pose, `code` the component
    /// `--ligand` names, if it does.
    pub fn new(
        query: &'a str,
        code: Option<&str>,
        smiles: Option<Smiles>,
        cache_dir: &'a Path,
    ) -> Self {
        Self {
            query,
            code: code.map(str::to_uppercase),
            smiles,
            cache_dir,
            entries: HashMap::new(),
        }
    }

    /// What `--smiles`, or else a residue's own SMILES in its file, states of
    /// the `ligand` residues of `protein`, with a warning for each string that
    /// does not read or pair.
    fn stated(
        &self,
        protein: &Protein,
        ligand: &AtomSet,
        pose: &str,
        warnings: &mut Warnings,
    ) -> Vec<Stated> {
        let stated = checks::stated(protein, ligand, self.smiles.as_ref());
        for code in &stated.unpaired {
            warnings.warn(match self.smiles {
                Some(_) => format!(
                    "--smiles pairs with no ligand residue of {pose} by bond graph, even through \
                     bonds past bonding distance: another molecule, so it states nothing of them"
                ),
                None => format!(
                    "the SMILES {pose} gives {code} does not pair with its atoms by bond graph, \
                     so it states nothing of them"
                ),
            });
        }
        for (code, error) in &stated.unreadable {
            warnings.warn(format!(
                "the SMILES {pose} gives {code} does not read ({error})"
            ));
        }
        if stated.found.iter().any(Stated::loose) {
            warnings.warn(format!(
                "the SMILES pairs with the ligand in {pose} only through bonds past bonding \
                 distance: a bond stretched or atoms clashing, which its geometry checks show"
            ));
        }

        stated.found
    }

    /// The dictionary's definition of `code`, whatever its case, read once; a
    /// failed download warns and counts as none.
    fn entry(&mut self, code: &str, warnings: &mut Warnings) -> Option<Component> {
        let code = code.to_uppercase();
        let cache_dir = self.cache_dir;
        self.entries
            .entry(code.clone())
            .or_insert_with(|| match load::chemical_component(&code, cache_dir) {
                Ok(entry) => entry.and_then(|entry| entry.get(&code).cloned()),
                Err(error) => {
                    warnings.warn(format!(
                        "no dictionary entry for {code} ({error:#}), so its stereochemistry \
                         goes unchecked"
                    ));
                    None
                }
            })
            .clone()
    }

    /// The dictionary entries of the residue names `codes` that have one,
    /// each under the name as the pose writes it.
    fn entries(&mut self, codes: &BTreeSet<String>, warnings: &mut Warnings) -> Dictionary {
        let mut found = Dictionary::default();
        for code in codes {
            if let Some(component) = self.entry(code, warnings) {
                found = found.with(code, component);
            }
        }
        found
    }

    /// Check model `model` of the pose file `pose`. A ligand residue the file
    /// does not define borrows the dictionary entry of its component when its
    /// atoms match it.
    pub fn check(
        &mut self,
        protein: Protein,
        links: &[Link],
        pose: &str,
        model: usize,
        warnings: &mut Warnings,
    ) -> Result<CheckRecord> {
        let ligand = load::chosen(&protein, Some(self.query))?;
        let codes: BTreeSet<String> = ligand
            .iter()
            .map(|atom| protein.atoms[atom].residue_name.to_string())
            .collect();
        let entries = self.entries(&codes, warnings);
        let (matched, refused) = checks::matching(&protein, &ligand, &entries);
        for code in refused {
            warnings.warn(format!(
                "{code} in {pose} is not the dictionary's {code} by its atoms, so it borrows \
                 nothing from that entry"
            ));
        }
        let protein = Protein {
            components: protein.components.with_fallback(&matched),
            ..protein
        };
        let stated = self.stated(&protein, &ligand, pose, warnings);
        let protein = checks::with_stated(protein, &stated);
        let confidence = ligand_chain(&protein, &ligand)
            .map(|chain| confidence::beside(Path::new(pose), &chain, &chain_order(&protein)))
            .transpose()?
            .flatten();

        let laid = Laid {
            attached: pose::linked(&protein, &ligand, links),
            protein: &protein,
            ligand: &ligand,
            ideals: &matched,
            stated: &stated,
        };

        Ok(scored(&laid, pose, model, confidence, warnings))
    }

    /// Check record `model` of the SD file `pose`, its ligand laid into
    /// `receptor`: named after the dictionary entry `--ligand` names where its
    /// bond graph pairs with it, and defined by its own bond table.
    pub fn check_docked(
        &mut self,
        record: &sdf::Record,
        receptor: &Protein,
        pose: &str,
        model: usize,
        warnings: &mut Warnings,
    ) -> CheckRecord {
        let Some(code) = self.code.clone() else {
            if self.smiles.is_none() {
                warnings.warn(
                    "the .sdf poses are checked on their own bonds alone; pass --ligand CODE or \
                     --smiles to check their stereochemistry"
                        .to_string(),
                );
            }
            let docked = record.defined(receptor, None);
            let stated = self.given_stated(&docked.protein, &docked.ligand, pose, warnings);
            let laid = Laid {
                protein: &docked.protein,
                ligand: &docked.ligand,
                ideals: &Dictionary::default(),
                attached: Attachments::default(),
                stated: &stated,
            };
            return scored(&laid, pose, model, None, warnings);
        };
        let entry = self.entry(&code, warnings);
        let reference = entry
            .as_ref()
            .and_then(|component| checks::ideal_reference(&code, component));
        let docked = record.defined(receptor, reference.as_ref());
        let ideals = match (docked.naming, entry) {
            (sdf::Naming::Unpaired, _) | (_, None) => {
                if reference.is_some() {
                    warnings.warn(format!(
                        "a ligand in {pose} is not {code} by its bond graph, so it is checked \
                         on its own bonds alone"
                    ));
                }
                Dictionary::default()
            }
            (naming, Some(component)) => {
                if naming == sdf::Naming::BondOrdersDiffer {
                    warnings.warn(format!(
                        "a ligand in {pose} has {code}'s bonds but not its bond orders, as \
                         another tautomer would; its geometry is checked on its own"
                    ));
                }
                Dictionary::default().with(&code, component)
            }
        };

        let stated = self.given_stated(&docked.protein, &docked.ligand, pose, warnings);
        let laid = Laid {
            protein: &docked.protein,
            ligand: &docked.ligand,
            ideals: &ideals,
            attached: Attachments::default(),
            stated: &stated,
        };
        scored(&laid, pose, model, None, warnings)
    }

    /// What `--smiles` states of a docked ligand, which carries its own
    /// definition: nothing without the flag.
    fn given_stated(
        &self,
        protein: &Protein,
        ligand: &AtomSet,
        pose: &str,
        warnings: &mut Warnings,
    ) -> Vec<Stated> {
        match self.smiles {
            Some(_) => self.stated(protein, ligand, pose, warnings),
            None => Vec::new(),
        }
    }

    /// `structure`, given as the protein .sdf poses were docked into, without
    /// its own copies of the component `--ligand` names, whose place the
    /// poses take.
    pub fn receptor(&self, structure: Protein, warnings: &mut Warnings) -> Protein {
        let Some(code) = &self.code else {
            return structure;
        };
        let own = |atom: &Atom| {
            !atom.is_polymer() && atom.residue_name.as_str().eq_ignore_ascii_case(code)
        };
        if structure.atoms.iter().any(own) {
            warnings.warn(format!(
                "the receptor's own {code} is left out, as the poses take its place"
            ));
        }

        structure.without(|_, atom| own(atom))
    }
}

/// The chain of the ligand's first atom, which a co-folding model's
/// confidence names it by.
fn ligand_chain(protein: &Protein, ligand: &AtomSet) -> Option<String> {
    ligand
        .iter()
        .next()
        .map(|atom| protein.atoms[atom].chain_id.to_string())
}

/// The model's chains in the order its atoms hold them.
fn chain_order(protein: &Protein) -> Vec<String> {
    let mut chains: Vec<String> = Vec::new();
    for atom in &protein.atoms {
        let chain = atom.chain_id.to_string();
        if chains.last() != Some(&chain) && !chains.contains(&chain) {
            chains.push(chain);
        }
    }
    chains
}

/// A pose's ligand laid out for checking: the structure it is in, its atoms,
/// the ideal conformers its stereochemistry is held to, and the bonds to the
/// protein its file records.
struct Laid<'a> {
    protein: &'a Protein,
    ligand: &'a AtomSet,
    ideals: &'a Dictionary,
    attached: Attachments,
    /// What a SMILES string states of the ligand's residues.
    stated: &'a [Stated],
}

/// The ligand `laid` out checked.
fn scored(
    laid: &Laid,
    pose: &str,
    model: usize,
    confidence: Option<Confidence>,
    warnings: &mut Warnings,
) -> CheckRecord {
    let Laid {
        protein,
        ligand,
        ideals,
        attached,
        stated,
    } = laid;
    let (protein, ligand, ideals) = (*protein, *ligand, *ideals);
    let geometry = match checks::geometry(protein, ligand) {
        Ok(geometry) => Some(GeometryRecord::from(&geometry)),
        Err(residue) => {
            warnings.warn(format!(
                "{residue} has no chemical definition here, and no dictionary entry its atoms \
                 match, so its geometry goes unchecked"
            ));
            None
        }
    };
    let stereo = checks::stereo(protein, ligand, ideals, stated);
    if !stereo.unchecked.is_empty() {
        let codes: Vec<String> = stereo.unchecked.iter().map(ToString::to_string).collect();
        warnings.warn(format!(
            "the stereochemistry of {} goes unchecked: no dictionary entry its atoms match",
            codes.join(", ")
        ));
    }
    let clashes = pose::clashes(protein, ligand, attached);
    let labels: Vec<String> = stereo
        .inverted
        .iter()
        .map(atom_label)
        .chain(stereo.flipped.iter().map(bond_label))
        .collect();

    let mut failed = Vec::new();
    if clashes > 0 {
        failed.push("clashes");
    }
    if let Some(geometry) = &geometry {
        for (check, count) in [
            ("bonds", geometry.bonds),
            ("angles", geometry.angles),
            ("internal_clashes", geometry.internal_clashes),
            ("flatness", geometry.unflat.unwrap_or(0)),
            ("flat_rings", geometry.flat_rings.unwrap_or(0)),
        ] {
            if count > 0 {
                failed.push(check);
            }
        }
    }
    if !labels.is_empty() {
        failed.push("stereo");
    }
    let mut unchecked = Vec::new();
    if geometry.is_none() {
        unchecked.push("geometry");
    }
    if !stereo.unchecked.is_empty() {
        unchecked.push("stereo");
    }

    CheckRecord {
        pose: pose.to_string(),
        model,
        failed,
        unchecked,
        contacts: pose::fingerprint(protein, &detect(protein), ligand).len(),
        clashes,
        geometry,
        stereo: (stereo.checked > 0).then_some(labels),
        confidence,
    }
}

/// The records best first: the fewest checks failed, a check that could not
/// be made counting as one (a pose distorted past matching its dictionary
/// entry goes unchecked), then the fewest actually failed, then the most
/// confident, then in the order given.
pub fn ranked(mut records: Vec<CheckRecord>) -> Vec<CheckRecord> {
    let score = |record: &CheckRecord| {
        record
            .confidence
            .as_ref()
            .map(|confidence| confidence.score)
    };
    let doubtful = |record: &CheckRecord| record.failed.len() + record.unchecked.len();
    records.sort_by(|a, b| {
        doubtful(a)
            .cmp(&doubtful(b))
            .then_with(|| a.failed.len().cmp(&b.failed.len()))
            .then_with(|| match (score(a), score(b)) {
                (Some(a), Some(b)) => b.partial_cmp(&a).unwrap_or(Ordering::Equal),
                (Some(_), None) => Ordering::Less,
                (None, Some(_)) => Ordering::Greater,
                (None, None) => Ordering::Equal,
            })
    });

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(
        model: usize,
        failed: Vec<&'static str>,
        unchecked: Vec<&'static str>,
        score: Option<f32>,
    ) -> CheckRecord {
        CheckRecord {
            pose: "pose.cif".to_string(),
            model,
            failed,
            unchecked,
            contacts: 0,
            clashes: 0,
            geometry: None,
            stereo: None,
            confidence: score.map(|score| Confidence {
                source: "boltz",
                score,
                ligand_iptm: None,
            }),
        }
    }

    #[test]
    fn a_check_left_unmade_ranks_as_one_failed() {
        let ranked = ranked(vec![
            record(1, vec![], vec!["geometry", "stereo"], Some(0.9)),
            record(2, vec!["bonds"], vec![], Some(0.1)),
            record(3, vec![], vec![], Some(0.2)),
            record(4, vec![], vec![], Some(0.5)),
            record(5, vec!["clashes"], vec!["stereo"], None),
        ]);

        let order: Vec<usize> = ranked.iter().map(|record| record.model).collect();
        assert_eq!(order, [4, 3, 2, 1, 5]);
    }
}
