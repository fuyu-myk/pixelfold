use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use glam::Vec3;
use mmcif_cif::{Block, BuildError, Document, DocumentBuilder, Entry};
use mmcif_model::{AtomRef, ChainRef, EntityKind, Model, PolymerType, Structure};
use std::{collections::BTreeSet, io::Read, path::Path};

use crate::{
    assembly,
    components::Dictionary,
    dssp,
    fixed_str::FixedStr,
    interactions::topology::{is_standard_nucleotide, is_standard_residue},
    structure::{AltlocPolicy, Atom, Molecule, Protein, SecondaryStructure, filter_altlocs},
};

/// A structure file as read from disk. An mmCIF document is kept whole: the
/// component definitions and the assembly details live in categories the
/// typed structure does not carry.
enum StructureFile {
    Mmcif(Document),
    Pdb(String),
}

impl StructureFile {
    /// Read `path`, gunzipping a `.gz` file first. The extension underneath
    /// picks the format: legacy PDB for `.pdb`, `.ent`, or a numbered assembly
    /// file such as `.pdb1`, and mmCIF for anything else.
    fn read(path: &Path) -> Result<Self> {
        let mut bytes = std::fs::read(path)?;
        let mut name = path;
        if extension_of(path).as_deref() == Some("gz") {
            let mut decompressed = Vec::new();
            MultiGzDecoder::new(bytes.as_slice()).read_to_end(&mut decompressed)?;
            bytes = decompressed;
            name = Path::new(path.file_stem().unwrap_or_default());
        }

        let file = match extension_of(name) {
            Some(extension) if is_pdb_extension(&extension) => Self::Pdb(String::from_utf8(bytes)?),
            _ => Self::Mmcif(Document::from_bytes(bytes)?),
        };

        Ok(file)
    }

    fn block(&self) -> Option<Block<'_>> {
        match self {
            Self::Mmcif(document) => document.block(0),
            Self::Pdb(_) => None,
        }
    }

    /// The typed structure. An mmCIF file whose coordinates read but whose
    /// assembly or entity records do not falls back to its coordinates alone,
    /// with a warning: none of those records is needed to show the structure.
    fn structure(&self) -> Result<Structure> {
        let document = match self {
            Self::Mmcif(document) => document,
            Self::Pdb(text) => return Ok(mmcif_model::parse_pdb(text)?),
        };

        let error = match Structure::from_document(document) {
            Ok(structure) => return Ok(structure),
            Err(error) => error,
        };
        let fallback = document
            .block(0)
            .and_then(|block| coordinates_only(&block).ok())
            .and_then(|reduced| Structure::from_document(&reduced).ok());
        match fallback {
            Some(structure) => {
                eprintln!("Warning: {error}; the file's assembly and entity records were ignored");
                Ok(structure)
            }
            None => Err(error.into()),
        }
    }

    fn title(&self, structure: &Structure) -> String {
        let block = self.block();
        structure
            .entry_id()
            .or_else(|| block.as_ref().map(Block::name))
            .unwrap_or("Unknown")
            .to_string()
    }
}

/// A structure file read and parsed once, from which any of its models can be
/// built. The viewer keeps one to step between models without re-reading.
pub struct ParsedStructure {
    file: StructureFile,
    structure: Structure,
}

impl ParsedStructure {
    pub fn read<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        StructureFile::read(path)
            .and_then(Self::from_file)
            .with_context(|| format!("failed to read {}", path.display()))
    }

    fn from_file(file: StructureFile) -> Result<Self> {
        let structure = file.structure()?;
        Ok(Self { file, structure })
    }

    /// The numbers of the file's models, in file order.
    pub fn model_numbers(&self) -> Vec<usize> {
        self.structure
            .models()
            .iter()
            .map(|model| model_number(model))
            .collect()
    }

    /// Model `number`, or the first when `None`, its alternate locations
    /// resolved by `altloc`, with the biological assembly and the component
    /// definitions the file declares.
    pub fn protein(&self, number: Option<usize>, altloc: AltlocPolicy) -> Result<Protein> {
        let model = self.model(number)?;
        let mut atoms = filter_altlocs(atoms_of(&self.structure, model, |_| true), altloc);
        let hbonds = dssp::assign(&mut atoms);

        let block = self.file.block();
        let assembly = match (&self.file, &block) {
            (StructureFile::Pdb(text), _) => assembly::detect_pdb(text),
            (StructureFile::Mmcif(_), Some(block)) => {
                assembly::detect_mmcif(&self.structure, block)
            }
            (StructureFile::Mmcif(_), None) => None,
        };
        let components = block
            .as_ref()
            .map(Dictionary::from_block)
            .unwrap_or_default();

        Ok(Protein {
            atoms,
            title: self.file.title(&self.structure),
            hbonds,
            assembly,
            components,
        })
    }

    /// The first model's atoms that `keep` accepts by name, every conformer
    /// included, without assembly or component definitions.
    fn subset(&self, keep: impl Fn(&str) -> bool) -> Result<Protein> {
        let mut atoms = atoms_of(&self.structure, first_model(&self.structure)?, keep);
        let hbonds = dssp::assign(&mut atoms);

        Ok(Protein {
            atoms,
            title: self.file.title(&self.structure),
            hbonds,
            assembly: None,
            components: Default::default(),
        })
    }

    fn model(&self, number: Option<usize>) -> Result<&Model> {
        let Some(number) = number else {
            return first_model(&self.structure);
        };

        self.structure
            .models()
            .iter()
            .map(|model| model.as_ref())
            .find(|model| model_number(model) == number)
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "there is no model {number}; the file has {}",
                    describe_models(&self.model_numbers())
                )
            })
    }
}

fn model_number(model: &Model) -> usize {
    usize::try_from(model.number().get()).unwrap_or_default()
}

/// `model 1`, `models 1-10`, or `models 1, 3, 4` where the numbers skip.
fn describe_models(numbers: &[usize]) -> String {
    match numbers {
        [] => "no models".to_string(),
        [only] => format!("model {only}"),
        [first, .., last] if numbers.windows(2).all(|pair| pair[1] == pair[0] + 1) => {
            format!("models {first}-{last}")
        }
        _ => {
            let list: Vec<String> = numbers.iter().map(ToString::to_string).collect();
            format!("models {}", list.join(", "))
        }
    }
}

fn extension_of(path: &Path) -> Option<String> {
    path.extension()
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase)
}

fn is_pdb_extension(extension: &str) -> bool {
    extension == "ent"
        || extension
            .strip_prefix("pdb")
            .is_some_and(|number| number.bytes().all(|b| b.is_ascii_digit()))
}

/// A copy of `block` keeping only its coordinates and entry id.
fn coordinates_only(block: &Block<'_>) -> Result<Document, BuildError> {
    let kept = |tag: &str| {
        ["_atom_site.", "_entry."].iter().any(|prefix| {
            tag.get(..prefix.len())
                .is_some_and(|start| start.eq_ignore_ascii_case(prefix))
        })
    };

    let mut builder = DocumentBuilder::new();
    builder.start_block(block.name())?;
    for entry in block.entries() {
        match entry {
            Entry::Pair(pair) if kept(pair.tag()) => builder.pair(pair.tag(), pair.value())?,
            Entry::Loop(rows) if rows.tags().next().is_some_and(kept) => {
                builder.start_loop(rows.tags())?;
                for value in rows.values() {
                    builder.push(value)?;
                }
            }
            _ => {}
        }
    }

    builder.finish()
}

fn first_model(structure: &Structure) -> Result<&Model> {
    structure
        .first_model()
        .map(|model| model.as_ref())
        .ok_or_else(|| anyhow::anyhow!("the file contains no atoms"))
}

/// What a chain's entity says it is. A polymer whose type the file does not
/// state, as a legacy PDB file's never is, is judged by its residue names.
fn molecule_of(structure: &Structure, chain: ChainRef<'_>) -> Molecule {
    let Some(entity) = chain.entity().and_then(|id| structure.entity(id.as_str())) else {
        return Molecule::Unknown;
    };

    match entity.kind() {
        EntityKind::Polymer => match entity.polymer().map(|polymer| polymer.kind()) {
            Some(PolymerType::PeptideL | PolymerType::PeptideD) => Molecule::Protein,
            Some(PolymerType::Dna | PolymerType::Rna | PolymerType::DnaRnaHybrid) => {
                Molecule::NucleicAcid
            }
            Some(PolymerType::Other(_)) | None => polymer_by_residue_names(chain),
            Some(_) => Molecule::OtherPolymer,
        },
        EntityKind::Branched => Molecule::Branched,
        EntityKind::NonPolymer | EntityKind::Macrolide => Molecule::NonPolymer,
        EntityKind::Water => Molecule::Water,
        EntityKind::Other(_) => Molecule::Unknown,
    }
}

fn polymer_by_residue_names(chain: ChainRef<'_>) -> Molecule {
    let (amino, nucleic) = chain.residues().fold((0, 0), |(amino, nucleic), residue| {
        let name = residue.comp().as_str();
        (
            amino + usize::from(is_standard_residue(name)),
            nucleic + usize::from(is_standard_nucleotide(name)),
        )
    });

    match (amino, nucleic) {
        (0, 0) => Molecule::Unknown,
        _ if amino >= nucleic => Molecule::Protein,
        _ => Molecule::NucleicAcid,
    }
}

/// Build `Atom` from one atom of a model, under author numbering.
fn build_atom(model: &Model, atom: AtomRef<'_>, molecule: Molecule) -> Atom {
    let residue = atom.residue();
    let [x, y, z] = atom.position();

    Atom {
        serial: atom
            .serial()
            .and_then(|serial| u32::try_from(serial).ok())
            .unwrap_or_default(),
        name: FixedStr::new(atom.name().as_str()),
        element: FixedStr::new(atom.element().symbol()),
        residue_name: FixedStr::new(residue.comp().as_str()),
        residue_seq: residue.auth_seq().get(),
        chain_id: FixedStr::new(residue.chain().auth_asym().as_str()),
        is_hetatm: residue.is_hetero(),
        molecule,
        altloc: atom.altloc().get(),
        occupancy: atom.occupancy() as f32,
        insertion_code: residue.ins_code().get(),
        model_number: model_number(model),
        position: Vec3::new(x as f32, y as f32, z as f32),
        b_factor: atom.b_factor() as f32,
        secondary_structure: SecondaryStructure::Coil,
    }
}

/// Every atom of `model` whose name `keep` accepts, warning once when an
/// identifier overflowed its fixed field.
fn atoms_of(structure: &Structure, model: &Model, keep: impl Fn(&str) -> bool) -> Vec<Atom> {
    let molecules: Vec<Molecule> = model
        .chains()
        .map(|chain| molecule_of(structure, chain))
        .collect();
    let mut truncated: BTreeSet<(&'static str, String)> = BTreeSet::new();
    let atoms = model
        .atoms()
        .filter(|source| keep(source.name().as_str()))
        .map(|source| {
            let molecule = molecules[source.residue().chain().index().get()];
            let atom = build_atom(model, source, molecule);
            let residue = source.residue();
            let fields = [
                ("atom name", atom.name.as_str(), source.name().as_str()),
                (
                    "residue name",
                    atom.residue_name.as_str(),
                    residue.comp().as_str(),
                ),
                (
                    "chain id",
                    atom.chain_id.as_str(),
                    residue.chain().auth_asym().as_str(),
                ),
            ];
            for (field, kept, read) in fields {
                if kept != read {
                    truncated.insert((field, read.to_string()));
                }
            }

            atom
        })
        .collect();

    if let Some(message) = truncation_warning(&truncated) {
        eprintln!("Warning: {message}");
    }

    atoms
}

/// Format a one-line warning for identifiers that overflowed their fixed field
/// capacity, or `None` when nothing was truncated.
fn truncation_warning(truncated: &BTreeSet<(&'static str, String)>) -> Option<String> {
    if truncated.is_empty() {
        return None;
    }
    let list: Vec<String> = truncated
        .iter()
        .map(|(field, value)| format!("{field} '{value}'"))
        .collect();
    Some(format!(
        "{} identifier(s) exceeded the fixed field capacity and were truncated: {}",
        truncated.len(),
        list.join(", ")
    ))
}

/// Load a protein structure from a PDB or mmCIF file
pub fn load_protein<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_with_options(path, AltlocPolicy::default())
}

/// Load a protein structure with additional options
pub fn load_protein_with_options<P: AsRef<Path>>(path: P, altloc: AltlocPolicy) -> Result<Protein> {
    let path = path.as_ref();
    ParsedStructure::read(path)?
        .protein(None, altloc)
        .with_context(|| format!("failed to read {}", path.display()))
}

/// Load only backbone atoms (CA, C, N, O) for simplified rendering
pub fn load_protein_backbone<P: AsRef<Path>>(path: P) -> Result<Protein> {
    let path = path.as_ref();
    ParsedStructure::read(path)?
        .subset(|name| matches!(name, "CA" | "C" | "N" | "O"))
        .with_context(|| format!("failed to read {}", path.display()))
}

/// Load only CA (alpha carbon) atoms for minimal rendering of large proteins
pub fn load_protein_ca_only<P: AsRef<Path>>(path: P) -> Result<Protein> {
    let path = path.as_ref();
    ParsedStructure::read(path)?
        .subset(|name| name == "CA")
        .with_context(|| format!("failed to read {}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_atom_reads_element_and_hetatm() {
        // Minimal mmCIF: an ATOM C-alpha carbon (name "CA", element "C") and a
        // HETATM calcium ion (name "CA", element "Ca"). The element column, not
        // the name, must distinguish them.
        let cif = "\
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A 1 0.0 0.0 0.0 1.0 0.0
HETATM 2 Ca CA CA B 2 5.0 0.0 0.0 1.0 0.0
";
        let file = StructureFile::Mmcif(Document::parse(cif).expect("minimal mmCIF should parse"));
        let atoms = ParsedStructure::from_file(file)
            .unwrap()
            .protein(None, AltlocPolicy::All)
            .unwrap()
            .atoms;
        assert_eq!(atoms.len(), 2);

        let carbon = atoms
            .iter()
            .find(|a| !a.is_hetatm)
            .expect("the ATOM record");
        assert_eq!(carbon.name, "CA"); // C-alpha, named CA
        assert_eq!(carbon.element.as_str().to_uppercase(), "C"); // but its element is carbon

        let calcium = atoms
            .iter()
            .find(|a| a.is_hetatm)
            .expect("the HETATM record");
        assert_eq!(calcium.element.as_str().to_uppercase(), "CA"); // calcium, not carbon
    }

    #[test]
    fn truncation_warning_lists_offenders() {
        let mut set = BTreeSet::new();
        assert_eq!(truncation_warning(&set), None);

        set.insert(("chain id", "AAAAA".to_string()));
        let message = truncation_warning(&set).expect("non-empty set warns");
        assert!(message.contains("chain id 'AAAAA'"), "{message}");
        assert!(message.contains("1 identifier"), "{message}");
    }

    #[test]
    fn loader_truncates_overlong_chain_id_without_crashing() {
        // A 9-character chain id exceeds the 8-byte field and must truncate
        // rather than panic; the loader still returns both atoms.
        let cif = "\
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA AAAAAAAAA AAAAAAAAA 1 0.0 0.0 0.0 1.0 0.0
ATOM 2 C CA ALA AAAAAAAAA AAAAAAAAA 2 3.8 0.0 0.0 1.0 0.0
";
        let path = std::env::temp_dir().join("pf_truncation_test.cif");
        std::fs::write(&path, cif).unwrap();
        let protein = load_protein_with_options(&path, AltlocPolicy::All).unwrap();

        assert_eq!(protein.atoms.len(), 2);
        assert_eq!(protein.atoms[0].chain_id.as_str(), "AAAAAAAA"); // truncated from 9
        let _ = std::fs::remove_file(&path);
    }

    // 6ZPA names one space group two inconsistent ways ('H 3 2' vs 'R 3 2 :H').
    // Crystal symmetry is never used, so the coordinates must still load.
    #[test]
    fn conflicting_space_group_records_still_load() {
        let cif = "\
data_test
_symmetry.Int_Tables_number 155
_symmetry.space_group_name_Hall
;R 3 2\"
;
_symmetry.space_group_name_H-M 'H 3 2'
_space_group.name_H-M_alt 'R 3 2 :H'
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A A 1 0.0 0.0 0.0 1.0 0.0
";
        let path = std::env::temp_dir().join("pf_conflicting_spacegroup_test.cif");
        std::fs::write(&path, cif).unwrap();
        let protein = load_protein_with_options(&path, AltlocPolicy::All).unwrap();
        assert_eq!(protein.atoms.len(), 1);
        assert_eq!(protein.atoms[0].chain_id.as_str(), "A");
        let _ = std::fs::remove_file(&path);
    }

    fn load_cif(name: &str, cif: &str, altloc: AltlocPolicy) -> Protein {
        let path = std::env::temp_dir().join(format!("pf_{name}.cif"));
        std::fs::write(&path, cif).unwrap();
        let protein = load_protein_with_options(&path, altloc).unwrap();
        let _ = std::fs::remove_file(&path);
        protein
    }

    const ATOM_SITE_HEADER: &str = "\
data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
";

    // The backbone is shared (no altloc) while the side chain is split, as in
    // 1GCI SER 3. The shared atoms must stay unlabelled at full occupancy, so
    // the side chain's occupancy decides which conformer is kept.
    const SPLIT_SIDE_CHAIN: &str = "\
ATOM 1 N N . SER A 1 0.0 0.0 0.0 1.00 5.0 3 A
ATOM 2 C CA . SER A 1 1.5 0.0 0.0 1.00 5.0 3 A
ATOM 3 C C . SER A 1 2.0 1.4 0.0 1.00 5.0 3 A
ATOM 4 O O . SER A 1 1.3 2.4 0.0 1.00 5.0 3 A
ATOM 5 C CB A SER A 1 2.0 -0.8 1.2 0.47 6.0 3 A
ATOM 6 C CB B SER A 1 2.0 -0.8 -1.2 0.53 6.0 3 A
ATOM 7 O OG A SER A 1 3.4 -0.8 1.2 0.47 6.0 3 A
ATOM 8 O OG B SER A 1 3.4 -0.8 -1.2 0.53 6.0 3 A
";

    #[test]
    fn shared_atoms_of_a_split_residue_stay_unlabelled_at_full_occupancy() {
        let cif = format!("{ATOM_SITE_HEADER}{SPLIT_SIDE_CHAIN}");
        let protein = load_cif("shared_atoms", &cif, AltlocPolicy::All);

        let backbone: Vec<&Atom> = protein
            .atoms
            .iter()
            .filter(|a| matches!(a.name.as_str(), "N" | "CA" | "C" | "O"))
            .collect();
        assert_eq!(backbone.len(), 4);
        assert!(backbone.iter().all(|a| a.altloc.is_none()));
        assert!(backbone.iter().all(|a| a.occupancy == 1.0));
    }

    #[test]
    fn default_policy_keeps_the_more_occupied_side_chain() {
        let cif = format!("{ATOM_SITE_HEADER}{SPLIT_SIDE_CHAIN}");
        let protein = load_cif("split_side_chain", &cif, AltlocPolicy::Occupancy);

        assert_eq!(protein.atoms.len(), 6);
        let cb = protein.atoms.iter().find(|a| a.name == "CB").unwrap();
        assert_eq!(cb.altloc, Some('B'));
    }

    #[test]
    fn policy_b_keeps_the_shared_atoms_of_a_residue_modelled_only_in_a() {
        let cif = format!(
            "{ATOM_SITE_HEADER}\
ATOM 1 N N . ALA A 1 0.0 0.0 0.0 1.00 5.0 1 A
ATOM 2 C CA . ALA A 1 1.5 0.0 0.0 1.00 5.0 1 A
ATOM 3 C CB A ALA A 1 2.0 -0.8 1.2 0.60 6.0 1 A
"
        );
        let protein = load_cif("a_only_residue", &cif, AltlocPolicy::B);

        let names: Vec<&str> = protein.atoms.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, vec!["N", "CA"]);
    }

    #[test]
    fn negative_author_numbers_are_kept() {
        let cif = format!(
            "{ATOM_SITE_HEADER}\
ATOM 1 C CA . GLY A 1 0.0 0.0 0.0 1.00 5.0 -3 A
ATOM 2 C CA . GLY A 2 3.8 0.0 0.0 1.00 5.0 -2 A
"
        );
        let protein = load_cif("negative_numbers", &cif, AltlocPolicy::Occupancy);

        let numbers: Vec<i64> = protein
            .atoms
            .iter()
            .map(|a| i64::from(a.residue_seq))
            .collect();
        assert_eq!(numbers, vec![-3, -2]);
    }

    #[test]
    fn an_unknown_atom_keeps_element_x() {
        let cif = format!(
            "{ATOM_SITE_HEADER}\
HETATM 1 X UNK . UNX B . 0.0 0.0 0.0 1.00 5.0 324 A
"
        );
        let protein = load_cif("unknown_atom", &cif, AltlocPolicy::Occupancy);

        assert!(protein.atoms[0].element.as_str().eq_ignore_ascii_case("X"));
    }

    #[test]
    fn a_numeric_looking_component_id_is_kept_verbatim() {
        let cif = format!(
            "{ATOM_SITE_HEADER}\
HETATM 1 C C1 . 4E2 B . 0.0 0.0 0.0 1.00 5.0 401 A
HETATM 2 Cl CL2 . 4E2 B . 1.7 0.0 0.0 1.00 5.0 401 A
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
4E2 C1 C
4E2 CL2 CL
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
4E2 C1 CL2
#
"
        );
        let protein = load_cif("numeric_component", &cif, AltlocPolicy::Occupancy);

        let residue = protein.atoms[0].residue_name;
        assert_eq!(residue.as_str(), "4E2");
        assert!(protein.components.get(residue.as_str()).is_some());
    }

    const OPERATORS: &str = "\
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0
2 'crystal symmetry operation' -1 0 0 10 0 -1 0 0 0 0 1 0
#
";

    const TWO_CHAINS: &str = "\
ATOM 1 C CA . GLY A 1 0.0 0.0 0.0 1.00 5.0 1 A
ATOM 2 C CA . GLY B 1 9.0 0.0 0.0 1.00 5.0 1 B
#
";

    // 7AEW: the whole asymmetric unit is listed twice, once per operator, in
    // separate generation records.
    #[test]
    fn an_assembly_copying_the_whole_unit_across_records_is_partial() {
        let cif = format!(
            "{ATOM_SITE_HEADER}{TWO_CHAINS}{OPERATORS}\
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A,B
1 2 A,B
#
"
        );
        let protein = load_cif("assembly_copies", &cif, AltlocPolicy::Occupancy);

        assert_eq!(protein.assembly.map(|a| a.operator_count), Some(2));
    }

    // 3BOG: each chain appears once, but one is placed by a symmetry operator,
    // so the deposited coordinates are not the biological arrangement.
    #[test]
    fn an_assembly_moving_one_chain_is_partial() {
        let cif = format!(
            "{ATOM_SITE_HEADER}{TWO_CHAINS}{OPERATORS}\
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A
1 2 B
#
"
        );
        let protein = load_cif("assembly_moves", &cif, AltlocPolicy::Occupancy);

        assert!(protein.assembly.is_some());
    }

    #[test]
    fn an_identity_assembly_split_across_records_is_complete() {
        let cif = format!(
            "{ATOM_SITE_HEADER}{TWO_CHAINS}{OPERATORS}\
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A
1 1 B
#
"
        );
        let protein = load_cif("assembly_identity", &cif, AltlocPolicy::Occupancy);

        assert_eq!(protein.assembly, None);
    }

    // A non-identity operator applied once to every chain moves the whole unit
    // rigidly; the arrangement on screen is still the complete assembly.
    #[test]
    fn a_rigid_move_of_the_whole_unit_is_complete() {
        let cif = format!(
            "{ATOM_SITE_HEADER}{TWO_CHAINS}{OPERATORS}\
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 2 A,B
#
"
        );
        let protein = load_cif("assembly_rigid", &cif, AltlocPolicy::Occupancy);

        assert_eq!(protein.assembly, None);
    }

    // Crystal symmetry and assemblies are never needed to show coordinates, so
    // a record that cannot be read costs the assembly note, not the load.
    #[test]
    fn an_unreadable_assembly_record_does_not_stop_the_load() {
        let cif = format!(
            "{ATOM_SITE_HEADER}{TWO_CHAINS}\
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 ? A,B
#
"
        );
        let protein = load_cif("assembly_unreadable", &cif, AltlocPolicy::Occupancy);

        assert_eq!(protein.atoms.len(), 2);
        assert_eq!(protein.assembly, None);
    }

    #[test]
    fn a_gzipped_file_loads() {
        use std::io::Write;

        let cif = format!("{ATOM_SITE_HEADER}{TWO_CHAINS}");
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(cif.as_bytes()).unwrap();
        let path = std::env::temp_dir().join("pf_gzipped.cif.gz");
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();

        let protein = load_protein_with_options(&path, AltlocPolicy::Occupancy).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(protein.atoms.len(), 2);
    }

    // RCSB serves each biological assembly as `<id>.pdb1`, `<id>.pdb2`, ...
    #[test]
    fn a_numbered_pdb_assembly_file_reads_as_pdb() {
        let pdb =
            "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  5.00           C\n";
        let path = std::env::temp_dir().join("pf_assembly.pdb1");
        std::fs::write(&path, pdb).unwrap();

        let protein = load_protein_with_options(&path, AltlocPolicy::Occupancy).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(protein.atoms.len(), 1);
        assert_eq!(protein.atoms[0].name, "CA");
    }

    fn selected_residues(protein: &Protein, query: &str) -> Vec<String> {
        let set = crate::select::select(query, protein).unwrap();
        let names: BTreeSet<String> = set
            .iter()
            .map(|i| protein.atoms[i].residue_name.as_str().to_string())
            .collect();
        names.into_iter().collect()
    }

    // Selenomethionine is a HETATM record inside the polymer chain: the file's
    // entities say it is protein, its record type and name would say ligand.
    #[test]
    fn entities_classify_a_modified_residue_as_protein() {
        let cif = "\
data_TEST
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
#
_entity_poly.entity_id 1
_entity_poly.type 'polypeptide(L)'
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
ATOM 1 C CA ALA A 1 1 0.0 0.0 0.0 1 A
HETATM 2 C CA MSE A 1 2 3.8 0.0 0.0 2 A
ATOM 3 C CA GLY A 1 3 7.6 0.0 0.0 3 A
HETATM 4 C C1 LIG B 2 . 0.0 9.0 0.0 101 A
HETATM 5 O O HOH C 3 . 0.0 0.0 9.0 201 A
";
        let protein = load_cif("entities", cif, AltlocPolicy::Occupancy);

        assert_eq!(
            selected_residues(&protein, "protein"),
            ["ALA", "GLY", "MSE"]
        );
        assert_eq!(selected_residues(&protein, "ligand"), ["LIG"]);
        assert_eq!(selected_residues(&protein, "water"), ["HOH"]);
    }

    fn pdb_line(record: &str, serial: u32, name: &str, residue: &str, seq: i32, x: f32) -> String {
        let element = &name.trim()[..1];
        format!(
            "{record:<6}{serial:>5} {name:<4} {residue:>3} A{seq:>4}    {x:>8.3}{:>8.3}{:>8.3}  1.00  5.00          {element:>2}\n",
            0.0, 0.0
        )
    }

    // A legacy PDB file names no entities; everything before a chain's last ATOM
    // record is its polymer, and the residue names decide what kind.
    #[test]
    fn a_pdb_files_polymer_keeps_its_modified_residues() {
        let pdb = [
            pdb_line("ATOM", 1, " CA ", "ALA", 1, 0.0),
            pdb_line("HETATM", 2, " CA ", "MSE", 2, 3.8),
            pdb_line("ATOM", 3, " CA ", "GLY", 3, 7.6),
            "TER\n".to_string(),
            pdb_line("HETATM", 4, " C1 ", "LIG", 101, 20.0),
            pdb_line("HETATM", 5, " O  ", "HOH", 201, 30.0),
        ]
        .concat();
        let path = std::env::temp_dir().join("pf_entities.pdb");
        std::fs::write(&path, pdb).unwrap();
        let protein = load_protein_with_options(&path, AltlocPolicy::Occupancy).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(
            selected_residues(&protein, "protein"),
            ["ALA", "GLY", "MSE"]
        );
        assert_eq!(selected_residues(&protein, "ligand"), ["LIG"]);
        assert_eq!(selected_residues(&protein, "water"), ["HOH"]);
    }

    // Without entity records the residue names and record types still decide.
    #[test]
    fn without_entities_a_hetatm_residue_is_a_ligand() {
        let cif = format!(
            "{ATOM_SITE_HEADER}\
ATOM 1 C CA . ALA A 1 0.0 0.0 0.0 1.00 5.0 1 A
HETATM 2 C CA . MSE A 2 3.8 0.0 0.0 1.00 5.0 2 A
"
        );
        let protein = load_cif("no_entities", &cif, AltlocPolicy::Occupancy);

        assert_eq!(selected_residues(&protein, "protein"), ["ALA"]);
        assert_eq!(selected_residues(&protein, "ligand"), ["MSE"]);
    }

    const TWO_MODELS: &str = "\
data_NMR
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA GLY A 1 0.0 0.0 0.0 1
ATOM 2 C CA GLY A 2 3.8 0.0 0.0 1
ATOM 3 C CA GLY A 1 0.5 0.0 0.0 2
ATOM 4 C CA GLY A 2 4.3 0.0 0.0 2
";

    fn parsed(name: &str, cif: &str) -> ParsedStructure {
        let path = std::env::temp_dir().join(format!("pf_{name}.cif"));
        std::fs::write(&path, cif).unwrap();
        let parsed = ParsedStructure::read(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        parsed
    }

    #[test]
    fn any_model_of_a_multi_model_file_loads() {
        let file = parsed("two_models", TWO_MODELS);
        assert_eq!(file.model_numbers(), [1, 2]);

        let first = file.protein(None, AltlocPolicy::Occupancy).unwrap();
        let second = file.protein(Some(2), AltlocPolicy::Occupancy).unwrap();

        assert!(first.atoms.iter().all(|a| a.model_number == 1));
        assert!(second.atoms.iter().all(|a| a.model_number == 2));
        assert_eq!(second.atoms.len(), 2);
        assert_eq!(second.atoms[0].position.x, 0.5);
    }

    #[test]
    fn an_absent_model_is_an_error_naming_the_models_there() {
        let file = parsed("absent_model", TWO_MODELS);

        let Err(error) = file.protein(Some(5), AltlocPolicy::Occupancy) else {
            panic!("model 5 is not in the file");
        };

        assert_eq!(
            error.to_string(),
            "there is no model 5; the file has models 1-2"
        );
    }

    #[test]
    fn model_lists_read_as_ranges_where_they_run_on() {
        assert_eq!(describe_models(&[1]), "model 1");
        assert_eq!(describe_models(&[1, 2, 3]), "models 1-3");
        assert_eq!(describe_models(&[1, 3, 4]), "models 1, 3, 4");
    }
}
