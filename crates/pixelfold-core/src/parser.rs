use anyhow::Result;
use glam::Vec3;
use pdbtbx::*;
use std::{
    collections::BTreeSet,
    io::{BufReader, Cursor},
    path::Path,
};

use crate::{
    dssp,
    fixed_str::FixedStr,
    structure::{AltlocPolicy, Atom, Protein, SecondaryStructure, filter_altlocs},
};

/// Build `Atom` from a pdbtbx hierarchy item.
fn build_atom(hierarchy: &impl ContainsAtomConformerResidueChainModel) -> Atom {
    let atom = hierarchy.atom();
    let residue = hierarchy.residue();
    let conformer = hierarchy.conformer();
    let chain = hierarchy.chain();
    let element = atom
        .element()
        .map(|e| e.symbol().to_string())
        .unwrap_or_else(|| infer_element(atom.name()));
    let altloc = conformer
        .alternative_location()
        .and_then(|s| s.chars().next());
    let insertion_code = residue.insertion_code().and_then(|s| s.chars().next());
    let position = Vec3::new(atom.x() as f32, atom.y() as f32, atom.z() as f32);

    Atom {
        serial: atom.serial_number() as u32,
        name: FixedStr::new(atom.name()),
        element: FixedStr::new(&element),
        residue_name: FixedStr::new(conformer.name()),
        residue_seq: residue.serial_number() as u32,
        chain_id: FixedStr::new(chain.id()),
        is_hetatm: atom.hetero(),
        altloc,
        occupancy: atom.occupancy() as f32,
        insertion_code,
        model_number: hierarchy.model().serial_number(),
        position,
        b_factor: atom.b_factor() as f32,
        secondary_structure: SecondaryStructure::Coil,
    }
}

/// Best-effort element inference from an atom name when the element column is
/// absent. Uses the leading letter, which is correct for standard protein
/// atoms but cannot disambiguate metals (for example a calcium named "CA")
/// without the column; those are expected to carry an explicit element symbol.
fn infer_element(name: &str) -> String {
    name.chars()
        .find(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase().to_string())
        .unwrap_or_default()
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

/// Read a structure with pdbtbx (Loose strictness, first model only), retrying
/// past one specific spurious failure.
///
/// pdbtbx raises an `InvalidatingError` when a file names one crystal space
/// group two non-identical ways. Crystal symmetry is never used, so retry with
/// the symmetry records stripped.
fn read_structure(path_str: &str) -> Result<(PDB, Vec<PDBError>)> {
    let first = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .read(path_str);

    match first {
        Ok(loaded) => Ok(loaded),
        Err(errors) if errors.iter().any(is_space_group_mismatch) => {
            let text = std::fs::read_to_string(path_str)
                .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;
            ReadOptions::default()
                .set_level(StrictnessLevel::Loose)
                .set_only_first_model(true)
                .set_format(Format::Mmcif)
                .read_raw(BufReader::new(Cursor::new(strip_symmetry_records(&text))))
                .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))
        }
        Err(errors) => Err(anyhow::anyhow!("Failed to open file: {:?}", errors)),
    }
}

fn is_space_group_mismatch(error: &PDBError) -> bool {
    error.short_description() == "Space group does not match"
}

/// Drop the single-value crystal space-group records that pdbtbx cross-checks
/// against one another. Leaves the `_space_group_symop` loop and the unit cell
/// intact; pixelfold reads none of them.
fn strip_symmetry_records(cif: &str) -> String {
    const DROP: [&str; 6] = [
        "_symmetry.Int_Tables_number",
        "_symmetry.space_group_name_H-M",
        "_symmetry.space_group_name_Hall",
        "_space_group.name_H-M_alt",
        "_space_group.name_Hall",
        "_space_group.IT_number",
    ];
    let mut out = String::with_capacity(cif.len());
    let mut lines = cif.lines().peekable();
    while let Some(line) = lines.next() {
        let trimmed = line.trim_start();
        if !DROP.iter().any(|tag| trimmed.starts_with(tag)) {
            out.push_str(line);
            out.push('\n');
            continue;
        }

        // Drop the key line.
        let inline_value = trimmed
            .split_once(char::is_whitespace)
            .is_some_and(|(_, rest)| !rest.trim().is_empty());
        if inline_value {
            continue;
        }

        match lines.peek() {
            Some(next) if next.trim_start().starts_with(';') => {
                lines.next(); // opening ';' line
                for block in lines.by_ref() {
                    if block.trim_start().starts_with(';') {
                        break; // closing ';' line
                    }
                }
            }
            Some(_) => {
                lines.next(); // value on its own line
            }
            None => {}
        }
    }
    out
}

/// Load a protein structure from a PDB or mmCIF file
pub fn load_protein<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_with_options(path, AltlocPolicy::default())
}

/// Load a protein structure with additional options
pub fn load_protein_with_options<P: AsRef<Path>>(path: P, altloc: AltlocPolicy) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = read_structure(path_str)?;

    if !errors.is_empty() {
        eprintln!("Warning: {} errors while parsing file", errors.len());
    }

    let title = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    // Build atoms, noting any identifier that overflowed its fixed-capacity
    // field.
    let mut truncated: BTreeSet<(&'static str, String)> = BTreeSet::new();
    let atoms: Vec<Atom> = pdb
        .atoms_with_hierarchy()
        .map(|hierarchy| {
            let atom = build_atom(&hierarchy);
            if atom.name.as_str() != hierarchy.atom().name() {
                truncated.insert(("atom name", hierarchy.atom().name().to_string()));
            }
            if atom.residue_name.as_str() != hierarchy.conformer().name() {
                truncated.insert(("residue name", hierarchy.conformer().name().to_string()));
            }
            if atom.chain_id.as_str() != hierarchy.chain().id() {
                truncated.insert(("chain id", hierarchy.chain().id().to_string()));
            }

            atom
        })
        .collect();

    if let Some(message) = truncation_warning(&truncated) {
        eprintln!("Warning: {message}");
    }

    let mut atoms = filter_altlocs(atoms, altloc);

    let hbonds = dssp::assign(&mut atoms);

    // Two things live only in the raw file, past what pdbtbx returns: a
    // biological unit that differs from the deposited coordinates, and the
    // chemical definitions of the components used.
    let file_text = std::fs::read_to_string(path).ok();
    let assembly = file_text
        .as_deref()
        .and_then(crate::assembly::detect_partial_assembly);
    let components = file_text
        .as_deref()
        .map(crate::components::Dictionary::parse)
        .unwrap_or_default();

    Ok(Protein {
        atoms,
        title,
        hbonds,
        assembly,
        components,
    })
}

/// Load only backbone atoms (CA, C, N, O) for simplified rendering
pub fn load_protein_backbone<P: AsRef<Path>>(path: P) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = read_structure(path_str)?;

    if !errors.is_empty() {
        eprintln!("Warning: {} errors while parsing file", errors.len());
    }

    let title = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    let mut atoms: Vec<Atom> = pdb
        .atoms_with_hierarchy()
        .filter(|hierarchy| {
            let name = hierarchy.atom().name();
            name == "CA" || name == "C" || name == "N" || name == "O"
        })
        .map(|hierarchy| build_atom(&hierarchy))
        .collect();

    let hbonds = dssp::assign(&mut atoms);

    Ok(Protein {
        atoms,
        title,
        hbonds,
        assembly: None,
        components: Default::default(),
    })
}

/// Load only CA (alpha carbon) atoms for minimal rendering of large proteins
pub fn load_protein_ca_only<P: AsRef<Path>>(path: P) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = read_structure(path_str)?;

    if !errors.is_empty() {
        eprintln!("Warning: {} errors while parsing file", errors.len());
    }

    let title = pdb.identifier.as_deref().unwrap_or("Unknown").to_string();

    let mut atoms: Vec<Atom> = pdb
        .atoms_with_hierarchy()
        .filter(|hierarchy| hierarchy.atom().name() == "CA")
        .map(|hierarchy| build_atom(&hierarchy))
        .collect();

    let hbonds = dssp::assign(&mut atoms);

    Ok(Protein {
        atoms,
        title,
        hbonds,
        assembly: None,
        components: Default::default(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

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
        let (pdb, _errors) = ReadOptions::default()
            .set_level(StrictnessLevel::Loose)
            .set_format(Format::Mmcif)
            .read_raw(BufReader::new(Cursor::new(cif)))
            .expect("minimal mmCIF should parse");

        let atoms: Vec<Atom> = pdb.atoms_with_hierarchy().map(|h| build_atom(&h)).collect();
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
        // A 5-character chain id exceeds the 4-byte field and must truncate
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
ATOM 1 C CA ALA AAAAA AAAAA 1 0.0 0.0 0.0 1.0 0.0
ATOM 2 C CA ALA AAAAA AAAAA 2 3.8 0.0 0.0 1.0 0.0
";
        let path = std::env::temp_dir().join("pf_truncation_test.cif");
        std::fs::write(&path, cif).unwrap();
        let protein = load_protein_with_options(&path, AltlocPolicy::All).unwrap();

        assert_eq!(protein.atoms.len(), 2);
        assert_eq!(protein.atoms[0].chain_id.as_str(), "AAAA"); // truncated from "AAAAA"
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn strip_symmetry_records_drops_names_keeps_symop_and_cell() {
        // The Hall symbol is a multi-line `;`-delimited value, exactly as in 6ZPA.
        let cif = "\
_symmetry.Int_Tables_number 155
_symmetry.space_group_name_Hall
;R 3 2\"
;
_symmetry.space_group_name_H-M 'H 3 2'
_space_group.name_H-M_alt 'R 3 2 :H'
_space_group.IT_number 155
_cell.length_a 141.489
loop_
_space_group_symop.id
_space_group_symop.operation_xyz
1 x,y,z
";
        let cleaned = strip_symmetry_records(cif);
        assert!(!cleaned.contains("space_group_name_H-M"));
        assert!(!cleaned.contains("name_H-M_alt"));
        assert!(!cleaned.contains("Int_Tables_number"));
        assert!(!cleaned.contains("_space_group.IT_number"));
        // The multi-line Hall value is removed with its key, not orphaned.
        assert!(!cleaned.contains("space_group_name_Hall"));
        assert!(!cleaned.contains(";R 3 2"));
        // The symop loop and unit cell are untouched.
        assert!(cleaned.contains("_space_group_symop.operation_xyz"));
        assert!(cleaned.contains("1 x,y,z"));
        assert!(cleaned.contains("_cell.length_a 141.489"));
    }

    // A file that names one space group two inconsistent ways (as 6ZPA does with
    // 'H 3 2' vs 'R 3 2 :H') makes pdbtbx fail even at Loose strictness; the
    // loader must recover by stripping the symmetry records.
    #[test]
    fn loader_recovers_from_conflicting_space_group() {
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
        // pdbtbx rejects the conflicting records outright...
        let raw = ReadOptions::default()
            .set_level(StrictnessLevel::Loose)
            .set_format(Format::Mmcif)
            .read_raw(BufReader::new(Cursor::new(cif.to_string())));
        assert!(
            raw.is_err(),
            "expected pdbtbx to reject the conflicting space group"
        );

        // ...but the loader strips symmetry and still returns the atom.
        let path = std::env::temp_dir().join("pf_conflicting_spacegroup_test.cif");
        std::fs::write(&path, cif).unwrap();
        let protein = load_protein_with_options(&path, AltlocPolicy::All).unwrap();
        assert_eq!(protein.atoms.len(), 1);
        assert_eq!(protein.atoms[0].chain_id.as_str(), "A");
        let _ = std::fs::remove_file(&path);
    }
}
