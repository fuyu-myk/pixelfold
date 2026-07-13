use anyhow::Result;
use glam::Vec3;
use pdbtbx::*;
use std::{collections::BTreeSet, path::Path};

use crate::{
    dssp,
    fixed_str::FixedStr,
    sasa::SurfaceCalculator,
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

/// Load a protein structure from a PDB or mmCIF file
pub fn load_protein<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_with_options(path, false, AltlocPolicy::default())
}

/// Load a protein structure with additional options
pub fn load_protein_with_options<P: AsRef<Path>>(
    path: P,
    skip_surface: bool,
    altloc: AltlocPolicy,
) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .read(path_str)
        .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;

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

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    // A structure whose biological unit differs from the deposited asymmetric
    // unit is only detectable from the raw file; a gzip or read error just
    // means no warning.
    let assembly = std::fs::read_to_string(path)
        .ok()
        .and_then(|text| crate::assembly::detect_partial_assembly(&text));

    Ok(Protein {
        atoms,
        title,
        surface_points,
        hbonds,
        assembly,
    })
}

/// Load only backbone atoms (CA, C, N, O) for simplified rendering
pub fn load_protein_backbone<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_backbone_with_options(path, false)
}

/// Load backbone atoms with additional options
pub fn load_protein_backbone_with_options<P: AsRef<Path>>(
    path: P,
    skip_surface: bool,
) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .read(path_str)
        .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;

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

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    Ok(Protein {
        atoms,
        title,
        surface_points,
        hbonds,
        assembly: None,
    })
}

/// Load only CA (alpha carbon) atoms for minimal rendering of large proteins
pub fn load_protein_ca_only<P: AsRef<Path>>(path: P) -> Result<Protein> {
    load_protein_ca_only_with_options(path, false)
}

/// Load CA atoms with additional options
pub fn load_protein_ca_only_with_options<P: AsRef<Path>>(
    path: P,
    skip_surface: bool,
) -> Result<Protein> {
    let path = path.as_ref();
    let path_str = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Invalid path: {}", path.display()))?;

    let (pdb, errors) = ReadOptions::default()
        .set_level(StrictnessLevel::Loose)
        .set_only_first_model(true)
        .read(path_str)
        .map_err(|e| anyhow::anyhow!("Failed to open file: {:?}", e))?;

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

    // Compute solvent-accessible surface (unless skipped for performance)
    let surface_points = if skip_surface {
        Vec::new()
    } else {
        let surface_calculator = SurfaceCalculator::default();
        surface_calculator.calculate_surface(&atoms)
    };

    Ok(Protein {
        atoms,
        title,
        surface_points,
        hbonds,
        assembly: None,
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
        let protein = load_protein_with_options(&path, true, AltlocPolicy::All).unwrap();

        assert_eq!(protein.atoms.len(), 2);
        assert_eq!(protein.atoms[0].chain_id.as_str(), "AAAA"); // truncated from "AAAAA"
        let _ = std::fs::remove_file(&path);
    }
}
