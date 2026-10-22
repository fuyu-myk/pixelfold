//! Detects whether a structure file declares a biological assembly that differs
//! from the deposited coordinates.
//!
//! An mmCIF file's assemblies are read through `mmcif-model`, operators and
//! all; a legacy PDB file's come from its `REMARK 350` `BIOMT` records. Either
//! way this reports only whether the coordinates on screen are a partial view
//! of the primary assembly; it does not generate the assembly.

use std::collections::{HashMap, HashSet};

use mmcif_cif::Block;
use mmcif_model::{AssemblyId, OperatorId, Structure, Transform};

use crate::structure::BiologicalAssembly;

/// Matrices are written to a few decimal places; two genuinely different
/// operators differ by far more than this.
const TRANSFORM_TOLERANCE: f64 = 1e-3;

/// Detect a primary assembly the deposited coordinates do not already show:
/// one that places a chain more than once, or places its chains by different
/// transforms. One transform for every chain, identity or not, moves the whole
/// unit rigidly and leaves its arrangement as deposited. The primary assembly
/// is "1", the id RCSB gives the preferred one, or else the first declared.
/// `None` when the coordinates are the biological unit or no assembly is
/// declared.
pub fn detect_mmcif(structure: &Structure, block: &Block<'_>) -> Option<BiologicalAssembly> {
    let assemblies = structure.assemblies();
    let primary = assemblies
        .iter()
        .find(|assembly| assembly.id().as_str() == "1")
        .or_else(|| assemblies.first())?;

    let mut placements: HashMap<&str, usize> = HashMap::new();
    let mut combinations: HashSet<Vec<OperatorId>> = HashSet::new();
    let mut first: Option<Transform> = None;
    let mut rearranged = false;
    for generator in primary.generators() {
        for combination in generator.operators().combinations() {
            let transform = compose(structure, &combination);
            match first {
                Some(first) => rearranged |= !same(&first, &transform),
                None => first = Some(transform),
            }
            for asym in generator.asym_ids() {
                *placements.entry(asym.as_str()).or_default() += 1;
            }
            combinations.insert(combination);
        }
    }

    let copied = placements.values().any(|&count| count > 1);
    if !copied && !rearranged {
        return None;
    }

    Some(BiologicalAssembly {
        operator_count: combinations.len(),
        oligomer: oligomeric_details(block, primary.id()),
    })
}

/// The transform a combination applies, its last operator first. An operator
/// the file does not define is taken as the identity, so a file without
/// `_pdbx_struct_oper_list` is judged by its copies alone.
fn compose(structure: &Structure, combination: &[OperatorId]) -> Transform {
    combination.iter().fold(Transform::IDENTITY, |applied, id| {
        structure
            .operator(id.as_str())
            .map_or(applied, |operator| applied.compose(operator.transform()))
    })
}

fn same(a: &Transform, b: &Transform) -> bool {
    let close = |x: f64, y: f64| (x - y).abs() <= TRANSFORM_TOLERANCE;
    (0..3).all(|i| {
        close(a.translation[i], b.translation[i])
            && (0..3).all(|j| close(a.rotation[i][j], b.rotation[i][j]))
    })
}

/// `_pdbx_struct_assembly.oligomeric_details` of the given assembly, which the
/// typed structure does not carry.
fn oligomeric_details(block: &Block<'_>, assembly: &AssemblyId) -> Option<String> {
    let category = block.category("pdbx_struct_assembly")?;
    let ids = category.column("id")?;
    let details = category.column("oligomeric_details")?;

    (0..category.row_count())
        .find(|&row| ids.get(row).and_then(|value| value.as_str()) == Some(assembly.as_str()))
        .and_then(|row| details.get(row)?.as_str())
        .filter(|value| !value.is_empty())
        .map(str::to_owned)
}

/// Detect a biological assembly that needs symmetry expansion from the
/// `REMARK 350` records of a legacy PDB file. Returns `Some` only when the first
/// biomolecule applies more than one operator.
pub fn detect_pdb(file_text: &str) -> Option<BiologicalAssembly> {
    let mut in_biomolecule_1 = false;
    let mut operator_count = 0usize;
    let mut oligomer: Option<String> = None;

    for line in file_text.lines() {
        let Some(rest) = line.strip_prefix("REMARK 350") else {
            continue;
        };
        let rest = rest.trim();

        if let Some(id) = rest.strip_prefix("BIOMOLECULE:") {
            if in_biomolecule_1 {
                break; // reached the next biomolecule; only the first is inspected
            }

            in_biomolecule_1 = id.trim() == "1";
            continue;
        }
        if !in_biomolecule_1 {
            continue;
        }

        if let Some(desc) = rest.strip_prefix("AUTHOR DETERMINED BIOLOGICAL UNIT:") {
            oligomer = Some(desc.trim().to_string()).filter(|s| !s.is_empty());
        } else if oligomer.is_none()
            && let Some(desc) = rest.strip_prefix("SOFTWARE DETERMINED QUATERNARY STRUCTURE:")
        {
            oligomer = Some(desc.trim().to_string()).filter(|s| !s.is_empty());
        }

        // Operator serials restart per chain group, so the largest serial is the
        // operator count of the biggest group. Greater than one means expansion.
        if rest.starts_with("BIOMT1")
            && let Some(serial) = rest.split_whitespace().nth(1).and_then(|s| s.parse().ok())
        {
            operator_count = operator_count.max(serial);
        }
    }

    if operator_count > 1 {
        Some(BiologicalAssembly {
            operator_count,
            oligomer,
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mmcif_cif::Document;

    const ATOMS: &str = "\
data_TEST
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
ATOM 1 C CA GLY A 1 0.0 0.0 0.0
ATOM 2 C CA GLY B 1 5.0 0.0 0.0
ATOM 3 C CA GLY C 1 10.0 0.0 0.0
ATOM 4 C CA GLY D 1 15.0 0.0 0.0
#
";

    fn detect(records: &str) -> Option<BiologicalAssembly> {
        let document = Document::parse(format!("{ATOMS}{records}")).unwrap();
        let block = document.block(0).unwrap();
        let structure = Structure::from_block(&block).unwrap();
        detect_mmcif(&structure, &block)
    }

    const MMCIF_DIMER: &str = "\
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
1 author_defined_assembly dimeric 2
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1,2 A
#
";

    const MMCIF_MONOMER: &str = "\
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
1 author_defined_assembly monomeric 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A
#
";

    // A crystal whose asymmetric unit already holds the full complex: the
    // assembly applies identity to two chain groups in separate records, so
    // this must not warn.
    const MMCIF_TWO_IDENTITY_RECORDS: &str = "\
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A,B
1 1 C,D
#
";

    const PDB_DIMER: &str = "\
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
";

    const PDB_MONOMER: &str = "\
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
";

    // RCSB mmCIF commonly writes single-value records (no `loop_`) with padded
    // values, as every bundled sample file does. The primary assembly expands
    // one chain into a dimer.
    const MMCIF_SINGLE_RECORD_DIMER: &str = "\
_pdbx_struct_assembly.id                   1
_pdbx_struct_assembly.details              author_defined_assembly
_pdbx_struct_assembly.oligomeric_details   dimeric
_pdbx_struct_assembly.oligomeric_count     2
#
_pdbx_struct_assembly_gen.assembly_id       1
_pdbx_struct_assembly_gen.oper_expression   1,2
_pdbx_struct_assembly_gen.asym_id_list      A
#
";

    // The 1UBQ case: the deposited coordinates already hold the stated oligomer
    // (chains A,B), so the generation operator is the identity and nothing is
    // partial. The oligomeric label must not, by itself, trigger a warning.
    const MMCIF_SINGLE_RECORD_COMPLETE: &str = "\
_pdbx_struct_assembly.id                   1
_pdbx_struct_assembly.details              author_defined_assembly
_pdbx_struct_assembly.oligomeric_details   dimeric
_pdbx_struct_assembly.oligomeric_count     2
#
_pdbx_struct_assembly_gen.assembly_id       1
_pdbx_struct_assembly_gen.oper_expression   1
_pdbx_struct_assembly_gen.asym_id_list      A,B
#
";

    // A capsid: one chain group expanded by sixty operators.
    const MMCIF_ICOSAHEDRAL: &str = "\
_pdbx_struct_assembly_gen.assembly_id       1
_pdbx_struct_assembly_gen.oper_expression   (1-60)
_pdbx_struct_assembly_gen.asym_id_list      A,B
#
";

    #[test]
    fn mmcif_multimeric_assembly_warns_with_oligomer() {
        let assembly = detect(MMCIF_DIMER).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        assert_eq!(assembly.oligomer.as_deref(), Some("dimeric"));
    }

    #[test]
    fn mmcif_monomeric_assembly_does_not_warn() {
        assert_eq!(detect(MMCIF_MONOMER), None);
    }

    #[test]
    fn mmcif_identity_split_across_records_does_not_warn() {
        assert_eq!(detect(MMCIF_TWO_IDENTITY_RECORDS), None);
    }

    #[test]
    fn mmcif_single_record_multimeric_warns() {
        let assembly = detect(MMCIF_SINGLE_RECORD_DIMER).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        assert_eq!(assembly.oligomer.as_deref(), Some("dimeric"));
    }

    #[test]
    fn mmcif_single_record_complete_asu_does_not_warn() {
        assert_eq!(detect(MMCIF_SINGLE_RECORD_COMPLETE), None);
    }

    #[test]
    fn mmcif_operator_ranges_count_every_combination() {
        let assembly = detect(MMCIF_ICOSAHEDRAL).expect("capsid should warn");
        assert_eq!(assembly.operator_count, 60);
        assert_eq!(assembly.oligomer, None);
    }

    #[test]
    fn no_mmcif_assembly_returns_none() {
        assert_eq!(detect(""), None);
    }

    #[test]
    fn pdb_dimer_warns_with_author_unit() {
        let assembly = detect_pdb(PDB_DIMER).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        assert_eq!(assembly.oligomer.as_deref(), Some("DIMERIC"));
    }

    #[test]
    fn pdb_monomer_does_not_warn() {
        assert_eq!(detect_pdb(PDB_MONOMER), None);
    }

    #[test]
    fn no_pdb_assembly_records_returns_none() {
        assert_eq!(detect_pdb("ATOM      1  N   ALA A   1\n"), None);
        assert_eq!(detect_pdb(""), None);
    }
}
