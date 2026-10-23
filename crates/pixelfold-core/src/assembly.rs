//! Detects whether a structure file declares a biological assembly that the
//! deposited coordinates do not show. The assemblies are `mmcif-model`'s, read
//! from mmCIF's assembly categories or a PDB file's `REMARK 350`; building one is
//! the loader's.

use std::collections::HashSet;

use mmcif_model::{Assembly, Structure};

use crate::structure::BiologicalAssembly;

/// Matrices are written to a few decimal places; two genuinely different
/// operators differ by far more than this.
const TRANSFORM_TOLERANCE: f64 = 1e-3;

/// Detect a [`primary`] assembly the deposited coordinates do not already show:
/// one that places a chain more than once, or places its chains by different
/// transforms. One transform for every chain, identity or not, moves the whole
/// unit rigidly and leaves its arrangement as deposited. `None` when the
/// coordinates are the biological unit, when no assembly is declared, and when
/// the assembly names an operator the file does not define.
pub fn detect(structure: &Structure) -> Option<BiologicalAssembly> {
    let primary = primary(structure)?;
    let deposited = structure
        .assembly_is_deposited(primary.id().as_str(), TRANSFORM_TOLERANCE)
        .unwrap_or(true);
    if deposited {
        return None;
    }

    let combinations: HashSet<_> = primary
        .generators()
        .iter()
        .flat_map(|generator| generator.operators().combinations())
        .collect();
    Some(BiologicalAssembly {
        operator_count: combinations.len(),
        oligomer: primary.oligomeric_details().map(str::to_owned),
    })
}

/// The primary assembly: "1", the id RCSB gives the preferred one, or else the
/// first declared.
pub(crate) fn primary(structure: &Structure) -> Option<&Assembly> {
    let assemblies = structure.assemblies();
    assemblies
        .iter()
        .find(|assembly| *assembly.id() == "1")
        .or_else(|| assemblies.first())
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// Operators 1 to `count`: operator `k` translates by `k - 1` along x, so
    /// operator 1 is the identity.
    fn operators(count: usize) -> String {
        let mut text = String::from(
            "loop_
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
",
        );
        for k in 1..=count {
            text += &format!("{k} op 1 0 0 {} 0 1 0 0 0 0 1 0\n", k - 1);
        }
        text + "#\n"
    }

    fn detect_mmcif(records: &str, operator_count: usize) -> Option<BiologicalAssembly> {
        let text = format!("{ATOMS}{}{records}", operators(operator_count));
        detect(&mmcif_model::parse_mmcif(&text).unwrap())
    }

    fn detect_pdb(remarks: &str) -> Option<BiologicalAssembly> {
        let atom =
            "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  5.00           C\n";
        detect(&mmcif_model::parse_pdb(&format!("{remarks}{atom}")).unwrap())
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

    // An assembly naming an operator the file does not define cannot be judged.
    const MMCIF_UNDEFINED_OPERATOR: &str = "\
_pdbx_struct_assembly_gen.assembly_id       1
_pdbx_struct_assembly_gen.oper_expression   1,7
_pdbx_struct_assembly_gen.asym_id_list      A
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

    #[test]
    fn mmcif_multimeric_assembly_warns_with_oligomer() {
        let assembly = detect_mmcif(MMCIF_DIMER, 2).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        assert_eq!(assembly.oligomer.as_deref(), Some("dimeric"));
    }

    #[test]
    fn mmcif_monomeric_assembly_does_not_warn() {
        assert_eq!(detect_mmcif(MMCIF_MONOMER, 2), None);
    }

    #[test]
    fn mmcif_identity_split_across_records_does_not_warn() {
        assert_eq!(detect_mmcif(MMCIF_TWO_IDENTITY_RECORDS, 2), None);
    }

    #[test]
    fn mmcif_single_record_multimeric_warns() {
        let assembly = detect_mmcif(MMCIF_SINGLE_RECORD_DIMER, 2).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        assert_eq!(assembly.oligomer.as_deref(), Some("dimeric"));
    }

    #[test]
    fn mmcif_single_record_complete_asu_does_not_warn() {
        assert_eq!(detect_mmcif(MMCIF_SINGLE_RECORD_COMPLETE, 2), None);
    }

    #[test]
    fn mmcif_operator_ranges_count_every_combination() {
        let assembly = detect_mmcif(MMCIF_ICOSAHEDRAL, 60).expect("capsid should warn");
        assert_eq!(assembly.operator_count, 60);
        assert_eq!(assembly.oligomer, None);
    }

    #[test]
    fn an_assembly_with_an_undefined_operator_is_not_judged() {
        assert_eq!(detect_mmcif(MMCIF_UNDEFINED_OPERATOR, 2), None);
    }

    #[test]
    fn no_mmcif_assembly_returns_none() {
        assert_eq!(detect_mmcif("", 2), None);
    }

    #[test]
    fn pdb_dimer_warns_with_author_unit() {
        let assembly = detect_pdb(PDB_DIMER).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        // Lower-cased as mmCIF writes it, so both formats of one entry agree.
        assert_eq!(assembly.oligomer.as_deref(), Some("dimeric"));
    }

    #[test]
    fn pdb_monomer_does_not_warn() {
        assert_eq!(detect_pdb(PDB_MONOMER), None);
    }

    #[test]
    fn no_pdb_assembly_records_returns_none() {
        assert_eq!(detect_pdb(""), None);
    }
}
