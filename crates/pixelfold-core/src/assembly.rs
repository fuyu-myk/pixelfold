//! Detects whether a structure file declares a biological assembly that differs
//! from the deposited asymmetric unit.
//!
//! pdbtbx does not parse assembly records, so this reads them straight from the
//! raw file text: mmCIF `_pdbx_struct_assembly` / `_pdbx_struct_assembly_gen`,
//! or legacy PDB `REMARK 350` `BIOMT` records. It reports only whether the
//! primary assembly applies symmetry beyond the identity (so the coordinates on
//! screen are partial); it does not generate the assembly.

use crate::mmcif::{read_loop, read_value};
use crate::structure::BiologicalAssembly;

/// Detect a biological assembly that requires symmetry expansion beyond the
/// deposited coordinates. Returns `Some` only when the primary assembly applies
/// more than one operator; `None` when the deposited coordinates already are the
/// biological unit, no assembly is declared, or the text cannot be interpreted.
pub fn detect_partial_assembly(file_text: &str) -> Option<BiologicalAssembly> {
    let lines: Vec<&str> = file_text.lines().collect();

    if lines
        .iter()
        .any(|l| l.trim_start().starts_with("_pdbx_struct_assembly"))
    {
        return detect_mmcif(&lines);
    }

    if lines.iter().any(|l| l.starts_with("REMARK 350")) {
        return detect_pdb(&lines);
    }

    None
}

fn detect_mmcif(lines: &[&str]) -> Option<BiologicalAssembly> {
    let gens = assembly_gen_records(lines);
    if gens.is_empty() {
        return None;
    }

    // RCSB numbers the preferred assembly "1"; fall back to the first declared.
    let primary_id = if gens.iter().any(|(id, _)| id == "1") {
        "1".to_string()
    } else {
        gens[0].0.clone()
    };

    // Warn when any single generation record expands the coordinates. Summing
    // across records would misfire on structures whose asymmetric unit already
    // holds the full complex as several identity-mapped chain groups.
    let operator_count = gens
        .iter()
        .filter(|(id, _)| *id == primary_id)
        .map(|(_, expr)| count_operators(expr))
        .max()
        .unwrap_or(0);

    if operator_count <= 1 {
        return None;
    }

    let oligomer = oligomeric_details(lines, &primary_id);
    Some(BiologicalAssembly {
        operator_count,
        oligomer,
    })
}

/// `(assembly_id, oper_expression)` pairs from `_pdbx_struct_assembly_gen`,
/// handling both the looped and single-record mmCIF layouts.
fn assembly_gen_records(lines: &[&str]) -> Vec<(String, String)> {
    const PREFIX: &str = "_pdbx_struct_assembly_gen.";

    if let Some((cols, rows)) = read_loop(lines, PREFIX) {
        let id = cols.iter().position(|c| c == "assembly_id");
        let expr = cols.iter().position(|c| c == "oper_expression");
        if let (Some(id), Some(expr)) = (id, expr) {
            return rows
                .iter()
                .filter_map(|r| Some((r.get(id)?.clone(), r.get(expr)?.clone())))
                .collect();
        }

        return Vec::new();
    }

    let id = read_value(lines, "_pdbx_struct_assembly_gen.assembly_id");
    let expr = read_value(lines, "_pdbx_struct_assembly_gen.oper_expression");
    match (id, expr) {
        (Some(id), Some(expr)) => vec![(id, expr)],
        _ => Vec::new(),
    }
}

/// `oligomeric_details` for the given assembly id, from either mmCIF layout.
fn oligomeric_details(lines: &[&str], assembly_id: &str) -> Option<String> {
    const PREFIX: &str = "_pdbx_struct_assembly.";

    let value = if let Some((cols, rows)) = read_loop(lines, PREFIX) {
        let id = cols.iter().position(|c| c == "id")?;
        let details = cols.iter().position(|c| c == "oligomeric_details")?;
        rows.iter()
            .find(|r| r.get(id).map(String::as_str) == Some(assembly_id))
            .and_then(|r| r.get(details).cloned())
    } else {
        match read_value(lines, "_pdbx_struct_assembly.id") {
            Some(id) if id == assembly_id => {
                read_value(lines, "_pdbx_struct_assembly.oligomeric_details")
            }
            _ => None,
        }
    };

    value.filter(|s| !s.is_empty() && s != "?" && s != ".")
}

fn detect_pdb(lines: &[&str]) -> Option<BiologicalAssembly> {
    let mut in_biomolecule_1 = false;
    let mut operator_count = 0usize;
    let mut oligomer: Option<String> = None;

    for line in lines {
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

/// Count the distinct transforms an `oper_expression` produces: comma-separated
/// terms, `a-b` ranges, and the Cartesian product of parenthesised groups.
fn count_operators(expr: &str) -> usize {
    let expr = expr.trim();

    let groups: Vec<String> = if expr.contains('(') {
        let mut groups = Vec::new();
        let mut depth = 0usize;
        let mut current = String::new();
        for c in expr.chars() {
            match c {
                '(' => {
                    depth += 1;
                    if depth == 1 {
                        current.clear();
                    } else {
                        current.push(c);
                    }
                }
                ')' => {
                    if depth == 1 {
                        groups.push(std::mem::take(&mut current));
                    } else if depth > 1 {
                        current.push(c);
                    }

                    depth = depth.saturating_sub(1);
                }
                _ if depth >= 1 => current.push(c),
                _ => {}
            }
        }

        groups
    } else {
        vec![expr.to_string()]
    };

    if groups.is_empty() {
        return 0;
    }

    groups.iter().map(|g| count_group(g)).product()
}

fn count_group(group: &str) -> usize {
    group.split(',').map(count_term).sum()
}

fn count_term(term: &str) -> usize {
    let term = term.trim();
    if term.is_empty() {
        return 0;
    }

    if let Some((start, end)) = term.split_once('-')
        && let (Ok(start), Ok(end)) = (start.trim().parse::<i64>(), end.trim().parse::<i64>())
        && end >= start
    {
        return (end - start + 1) as usize;
    }

    1
}

#[cfg(test)]
mod tests {
    use super::*;

    const MMCIF_DIMER: &str = "\
data_TEST
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
data_TEST
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
    // assembly applies identity to two chain groups in separate records. The
    // per-record operator count is one, so this must not warn.
    const MMCIF_TWO_IDENTITY_RECORDS: &str = "\
data_TEST
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
data_TEST
#
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
data_TEST
#
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

    #[test]
    fn count_operators_handles_ranges_lists_and_products() {
        assert_eq!(count_operators("1"), 1);
        assert_eq!(count_operators("1,2"), 2);
        assert_eq!(count_operators("1,2,3,4"), 4);
        assert_eq!(count_operators("1-24"), 24);
        assert_eq!(count_operators("(1-60)"), 60);
        assert_eq!(count_operators("(1-60)(P)"), 60);
        assert_eq!(count_operators("(1,2,3)(4,5)"), 6);
    }

    #[test]
    fn mmcif_multimeric_assembly_warns_with_oligomer() {
        let assembly = detect_partial_assembly(MMCIF_DIMER).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        assert_eq!(assembly.oligomer.as_deref(), Some("dimeric"));
    }

    #[test]
    fn mmcif_monomeric_assembly_does_not_warn() {
        assert_eq!(detect_partial_assembly(MMCIF_MONOMER), None);
    }

    #[test]
    fn mmcif_identity_split_across_records_does_not_warn() {
        assert_eq!(detect_partial_assembly(MMCIF_TWO_IDENTITY_RECORDS), None);
    }

    #[test]
    fn mmcif_single_record_multimeric_warns() {
        let assembly =
            detect_partial_assembly(MMCIF_SINGLE_RECORD_DIMER).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        assert_eq!(assembly.oligomer.as_deref(), Some("dimeric"));
    }

    #[test]
    fn mmcif_single_record_complete_asu_does_not_warn() {
        assert_eq!(detect_partial_assembly(MMCIF_SINGLE_RECORD_COMPLETE), None);
    }

    #[test]
    fn pdb_dimer_warns_with_author_unit() {
        let assembly = detect_partial_assembly(PDB_DIMER).expect("dimer should warn");
        assert_eq!(assembly.operator_count, 2);
        assert_eq!(assembly.oligomer.as_deref(), Some("DIMERIC"));
    }

    #[test]
    fn pdb_monomer_does_not_warn() {
        assert_eq!(detect_partial_assembly(PDB_MONOMER), None);
    }

    #[test]
    fn no_assembly_records_returns_none() {
        assert_eq!(
            detect_partial_assembly("ATOM      1  N   ALA A   1\n"),
            None
        );
        assert_eq!(detect_partial_assembly(""), None);
    }
}
