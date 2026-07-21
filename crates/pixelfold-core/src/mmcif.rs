//! Small mmCIF readers for the categories pdbtbx does not expose.
//!
//! pdbtbx returns coordinates; the biological assembly and the chemical
//! component definitions live in categories it drops, so those are read from the
//! raw file text. Only what those readers need is implemented here: looped and
//! single-record categories, with quoted values honoured.

/// Read a looped mmCIF category into `(column names without prefix, data rows)`.
/// Returns the first loop whose columns belong to `prefix`.
pub(crate) fn read_loop(lines: &[&str], prefix: &str) -> Option<(Vec<String>, Vec<Vec<String>>)> {
    let mut i = 0;
    while i < lines.len() {
        if lines[i].trim() != "loop_" {
            i += 1;
            continue;
        }

        let mut columns = Vec::new();
        let mut j = i + 1;
        while j < lines.len() && lines[j].trim_start().starts_with('_') {
            columns.push(lines[j].trim().to_string());
            j += 1;
        }

        if columns.first().is_some_and(|c| c.starts_with(prefix)) {
            let mut rows = Vec::new();
            while j < lines.len() {
                let line = lines[j].trim();
                if line.is_empty()
                    || line == "loop_"
                    || line.starts_with('_')
                    || line.starts_with('#')
                    || line.starts_with("data_")
                {
                    break;
                }

                rows.push(tokenize(lines[j]));
                j += 1;
            }
            let columns = columns
                .iter()
                .map(|c| c.trim_start_matches(prefix).to_string())
                .collect();
            return Some((columns, rows));
        }

        i = j;
    }

    None
}

/// Read a single-record mmCIF item (`_category.item value`).
pub(crate) fn read_value(lines: &[&str], item: &str) -> Option<String> {
    lines.iter().find_map(|line| {
        let rest = line.trim().strip_prefix(item)?;
        if !rest.starts_with(char::is_whitespace) {
            return None; // guard against a longer item sharing this prefix
        }

        tokenize(rest).into_iter().next()
    })
}

/// Split an mmCIF line into tokens, honouring single/double quoted values.
pub(crate) fn tokenize(line: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut chars = line.chars().peekable();
    while let Some(&c) = chars.peek() {
        if c.is_whitespace() {
            chars.next();
        } else if c == '\'' || c == '"' {
            chars.next();
            let mut value = String::new();
            for d in chars.by_ref() {
                if d == c {
                    break;
                }

                value.push(d);
            }

            tokens.push(value);
        } else {
            let mut value = String::new();
            while let Some(&d) = chars.peek() {
                if d.is_whitespace() {
                    break;
                }

                value.push(d);
                chars.next();
            }

            tokens.push(value);
        }
    }

    tokens
}

/// The index of `column` in a loop's column list.
pub(crate) fn column_of(columns: &[String], column: &str) -> Option<usize> {
    columns.iter().position(|c| c == column)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quoted_values_survive_their_spaces() {
        assert_eq!(
            tokenize("ALA 'L-peptide linking' y ALANINE"),
            vec!["ALA", "L-peptide linking", "y", "ALANINE"]
        );
        // An atom name carrying a prime needs no quoting and must not lose it.
        assert_eq!(tokenize("DA C1' C2'"), vec!["DA", "C1'", "C2'"]);
    }

    #[test]
    fn a_loop_reads_its_columns_and_rows() {
        let text = ["loop_", "_thing.a ", "_thing.b ", "x 1 ", "y 2 ", "# "];
        let (columns, rows) = read_loop(&text, "_thing.").unwrap();

        assert_eq!(columns, vec!["a", "b"]);
        assert_eq!(rows, vec![vec!["x", "1"], vec!["y", "2"]]);
        assert_eq!(column_of(&columns, "b"), Some(1));
        assert!(column_of(&columns, "missing").is_none());
    }

    #[test]
    fn a_loop_of_another_category_is_skipped() {
        let text = [
            "loop_",
            "_other.a ",
            "junk ",
            "# ",
            "loop_",
            "_thing.a ",
            "wanted ",
            "# ",
        ];
        let (_, rows) = read_loop(&text, "_thing.").unwrap();

        assert_eq!(rows, vec![vec!["wanted"]]);
        assert!(read_loop(&text, "_absent.").is_none());
    }
}
