//! How an analysis reaches the terminal, or a pipe.
//!
//! Every command produces a list of records and hands it here. A record knows
//! its own columns, so the three formats are written. This includes an aligned
//! table to read, tab-separated values to pipe, and JSON to parse.

use std::io::Write;

use anyhow::Result;
use clap::ValueEnum;
use serde::Serialize;

/// How to write a command's records.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum Format {
    /// Aligned columns, for reading (default)
    #[default]
    Table,
    /// Tab-separated values, for piping into cut, awk, or a dataframe
    Tsv,
    /// JSON array of objects
    Json,
}

/// A record that can be written as a row of cells.
///
/// The column names are fixed per record type rather than per value, so an empty
/// result still prints its header and a consumer can tell "nothing found" from
/// "wrong command".
pub trait Row {
    fn columns() -> &'static [&'static str];
    fn cells(&self) -> Vec<String>;
}

/// Write `records` to `out` in `format`.
pub fn write<T, W>(out: &mut W, records: &[T], format: Format) -> Result<()>
where
    T: Row + Serialize,
    W: Write,
{
    match format {
        Format::Json => {
            serde_json::to_writer_pretty(&mut *out, records)?;
            writeln!(out)?;
        }
        Format::Tsv => {
            writeln!(out, "{}", T::columns().join("\t"))?;
            for record in records {
                writeln!(out, "{}", record.cells().join("\t"))?;
            }
        }
        Format::Table => write_table(out, records)?,
    }

    Ok(())
}

/// Write the records as columns padded to their widest cell.
///
/// Numeric-looking columns are right-aligned and empty columns
/// are left out. Either way, the tables are always formatted
/// appropriately.
fn write_table<T, W>(out: &mut W, records: &[T]) -> Result<()>
where
    T: Row,
    W: Write,
{
    let all = T::columns();
    let rows: Vec<Vec<String>> = records.iter().map(Row::cells).collect();

    let keep: Vec<usize> = (0..all.len())
        .filter(|&i| {
            rows.is_empty()
                || rows
                    .iter()
                    .any(|row| row.get(i).is_some_and(|cell| !cell.is_empty()))
        })
        .collect();
    let columns: Vec<&str> = keep.iter().map(|&i| all[i]).collect();
    let rows: Vec<Vec<String>> = rows
        .iter()
        .map(|row| {
            keep.iter()
                .map(|&i| row.get(i).cloned().unwrap_or_default())
                .collect()
        })
        .collect();

    let widths: Vec<usize> = columns
        .iter()
        .enumerate()
        .map(|(i, name)| {
            rows.iter()
                .filter_map(|row| row.get(i))
                .map(|cell| cell.chars().count())
                .chain(std::iter::once(name.chars().count()))
                .max()
                .unwrap_or(0)
        })
        .collect();
    let numeric: Vec<bool> = (0..columns.len())
        .map(|i| {
            rows.iter()
                .filter_map(|row| row.get(i))
                .all(|cell| cell.parse::<f64>().is_ok())
        })
        .collect();

    let line = |cells: &[String], out: &mut W| -> Result<()> {
        let mut rendered = Vec::with_capacity(cells.len());
        for (i, cell) in cells.iter().enumerate() {
            let width = widths[i];
            rendered.push(if numeric[i] {
                format!("{cell:>width$}")
            } else {
                format!("{cell:<width$}")
            });
        }
        writeln!(out, "{}", rendered.join("  ").trim_end())?;

        Ok(())
    };

    line(
        &columns.iter().map(|c| (*c).to_string()).collect::<Vec<_>>(),
        out,
    )?;
    writeln!(
        out,
        "{}",
        widths
            .iter()
            .map(|&w| "-".repeat(w))
            .collect::<Vec<_>>()
            .join("  ")
    )?;
    for row in &rows {
        line(row, out)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Serialize)]
    struct Sample {
        chain: String,
        distance: f32,
    }

    impl Row for Sample {
        fn columns() -> &'static [&'static str] {
            &["chain", "distance"]
        }

        fn cells(&self) -> Vec<String> {
            vec![self.chain.clone(), format!("{:.2}", self.distance)]
        }
    }

    fn rendered(records: &[Sample], format: Format) -> String {
        let mut out = Vec::new();
        write(&mut out, records, format).expect("writing to a vec cannot fail");

        String::from_utf8(out).expect("output is utf8")
    }

    fn sample() -> Vec<Sample> {
        vec![
            Sample {
                chain: "A".into(),
                distance: 2.5,
            },
            Sample {
                chain: "LONGER".into(),
                distance: 12.75,
            },
        ]
    }

    #[test]
    fn tsv_is_one_header_and_one_line_per_record() {
        let text = rendered(&sample(), Format::Tsv);
        let lines: Vec<&str> = text.lines().collect();

        assert_eq!(lines[0], "chain\tdistance");
        assert_eq!(lines[1], "A\t2.50");
        assert_eq!(lines[2], "LONGER\t12.75");
    }

    #[test]
    fn a_table_pads_text_left_and_numbers_right() {
        let text = rendered(&sample(), Format::Table);
        let lines: Vec<&str> = text.lines().collect();

        assert_eq!(lines[0], "chain   distance");
        assert_eq!(lines[1], "------  --------");
        // The chain pads out to the width of longer; the distances line up on
        // their last digit.
        assert_eq!(lines[2], "A           2.50");
        assert_eq!(lines[3], "LONGER     12.75");
    }

    #[test]
    fn json_is_an_array_of_objects() {
        let text = rendered(&sample(), Format::Json);
        let parsed: serde_json::Value = serde_json::from_str(&text).expect("valid json");

        assert_eq!(parsed[0]["chain"], "A");
        assert_eq!(parsed[1]["distance"], 12.75);
    }

    /// An empty column should not affect machine outputs.
    #[test]
    fn the_table_drops_an_empty_column_and_the_machine_formats_keep_it() {
        #[derive(Serialize)]
        struct Sparse {
            name: String,
            angle: String,
        }

        impl Row for Sparse {
            fn columns() -> &'static [&'static str] {
                &["name", "angle"]
            }

            fn cells(&self) -> Vec<String> {
                vec![self.name.clone(), self.angle.clone()]
            }
        }

        let records = vec![Sparse {
            name: "one".into(),
            angle: String::new(),
        }];
        let render = |format| {
            let mut out = Vec::new();
            write(&mut out, &records, format).expect("writing to a vec cannot fail");
            String::from_utf8(out).expect("utf8")
        };

        assert_eq!(render(Format::Table).lines().next(), Some("name"));
        assert_eq!(render(Format::Tsv).lines().next(), Some("name\tangle"));
        assert!(render(Format::Json).contains("\"angle\""));
    }

    /// An empty result still says what it would have contained, so a consumer
    /// can tell nothing-found from wrong-command.
    #[test]
    fn an_empty_result_still_carries_its_columns() {
        let empty: Vec<Sample> = Vec::new();

        assert_eq!(rendered(&empty, Format::Tsv).trim(), "chain\tdistance");
        assert!(rendered(&empty, Format::Table).starts_with("chain  distance"));
        assert_eq!(rendered(&empty, Format::Json).trim(), "[]");
    }
}
