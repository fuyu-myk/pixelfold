//! SD files: the ligand poses a docking program writes, one record per pose.
//!
//! A record holds a ligand without its protein, and its atoms carry no names.
//! A pose is therefore laid into a receptor, and its atoms are named after the
//! reference ligand's by pairing the two bond graphs, so the reference's
//! chemical definitions type it as they would a pose file's own ligand. Only
//! the V2000 layout is read, the one docking programs write for a small
//! molecule, gzipped or not.

use std::collections::HashSet;
use std::io::Read;
use std::path::Path;

use anyhow::{Context, Result, bail};
use glam::Vec3;

use crate::components::BondOrder;
use crate::fixed_str::FixedStr;
use crate::interactions::connectivity;
use crate::rmsd::RmsdReference;
use crate::select::AtomSet;
use crate::structure::{Atom, ChainId, Molecule, Protein, SecondaryStructure};

/// The file extensions an SD file goes by.
const EXTENSIONS: [&str; 3] = ["sdf", "sd", "mol"];

/// The residue name of a ligand that is not the reference's: the PDB's own
/// for an unknown ligand.
const UNKNOWN_LIGAND: &str = "UNL";

/// One record of an SD file: its heavy atoms and the bonds between them.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    elements: Vec<FixedStr<2>>,
    positions: Vec<Vec3>,
    /// Bonds between heavy atoms, as indices into `elements`, each with its
    /// order, or `None` for a bond type that names none (a query bond).
    bonds: Vec<(usize, usize, Option<BondOrder>)>,
}

/// A record's ligand laid into a receptor.
pub struct Docked {
    pub protein: Protein,
    /// The ligand's atoms in `protein`.
    pub ligand: AtomSet,
    pub naming: Naming,
}

/// How a docked ligand's atoms were named.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Naming {
    /// After the reference ligand's atoms, bond orders and all where the
    /// reference's definitions give them.
    Named,
    /// After the reference ligand's atoms by their bonds alone: no pairing
    /// keeps the bond orders, as for another tautomer, or another molecule of
    /// the same skeleton.
    BondOrdersDiffer,
    /// As residue `UNL`, its atoms by element and number: the ligand pairs
    /// with the reference's in no way.
    Unpaired,
}

/// Whether `path` names an SD file, by its extension, gzipped (`.sdf.gz`) or
/// not.
pub fn is_sdf(path: &Path) -> bool {
    let named = |path: &Path| {
        path.extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| {
                EXTENSIONS
                    .iter()
                    .any(|known| extension.eq_ignore_ascii_case(known))
            })
    };

    if is_gzipped(path) {
        path.file_stem().is_some_and(|stem| named(Path::new(stem)))
    } else {
        named(path)
    }
}

fn is_gzipped(path: &Path) -> bool {
    path.extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
}

/// Every record of the SD file at `path`.
pub fn read(path: &Path) -> Result<Vec<Record>> {
    contents(path)
        .and_then(|text| parse(&text))
        .with_context(|| format!("failed to read {}", path.display()))
}

/// The text of the file at `path`, gunzipped when it is gzipped. A byte that
/// is not UTF-8, as a property value may hold, is replaced rather than
/// refused, since only the atoms and bonds are read.
fn contents(path: &Path) -> Result<String> {
    let bytes = std::fs::read(path)?;
    let bytes = if is_gzipped(path) {
        let mut text = Vec::new();
        flate2::read::GzDecoder::new(bytes.as_slice()).read_to_end(&mut text)?;
        text
    } else {
        bytes
    };

    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

/// Every record of SD file text, in file order.
pub fn parse(text: &str) -> Result<Vec<Record>> {
    let lines: Vec<&str> = text.lines().collect();
    let mut records = Vec::new();
    let mut start = 0;
    while start < lines.len() {
        let end = lines[start..]
            .iter()
            .position(|line| line.trim_end() == "$$$$")
            .map_or(lines.len(), |at| start + at);
        if lines[start..end].iter().any(|line| !line.trim().is_empty()) {
            let record = parse_record(&lines[start..end], start)
                .with_context(|| format!("record {}", records.len() + 1))?;
            records.push(record);
        }
        start = end + 1;
    }
    if records.is_empty() {
        bail!("no records");
    }

    Ok(records)
}

/// One record from its `lines`, the first of which is the file's line
/// `offset + 1`.
fn parse_record(lines: &[&str], offset: usize) -> Result<Record> {
    let line = |at: usize| {
        lines
            .get(at)
            .copied()
            .with_context(|| format!("line {}: the record ends early", offset + at + 1))
    };
    let counts = line(3)?;
    if counts.contains("V3000") {
        bail!(
            "line {}: V3000 records are not read; write the file as V2000",
            offset + 4
        );
    }
    let count = |start: usize| {
        columns(counts, start, 3)
            .parse::<usize>()
            .with_context(|| format!("line {}: unreadable atom and bond counts", offset + 4))
    };
    let (atoms, bonds) = (count(0)?, count(3)?);

    // Each atom's index among the heavy atoms, `None` for a hydrogen.
    let mut heavy = Vec::with_capacity(atoms);
    let mut elements = Vec::new();
    let mut positions = Vec::new();
    for at in 4..4 + atoms {
        let text = line(at)?;
        let coordinate = |start: usize| columns(text, start, 10).parse::<f32>().ok();
        let symbol = columns(text, 31, 3);
        let (Some(x), Some(y), Some(z)) = (coordinate(0), coordinate(10), coordinate(20)) else {
            bail!("line {}: unreadable atom", offset + at + 1);
        };
        if symbol.is_empty() {
            bail!("line {}: an atom without an element", offset + at + 1);
        }
        if matches!(symbol, "H" | "D" | "T") {
            heavy.push(None);
            continue;
        }
        heavy.push(Some(elements.len()));
        elements.push(FixedStr::new(symbol));
        positions.push(Vec3::new(x, y, z));
    }
    if elements.is_empty() {
        bail!("line {}: no heavy atoms", offset + 4);
    }

    let mut heavy_bonds = Vec::with_capacity(bonds);
    for at in 4 + atoms..4 + atoms + bonds {
        let text = line(at)?;
        let end = |start: usize| {
            columns(text, start, 3)
                .parse::<usize>()
                .ok()
                .and_then(|number| heavy.get(number.checked_sub(1)?).copied())
        };
        let (Some(first), Some(second)) = (end(0), end(3)) else {
            bail!(
                "line {}: a bond to an atom the record lacks",
                offset + at + 1
            );
        };
        let order = match columns(text, 6, 3) {
            "1" => Some(BondOrder::Single),
            "2" => Some(BondOrder::Double),
            "3" => Some(BondOrder::Triple),
            "4" => Some(BondOrder::Aromatic),
            _ => None,
        };
        if let (Some(first), Some(second)) = (first, second) {
            heavy_bonds.push((first, second, order));
        }
    }

    Ok(Record {
        elements,
        positions,
        bonds: heavy_bonds,
    })
}

/// The `width` characters of `line` from column `start`, trimmed, or as many
/// as the line holds.
fn columns(line: &str, start: usize, width: usize) -> &str {
    line.get(start..)
        .map_or("", |rest| rest.get(..width).unwrap_or(rest))
        .trim()
}

impl Record {
    /// This record's ligand added to `receptor` in a chain of its own, its
    /// atoms named after the `reference` ligand atoms they pair with by bond
    /// graph, in the reference's order.
    ///
    /// The pairing is the nearest in place among those the bonds, and their
    /// orders where both sides give them, allow, with the record's coordinates
    /// taken as they are rather than laid on the reference's pocket: any such
    /// pairing names each atom after one of the same chemistry, and the RMSD
    /// pairs the two again once the pose is superposed. Without a pairing
    /// that keeps the bond orders, the bonds alone pair the atoms, since a
    /// docking program may write another tautomer; without one by the bond
    /// table, the bonds the record's distances imply, as the reference's own
    /// are read, which differ for a metal the distances bond to nothing.
    pub fn docked(&self, receptor: &Protein, reference: &RmsdReference) -> Docked {
        let skeleton: Vec<(usize, usize)> = self.bonds.iter().map(|&(a, b, _)| (a, b)).collect();
        let pair = |bonds: &[(usize, usize)], sums: Option<&[u8]>| {
            reference.counterparts(&self.elements, &self.positions, bonds, sums)
        };
        let sums = reference
            .has_bond_orders()
            .then(|| self.order_sums())
            .flatten();
        let (counterparts, naming) =
            if let Some(pairs) = sums.as_deref().and_then(|sums| pair(&skeleton, Some(sums))) {
                (Some(pairs), Naming::Named)
            } else if let Some(pairs) = pair(&skeleton, None) {
                let naming = match sums {
                    Some(_) => Naming::BondOrdersDiffer,
                    None => Naming::Named,
                };
                (Some(pairs), naming)
            } else {
                let pairs = pair(&self.implied_bonds(), None);
                let naming = match pairs {
                    Some(_) => Naming::Named,
                    None => Naming::Unpaired,
                };
                (pairs, naming)
            };

        // Each ligand atom as its residue, number, insertion code and name,
        // with the record atom it is.
        let identities: Vec<(FixedStr<6>, i32, Option<char>, FixedStr<4>, usize)> =
            match counterparts {
                Some(pairs) => pairs
                    .into_iter()
                    .map(|(((_, seq, icode), residue, name), at)| (residue, seq, icode, name, at))
                    .collect(),
                None => (0..self.elements.len())
                    .map(|at| {
                        let name = FixedStr::new(&format!("{}{}", self.elements[at], at + 1));
                        (FixedStr::new(UNKNOWN_LIGAND), 1, None, name, at)
                    })
                    .collect(),
            };

        let chain = free_chain(receptor);
        let first_serial = receptor
            .atoms
            .iter()
            .map(|atom| atom.serial)
            .max()
            .unwrap_or(0)
            + 1;
        let model_number = receptor.atoms.first().map_or(1, |atom| atom.model_number);
        let start = receptor.atoms.len();
        let mut protein = receptor.clone();
        protein.atoms.extend(identities.into_iter().enumerate().map(
            |(order, (residue_name, residue_seq, insertion_code, name, at))| Atom {
                serial: first_serial + order as u32,
                name,
                residue_name,
                residue_seq,
                chain_id: chain,
                insertion_code,
                model_number,
                ..self.atom(at)
            },
        ));
        let mut ligand = AtomSet::empty(protein.atoms.len());
        for atom in start..protein.atoms.len() {
            ligand.insert(atom);
        }

        Docked {
            protein,
            ligand,
            naming,
        }
    }

    /// Each heavy atom's bond orders summed, in half bonds, or `None` when a
    /// bond names no order.
    fn order_sums(&self) -> Option<Vec<u8>> {
        let mut sums = vec![0u8; self.elements.len()];
        for &(a, b, order) in &self.bonds {
            let halves = order?.halves();
            sums[a] = sums[a].saturating_add(halves);
            sums[b] = sums[b].saturating_add(halves);
        }

        Some(sums)
    }

    /// The bonds the heavy atoms' distances imply.
    fn implied_bonds(&self) -> Vec<(usize, usize)> {
        let atoms: Vec<Atom> = (0..self.elements.len()).map(|at| self.atom(at)).collect();
        let mut bonds = Vec::new();
        for (a, first) in atoms.iter().enumerate() {
            for (b, second) in atoms.iter().enumerate().skip(a + 1) {
                if connectivity::at_bonding_distance(first, second) {
                    bonds.push((a, b));
                }
            }
        }

        bonds
    }

    /// Heavy atom `at` as an atom of an unknown ligand, in no chain.
    fn atom(&self, at: usize) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(""),
            element: self.elements[at],
            residue_name: FixedStr::new(UNKNOWN_LIGAND),
            residue_seq: 1,
            chain_id: ChainId::new(""),
            is_hetatm: true,
            molecule: Molecule::NonPolymer,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: self.positions[at],
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }
}

/// A chain id `protein` does not use: `L`, else `L1`, `L2` and on.
fn free_chain(protein: &Protein) -> ChainId {
    let used: HashSet<ChainId> = protein.atoms.iter().map(|atom| atom.chain_id).collect();
    std::iter::once(ChainId::new("L"))
        .chain((1..).map(|number| ChainId::new(&format!("L{number}"))))
        .find(|chain| !used.contains(chain))
        .expect("a structure uses finitely many chain ids")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::components::Dictionary;
    use crate::select;

    /// Acetate as a docking program writes it: atoms in an order of its own,
    /// a hydrogen on the methyl carbon, and a score after the bonds.
    const ACETATE: &str = "\
acetate
  docked

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.6200    1.0700    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6200   -1.0700    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9000    1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  2  4  1  0
  4  5  1  0
M  END
> <minimizedAffinity>
-5.2

$$$$
";

    /// Ethanol, a molecule other than acetate.
    const ETHANOL: &str = "\
ethanol
  docked

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5000    0.0000    0.0000 C   0  0
    2.2000    1.2000    0.0000 O   0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
";

    /// Glycolaldehyde's definition.
    const GLYCOLALDEHYDE: &str = "\
data_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
LIG O1 O N
LIG C1 C N
LIG C2 C N
LIG O2 O N
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG O1 C1 DOUB
LIG C1 C2 SING
LIG C2 O2 SING
";

    /// Glycolaldehyde turned end for end, its carbonyl oxygen first.
    const TURNED_ABOUT: &str = "\
turned
  docked

  4  3  0  0  0  0  0  0  0  0999 V2000
    2.2000   -1.2000    0.0000 O   0  0
    1.5000    0.0000    0.0000 C   0  0
    0.0000    0.0000    0.0000 C   0  0
    0.0000    1.2000    0.0000 O   0  0
  1  2  2  0
  2  3  1  0
  3  4  1  0
M  END
";

    /// A ferrous compound, Fe-O-C, as a docking program bonds it.
    const IRON: &str = "\
iron
  docked

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Fe  0  0
    1.9000    0.0000    0.0000 O   0  0
    3.3000    0.0000    0.0000 C   0  0
  1  2  1  0
  2  3  1  0
M  END
";

    fn atom(chain: &str, seq: i32, residue: &str, name: &str, element: &str, at: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new(chain),
            is_hetatm: residue == "LIG",
            molecule: Molecule::Unknown,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: at,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn structure(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    /// The reference made of the ligand `atoms`, defined by `definitions`.
    fn reference_of(atoms: Vec<Atom>, definitions: &str) -> RmsdReference {
        let crystal = Protein {
            components: Dictionary::parse(definitions),
            ..structure(atoms)
        };

        RmsdReference::new(
            &crystal,
            &select("resn LIG", &crystal).expect("a valid query"),
        )
    }

    /// Acetate as the reference holds it, named as its definition names it.
    fn reference() -> RmsdReference {
        reference_of(
            vec![
                atom("A", 900, "LIG", "C2", "C", Vec3::new(-1.5, 0.0, 0.0)),
                atom("A", 900, "LIG", "C1", "C", Vec3::ZERO),
                atom("A", 900, "LIG", "O1", "O", Vec3::new(0.62, 1.07, 0.0)),
                atom("A", 900, "LIG", "O2", "O", Vec3::new(0.62, -1.07, 0.0)),
            ],
            "",
        )
    }

    /// Glycolaldehyde, O1=C1-C2-O2, whose bond graph without its orders reads
    /// the same from either end.
    fn glycolaldehyde() -> RmsdReference {
        reference_of(
            vec![
                atom("A", 900, "LIG", "O1", "O", Vec3::new(0.0, 1.2, 0.0)),
                atom("A", 900, "LIG", "C1", "C", Vec3::ZERO),
                atom("A", 900, "LIG", "C2", "C", Vec3::new(1.5, 0.0, 0.0)),
                atom("A", 900, "LIG", "O2", "O", Vec3::new(2.2, -1.2, 0.0)),
            ],
            GLYCOLALDEHYDE,
        )
    }

    /// One record of glycolaldehyde's skeleton, its carbonyl, when `carbonyl`,
    /// at the end where the reference's hydroxyl sits.
    fn turned_about(carbonyl: bool) -> Record {
        let order = if carbonyl { 2 } else { 1 };
        let text = TURNED_ABOUT.replace("  1  2  2  0", &format!("  1  2  {order}  0"));

        parse(&text).expect("a record").remove(0)
    }

    fn name_at(docked: &Docked, at: Vec3) -> String {
        let atom = docked
            .ligand
            .iter()
            .map(|atom| &docked.protein.atoms[atom])
            .find(|atom| atom.position == at)
            .expect("an atom there");
        format!("{} {}", atom.residue_name, atom.name)
    }

    /// A receptor whose chain L leaves the ligand chain L1.
    fn receptor() -> Protein {
        structure(vec![
            atom("A", 1, "GLY", "CA", "C", Vec3::new(4.0, 0.0, 0.0)),
            atom("L", 2, "GLY", "CA", "C", Vec3::new(0.0, 5.0, 1.0)),
        ])
    }

    #[test]
    fn each_record_is_a_pose_of_its_heavy_atoms_and_their_bonds() {
        let records = parse(&format!("{ACETATE}{ETHANOL}")).expect("two records");

        assert_eq!(records.len(), 2);
        let acetate = &records[0];
        let elements: Vec<&str> = acetate.elements.iter().map(FixedStr::as_str).collect();
        assert_eq!(elements, ["O", "C", "O", "C"]);
        assert_eq!(acetate.positions[3], Vec3::new(-1.5, 0.0, 0.0));
        assert_eq!(
            acetate.bonds,
            [
                (0, 1, Some(BondOrder::Double)),
                (1, 2, Some(BondOrder::Single)),
                (1, 3, Some(BondOrder::Single))
            ]
        );
        assert_eq!(records[1].elements.len(), 3);
    }

    #[test]
    fn an_sd_file_is_known_by_its_extension() {
        for path in ["poses.sdf", "rank1.SDF", "out.sd", "ligand.mol"] {
            assert!(is_sdf(Path::new(path)), "{path}");
        }
        for path in ["pose.cif", "pose.pdb", "sdf"] {
            assert!(!is_sdf(Path::new(path)), "{path}");
        }
    }

    #[test]
    fn a_v3000_record_is_refused() {
        let text = ACETATE.replace("0999 V2000", "0999 V3000");

        let error = parse(&text).expect_err("V3000 is not read");
        assert!(format!("{error:#}").contains("V3000"), "{error:#}");
    }

    #[test]
    fn a_record_cut_short_names_the_line_it_ends_on() {
        let cut: String = ACETATE
            .lines()
            .take(6)
            .map(|line| format!("{line}\n"))
            .collect();

        let error = parse(&cut).expect_err("two atoms are missing");
        assert!(format!("{error:#}").contains("line 7"), "{error:#}");
    }

    #[test]
    fn a_bond_to_an_atom_the_record_lacks_is_refused() {
        let text = ACETATE.replace("  4  5  1  0", "  4  9  1  0");

        let error = parse(&text).expect_err("atom 9 does not exist");
        assert!(format!("{error:#}").contains("line 13"), "{error:#}");
    }

    #[test]
    fn a_docked_ligand_takes_the_reference_names_in_a_chain_of_its_own() {
        let record = &parse(ACETATE).expect("a record")[0];

        let docked = record.docked(&receptor(), &reference());

        assert_eq!(docked.naming, Naming::Named);
        assert_eq!(docked.protein.atoms.len(), 6);
        let ligand: Vec<&Atom> = docked
            .ligand
            .iter()
            .map(|atom| &docked.protein.atoms[atom])
            .collect();
        let names: Vec<&str> = ligand.iter().map(|atom| atom.name.as_str()).collect();
        assert_eq!(names, ["C2", "C1", "O1", "O2"]);
        assert!(ligand.iter().all(|atom| atom.residue_name == "LIG"
            && atom.residue_seq == 900
            && atom.chain_id == "L1"
            && atom.is_hetatm
            && atom.molecule == Molecule::NonPolymer));
        assert_eq!(ligand[2].position, Vec3::new(0.62, 1.07, 0.0));
    }

    #[test]
    fn a_molecule_other_than_the_reference_ligand_keeps_names_of_its_own() {
        let record = &parse(ETHANOL).expect("a record")[0];

        let docked = record.docked(&receptor(), &reference());

        assert_eq!(docked.naming, Naming::Unpaired);
        let names: Vec<String> = docked
            .ligand
            .iter()
            .map(|atom| {
                let atom = &docked.protein.atoms[atom];
                format!("{} {}", atom.residue_name, atom.name)
            })
            .collect();
        assert_eq!(names, ["UNL C1", "UNL C2", "UNL O3"]);
    }

    #[test]
    fn a_pose_turned_about_takes_the_names_its_bond_orders_allow() {
        // Placed as it is, the carbonyl oxygen sits on the reference's
        // hydroxyl; the bond orders still name it the carbonyl.
        let docked = turned_about(true).docked(&receptor(), &glycolaldehyde());

        assert_eq!(docked.naming, Naming::Named);
        assert_eq!(name_at(&docked, Vec3::new(2.2, -1.2, 0.0)), "LIG O1");
        assert_eq!(name_at(&docked, Vec3::new(0.0, 1.2, 0.0)), "LIG O2");
    }

    #[test]
    fn a_skeleton_whose_bond_orders_differ_is_named_by_its_bonds_alone() {
        // Ethylene glycol: glycolaldehyde's skeleton, every bond single.
        let docked = turned_about(false).docked(&receptor(), &glycolaldehyde());

        assert_eq!(docked.naming, Naming::BondOrdersDiffer);
        assert_eq!(name_at(&docked, Vec3::new(2.2, -1.2, 0.0)), "LIG O2");
    }

    #[test]
    fn a_bond_the_reference_geometry_does_not_draw_still_pairs() {
        // No covalent radius is known for iron, so the reference's geometry
        // bonds it to nothing while the record's bond table bonds it to O1.
        let reference = reference_of(
            vec![
                atom("A", 900, "LIG", "FE", "Fe", Vec3::ZERO),
                atom("A", 900, "LIG", "O1", "O", Vec3::new(1.9, 0.0, 0.0)),
                atom("A", 900, "LIG", "C1", "C", Vec3::new(3.3, 0.0, 0.0)),
            ],
            "",
        );
        let record = parse(IRON).expect("a record").remove(0);

        let docked = record.docked(&receptor(), &reference);

        assert_eq!(docked.naming, Naming::Named);
        assert_eq!(name_at(&docked, Vec3::ZERO), "LIG FE");
    }

    #[test]
    fn a_gzipped_file_with_a_byte_that_is_not_utf8_is_read() {
        use std::io::Write;

        let mut text = ACETATE.replace("-5.2", "-5.2 ~").into_bytes();
        let stray = text
            .iter()
            .position(|&byte| byte == b'~')
            .expect("the marker");
        text[stray] = 0xB5;
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(&text).unwrap();
        let path = std::env::temp_dir().join("pf_sdf_stray_byte.sdf.gz");
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();

        assert!(is_sdf(&path));
        let records = read(&path);
        let _ = std::fs::remove_file(&path);
        assert_eq!(records.expect("one record").len(), 1);
    }
}
