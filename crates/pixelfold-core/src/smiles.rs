//! SMILES strings, the way a co-folding tool is given a ligand that has no
//! code in the Chemical Component Dictionary. AlphaFold 3 writes the string it
//! was given into its model file (`_chem_comp.pdbx_smiles`), so the string
//! says what the ligand's atoms are bonded to, by what order of bond, with
//! what charge and how many hydrogens, and which way its stereocentres turn
//! and its double bonds lie, where the pose itself says none of it.
//!
//! The reader follows OpenSMILES: the organic subset outside brackets with
//! its hydrogens implied by the lowest valence that fits (an aromatic atom
//! keeping one bond for the ring's pi system), bracket atoms with their
//! hydrogens, charges and `@`/`@@`, ring closures, branches, `/` and `\` bond
//! directions, and `.` between pieces. Hydrogens written as atoms become
//! counts on the atoms they hang from. Aromatic bonds stay aromatic, as a
//! definition file may write them.

use std::collections::HashMap;
use std::fmt;

use crate::components::BondOrder;
use crate::fixed_str::FixedStr;

/// The most charge, either way, and the most hydrogens a bracket atom may
/// give: OpenSMILES' bound on charge, and past any real atom's hydrogens.
const MOST: i8 = 15;

/// Which way a stereocentre's later neighbours turn, looking from its first.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Turn {
    /// `@`: anticlockwise.
    Anticlockwise,
    /// `@@`: clockwise.
    Clockwise,
}

/// A neighbour of a stereocentre, in the order the string gives them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Neighbour {
    Atom(usize),
    /// The atom's hydrogen, or its lone pair for a centre of three bonds.
    Hydrogen,
}

/// One heavy atom of a SMILES string.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmilesAtom {
    /// The element, in capitals.
    pub element: FixedStr<2>,
    pub aromatic: bool,
    pub charge: i8,
    pub hydrogens: u8,
    pub turn: Option<Turn>,
    /// The atom's neighbours in the order the string gives them, which its
    /// `turn` is read in.
    pub order: Vec<Neighbour>,
    /// Whether a `/` or `\\` on the bond to a hydrogen written as an atom
    /// puts that hydrogen above a double bond from this atom, which sets its
    /// other substituent below.
    pub hydrogen_above: Option<bool>,
}

/// Which side of a double bond a single bond's far atom sits, as `/` and `\`
/// write it: the bond from `from` to `to` goes up or down.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Direction {
    pub from: usize,
    pub to: usize,
    pub up: bool,
}

/// One bond between heavy atoms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmilesBond {
    pub atoms: (usize, usize),
    pub order: BondOrder,
    pub direction: Option<Direction>,
}

/// A molecule read from a SMILES string.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Smiles {
    pub atoms: Vec<SmilesAtom>,
    pub bonds: Vec<SmilesBond>,
}

/// Why a string is not read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmilesError {
    pub at: usize,
    pub reason: String,
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "at character {}: {}", self.at + 1, self.reason)
    }
}

impl std::error::Error for SmilesError {}

/// An atom as the parser keeps it, hydrogens written as atoms included.
struct Parsed {
    element: FixedStr<2>,
    aromatic: bool,
    charge: i8,
    /// The hydrogens a bracket gives, or `None` for an organic-subset atom,
    /// whose hydrogens are implied.
    hydrogens: Option<u8>,
    turn: Option<Turn>,
    order: Vec<Slot>,
    /// Whether an atom comes before this one in its piece, which puts a lone
    /// pair or a hydrogen second in its order rather than first.
    preceded: bool,
}

#[derive(Clone, Copy, Debug)]
enum Slot {
    Atom(usize),
    Hydrogen,
    /// A ring bond opened here and not yet closed.
    Open(u32),
}

/// A bond symbol as written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Symbol {
    Order(BondOrder),
    Up,
    Down,
}

/// The molecule `text` writes.
pub fn parse(text: &str) -> Result<Smiles, SmilesError> {
    let chars: Vec<char> = text.trim().chars().collect();
    let fail = |at: usize, reason: &str| SmilesError {
        at,
        reason: reason.to_string(),
    };
    let mut atoms: Vec<Parsed> = Vec::new();
    // Each bond: its atoms, its symbol if written, as written from first to
    // second.
    let mut bonds: Vec<(usize, usize, Option<Symbol>)> = Vec::new();
    let mut branches: Vec<usize> = Vec::new();
    let mut rings: HashMap<u32, (usize, Option<Symbol>)> = HashMap::new();
    let mut previous: Option<usize> = None;
    let mut symbol: Option<Symbol> = None;
    let mut broken = false;
    // Whether a branch has opened with no atom in it yet.
    let mut branch_opened = false;
    let mut at = 0;

    while at < chars.len() {
        let c = chars[at];
        match c {
            '(' => {
                let Some(from) = previous else {
                    return Err(fail(at, "a branch with nothing before it"));
                };
                branches.push(from);
                branch_opened = true;
                at += 1;
            }
            ')' => {
                previous = Some(
                    branches
                        .pop()
                        .ok_or_else(|| fail(at, "a `)` with no `(`"))?,
                );
                symbol = None;
                at += 1;
            }
            '-' => (symbol, at) = (Some(Symbol::Order(BondOrder::Single)), at + 1),
            '=' => (symbol, at) = (Some(Symbol::Order(BondOrder::Double)), at + 1),
            '#' => (symbol, at) = (Some(Symbol::Order(BondOrder::Triple)), at + 1),
            ':' => (symbol, at) = (Some(Symbol::Order(BondOrder::Aromatic)), at + 1),
            '/' => (symbol, at) = (Some(Symbol::Up), at + 1),
            '\\' => (symbol, at) = (Some(Symbol::Down), at + 1),
            '$' => return Err(fail(at, "quadruple bonds are not read")),
            '.' => {
                broken = true;
                symbol = None;
                at += 1;
            }
            '%' | '0'..='9' => {
                let Some(here) = previous.filter(|_| !branch_opened) else {
                    return Err(fail(at, "a ring bond with no atom before it"));
                };
                let (number, next) = if c == '%' {
                    let digits: String = chars.iter().skip(at + 1).take(2).collect();
                    let number = digits
                        .parse()
                        .ok()
                        .filter(|_| digits.len() == 2 && digits.chars().all(|c| c.is_ascii_digit()))
                        .ok_or_else(|| fail(at, "`%` takes two digits"))?;
                    (number, at + 3)
                } else {
                    (c.to_digit(10).expect("a digit"), at + 1)
                };
                match rings.remove(&number) {
                    Some((other, opened)) => {
                        if other == here {
                            return Err(fail(at, "a ring bond from an atom to itself"));
                        }
                        let slot = atoms[other]
                            .order
                            .iter_mut()
                            .find(|slot| matches!(slot, Slot::Open(open) if *open == number))
                            .expect("an open ring bond keeps its slot");
                        *slot = Slot::Atom(here);
                        atoms[here].order.push(Slot::Atom(other));
                        // A direction written where the ring bond opens reads
                        // from the opening atom, one written where it closes
                        // from the closing atom.
                        match (opened, symbol) {
                            (Some(written), _) => bonds.push((other, here, Some(written))),
                            (None, written) => bonds.push((here, other, written)),
                        }
                    }
                    None => {
                        atoms[here].order.push(Slot::Open(number));
                        rings.insert(number, (here, symbol));
                    }
                }
                symbol = None;
                at = next;
            }
            _ => {
                let (atom, next) = read_atom(&chars, at)?;
                let index = atoms.len();
                atoms.push(atom);
                if let Some(from) = previous.filter(|_| !broken) {
                    atoms[index].preceded = true;
                    atoms[index].order.push(Slot::Atom(from));
                    atoms[from].order.push(Slot::Atom(index));
                    bonds.push((from, index, symbol));
                }
                if atoms[index].hydrogens.is_some_and(|count| count > 0) {
                    atoms[index].order.push(Slot::Hydrogen);
                }
                previous = Some(index);
                symbol = None;
                broken = false;
                branch_opened = false;
                at = next;
            }
        }
    }
    if let Some((_, _)) = rings.values().next() {
        let open: Vec<String> = rings.keys().map(ToString::to_string).collect();
        return Err(fail(
            chars.len(),
            &format!("ring bond {} left open", open.join(", ")),
        ));
    }
    if !branches.is_empty() {
        return Err(fail(chars.len(), "a `(` with no `)`"));
    }
    if atoms.is_empty() {
        return Err(fail(0, "no atoms"));
    }

    Ok(assemble(atoms, bonds))
}

/// One atom from `at`, and where the text after it starts.
fn read_atom(chars: &[char], at: usize) -> Result<(Parsed, usize), SmilesError> {
    let fail = |reason: &str| SmilesError {
        at,
        reason: reason.to_string(),
    };
    if chars[at] != '[' {
        let two: String = chars.iter().skip(at).take(2).collect();
        let (symbol, length) = match two.as_str() {
            "Cl" | "Br" => (two.clone(), 2),
            _ => (chars[at].to_string(), 1),
        };
        if !matches!(
            symbol.as_str(),
            "B" | "C"
                | "N"
                | "O"
                | "P"
                | "S"
                | "F"
                | "Cl"
                | "Br"
                | "I"
                | "b"
                | "c"
                | "n"
                | "o"
                | "p"
                | "s"
        ) {
            return Err(fail(&format!("`{symbol}` is no atom outside brackets")));
        }
        let aromatic = symbol.chars().all(|c| c.is_ascii_lowercase());
        return Ok((
            Parsed {
                element: FixedStr::uppercased(&symbol),
                aromatic,
                charge: 0,
                hydrogens: None,
                turn: None,
                order: Vec::new(),
                preceded: false,
            },
            at + length,
        ));
    }

    let end = chars[at..]
        .iter()
        .position(|&c| c == ']')
        .map(|offset| at + offset)
        .ok_or_else(|| fail("a `[` with no `]`"))?;
    let inside: Vec<char> = chars[at + 1..end].to_vec();
    let mut i = 0;
    while i < inside.len() && inside[i].is_ascii_digit() {
        i += 1; // an isotope, which changes nothing here
    }
    let start = i;
    if i < inside.len() && inside[i].is_ascii_alphabetic() {
        i += 1;
        let first = inside[start];
        let pair: String = inside[start..(start + 2).min(inside.len())]
            .iter()
            .collect();
        let two_letters = inside.get(i).is_some_and(|c| c.is_ascii_lowercase())
            && (first.is_ascii_uppercase() || matches!(pair.as_str(), "se" | "as"));
        if two_letters {
            i += 1;
        }
    }
    let symbol: String = inside[start..i].iter().collect();
    if symbol.is_empty() || symbol == "*" {
        return Err(fail("an atom without an element"));
    }
    let aromatic = symbol
        .chars()
        .next()
        .is_some_and(|c| c.is_ascii_lowercase());
    let element = FixedStr::uppercased(&symbol);
    if crate::molecule::Element::of(element.as_str()).number == 0 && element.as_str() != "H" {
        return Err(fail(&format!("`{symbol}` is no element read here")));
    }

    let mut turn = None;
    if inside.get(i) == Some(&'@') {
        i += 1;
        turn = Some(Turn::Anticlockwise);
        if inside.get(i) == Some(&'@') {
            i += 1;
            turn = Some(Turn::Clockwise);
        } else if inside[i..].starts_with(&['T', 'H']) {
            turn = match inside.get(i + 2) {
                Some('1') => Some(Turn::Anticlockwise),
                Some('2') => Some(Turn::Clockwise),
                _ => return Err(fail("an unknown tetrahedral class")),
            };
            i += 3;
        } else if [['S', 'P'], ['T', 'B'], ['O', 'H']]
            .iter()
            .any(|class| inside[i..].starts_with(class))
        {
            // A square planar, bipyramidal or octahedral class: not read.
            turn = None;
            while inside.get(i).is_some_and(|c| c.is_ascii_alphanumeric()) {
                i += 1;
            }
        }
    }
    let mut hydrogens = 0u8;
    if inside.get(i) == Some(&'H') {
        i += 1;
        hydrogens = 1;
        let digits: String = inside[i..]
            .iter()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        if !digits.is_empty() {
            hydrogens = digits
                .parse()
                .ok()
                .filter(|&count: &u8| count <= MOST as u8)
                .ok_or_else(|| fail("a hydrogen count past any atom's"))?;
            i += digits.len();
        }
    }
    let mut charge = 0i8;
    while let Some(&sign) = inside.get(i).filter(|c| matches!(c, '+' | '-')) {
        i += 1;
        let digits: String = inside[i..]
            .iter()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        let size: i8 = if digits.is_empty() {
            1
        } else {
            i += digits.len();
            digits.parse().map_err(|_| fail("an unreadable charge"))?
        };
        charge = charge
            .checked_add(if sign == '+' { size } else { -size })
            .filter(|charge| charge.abs() <= MOST)
            .ok_or_else(|| fail("a charge past fifteen"))?;
    }
    if inside.get(i) == Some(&':') {
        i = inside.len(); // an atom class, which changes nothing here
    }
    if i != inside.len() {
        return Err(fail("an unreadable bracket atom"));
    }

    Ok((
        Parsed {
            element,
            aromatic,
            charge,
            hydrogens: Some(hydrogens),
            turn,
            order: Vec::new(),
            preceded: false,
        },
        end + 1,
    ))
}

/// The molecule the parsed atoms and bonds make: hydrogens written as atoms
/// folded into counts, implied hydrogens counted, bonds without a symbol made
/// aromatic between aromatic atoms and single otherwise.
fn assemble(atoms: Vec<Parsed>, bonds: Vec<(usize, usize, Option<Symbol>)>) -> Smiles {
    let is_hydrogen = |atom: usize| atoms[atom].element.as_str() == "H";
    let mut index = vec![usize::MAX; atoms.len()];
    let mut next = 0;
    for (atom, slot) in index.iter_mut().enumerate() {
        if !is_hydrogen(atom) {
            *slot = next;
            next += 1;
        }
    }

    let mut written_hydrogens = vec![0u8; atoms.len()];
    let mut hydrogen_above: Vec<Option<bool>> = vec![None; atoms.len()];
    let mut out_bonds = Vec::new();
    // The bonds no symbol was written for between aromatic atoms, which are
    // aromatic only inside a ring.
    let mut unwritten_aromatic = Vec::new();
    for &(a, b, symbol) in &bonds {
        match (is_hydrogen(a), is_hydrogen(b)) {
            (false, false) => {
                if symbol.is_none() && atoms[a].aromatic && atoms[b].aromatic {
                    unwritten_aromatic.push(out_bonds.len());
                }
                let order = match symbol {
                    Some(Symbol::Order(order)) => order,
                    Some(Symbol::Up | Symbol::Down) => BondOrder::Single,
                    None if atoms[a].aromatic && atoms[b].aromatic => BondOrder::Aromatic,
                    None => BondOrder::Single,
                };
                let direction = match symbol {
                    Some(Symbol::Up) => Some(true),
                    Some(Symbol::Down) => Some(false),
                    _ => None,
                }
                .map(|up| Direction {
                    from: index[a],
                    to: index[b],
                    up,
                });
                out_bonds.push(SmilesBond {
                    atoms: (index[a], index[b]),
                    order,
                    direction,
                });
            }
            (true, false) | (false, true) => {
                let (hydrogen, atom) = if is_hydrogen(a) { (a, b) } else { (b, a) };
                written_hydrogens[atom] = written_hydrogens[atom].saturating_add(1);
                // Up from the hydrogen puts it below; up to it, above.
                let up = match symbol {
                    Some(Symbol::Up) => Some(true),
                    Some(Symbol::Down) => Some(false),
                    _ => None,
                };
                if let Some(up) = up {
                    hydrogen_above[atom] = Some(if hydrogen == a { !up } else { up });
                }
            }
            (true, true) => {}
        }
    }
    let mut neighbours = vec![Vec::new(); next];
    for bond in &out_bonds {
        neighbours[bond.atoms.0].push(bond.atoms.1);
        neighbours[bond.atoms.1].push(bond.atoms.0);
    }
    let rings = crate::rings::Rings::of(&neighbours);
    for at in unwritten_aromatic {
        let (a, b) = out_bonds[at].atoms;
        if rings.of_bond(a, b).next().is_none() {
            out_bonds[at].order = BondOrder::Single;
        }
    }

    let mut out_atoms: Vec<SmilesAtom> = Vec::new();
    for (atom, parsed) in atoms.iter().enumerate() {
        if is_hydrogen(atom) {
            continue;
        }
        let valence: u8 = out_bonds
            .iter()
            .filter(|bond| bond.atoms.0 == index[atom] || bond.atoms.1 == index[atom])
            .map(|bond| match bond.order {
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                _ => 1,
            })
            .fold(written_hydrogens[atom], u8::saturating_add);
        let hydrogens = written_hydrogens[atom].saturating_add(
            parsed
                .hydrogens
                .unwrap_or_else(|| implied_hydrogens(parsed.element, parsed.aromatic, valence)),
        );
        let mut order: Vec<Neighbour> = parsed
            .order
            .iter()
            .filter_map(|slot| match *slot {
                Slot::Atom(other) if is_hydrogen(other) => Some(Neighbour::Hydrogen),
                Slot::Atom(other) => Some(Neighbour::Atom(index[other])),
                Slot::Hydrogen => Some(Neighbour::Hydrogen),
                Slot::Open(_) => None,
            })
            .collect();
        // A stereocentre of three bonds turns about its lone pair, which
        // stands where a hydrogen would.
        if parsed.turn.is_some() && order.len() == 3 && !order.contains(&Neighbour::Hydrogen) {
            order.insert(
                usize::from(parsed.preceded).min(order.len()),
                Neighbour::Hydrogen,
            );
        }
        out_atoms.push(SmilesAtom {
            element: parsed.element,
            aromatic: parsed.aromatic,
            charge: parsed.charge,
            hydrogens,
            turn: parsed.turn,
            order,
            hydrogen_above: hydrogen_above[atom],
        });
    }

    Smiles {
        atoms: out_atoms,
        bonds: out_bonds,
    }
}

/// The hydrogens an organic-subset atom implies: up to the lowest valence it
/// allows at or above its bonds' (aromatic bonds counting one), an aromatic
/// carbon, nitrogen, phosphorus or boron keeping one more for the ring.
fn implied_hydrogens(element: FixedStr<2>, aromatic: bool, valence: u8) -> u8 {
    let allowed: &[u8] = match element.as_str() {
        "B" => &[3],
        "C" => &[4],
        "N" => &[3, 5],
        "O" => &[2],
        "P" => &[3, 5],
        "S" => &[2, 4, 6],
        _ => &[1],
    };
    if aromatic {
        // An aromatic atom takes its lowest valence, one bond's worth kept
        // for the ring where it can hold a double bond: a pyrrole's
        // N-methyl nitrogen keeps none, its NH is written in brackets.
        let ring = u8::from(matches!(element.as_str(), "B" | "C" | "N" | "P"));
        return allowed[0].saturating_sub(valence.saturating_add(ring));
    }

    allowed
        .iter()
        .find(|&&allowed| allowed >= valence)
        .map_or(0, |&allowed| allowed - valence)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hydrogens(smiles: &Smiles) -> Vec<u8> {
        smiles.atoms.iter().map(|atom| atom.hydrogens).collect()
    }

    #[test]
    fn the_organic_subset_implies_its_hydrogens() {
        let ethanol = parse("CCO").unwrap();
        assert_eq!(hydrogens(&ethanol), [3, 2, 1]);
        assert!(
            ethanol
                .bonds
                .iter()
                .all(|bond| bond.order == BondOrder::Single)
        );

        let acetate = parse("CC(=O)[O-]").unwrap();
        assert_eq!(hydrogens(&acetate), [3, 0, 0, 0]);
        assert_eq!(acetate.atoms[3].charge, -1);
        assert_eq!(acetate.bonds[1].order, BondOrder::Double);

        let ammonium = parse("[NH4+]").unwrap();
        assert_eq!(
            (ammonium.atoms[0].hydrogens, ammonium.atoms[0].charge),
            (4, 1)
        );
    }

    #[test]
    fn aromatic_rings_close_and_keep_one_bond_for_the_ring() {
        let benzene = parse("c1ccccc1").unwrap();
        assert_eq!(hydrogens(&benzene), [1; 6]);
        assert_eq!(benzene.bonds.len(), 6);
        assert!(
            benzene
                .bonds
                .iter()
                .all(|bond| bond.order == BondOrder::Aromatic)
        );

        assert_eq!(hydrogens(&parse("n1ccccc1").unwrap())[0], 0);
        assert_eq!(hydrogens(&parse("[nH]1cccc1").unwrap())[0], 1);
        assert_eq!(hydrogens(&parse("o1cccc1").unwrap())[0], 0);
        assert_eq!(hydrogens(&parse("Cn1cccc1").unwrap())[1], 0);
        // 2-pyridone: the carbonyl carbon keeps nothing for hydrogens.
        assert_eq!(hydrogens(&parse("O=c1cccc[nH]1").unwrap())[1], 0);
        // Two-digit ring bonds, and a biphenyl's single bond between rings.
        let biphenyl = parse("c%10ccccc%10-c1ccccc1").unwrap();
        assert_eq!(
            biphenyl
                .bonds
                .iter()
                .filter(|bond| bond.order == BondOrder::Single)
                .count(),
            1
        );
    }

    #[test]
    fn a_stereocentre_lists_its_neighbours_in_the_order_written() {
        // L-alanine: from N, the hydrogen, methyl and carboxyl turn clockwise.
        let alanine = parse("N[C@@H](C)C(=O)O").unwrap();
        let centre = &alanine.atoms[1];
        assert_eq!(centre.turn, Some(Turn::Clockwise));
        assert_eq!(
            centre.order,
            [
                Neighbour::Atom(0),
                Neighbour::Hydrogen,
                Neighbour::Atom(2),
                Neighbour::Atom(3)
            ]
        );

        // A ring bond takes its place where its digit is written.
        let ring = parse("C[C@H]1CCCN1").unwrap();
        assert_eq!(
            ring.atoms[1].order,
            [
                Neighbour::Atom(0),
                Neighbour::Hydrogen,
                Neighbour::Atom(5),
                Neighbour::Atom(2)
            ]
        );

        // A hydrogen written as an atom becomes the atom's count, in place.
        let written = parse("[H][C@@](F)(Cl)Br").unwrap();
        assert_eq!(written.atoms.len(), 4);
        assert_eq!(written.atoms[0].hydrogens, 1);
        assert_eq!(written.atoms[0].order[0], Neighbour::Hydrogen);
    }

    #[test]
    fn bond_directions_say_which_side_of_a_double_bond_each_end_sits() {
        let trans = parse("F/C=C/F").unwrap();
        let directions: Vec<Direction> = trans
            .bonds
            .iter()
            .filter_map(|bond| bond.direction)
            .collect();
        assert_eq!(
            directions,
            [
                Direction {
                    from: 0,
                    to: 1,
                    up: true
                },
                Direction {
                    from: 2,
                    to: 3,
                    up: true
                },
            ]
        );
        let cis = parse("C(/F)=C/F").unwrap();
        assert_eq!(
            cis.bonds[0].direction,
            Some(Direction {
                from: 0,
                to: 1,
                up: true
            })
        );
    }

    #[test]
    fn a_string_that_does_not_read_is_refused_with_where() {
        for (text, at) in [
            ("C1CC", 4),
            ("C(C", 3),
            ("CC)", 2),
            ("[Xx]", 0),
            ("C*", 1),
            ("C$C", 1),
        ] {
            let error = parse(text).unwrap_err();
            assert_eq!(error.at, at, "{text}: {error}");
        }
    }

    #[test]
    fn counts_out_of_reach_are_refused_rather_than_wrapped() {
        for text in ["[C+99+99]", "[CH255][H]", "C%+1CC%+1", "C(1CC1)"] {
            assert!(parse(text).is_err(), "{text} read");
        }
    }

    #[test]
    fn an_unwritten_bond_between_aromatic_rings_is_single() {
        let biphenyl = parse("c1ccccc1c1ccccc1").unwrap();
        let between = biphenyl
            .bonds
            .iter()
            .find(|bond| bond.atoms == (5, 6))
            .unwrap();
        assert_eq!(between.order, BondOrder::Single);
        assert_eq!(
            biphenyl
                .bonds
                .iter()
                .filter(|bond| bond.order == BondOrder::Aromatic)
                .count(),
            12
        );
    }

    #[test]
    fn a_direction_on_a_written_hydrogen_is_kept_for_its_atom() {
        let marked = parse("[H]/C(F)=C/F").unwrap();
        // Up from the hydrogen to C: the hydrogen below, so not above.
        assert_eq!(marked.atoms[0].hydrogen_above, Some(false));
        assert_eq!(marked.atoms[1].hydrogen_above, None);
    }

    #[test]
    fn pieces_apart_are_not_bonded() {
        let salt = parse("[Na+].[Cl-]").unwrap();
        assert_eq!(salt.atoms.len(), 2);
        assert!(salt.bonds.is_empty());
    }
}
