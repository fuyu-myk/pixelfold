//! What a SMILES string states of one residue of a pose: its chemistry, as a
//! definition named after the pose's atoms, and whether the pose's
//! stereocentres turn and its double bonds lie as the string says.
//!
//! A pose's atoms are paired with the string's by bond graph, the pose's bonds
//! taken from its own definition where it has one and from distances where it
//! has none, as a co-folding model's ligand from a SMILES input has none. A
//! graph with symmetries pairs in several ways, which can differ in where they
//! put the string's double bonds, so the pairing kept is the one whose bond
//! lengths best fit the string's bond orders, and of those, the one that
//! leaves the fewest stereocentres inverted (a meso compound turned about is
//! the same molecule).

use std::collections::HashMap;

use glam::Vec3;
use petgraph::algo::isomorphism::subgraph_isomorphisms_iter;
use petgraph::graph::{NodeIndex, UnGraph};

use crate::components::{BondOrder, Component, DefinedAtom};
use crate::fixed_str::FixedStr;
use crate::interactions::connectivity;
use crate::numbering::{ResidueKey, residue_key};
use crate::pose::is_hydrogen;
use crate::rings::Rings;
use crate::rmsd::AtomId;
use crate::smiles::{Neighbour, Smiles, Turn};
use crate::stereo::{self, PLANAR, RING_WITH_TRANS, TETRAHEDRAL};
use crate::structure::Protein;

/// The most pairings weighed: a ligand of many equivalent groups pairs in
/// more ways than are worth trying.
const MAX_PAIRINGS: usize = 256;

/// How much further than their covalent radii summed two pose atoms may lie
/// and still stand for a bond of the string, when the pose's own bonds do not
/// pair, in angstroms: a bond stretched by half again still pairs.
const LOOSE_REACH: f32 = 0.9;

/// The most candidates the loose search tries before keeping what it found.
const MAX_STEPS: usize = 2_000_000;

/// The most atoms the loose search pairs; a larger ligand pairs exactly or
/// not at all.
const MAX_LOOSE_ATOMS: usize = 300;

/// A SMILES string's statement about one residue of a pose.
#[derive(Clone, Debug)]
pub struct Stated {
    residue: ResidueKey,
    component: Component,
    inverted: Vec<AtomId>,
    flipped: Vec<[AtomId; 2]>,
    loose: bool,
}

impl Stated {
    /// The residue stated of.
    pub fn residue(&self) -> ResidueKey {
        self.residue
    }

    /// The residue's definition as the string gives it, named after the
    /// pose's atoms.
    pub fn component(&self) -> &Component {
        &self.component
    }

    /// The string's stereocentres the pose turns the other way.
    pub fn inverted(&self) -> &[AtomId] {
        &self.inverted
    }

    /// The string's E or Z double bonds the pose turns the other way.
    pub fn flipped(&self) -> &[[AtomId; 2]] {
        &self.flipped
    }

    /// Whether the pose's atoms paired with the string's only once bonds past
    /// bonding distance were allowed: a bond stretched, or atoms clashing.
    pub fn loose(&self) -> bool {
        self.loose
    }
}

/// One pairing of the string's atoms with the pose's, and what it makes of
/// them.
struct Candidate {
    /// How far the pose's bond lengths stray from those the string's bond
    /// orders give, squared and summed.
    fit: f32,
    pairing: Vec<usize>,
    inverted: Vec<usize>,
    flipped: Vec<usize>,
}

impl Candidate {
    /// Whether this pairing fits the string's bond orders better than
    /// `other`, or as well with fewer centres or bonds turned.
    fn beats(&self, other: &Candidate) -> bool {
        let wrong = |candidate: &Candidate| candidate.inverted.len() + candidate.flipped.len();
        self.fit < other.fit - 1e-3 || (self.fit <= other.fit + 1e-3 && wrong(self) < wrong(other))
    }
}

/// A stereocentre the string states, in the string's atoms: the centre, its
/// neighbours in the string's order (`None` for a hydrogen or lone pair) and
/// their turn.
struct Centre {
    atom: usize,
    order: Vec<Option<usize>>,
    turn: Turn,
}

/// A double bond the string states, in the string's atoms: one end, a
/// substituent on it, the other end and a substituent on that, as
/// `stereo::dihedral_cosine` takes them, and whether the two substituents
/// are held cis.
struct DoubleBond {
    atoms: [usize; 4],
    cis: bool,
}

/// What `smiles` states of the residue of `pose` whose atoms are `atoms`, or
/// `None` when the two bond graphs do not pair.
pub fn state(pose: &Protein, atoms: &[usize], smiles: &Smiles) -> Option<Stated> {
    let heavy: Vec<usize> = atoms
        .iter()
        .copied()
        .filter(|&atom| !is_hydrogen(&pose.atoms[atom]))
        .collect();
    if heavy.len() != smiles.atoms.len() || heavy.is_empty() {
        return None;
    }
    let pose_bonds = bonds_of(pose, &heavy);
    let positions: Vec<Vec3> = heavy
        .iter()
        .map(|&atom| pose.atoms[atom].position)
        .collect();
    let elements: Vec<FixedStr<2>> = heavy
        .iter()
        .map(|&atom| FixedStr::uppercased(pose.atoms[atom].element.as_str()))
        .collect();
    let mut pairings = if pose_bonds.len() == smiles.bonds.len() {
        exact_pairings(smiles, &elements, &pose_bonds)
    } else {
        Vec::new()
    };
    let loose = pairings.is_empty();
    if loose {
        pairings = loose_pairings(smiles, &elements, &positions);
    }

    let (centres, double_bonds) = statements(smiles);
    let mut best: Option<Candidate> = None;
    for pairing in pairings {
        let candidate = Candidate {
            fit: length_fit(smiles, &pairing, &positions),
            inverted: inverted_centres(&centres, &pairing, &positions),
            flipped: flipped_bonds(&double_bonds, &pairing, &positions),
            pairing,
        };
        if best.as_ref().is_none_or(|best| candidate.beats(best)) {
            best = Some(candidate);
        }
    }
    let Candidate {
        pairing,
        inverted,
        flipped,
        ..
    } = best?;

    let id = |at: usize| {
        let atom = &pose.atoms[heavy[at]];
        (residue_key(atom), atom.residue_name, atom.name)
    };
    Some(Stated {
        residue: residue_key(&pose.atoms[heavy[0]]),
        component: definition(pose, &heavy, smiles, &pairing),
        inverted: inverted
            .into_iter()
            .map(|at| id(pairing[centres[at].atom]))
            .collect(),
        flipped: flipped
            .into_iter()
            .map(|at| {
                let [b, _, c, _] = double_bonds[at].atoms;
                [id(pairing[b]), id(pairing[c])]
            })
            .collect(),
        loose,
    })
}

/// The pairings, up to [`MAX_PAIRINGS`], of the string's atoms with the
/// pose's whose bonds match exactly, by VF2.
fn exact_pairings(
    smiles: &Smiles,
    elements: &[FixedStr<2>],
    pose_bonds: &[(usize, usize)],
) -> Vec<Vec<usize>> {
    let mut ours: UnGraph<FixedStr<2>, ()> = UnGraph::new_undirected();
    for atom in &smiles.atoms {
        ours.add_node(atom.element);
    }
    for bond in &smiles.bonds {
        ours.add_edge(
            NodeIndex::new(bond.atoms.0),
            NodeIndex::new(bond.atoms.1),
            (),
        );
    }
    let mut theirs: UnGraph<FixedStr<2>, ()> = UnGraph::new_undirected();
    for &element in elements {
        theirs.add_node(element);
    }
    for &(a, b) in pose_bonds {
        theirs.add_edge(NodeIndex::new(a), NodeIndex::new(b), ());
    }

    let mut node_match = |a: &FixedStr<2>, b: &FixedStr<2>| a == b;
    let mut edge_match = |_: &(), _: &()| true;
    let (ours, theirs) = (&ours, &theirs);
    let Some(found) = subgraph_isomorphisms_iter(&ours, &theirs, &mut node_match, &mut edge_match)
    else {
        return Vec::new();
    };
    // petgraph 0.6.5's size hint indexes past its table at 21 nodes, so the
    // pairings are pushed one by one rather than collected.
    let mut pairings = Vec::new();
    for pairing in found.take(MAX_PAIRINGS) {
        pairings.push(pairing);
    }
    pairings
}

/// The pairings, up to [`MAX_PAIRINGS`] and within [`MAX_STEPS`] candidates
/// tried, of the string's atoms with the pose's that keep elements and put
/// every bond of the string between pose atoms within [`LOOSE_REACH`] of
/// covalent distance, the pose free to hold more such pairs. A pose distorted
/// past reading its bonds from distances, a bond stretched or two atoms
/// clashing, so still pairs, and fails its geometry rather than going
/// unchecked.
fn loose_pairings(
    smiles: &Smiles,
    elements: &[FixedStr<2>],
    positions: &[Vec3],
) -> Vec<Vec<usize>> {
    let count = smiles.atoms.len();
    if count > MAX_LOOSE_ATOMS {
        return Vec::new();
    }
    let mut ours = vec![Vec::new(); count];
    for bond in &smiles.bonds {
        ours[bond.atoms.0].push(bond.atoms.1);
        ours[bond.atoms.1].push(bond.atoms.0);
    }
    let near: Vec<Vec<bool>> = (0..count)
        .map(|a| {
            (0..count)
                .map(|b| {
                    a != b
                        && positions[a].distance(positions[b])
                            < covalent(elements[a]) + covalent(elements[b]) + LOOSE_REACH
                })
                .collect()
        })
        .collect();

    // The string's atoms in breadth-first order, piece by piece, so each but
    // a piece's first has a neighbour paired before it.
    let mut order = Vec::with_capacity(count);
    let mut seen = vec![false; count];
    for start in 0..count {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        order.push(start);
        let mut at = order.len() - 1;
        while at < order.len() {
            for &next in &ours[order[at]] {
                if !seen[next] {
                    seen[next] = true;
                    order.push(next);
                }
            }
            at += 1;
        }
    }

    let mut search = Search {
        smiles,
        ours: &ours,
        elements,
        near: &near,
        order: &order,
        pairing: vec![usize::MAX; count],
        used: vec![false; count],
        found: Vec::new(),
        steps: 0,
    };
    search.extend(0);
    search.found
}

/// A depth-first search for loose pairings.
struct Search<'a> {
    smiles: &'a Smiles,
    ours: &'a [Vec<usize>],
    elements: &'a [FixedStr<2>],
    near: &'a [Vec<bool>],
    order: &'a [usize],
    pairing: Vec<usize>,
    used: Vec<bool>,
    found: Vec<Vec<usize>>,
    steps: usize,
}

impl Search<'_> {
    fn extend(&mut self, depth: usize) {
        if self.found.len() >= MAX_PAIRINGS || self.steps >= MAX_STEPS {
            return;
        }
        let Some(&atom) = self.order.get(depth) else {
            self.found.push(self.pairing.clone());
            return;
        };
        let element = self.smiles.atoms[atom].element;
        for candidate in 0..self.pairing.len() {
            self.steps += 1;
            let fits = !self.used[candidate]
                && self.elements[candidate] == element
                && self.ours[atom].iter().all(|&neighbour| {
                    let paired = self.pairing[neighbour];
                    paired == usize::MAX || self.near[candidate][paired]
                });
            if !fits {
                continue;
            }
            self.pairing[atom] = candidate;
            self.used[candidate] = true;
            self.extend(depth + 1);
            self.used[candidate] = false;
            self.pairing[atom] = usize::MAX;
            if self.found.len() >= MAX_PAIRINGS || self.steps >= MAX_STEPS {
                return;
            }
        }
    }
}

/// The bonds between the `heavy` atoms of `pose`, as indices into them: its
/// definition's where it defines the residue, else those its distances imply.
fn bonds_of(pose: &Protein, heavy: &[usize]) -> Vec<(usize, usize)> {
    let atom = |at: usize| &pose.atoms[heavy[at]];
    let definition = pose
        .components
        .get(atom(0).residue_name.as_str())
        .filter(|component| {
            (0..heavy.len()).all(|at| component.element_of(atom(at).name.as_str()).is_some())
        });
    let mut bonds = Vec::new();
    for a in 0..heavy.len() {
        for b in a + 1..heavy.len() {
            let bonded = match definition {
                Some(component) => component
                    .neighbours(atom(a).name.as_str())
                    .any(|name| name == atom(b).name.as_str()),
                None => connectivity::at_bonding_distance(atom(a), atom(b)),
            };
            if bonded {
                bonds.push((a, b));
            }
        }
    }
    bonds
}

/// The stereocentres and double bonds the string states that its own bond
/// graph lets be told apart: a centre whose heavy neighbours are each of
/// their own class under colour refinement (a phosphate's three terminal
/// oxygens, which the string may mark, are not), not a nitrogen of three
/// bonds, which inverts; a double bond outside any ring of fewer than eight
/// atoms with a direction written on a substituent at each end.
fn statements(smiles: &Smiles) -> (Vec<Centre>, Vec<DoubleBond>) {
    let elements: Vec<FixedStr<2>> = smiles.atoms.iter().map(|atom| atom.element).collect();
    let mut neighbours = vec![Vec::new(); elements.len()];
    for bond in &smiles.bonds {
        neighbours[bond.atoms.0].push(bond.atoms.1);
        neighbours[bond.atoms.1].push(bond.atoms.0);
    }
    let classes = stereo::classes(&elements, &neighbours);

    let centres = smiles
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(at, atom)| {
            let turn = atom.turn?;
            let heavy: Vec<usize> = atom
                .order
                .iter()
                .filter_map(|neighbour| match neighbour {
                    Neighbour::Atom(other) => Some(*other),
                    Neighbour::Hydrogen => None,
                })
                .collect();
            let mut seen: Vec<usize> = heavy.iter().map(|&other| classes[other]).collect();
            seen.sort_unstable();
            seen.dedup();
            let inverts = atom.element.as_str() == "N"
                && neighbours[at].len() + usize::from(atom.hydrogens) == 3;
            (atom.order.len() == 4 && heavy.len() >= 3 && seen.len() == heavy.len() && !inverts)
                .then(|| Centre {
                    atom: at,
                    order: atom
                        .order
                        .iter()
                        .map(|neighbour| match neighbour {
                            Neighbour::Atom(other) => Some(*other),
                            Neighbour::Hydrogen => None,
                        })
                        .collect(),
                    turn,
                })
        })
        .collect();

    let rings = Rings::of(&neighbours);
    // Each atom's substituents with the side of the double bond `/` and `\`
    // put them on: `true` above.
    let mut sides: HashMap<(usize, usize), bool> = HashMap::new();
    for bond in &smiles.bonds {
        if let Some(direction) = bond.direction {
            // Written from the substituent to the double bond's atom, up
            // puts the substituent below; from the atom, above.
            sides.insert((direction.to, direction.from), !direction.up);
            sides.insert((direction.from, direction.to), direction.up);
        }
    }
    let double_bonds = smiles
        .bonds
        .iter()
        .filter(|bond| bond.order == BondOrder::Double)
        .filter(|bond| {
            let (b, c) = bond.atoms;
            !rings
                .of_bond(b, c)
                .any(|ring| rings.all()[ring].len() < RING_WITH_TRANS)
        })
        .filter_map(|bond| {
            let (b, c) = bond.atoms;
            let marked = |end: usize, other: usize| {
                let mut substituents = neighbours[end]
                    .iter()
                    .copied()
                    .filter(|&substituent| substituent != other);
                let written = substituents.clone().find_map(|substituent| {
                    sides.get(&(end, substituent)).map(|&up| (substituent, up))
                });
                // A direction on a hydrogen written as an atom puts the
                // end's other substituent on the other side.
                written.or_else(|| {
                    let above = smiles.atoms[end].hydrogen_above?;
                    let substituent = substituents.next()?;
                    substituents
                        .next()
                        .is_none()
                        .then_some((substituent, !above))
                })
            };
            let ((a, a_up), (d, d_up)) = (marked(b, c)?, marked(c, b)?);
            Some(DoubleBond {
                atoms: [b, a, c, d],
                cis: a_up == d_up,
            })
        })
        .collect();

    (centres, double_bonds)
}

/// How far the pose's bond lengths, under `pairing`, stray from the lengths
/// the string's bond orders give, squared and summed.
fn length_fit(smiles: &Smiles, pairing: &[usize], positions: &[Vec3]) -> f32 {
    smiles
        .bonds
        .iter()
        .map(|bond| {
            let (a, b) = bond.atoms;
            let radii = covalent(smiles.atoms[a].element) + covalent(smiles.atoms[b].element);
            let scale = match bond.order {
                BondOrder::Single => 1.0,
                BondOrder::Aromatic => 0.92,
                BondOrder::Double => 0.87,
                BondOrder::Triple => 0.79,
            };
            let length = positions[pairing[a]].distance(positions[pairing[b]]);
            (length - radii * scale).powi(2)
        })
        .sum()
}

/// A covalent radius for the bond-length fit, in angstroms.
fn covalent(element: FixedStr<2>) -> f32 {
    match element.as_str() {
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "S" => 1.05,
        "P" => 1.07,
        "CL" => 1.02,
        "BR" => 1.20,
        "I" => 1.39,
        "B" => 0.84,
        "SI" => 1.11,
        "SE" => 1.20,
        _ => 0.8,
    }
}

/// The `centres` the pose's atoms, under `pairing`, turn the other way, by
/// index: a centre left flat counts as neither.
fn inverted_centres(centres: &[Centre], pairing: &[usize], positions: &[Vec3]) -> Vec<usize> {
    centres
        .iter()
        .enumerate()
        .filter(|(_, centre)| {
            // Looking from the first neighbour, the other three turn
            // anticlockwise (`@`) when their unit bonds span a negative
            // volume; a hydrogen, moved to the front, flips the turn once for
            // each place it moves.
            let mut order = centre.order.clone();
            let mut anticlockwise = centre.turn == Turn::Anticlockwise;
            if let Some(place) = order.iter().position(Option::is_none) {
                order.remove(place);
                anticlockwise ^= place % 2 == 1;
            } else {
                order.remove(0);
            }
            let [a, b, c] =
                [order[0], order[1], order[2]].map(|at| at.expect("one hydrogen at most"));
            let spanned = stereo::volume([centre.atom, a, b, c].map(|at| positions[pairing[at]]));
            let expected = if anticlockwise { -1.0 } else { 1.0 };
            spanned * expected <= -TETRAHEDRAL
        })
        .map(|(at, _)| at)
        .collect()
}

/// The `bonds` the pose's atoms, under `pairing`, turn the other way, by
/// index: a bond left twisted across counts as neither.
fn flipped_bonds(bonds: &[DoubleBond], pairing: &[usize], positions: &[Vec3]) -> Vec<usize> {
    bonds
        .iter()
        .enumerate()
        .filter(|(_, bond)| {
            let cosine = stereo::dihedral_cosine(bond.atoms.map(|at| positions[pairing[at]]));
            let expected = if bond.cis { 1.0 } else { -1.0 };
            cosine * expected <= -PLANAR
        })
        .map(|(at, _)| at)
        .collect()
}

/// The residue's definition from the string, its atoms named after the pose
/// atoms `pairing` puts them on.
fn definition(pose: &Protein, heavy: &[usize], smiles: &Smiles, pairing: &[usize]) -> Component {
    let name = |at: usize| pose.atoms[heavy[pairing[at]]].name.to_string();
    let atoms = smiles
        .atoms
        .iter()
        .enumerate()
        .map(|(at, atom)| DefinedAtom {
            name: name(at),
            element: atom.element.to_string(),
            aromatic: atom.aromatic,
            hydrogens: usize::from(atom.hydrogens),
            charge: atom.charge,
        })
        .collect();
    let bonds = smiles
        .bonds
        .iter()
        .map(|bond| (name(bond.atoms.0), name(bond.atoms.1), bond.order))
        .collect();

    Component::of_parts(atoms, bonds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::components::Dictionary;
    use crate::smiles;
    use crate::structure::{Atom, ChainId, Molecule, SecondaryStructure};

    const NAG: &str = include_str!("../tests/data/NAG.cif");
    const GPP: &str = include_str!("../tests/data/GPP.cif");

    /// The dictionary entry `code` of `text` at its ideal coordinates, named
    /// as co-folding output names a SMILES ligand: `LIG_A`, its atoms by
    /// element and number, and no definition.
    fn smiles_pose(text: &str, code: &str) -> Protein {
        let dictionary = Dictionary::parse(text);
        let component = dictionary.get(code).unwrap();
        let atoms = component
            .atoms_in_order()
            .enumerate()
            .map(|(at, name)| {
                let element = component.element_of(name).unwrap();
                Atom {
                    serial: at as u32 + 1,
                    name: FixedStr::new(&format!("{element}{}", at + 1)),
                    element: FixedStr::new(element),
                    residue_name: FixedStr::new("LIG_A"),
                    residue_seq: 1,
                    chain_id: ChainId::new("A"),
                    is_hetatm: true,
                    molecule: Molecule::NonPolymer,
                    altloc: None,
                    occupancy: 1.0,
                    insertion_code: None,
                    model_number: 1,
                    position: component.ideal(name).unwrap(),
                    b_factor: 0.0,
                    secondary_structure: SecondaryStructure::Coil,
                }
            })
            .collect();

        Protein {
            atoms,
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Dictionary::default(),
        }
    }

    /// The OpenEye canonical SMILES the dictionary entry gives.
    fn canonical(text: &str) -> Smiles {
        let line = text
            .lines()
            .find(|line| line.contains("SMILES_CANONICAL") && line.contains("OpenEye"))
            .unwrap();
        let quoted = line.rsplit('"').nth(1).unwrap();
        smiles::parse(quoted).unwrap()
    }

    fn all(pose: &Protein) -> Vec<usize> {
        (0..pose.atoms.len()).collect()
    }

    fn names(ids: &[AtomId]) -> Vec<String> {
        let mut names: Vec<String> = ids.iter().map(|(_, _, name)| name.to_string()).collect();
        names.sort();
        names
    }

    #[test]
    fn a_ligand_of_twenty_one_heavy_atoms_pairs() {
        let atoms = (0..21)
            .map(|at| Atom {
                serial: at as u32 + 1,
                name: FixedStr::new(&format!("C{}", at + 1)),
                element: FixedStr::new("C"),
                residue_name: FixedStr::new("LIG_A"),
                residue_seq: 1,
                chain_id: ChainId::new("A"),
                is_hetatm: true,
                molecule: Molecule::NonPolymer,
                altloc: None,
                occupancy: 1.0,
                insertion_code: None,
                model_number: 1,
                position: glam::Vec3::new(at as f32 * 1.26, (at % 2) as f32 * 0.89, 0.0),
                b_factor: 0.0,
                secondary_structure: SecondaryStructure::Coil,
            })
            .collect();
        let pose = Protein {
            atoms,
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Dictionary::default(),
        };
        let chain = smiles::parse(&"C".repeat(21)).unwrap();

        let stated = state(&pose, &all(&pose), &chain).expect("the chains pair");

        assert!(!stated.loose());
    }

    #[test]
    fn a_ligand_as_its_smiles_states_it_has_nothing_turned() {
        for (text, code) in [(NAG, "NAG"), (GPP, "GPP")] {
            let pose = smiles_pose(text, code);
            let stated = state(&pose, &all(&pose), &canonical(text)).expect("the graphs pair");

            assert!(
                stated.inverted().is_empty(),
                "{code}: {:?}",
                stated.inverted()
            );
            assert!(stated.flipped().is_empty(), "{code}");
            assert!(!stated.loose(), "{code}");
            assert_eq!(stated.component().len(), pose.atoms.len());
        }
    }

    #[test]
    fn the_mirror_image_inverts_every_centre_the_string_states() {
        let mut pose = smiles_pose(NAG, "NAG");
        for atom in &mut pose.atoms {
            atom.position.x = -atom.position.x;
        }
        let stated = state(&pose, &all(&pose), &canonical(NAG)).unwrap();

        // NAG's C1 to C5, named here by element and order.
        assert_eq!(names(stated.inverted()), ["C1", "C2", "C3", "C4", "C5"]);
    }

    #[test]
    fn the_string_s_bond_orders_and_hydrogens_define_the_residue() {
        let pose = smiles_pose(GPP, "GPP");
        let stated = state(&pose, &all(&pose), &canonical(GPP)).unwrap();
        let component = stated.component();

        // C2=C3, the dictionary's third and fourth atoms.
        assert!(
            component
                .bonded("C3")
                .any(|(other, order)| other == "C4" && order == BondOrder::Double)
        );
        assert_eq!(component.hydrogen_count("C1"), 2);
    }

    #[test]
    fn geranyl_diphosphate_turned_about_its_e_bond_is_flipped() {
        let mut pose = smiles_pose(GPP, "GPP");
        // Every atom on C2's side of C2=C3 (the dictionary's C1, O1 and the
        // diphosphate: its first two atoms and its last nine) turned a half
        // turn about the bond.
        let (c2, c3) = (pose.atoms[2].position, pose.atoms[3].position);
        let turn = glam::Quat::from_axis_angle((c3 - c2).normalize(), std::f32::consts::PI);
        let side: Vec<usize> = (0..2)
            .chain(pose.atoms.len() - 9..pose.atoms.len())
            .collect();
        for at in side {
            pose.atoms[at].position = c2 + turn * (pose.atoms[at].position - c2);
        }
        // The turned diphosphate comes within bonding distance of the chain,
        // so the pose carries its bonds: the dictionary's, under its names.
        let dictionary = Dictionary::parse(GPP);
        let entry = dictionary.get("GPP").unwrap();
        for (atom, name) in pose.atoms.iter_mut().zip(entry.atoms_in_order()) {
            atom.name = FixedStr::new(name);
        }
        pose.components = Dictionary::default().with("LIG_A", entry.clone());
        let stated = state(&pose, &all(&pose), &canonical(GPP)).unwrap();

        let bonds: Vec<[String; 2]> = stated
            .flipped()
            .iter()
            .map(|ends| ends.map(|(_, _, name)| name.to_string()))
            .collect();
        assert_eq!(bonds.len(), 1);
        let mut ends = bonds[0].clone();
        ends.sort();
        assert_eq!(ends, ["C2", "C3"]);
    }

    #[test]
    fn a_pose_with_a_bond_stretched_past_reading_still_pairs_and_fails_it() {
        // NAG's acetyl methyl pulled 0.8 A along its bond: no longer at
        // bonding distance, so its distances alone give another graph.
        let mut pose = smiles_pose(NAG, "NAG");
        let (carbonyl, methyl) = (6, 7);
        let along = (pose.atoms[methyl].position - pose.atoms[carbonyl].position).normalize();
        pose.atoms[methyl].position += along * 0.8;
        let ligand = all(&pose);
        let stated = state(&pose, &ligand, &canonical(NAG)).expect("it pairs loosely");
        assert!(stated.loose());

        let mut set = crate::select::AtomSet::empty(pose.atoms.len());
        for atom in ligand {
            set.insert(atom);
        }
        let defined = crate::checks::with_stated(pose, std::slice::from_ref(&stated));
        assert!(crate::checks::geometry(&defined, &set).unwrap().bonds >= 1);
    }

    /// 1,2-difluoroethene, E: each fluorine on its own side of the double bond.
    fn difluoroethene() -> Protein {
        let mut pose = smiles_pose(NAG, "NAG");
        pose.atoms.truncate(4);
        for (atom, (element, at)) in pose.atoms.iter_mut().zip([
            ("C", Vec3::ZERO),
            ("C", Vec3::new(1.34, 0.0, 0.0)),
            ("F", Vec3::new(-0.68, 1.16, 0.0)),
            ("F", Vec3::new(2.02, -1.16, 0.0)),
        ]) {
            atom.element = FixedStr::new(element);
            atom.position = at;
        }
        pose
    }

    #[test]
    fn a_direction_written_on_a_hydrogen_states_the_double_bond() {
        let pose = difluoroethene();
        let flipped = |text: &str| {
            state(&pose, &all(&pose), &smiles::parse(text).unwrap())
                .unwrap()
                .flipped()
                .len()
        };

        assert_eq!(flipped("F/C=C/F"), 0);
        // The hydrogen trans to the second fluorine puts the two fluorines cis.
        assert_eq!(flipped("[H]/C(F)=C/F"), 1);
    }

    #[test]
    fn graphs_that_differ_do_not_pair() {
        let pose = smiles_pose(NAG, "NAG");
        assert!(state(&pose, &all(&pose), &canonical(GPP)).is_none());
        assert!(state(&pose, &all(&pose), &smiles::parse("CCO").unwrap()).is_none());
    }
}
