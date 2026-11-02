//! Ligand RMSD between a pose and a reference, the pose superposed on the
//! reference's binding pocket.
//!
//! A co-folding model and a crystal structure sit in unrelated frames, so the
//! pose is first superposed on the pocket: the backbone atoms (N, CA, C) within
//! [`POCKET_REACH`] of the reference ligand, on the protein chain holding the
//! most of them, as AlphaFold 3's evaluation does (Abramson et al., Nature
//! 630:493, 2024). The ligand's heavy atoms then pair with the reference's by
//! name where the pose keeps the reference's atom names, and by bond graph
//! where it does not, as output from a SMILES input does not. Either way the
//! pairing is the one closest in space among those the bond graph allows, so a
//! carboxylate or a ring that swaps equivalent atoms costs nothing, and nor do
//! copies of a ligand listed in another order. The same pairing says which of
//! the reference ligand's stereocentres the pose inverts, which of its double
//! bonds it turns from E to Z or back, and what of the pose's own geometry
//! strays from the reference ligand's.

use std::collections::{HashMap, HashSet};
use std::fmt;

use glam::Vec3;
use nalgebra::{Matrix3, Vector3};

use crate::components::BondOrder;
use crate::fixed_str::FixedStr;
use crate::geometry::{Chemistry, Template};
use crate::interactions::connectivity;
use crate::numbering::{ResidueKey, ResidueMap, residue_key};
use crate::pose::is_hydrogen;
use crate::select::AtomSet;
use crate::stereo::{self, Stereocentre};
use crate::structure::{ChainId, Protein};

/// How far from the reference ligand a backbone atom may lie and still belong
/// to the pocket, in angstroms.
pub const POCKET_REACH: f32 = 10.0;

/// The most candidate atoms weighed while pairing one ligand; past that the
/// closest pairing found stands. A ligand of many equivalent groups, or many
/// copies of one, could otherwise be searched without end.
const MAX_STEPS: usize = 1_000_000;

/// Why a pose has no RMSD.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NoRmsd {
    /// Fewer than three atoms of any pocket chain have a counterpart in the
    /// pose.
    Pocket,
    /// The pose's ligand has another number of heavy atoms than the
    /// reference's.
    Atoms { pose: usize, reference: usize },
    /// The ligand atoms pair with the reference's neither by name nor by bond
    /// graph.
    Unpaired,
}

impl fmt::Display for NoRmsd {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NoRmsd::Pocket => write!(
                f,
                "fewer than three of the reference's pocket atoms have a counterpart"
            ),
            NoRmsd::Atoms { pose, reference } => write!(
                f,
                "its ligand has {pose} heavy atoms and the reference's {reference}"
            ),
            NoRmsd::Unpaired => write!(
                f,
                "its ligand atoms pair with the reference's neither by name nor by bond graph"
            ),
        }
    }
}

/// A pose's ligand measured against the reference's.
#[derive(Clone, Debug, PartialEq)]
pub struct Fit {
    /// The heavy-atom RMSD in angstroms.
    pub rmsd: f32,
    /// The reference ligand's stereocentres the pose turns the other way, in
    /// the reference's order.
    pub inverted: Vec<AtomId>,
}

pub use crate::geometry::Geometry;

/// A pose's ligand measured against the reference's, its stereochemistry and
/// its own geometry included.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Assessment {
    /// The heavy-atom RMSD in angstroms.
    pub rmsd: f32,
    /// The reference ligand's stereocentres the pose turns the other way, in
    /// the reference's order.
    pub inverted: Vec<AtomId>,
    /// The reference ligand's E or Z double bonds the pose turns the other
    /// way, each by its two atoms.
    pub flipped: Vec<[AtomId; 2]>,
    /// What of the pose's bond lengths, angles, internal distances and flat
    /// groups strays from the reference ligand's.
    pub geometry: Geometry,
}

/// A reference atom by residue, component and name.
pub type AtomId = (ResidueKey, FixedStr<6>, FixedStr<4>);

/// A heavy atom's place in a ligand across structures: its component, which
/// residue of that component it is, counted in order, and its name. A pose
/// may number its residues otherwise, but lists a ligand's residues in the
/// same order.
type AtomKey = (FixedStr<6>, usize, FixedStr<4>);

/// A pocket atom: residue, atom name and position.
type PocketAtom = (ResidueKey, FixedStr<4>, Vec3);

/// What a pose's ligand RMSD is measured against: the reference ligand's heavy
/// atoms, their bond graph and stereocentres, and the pocket atoms a pose is
/// superposed on, by chain, the chain holding the most first.
#[derive(Clone, Debug)]
pub struct RmsdReference {
    keys: Vec<AtomKey>,
    atoms: Vec<AtomId>,
    graph: LigandGraph,
    /// Each atom's bond orders summed, when its definition gives them.
    order_sums: Option<Vec<u8>>,
    stereocentres: Vec<Stereocentre>,
    /// The E or Z double bonds, each with whether a symmetry of the bond graph
    /// may carry one of its ends onto an atom of other chemistry.
    double_bonds: Vec<(Stereocentre, bool)>,
    template: Template,
    pocket: Vec<Vec<PocketAtom>>,
}

impl RmsdReference {
    /// The reference made of the `ligand` atoms of `reference`.
    pub fn new(reference: &Protein, ligand: &AtomSet) -> Self {
        let heavy = heavy_atoms(reference, ligand);
        let graph = LigandGraph::of(reference, &heavy);
        let chemistry = Chemistry::of(reference, &heavy, &graph.neighbours);
        let doubles = chemistry
            .as_ref()
            .map_or(&[][..], |chemistry| &chemistry.doubles[..]);
        let classes = stereo::classes(&graph.elements, &graph.neighbours);
        let exchangeable = chemistry
            .as_ref()
            .map_or_else(Vec::new, |chemistry| chemistry.exchangeable(&classes));

        Self {
            keys: atom_keys(reference, &heavy),
            atoms: heavy
                .iter()
                .map(|&atom| {
                    let atom = &reference.atoms[atom];
                    (residue_key(atom), atom.residue_name, atom.name)
                })
                .collect(),
            stereocentres: stereo::stereocentres(
                &graph.elements,
                &graph.positions,
                &graph.neighbours,
            ),
            double_bonds: stereo::double_bonds(
                &graph.elements,
                &graph.positions,
                &graph.neighbours,
                doubles,
            )
            .into_iter()
            .map(|bond| {
                let ends = bond.bond().expect("a double bond has two ends");
                let traded = ends.iter().any(|&end| exchangeable[end]);
                (bond, traded)
            })
            .collect(),
            template: Template::new(
                &graph.elements,
                &graph.positions,
                &graph.neighbours,
                chemistry.as_ref().map(|chemistry| {
                    chemistry.flat(&graph.elements, &graph.neighbours, &exchangeable)
                }),
            ),
            order_sums: bond_order_sums(reference, &heavy, &graph),
            graph,
            pocket: pocket(reference, ligand, &heavy),
        }
    }

    /// True when some chain holds the three pocket atoms a superposition
    /// needs; without them no pose gets an RMSD.
    pub fn has_pocket(&self) -> bool {
        self.pocket.iter().any(|atoms| atoms.len() >= 3)
    }

    /// The `ligand` atoms of `pose` measured against the reference's once the
    /// pose is superposed on the reference's pocket through `map`.
    pub fn fit(&self, pose: &Protein, ligand: &AtomSet, map: &ResidueMap) -> Result<Fit, NoRmsd> {
        self.assess(pose, ligand, map).map(|assessment| Fit {
            rmsd: assessment.rmsd,
            inverted: assessment.inverted,
        })
    }

    /// The `ligand` atoms of `pose` measured as [`RmsdReference::fit`]
    /// measures them, with the reference's double bonds they turn from E to Z
    /// or back and the checks of their own geometry against the reference
    /// ligand's.
    pub fn assess(
        &self,
        pose: &Protein,
        ligand: &AtomSet,
        map: &ResidueMap,
    ) -> Result<Assessment, NoRmsd> {
        let heavy = heavy_atoms(pose, ligand);
        if heavy.len() != self.keys.len() {
            return Err(NoRmsd::Atoms {
                pose: heavy.len(),
                reference: self.keys.len(),
            });
        }
        let superposition = self.superposition(pose, map)?;
        let placed: Vec<Vec3> = heavy
            .iter()
            .map(|&atom| superposition.apply(pose.atoms[atom].position))
            .collect();

        let named = self.named(pose, &heavy, &placed);
        let by_name = named.is_some();
        let mobile = named.unwrap_or_else(|| LigandGraph::of(pose, &heavy).placed_at(placed));
        let (squared, pairing) =
            closest_pairing(&mobile, &self.graph, &[]).ok_or(NoRmsd::Unpaired)?;
        let mut paired = vec![Vec3::ZERO; pairing.len()];
        for (&position, &partner) in mobile.positions.iter().zip(&pairing) {
            paired[partner] = position;
        }

        // The closest pairing keeps bonds but not bond orders, so a symmetry
        // of the bond graph may carry a double bond onto a single one; atoms
        // named as the reference's are read by name instead.
        let read = if by_name { &mobile.positions } else { &paired };
        let flipped: Vec<[usize; 2]> = if by_name {
            self.double_bonds
                .iter()
                .map(|(bond, _)| bond)
                .filter(|bond| bond.inverted_at(bond.atoms().map(|atom| read[atom])))
                .filter_map(Stereocentre::bond)
                .collect()
        } else {
            let fixed: Vec<Stereocentre> = self
                .double_bonds
                .iter()
                .filter(|(_, traded)| !traded)
                .map(|(bond, _)| bond.clone())
                .collect();
            self.inverted(&mobile, &paired, &fixed)
                .into_iter()
                .filter_map(|at| fixed[at].bond())
                .collect()
        };

        Ok(Assessment {
            rmsd: (squared / heavy.len().max(1) as f32).sqrt(),
            inverted: self
                .inverted(&mobile, &paired, &self.stereocentres)
                .into_iter()
                .map(|at| self.atoms[self.stereocentres[at].centre])
                .collect(),
            flipped: flipped
                .into_iter()
                .map(|ends| ends.map(|end| self.atoms[end]))
                .collect(),
            geometry: self.template.check(read, by_name),
        })
    }

    /// Which of the reference's `elements`, stereocentres or else double
    /// bonds, the `mobile` atoms turn the other way, by index, `paired`
    /// holding their positions in the reference's order. A symmetry of the
    /// ligand may carry one onto another of its class that turns the other
    /// way, as in a meso compound, so before one such is reported, a pairing
    /// that keeps every one is sought: if there is one, the pose is the same
    /// molecule turned about.
    fn inverted(
        &self,
        mobile: &LigandGraph,
        paired: &[Vec3],
        elements: &[Stereocentre],
    ) -> Vec<usize> {
        let inverted: Vec<usize> = (0..elements.len())
            .filter(|&at| elements[at].inverted_at(elements[at].atoms().map(|atom| paired[atom])))
            .collect();

        let exchangeable = inverted.iter().any(|&at| {
            elements
                .iter()
                .filter(|other| other.class == elements[at].class)
                .count()
                > 1
        });
        if exchangeable && closest_pairing(mobile, &self.graph, elements).is_some() {
            return Vec::new();
        }

        inverted
    }

    /// The superposition of `pose` on the first pocket chain with at least
    /// three atoms `map` finds a counterpart for, whatever the pose's ligand.
    pub fn superposition(&self, pose: &Protein, map: &ResidueMap) -> Result<Superposition, NoRmsd> {
        let mut positions: HashMap<(ResidueKey, FixedStr<4>), Vec3> = HashMap::new();
        for atom in &pose.atoms {
            positions
                .entry((residue_key(atom), atom.name))
                .or_insert(atom.position);
        }

        self.pocket
            .iter()
            .find_map(|atoms| {
                let (mobile, target): (Vec<Vec3>, Vec<Vec3>) = atoms
                    .iter()
                    .filter_map(|&(residue, name, position)| {
                        let theirs = map.pose_residue(&residue)?;
                        Some((*positions.get(&(theirs, name))?, position))
                    })
                    .unzip();
                (mobile.len() >= 3).then(|| Superposition::of(&mobile, &target))
            })
            .ok_or(NoRmsd::Pocket)
    }

    /// The reference ligand's atoms in order, each with the one of a
    /// molecule's `elements` that stands for it: the pairing nearest in place,
    /// with `positions` taken as they are, among those that keep elements and
    /// the molecule's `bonds`. `None` when the two bond graphs differ.
    ///
    /// With `sums`, each atom's bond orders summed, a pairing must also keep
    /// those, where the reference's definitions give its own.
    pub(crate) fn counterparts(
        &self,
        elements: &[FixedStr<2>],
        positions: &[Vec3],
        bonds: &[(usize, usize)],
        sums: Option<&[u8]>,
    ) -> Option<Vec<(AtomId, usize)>> {
        let mut neighbours = vec![Vec::new(); elements.len()];
        for &(a, b) in bonds {
            if a != b && !neighbours[a].contains(&b) {
                neighbours[a].push(b);
                neighbours[b].push(a);
            }
        }
        let (labels, own) = match (sums, &self.order_sums) {
            (Some(sums), Some(own)) => (sums.to_vec(), own.clone()),
            _ => (vec![0; elements.len()], vec![0; self.graph.elements.len()]),
        };
        let molecule = LigandGraph {
            elements: elements
                .iter()
                .map(|element| FixedStr::uppercased(element.as_str()))
                .collect(),
            positions: positions.to_vec(),
            neighbours,
            labels,
        };
        let reference = LigandGraph {
            labels: own,
            ..self.graph.clone()
        };

        let (_, pairing) = closest_pairing(&molecule, &reference, &[])?;
        let mut standing = vec![0; pairing.len()];
        for (atom, &partner) in pairing.iter().enumerate() {
            standing[partner] = atom;
        }

        Some(self.atoms.iter().copied().zip(standing).collect())
    }

    /// Whether the reference's definitions give its atoms' bond orders.
    pub(crate) fn has_bond_orders(&self) -> bool {
        self.order_sums.is_some()
    }

    /// The reference's bond graph with each atom at the `placed` position of
    /// the pose atom of the same key and element, when every pose atom has
    /// one: the pose is then measured on the reference's bonds, so a distorted
    /// pose still pairs.
    fn named(&self, pose: &Protein, heavy: &[usize], placed: &[Vec3]) -> Option<LigandGraph> {
        let index: HashMap<AtomKey, usize> = self
            .keys
            .iter()
            .enumerate()
            .map(|(i, &key)| (key, i))
            .collect();
        let mut positions = vec![None; self.keys.len()];
        for ((key, &atom), &position) in atom_keys(pose, heavy).iter().zip(heavy).zip(placed) {
            let at = *index.get(key)?;
            let element = FixedStr::uppercased(pose.atoms[atom].element.as_str());
            if element != self.graph.elements[at] || positions[at].is_some() {
                return None;
            }
            positions[at] = Some(position);
        }

        let positions = positions.into_iter().collect::<Option<Vec<Vec3>>>()?;
        Some(self.graph.placed_at(positions))
    }
}

/// Each of the `heavy` atoms' bond orders summed, in half bonds, over its
/// bonds to the others: within a residue as the residue's definition orders
/// them, between residues (`graph`'s) as single bonds. `None` when a
/// definition is missing or lacks one of the atoms.
fn bond_order_sums(protein: &Protein, heavy: &[usize], graph: &LigandGraph) -> Option<Vec<u8>> {
    let present: HashSet<(ResidueKey, FixedStr<4>)> = heavy
        .iter()
        .map(|&atom| (residue_key(&protein.atoms[atom]), protein.atoms[atom].name))
        .collect();

    heavy
        .iter()
        .enumerate()
        .map(|(at, &index)| {
            let atom = &protein.atoms[index];
            let residue = residue_key(atom);
            let component = protein.components.get(atom.residue_name.as_str())?;
            component.element_of(atom.name.as_str())?;
            let within: u8 = component
                .bonded(atom.name.as_str())
                .filter(|(name, _)| present.contains(&(residue, FixedStr::new(name))))
                .map(|(_, order)| order.halves())
                .sum();
            let between = graph.neighbours[at]
                .iter()
                .filter(|&&other| residue_key(&protein.atoms[heavy[other]]) != residue)
                .count() as u8;

            Some(within + between * BondOrder::Single.halves())
        })
        .collect()
}

/// The ligand atoms that are not hydrogen, in order.
fn heavy_atoms(protein: &Protein, ligand: &AtomSet) -> Vec<usize> {
    ligand
        .iter()
        .filter(|&atom| !is_hydrogen(&protein.atoms[atom]))
        .collect()
}

/// The key of each of the `heavy` atoms.
fn atom_keys(protein: &Protein, heavy: &[usize]) -> Vec<AtomKey> {
    let mut seen: HashMap<FixedStr<6>, Vec<ResidueKey>> = HashMap::new();
    heavy
        .iter()
        .map(|&atom| {
            let atom = &protein.atoms[atom];
            let residues = seen.entry(atom.residue_name).or_default();
            let residue = residue_key(atom);
            let ordinal = residues
                .iter()
                .position(|&known| known == residue)
                .unwrap_or_else(|| {
                    residues.push(residue);
                    residues.len() - 1
                });
            (atom.residue_name, ordinal, atom.name)
        })
        .collect()
}

/// A ligand's heavy atoms as a graph: each atom's element, position and
/// bonded neighbours.
#[derive(Clone, Debug, Default)]
struct LigandGraph {
    elements: Vec<FixedStr<2>>,
    positions: Vec<Vec3>,
    neighbours: Vec<Vec<usize>>,
    /// What each atom must share with its partner beside its element: its
    /// bond orders summed where both graphs know them, else nothing (0).
    labels: Vec<u8>,
}

impl LigandGraph {
    /// The graph of `atoms`, with bonds read from the geometry: a pose may
    /// define no chemistry for its ligand, nor name its atoms as the reference
    /// does.
    fn of(protein: &Protein, atoms: &[usize]) -> LigandGraph {
        let mut neighbours = vec![Vec::new(); atoms.len()];
        for (a, &i) in atoms.iter().enumerate() {
            for (b, &j) in atoms.iter().enumerate().skip(a + 1) {
                if connectivity::at_bonding_distance(&protein.atoms[i], &protein.atoms[j]) {
                    neighbours[a].push(b);
                    neighbours[b].push(a);
                }
            }
        }

        LigandGraph {
            elements: atoms
                .iter()
                .map(|&atom| FixedStr::uppercased(protein.atoms[atom].element.as_str()))
                .collect(),
            positions: atoms
                .iter()
                .map(|&atom| protein.atoms[atom].position)
                .collect(),
            neighbours,
            labels: vec![0; atoms.len()],
        }
    }

    /// The same atoms and bonds at `positions`.
    fn placed_at(&self, positions: Vec<Vec3>) -> LigandGraph {
        LigandGraph {
            positions,
            ..self.clone()
        }
    }

    /// Each atom's element, number of bonds and label, sorted, which two
    /// graphs that can pair share.
    fn signature(&self) -> Vec<(FixedStr<2>, usize, u8)> {
        let mut signature: Vec<(FixedStr<2>, usize, u8)> = self
            .elements
            .iter()
            .zip(&self.neighbours)
            .zip(&self.labels)
            .map(|((&element, bonded), &label)| (element, bonded.len(), label))
            .collect();
        signature.sort();
        signature
    }
}

/// The least sum of squared distances over the pairings of `mobile`'s atoms
/// with `target`'s that keep elements and bonds, and turn none of `target`'s
/// stereocentres in `keep` the other way, with the pairing, each `mobile`
/// atom's partner, or `None` when there is no such pairing.
///
/// A branch-and-bound search: atoms go in bond order, so each meets the
/// neighbours of atoms already paired; each tries its nearest candidates
/// first; and a branch stops once its sum reaches the best so far. The first
/// pairing found is the nearest-first one, which for a pose that sits where the
/// reference does is the answer, and the rest are cut short. Past [`MAX_STEPS`]
/// the best found stands.
fn closest_pairing(
    mobile: &LigandGraph,
    target: &LigandGraph,
    keep: &[Stereocentre],
) -> Option<(f32, Vec<usize>)> {
    if mobile.signature() != target.signature() {
        return None;
    }
    let mut search = Pairing {
        mobile,
        target,
        keep,
        order: bond_order(mobile),
        paired: vec![None; mobile.elements.len()],
        owner: vec![None; target.elements.len()],
        best: None,
        steps: 0,
    };
    search.run(0, 0.0);

    search.best
}

/// Every atom of `graph`, each component in turn from its most bonded atom,
/// breadth first.
fn bond_order(graph: &LigandGraph) -> Vec<usize> {
    let mut order = Vec::with_capacity(graph.elements.len());
    let mut placed = vec![false; graph.elements.len()];
    while order.len() < graph.elements.len() {
        let start = (0..graph.elements.len())
            .filter(|&atom| !placed[atom])
            .max_by_key(|&atom| (graph.neighbours[atom].len(), std::cmp::Reverse(atom)))
            .expect("an atom is left to place");
        placed[start] = true;
        let mut next = order.len();
        order.push(start);
        while next < order.len() {
            for &bonded in &graph.neighbours[order[next]] {
                if !placed[bonded] {
                    placed[bonded] = true;
                    order.push(bonded);
                }
            }
            next += 1;
        }
    }

    order
}

/// The search [`closest_pairing`] runs.
struct Pairing<'a> {
    mobile: &'a LigandGraph,
    target: &'a LigandGraph,
    keep: &'a [Stereocentre],
    order: Vec<usize>,
    /// Each mobile atom's partner so far.
    paired: Vec<Option<usize>>,
    /// Each target atom's partner so far.
    owner: Vec<Option<usize>>,
    best: Option<(f32, Vec<usize>)>,
    steps: usize,
}

impl Pairing<'_> {
    fn run(&mut self, depth: usize, sum: f32) {
        let Some(&atom) = self.order.get(depth) else {
            if self.best.as_ref().is_none_or(|&(best, _)| sum < best) {
                let pairing = self
                    .paired
                    .iter()
                    .map(|partner| partner.expect("every atom is paired at a leaf"))
                    .collect();
                self.best = Some((sum, pairing));
            }
            return;
        };

        for (squared, candidate) in self.candidates(atom) {
            self.steps += 1;
            let hopeless = self
                .best
                .as_ref()
                .is_some_and(|&(best, _)| sum + squared >= best);
            if hopeless || self.steps > MAX_STEPS {
                break;
            }
            if !self.consistent(atom, candidate) || !self.keeps_turns(atom, candidate) {
                continue;
            }
            self.paired[atom] = Some(candidate);
            self.owner[candidate] = Some(atom);
            self.run(depth + 1, sum + squared);
            self.owner[candidate] = None;
            self.paired[atom] = None;
        }
    }

    /// The target atoms `atom` may pair with, nearest first: free, of its
    /// element and number of bonds, and bonded to the partner of an atom it is
    /// bonded to, if it has one paired yet.
    fn candidates(&self, atom: usize) -> Vec<(f32, usize)> {
        let (mobile, target) = (self.mobile, self.target);
        let anchor = mobile.neighbours[atom]
            .iter()
            .find_map(|&bonded| self.paired[bonded]);
        let pool: Vec<usize> = match anchor {
            Some(partner) => target.neighbours[partner].clone(),
            None => (0..target.elements.len()).collect(),
        };

        let mut candidates: Vec<(f32, usize)> = pool
            .into_iter()
            .filter(|&other| {
                self.owner[other].is_none()
                    && target.elements[other] == mobile.elements[atom]
                    && target.labels[other] == mobile.labels[atom]
                    && target.neighbours[other].len() == mobile.neighbours[atom].len()
            })
            .map(|other| {
                let squared = (mobile.positions[atom] - target.positions[other]).length_squared();
                (squared, other)
            })
            .collect();
        candidates.sort_by(|a, b| a.0.total_cmp(&b.0));
        candidates
    }

    /// True when pairing `atom` with `candidate` keeps every bond, and every
    /// absence of one, to the atoms already paired.
    fn consistent(&self, atom: usize, candidate: usize) -> bool {
        let (mobile, target) = (self.mobile, self.target);
        let mut bonded_paired = 0;
        for &bonded in &mobile.neighbours[atom] {
            if let Some(partner) = self.paired[bonded] {
                if !target.neighbours[candidate].contains(&partner) {
                    return false;
                }
                bonded_paired += 1;
            }
        }
        let partners_bonded = target.neighbours[candidate]
            .iter()
            .filter(|&&other| self.owner[other].is_some())
            .count();

        bonded_paired == partners_bonded
    }

    /// True when pairing `atom` with `candidate` turns none of the
    /// stereocentres to keep that it completes the other way.
    fn keeps_turns(&self, atom: usize, candidate: usize) -> bool {
        self.keep.iter().all(|stereocentre| {
            let atoms = stereocentre.atoms();
            if !atoms.contains(&candidate) {
                return true;
            }
            let position = |target: usize| {
                let partner = match target == candidate {
                    true => Some(atom),
                    false => self.owner[target],
                };
                partner.map(|partner| self.mobile.positions[partner])
            };
            match atoms.map(position) {
                [Some(centre), Some(a), Some(b), Some(c)] => {
                    !stereocentre.inverted_at([centre, a, b, c])
                }
                _ => true,
            }
        })
    }
}

/// The backbone atoms (N, CA, C) within [`POCKET_REACH`] of the `heavy` ligand
/// atoms, each once, by polymer chain, the chain holding the most first.
fn pocket(protein: &Protein, ligand: &AtomSet, heavy: &[usize]) -> Vec<Vec<PocketAtom>> {
    let reach = POCKET_REACH * POCKET_REACH;
    let near = |position: Vec3| {
        heavy
            .iter()
            .any(|&atom| (protein.atoms[atom].position - position).length_squared() <= reach)
    };

    let mut seen = HashSet::new();
    let mut chains: Vec<(ChainId, Vec<PocketAtom>)> = Vec::new();
    for (index, atom) in protein.atoms.iter().enumerate() {
        let backbone = matches!(atom.name.as_str(), "N" | "CA" | "C");
        if ligand.contains(index) || !atom.is_polymer() || !backbone || !near(atom.position) {
            continue;
        }
        if !seen.insert((residue_key(atom), atom.name)) {
            continue;
        }
        let pocket_atom = (residue_key(atom), atom.name, atom.position);
        match chains.iter_mut().find(|(chain, _)| *chain == atom.chain_id) {
            Some((_, atoms)) => atoms.push(pocket_atom),
            None => chains.push((atom.chain_id, vec![pocket_atom])),
        }
    }
    chains.sort_by_key(|(_, atoms)| std::cmp::Reverse(atoms.len()));

    chains.into_iter().map(|(_, atoms)| atoms).collect()
}

/// The rotation and translation that best superpose one set of points on
/// another, pair by pair, by least squares (Kabsch, Acta Cryst A32:922, 1976).
#[derive(Clone, Copy, Debug)]
pub struct Superposition {
    rotation: Matrix3<f64>,
    from: Vector3<f64>,
    to: Vector3<f64>,
}

impl Superposition {
    /// The superposition of `mobile` on `target`, which pair index by index;
    /// neither may be empty.
    pub fn of(mobile: &[Vec3], target: &[Vec3]) -> Superposition {
        let centre = |points: &[Vec3]| {
            points.iter().map(|&p| vector(p)).sum::<Vector3<f64>>() / points.len() as f64
        };
        let (from, to) = (centre(mobile), centre(target));
        let covariance: Matrix3<f64> = mobile
            .iter()
            .zip(target)
            .map(|(&m, &t)| (vector(m) - from) * (vector(t) - to).transpose())
            .sum();

        let svd = covariance.svd(true, true);
        let u = svd.u.expect("the left singular vectors were asked for");
        let v = svd
            .v_t
            .expect("the right singular vectors were asked for")
            .transpose();
        // A reflection fits as well as a rotation when the points are planar;
        // flipping the weakest axis keeps the transform proper.
        let handedness = if (v * u.transpose()).determinant() < 0.0 {
            -1.0
        } else {
            1.0
        };
        let rotation =
            v * Matrix3::from_diagonal(&Vector3::new(1.0, 1.0, handedness)) * u.transpose();

        Superposition { rotation, from, to }
    }

    /// Where `point` of the mobile set's frame lands in the target's.
    pub fn apply(&self, point: Vec3) -> Vec3 {
        let moved = self.rotation * (vector(point) - self.from) + self.to;
        Vec3::new(moved.x as f32, moved.y as f32, moved.z as f32)
    }
}

fn vector(point: Vec3) -> Vector3<f64> {
    Vector3::new(f64::from(point.x), f64::from(point.y), f64::from(point.z))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::components::Dictionary;
    use crate::numbering::{polymer_chains, residue_map};
    use crate::select;
    use crate::structure::{Atom, Molecule, SecondaryStructure};
    use glam::Quat;

    fn atom(seq: i32, residue: &str, name: &str, element: &str, position: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new("A"),
            is_hetatm: matches!(residue, "LIG" | "SO4"),
            molecule: Molecule::Unknown,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    /// Three glycines around the origin and, at it, a ligand shaped like acetate:
    /// a methyl carbon, a carboxyl carbon and its two equivalent oxygens.
    fn site() -> Protein {
        let mut atoms = Vec::new();
        for (seq, centre) in [
            (1, Vec3::new(4.0, 0.0, 0.0)),
            (2, Vec3::new(0.0, 5.0, 1.0)),
            (3, Vec3::new(-3.0, -2.0, 4.0)),
        ] {
            atoms.push(atom(
                seq,
                "GLY",
                "N",
                "N",
                centre + Vec3::new(-1.2, 0.5, 0.0),
            ));
            atoms.push(atom(seq, "GLY", "CA", "C", centre));
            atoms.push(atom(
                seq,
                "GLY",
                "C",
                "C",
                centre + Vec3::new(1.0, 0.8, 0.6),
            ));
        }
        atoms.extend([
            atom(900, "LIG", "C2", "C", Vec3::new(-1.5, 0.0, 0.0)),
            atom(900, "LIG", "C1", "C", Vec3::ZERO),
            atom(900, "LIG", "O1", "O", Vec3::new(0.62, 1.07, 0.0)),
            atom(900, "LIG", "O2", "O", Vec3::new(0.62, -1.07, 0.0)),
        ]);

        Protein {
            atoms,
            title: String::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    fn ligand_of(structure: &Protein) -> AtomSet {
        select("ligand", structure).expect("a valid query")
    }

    fn rmsd_of(pose: &Protein, reference: &Protein) -> Result<f32, NoRmsd> {
        fit_of(pose, reference).map(|fit| fit.rmsd)
    }

    fn fit_of(pose: &Protein, reference: &Protein) -> Result<Fit, NoRmsd> {
        assessment_of(pose, reference).map(|assessment| Fit {
            rmsd: assessment.rmsd,
            inverted: assessment.inverted,
        })
    }

    fn assessment_of(pose: &Protein, reference: &Protein) -> Result<Assessment, NoRmsd> {
        let map = residue_map(
            &polymer_chains(pose),
            &polymer_chains(reference),
            &[],
            |_| 0,
        );

        RmsdReference::new(reference, &ligand_of(reference)).assess(pose, &ligand_of(pose), &map)
    }

    const CHLOROFLUOROETHENE: &str = "\
data_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
LIG C1 C N
LIG C2 C N
LIG CL1 CL N
LIG F2 F N
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG C1 C2 DOUB
LIG C1 CL1 SING
LIG C2 F2 SING
";

    /// The glycines of [`site`] around 1-chloro-2-fluoroethene, defined by
    /// `definitions`: E, or Z when `cis`, its bond lengths and angles alike.
    fn alkene_site(cis: bool, definitions: &str) -> Protein {
        let mut structure = site();
        structure.atoms.retain(|atom| atom.residue_name != "LIG");
        let side = if cis { 1.0 } else { -1.0 };
        for (name, element, position) in [
            ("C1", "C", Vec3::new(-0.67, 0.0, 0.0)),
            ("C2", "C", Vec3::new(0.67, 0.0, 0.0)),
            ("CL1", "CL", Vec3::new(-1.53, 1.49, 0.0)),
            ("F2", "F", Vec3::new(1.35, side * 1.17, 0.0)),
        ] {
            structure
                .atoms
                .push(atom(900, "LIG", name, element, position));
        }
        structure.components = Dictionary::parse(definitions);

        structure
    }

    /// The glycines of [`site`] around a ligand with one stereocentre: C1 at
    /// the origin, bonded to a methyl carbon, a nitrogen, an oxygen and a
    /// fluorine at the corners of a tetrahedron.
    fn chiral_site() -> Protein {
        let mut structure = site();
        structure.atoms.retain(|atom| atom.residue_name != "LIG");
        let corners = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
        ];
        structure
            .atoms
            .push(atom(900, "LIG", "C1", "C", Vec3::ZERO));
        for ((name, element), corner) in [("C2", "C"), ("N1", "N"), ("O1", "O"), ("F1", "F")]
            .into_iter()
            .zip(corners)
        {
            let position = corner.normalize() * 1.5;
            structure
                .atoms
                .push(atom(900, "LIG", name, element, position));
        }

        structure
    }

    /// The glycines of [`site`] around meso-2,3-dichlorobutane, staggered so
    /// that C2 and C3, which turn opposite ways, face each other through the
    /// origin; each keeps its hydrogen implicit. `chlorine` places C2's
    /// chlorine, which the hydrogen's corner instead makes the molecule
    /// chiral.
    fn meso_site(chlorine: Vec3) -> Protein {
        let mut structure = site();
        structure.atoms.retain(|atom| atom.residue_name != "LIG");
        let c2 = Vec3::new(-0.77, 0.0, 0.0);
        let c1 = c2 + Vec3::new(-1.0 / 3.0, 0.9428, 0.0) * 1.53;
        let cl3 = -(c2 + Vec3::new(-1.0 / 3.0, -0.4714, 0.8165) * 1.78);
        for (name, element, position) in [
            ("C1", "C", c1),
            ("C2", "C", c2),
            ("CL2", "CL", c2 + chlorine * 1.78),
            ("C3", "C", -c2),
            ("C4", "C", -c1),
            ("CL3", "CL", cl3),
        ] {
            structure
                .atoms
                .push(atom(900, "LIG", name, element, position));
        }

        structure
    }

    const MESO_CHLORINE: Vec3 = Vec3::new(-1.0 / 3.0, -0.4714, 0.8165);
    const CHIRAL_CHLORINE: Vec3 = Vec3::new(-1.0 / 3.0, -0.4714, -0.8165);

    fn named(structure: &Protein, name: &str) -> usize {
        structure
            .atoms
            .iter()
            .position(|atom| atom.name == name)
            .unwrap()
    }

    fn moved(structure: &Protein, motion: impl Fn(Vec3) -> Vec3) -> Protein {
        let mut moved = structure.clone();
        for atom in &mut moved.atoms {
            atom.position = motion(atom.position);
        }
        moved
    }

    #[test]
    fn a_superposition_undoes_a_rigid_motion() {
        let points: Vec<Vec3> = site().atoms.iter().map(|atom| atom.position).collect();
        let turn = Quat::from_euler(glam::EulerRot::XYZ, 0.3, -1.1, 2.0);
        let mobile: Vec<Vec3> = points
            .iter()
            .map(|&p| turn * p + Vec3::new(7.0, -3.0, 12.0))
            .collect();

        let superposition = Superposition::of(&mobile, &points);

        for (&m, &p) in mobile.iter().zip(&points) {
            assert!((superposition.apply(m) - p).length() < 1e-4, "{m} {p}");
        }
    }

    #[test]
    fn a_pose_moved_rigidly_with_its_pocket_has_no_rmsd() {
        let reference = site();
        let turn = Quat::from_rotation_z(1.9);
        let pose = moved(&reference, |p| turn * p + Vec3::new(-20.0, 4.0, 9.0));

        assert!(rmsd_of(&pose, &reference).unwrap() < 1e-3);
    }

    #[test]
    fn a_pose_is_laid_on_the_reference_pocket_whatever_its_ligand() {
        let reference = site();
        let turn = Quat::from_rotation_z(1.9);
        let mut pose = moved(&reference, |p| turn * p + Vec3::new(-20.0, 4.0, 9.0));
        // A ligand of another size has no RMSD, but its pocket still fits.
        pose.atoms.pop();
        let map = residue_map(
            &polymer_chains(&pose),
            &polymer_chains(&reference),
            &[],
            |_| 0,
        );

        let superposition = RmsdReference::new(&reference, &ligand_of(&reference))
            .superposition(&pose, &map)
            .expect("the pocket has counterparts");

        for (placed, original) in pose.atoms.iter().zip(&reference.atoms) {
            let landed = superposition.apply(placed.position);
            assert!((landed - original.position).length() < 1e-3, "{landed}");
        }
    }

    #[test]
    fn a_displaced_ligand_has_the_rmsd_of_its_displacement() {
        let reference = site();
        let mut pose = reference.clone();
        for atom in pose
            .atoms
            .iter_mut()
            .filter(|atom| atom.residue_name == "LIG")
        {
            atom.position += Vec3::new(0.0, 0.0, 1.5);
        }

        assert!((rmsd_of(&pose, &reference).unwrap() - 1.5).abs() < 1e-4);
    }

    #[test]
    fn equivalent_atoms_trading_places_cost_nothing() {
        let reference = site();
        let mut pose = reference.clone();
        let [o1, o2] = ["O1", "O2"].map(|name| {
            pose.atoms
                .iter()
                .position(|atom| atom.name == name)
                .unwrap()
        });
        let (a, b) = (pose.atoms[o1].position, pose.atoms[o2].position);
        pose.atoms[o1].position = b;
        pose.atoms[o2].position = a;

        assert!(rmsd_of(&pose, &reference).unwrap() < 1e-4);
    }

    #[test]
    fn a_ligand_named_otherwise_pairs_by_its_bond_graph() {
        let reference = site();
        let mut pose = reference.clone();
        for atom in pose
            .atoms
            .iter_mut()
            .filter(|atom| atom.residue_name == "LIG")
        {
            atom.position += Vec3::new(0.0, 0.8, 0.0);
        }
        let named_alike = rmsd_of(&pose, &reference).unwrap();

        // Generated names in another order, as output from a SMILES input has.
        let ligand_start = pose.atoms.len() - 4;
        pose.atoms[ligand_start..].reverse();
        for (i, atom) in pose.atoms[ligand_start..].iter_mut().enumerate() {
            atom.name = FixedStr::new(&format!("X{i}"));
        }

        assert!((rmsd_of(&pose, &reference).unwrap() - named_alike).abs() < 1e-4);
        assert!((named_alike - 0.8).abs() < 1e-4);
    }

    #[test]
    fn a_pose_without_the_pocket_or_the_ligand_atoms_has_no_rmsd() {
        let reference = site();

        let mut unrelated = reference.clone();
        for atom in unrelated
            .atoms
            .iter_mut()
            .filter(|atom| atom.residue_name == "GLY")
        {
            atom.residue_name = FixedStr::new("TRP");
        }
        assert_eq!(rmsd_of(&unrelated, &reference), Err(NoRmsd::Pocket));

        let mut missing = reference.clone();
        missing.atoms.pop();
        for atom in missing
            .atoms
            .iter_mut()
            .filter(|atom| atom.residue_name == "LIG")
        {
            atom.name = FixedStr::new(&format!("Y{}", atom.name));
        }
        assert_eq!(
            rmsd_of(&missing, &reference),
            Err(NoRmsd::Atoms {
                pose: 3,
                reference: 4
            })
        );
    }

    /// The glycines of [`site`] and `copies` sulfates around them, each a
    /// sulfur and four oxygens 1.48 A out along the corners of a tetrahedron.
    fn sulfates(copies: usize) -> Protein {
        let mut structure = site();
        structure.atoms.retain(|atom| atom.residue_name == "GLY");
        let corners = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
        ];
        for copy in 0..copies {
            let centre = Vec3::new(4.5 * copy as f32 - 9.0, 0.0, 2.0 * (copy % 2) as f32);
            let seq = 901 + copy as i32;
            structure.atoms.push(atom(seq, "SO4", "S", "S", centre));
            for (corner, name) in corners.iter().zip(["O1", "O2", "O3", "O4"]) {
                let position = centre + corner.normalize() * 1.48;
                structure.atoms.push(atom(seq, "SO4", name, "O", position));
            }
        }
        structure
    }

    #[test]
    fn copies_in_another_order_with_equivalent_atoms_traded_cost_nothing() {
        let reference = sulfates(5);
        // The copies numbered in reverse and, in each, the oxygens named round
        // by one, as two co-folding samples order and name the same answer.
        let mut pose = reference.clone();
        for atom in pose
            .atoms
            .iter_mut()
            .filter(|atom| atom.residue_name == "SO4")
        {
            atom.residue_seq = 1806 - atom.residue_seq;
            let next = match atom.name.as_str() {
                "O1" => "O2",
                "O2" => "O3",
                "O3" => "O4",
                "O4" => "O1",
                other => other,
            };
            atom.name = FixedStr::new(next);
        }
        pose.atoms
            .sort_by_key(|atom| (atom.residue_name != "GLY", atom.residue_seq));

        assert!(rmsd_of(&pose, &reference).unwrap() < 1e-4);
    }

    #[test]
    fn a_broken_bond_among_many_copies_settles_quickly() {
        let reference = sulfates(6);
        let mut pose = reference.clone();
        let stretched = pose
            .atoms
            .iter()
            .position(|atom| atom.residue_seq == 903 && atom.name == "O2")
            .unwrap();
        let sulfur = pose.atoms[stretched - 2].position;
        let outward = (pose.atoms[stretched].position - sulfur).normalize();
        pose.atoms[stretched].position = sulfur + outward * 2.4;

        // Named as the reference, the pose is measured on the reference's bonds.
        let named = rmsd_of(&pose, &reference).unwrap();
        assert!(named > 0.0 && named < 0.5, "{named}");

        // Named otherwise, its bonds must match the reference's, and one is gone.
        for atom in pose
            .atoms
            .iter_mut()
            .filter(|atom| atom.residue_name == "SO4")
        {
            atom.name = FixedStr::new(&format!("X{}", atom.name));
        }
        let started = std::time::Instant::now();
        assert!(rmsd_of(&pose, &reference).is_err());
        assert!(started.elapsed().as_secs_f32() < 1.0);
    }

    #[test]
    fn atoms_named_alike_but_of_other_elements_pair_by_bond_graph() {
        let reference = site();
        // The carboxyl carbon and an oxygen trade names but not places.
        let mut pose = reference.clone();
        for atom in pose.atoms.iter_mut() {
            match atom.name.as_str() {
                "C1" if atom.residue_name == "LIG" => atom.name = FixedStr::new("O1"),
                "O1" => atom.name = FixedStr::new("C1"),
                _ => {}
            }
        }

        assert!(rmsd_of(&pose, &reference).unwrap() < 1e-4);
    }

    #[test]
    fn a_pocket_chain_without_a_counterpart_gives_way_to_another() {
        // Chain A lines most of the pocket, chain B the rest.
        let mut reference = site();
        for (seq, centre) in [
            (1, Vec3::new(0.0, -5.0, -2.0)),
            (2, Vec3::new(2.0, -4.0, 3.0)),
        ] {
            for (name, offset) in [
                ("N", Vec3::new(-1.2, 0.5, 0.0)),
                ("CA", Vec3::ZERO),
                ("C", Vec3::new(1.0, 0.8, 0.6)),
            ] {
                let element = &name[..1];
                reference.atoms.push(Atom {
                    chain_id: FixedStr::new("B"),
                    ..atom(seq, "ALA", name, element, centre + offset)
                });
            }
        }
        // A pose whose chain A is another protein altogether.
        let mut pose = reference.clone();
        for atom in pose
            .atoms
            .iter_mut()
            .filter(|atom| atom.residue_name == "GLY")
        {
            atom.residue_name = FixedStr::new("TRP");
        }

        assert!(rmsd_of(&pose, &reference).unwrap() < 1e-3);
    }

    #[test]
    fn a_pose_that_turns_a_double_bond_from_e_to_z_is_found_out() {
        let reference = alkene_site(false, CHLOROFLUOROETHENE);

        let same = assessment_of(&reference, &reference).unwrap();
        let turned = assessment_of(&alkene_site(true, CHLOROFLUOROETHENE), &reference).unwrap();

        assert!(same.flipped.is_empty());
        assert_eq!(
            same.geometry,
            Geometry {
                unflat: Some(0),
                ..Geometry::default()
            }
        );
        let ends = turned
            .flipped
            .iter()
            .map(|ends| ends.map(|(_, _, name)| name));
        assert_eq!(
            ends.collect::<Vec<_>>(),
            [[FixedStr::new("C1"), FixedStr::new("C2")]]
        );
        assert_eq!(
            turned.geometry, same.geometry,
            "its lengths and angles kept"
        );
    }

    #[test]
    fn without_definitions_neither_e_and_z_nor_flatness_is_known() {
        let reference = alkene_site(false, "");

        let turned = assessment_of(&alkene_site(true, ""), &reference).unwrap();

        assert!(turned.flipped.is_empty());
        assert_eq!(turned.geometry.unflat, None);
    }

    const PENTENEDIOL: &str = "\
data_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
LIG O1 O N
LIG C1 C N
LIG C2 C N
LIG C3 C N
LIG C4 C N
LIG C5 C N
LIG O5 O N
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG O1 C1 SING
LIG C1 C2 SING
LIG C2 C3 DOUB
LIG C3 C4 SING
LIG C4 C5 SING
LIG C5 O5 SING
";

    /// The glycines of [`site`] around (E)-pent-2-ene-1,5-diol laid flat,
    /// C5 and O5 turned `torsion` degrees about C3-C4, the whole turned end for
    /// end when `reversed`, its atoms named `names`. Its bond graph alone is
    /// symmetric end to end, so only bond orders tell C2=C3 from C3-C4.
    fn pentenediol(torsion: f32, reversed: bool, names: [&str; 7]) -> Protein {
        let mut structure = site();
        structure.atoms.retain(|atom| atom.residue_name != "LIG");
        let mut positions: Vec<Vec3> = (0..7)
            .map(|at| Vec3::new(at as f32 * 1.27 - 3.81, (at % 2) as f32 * 0.86 - 0.43, 0.0))
            .collect();
        let axis = (positions[4] - positions[3]).normalize();
        let turn = Quat::from_axis_angle(axis, torsion.to_radians());
        for at in 5..7 {
            positions[at] = positions[4] + turn * (positions[at] - positions[4]);
        }
        let centre = positions.iter().copied().sum::<Vec3>() / 7.0;
        let elements = ["O", "C", "C", "C", "C", "C", "O"];
        for ((name, element), position) in names.into_iter().zip(elements).zip(positions) {
            let offset = position - centre;
            let placed = match reversed {
                true => centre + Vec3::new(-offset.x, offset.y, -offset.z),
                false => position,
            };
            structure
                .atoms
                .push(atom(900, "LIG", name, element, placed));
        }
        structure.components = Dictionary::parse(PENTENEDIOL);

        structure
    }

    const NAMED: [&str; 7] = ["O1", "C1", "C2", "C3", "C4", "C5", "O5"];

    #[test]
    fn a_pose_turned_end_for_end_is_read_by_its_names() {
        let reference = pentenediol(0.0, false, NAMED);

        let turned = assessment_of(&pentenediol(150.0, true, NAMED), &reference).unwrap();

        assert!(turned.flipped.is_empty(), "{:?}", turned.flipped);
        assert_eq!(turned.geometry.unflat, Some(0));
    }

    #[test]
    fn a_double_bond_a_symmetry_could_trade_for_a_single_one_is_not_read_by_bond_graph() {
        let reference = pentenediol(0.0, false, NAMED);
        let renamed = ["X1", "X2", "X3", "X4", "X5", "X6", "X7"];

        let turned = assessment_of(&pentenediol(150.0, true, renamed), &reference).unwrap();

        assert!(turned.flipped.is_empty(), "{:?}", turned.flipped);
        assert_eq!(turned.geometry.unflat, Some(0));
    }

    /// The glycines of [`site`] around meso-2,5-dichlorohex-3-ene, C3 and C4
    /// joined by a bond of `order`: E, with a centre of inversion, or, when
    /// `z`, Z with a mirror plane, turned 180 degrees about z.
    fn dichlorohexene(z: bool, order: &str) -> Protein {
        let mut structure = site();
        structure.atoms.retain(|atom| atom.residue_name != "LIG");
        let c3 = Vec3::new(-0.67, 0.0, 0.0);
        let c2 = c3 + Vec3::new(-0.5, 0.866, 0.0) * 1.5;
        let corners = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
        ]
        .map(Vec3::normalize);
        let onto = Quat::from_rotation_arc(corners[0], (c3 - c2).normalize());
        let c1 = c2 + onto * corners[1] * 1.53;
        let cl2 = c2 + onto * corners[3] * 1.78;
        let other = |p: Vec3| if z { Vec3::new(-p.x, p.y, p.z) } else { -p };
        let turned = |p: Vec3| if z { Vec3::new(-p.x, -p.y, p.z) } else { p };
        for (name, element, position) in [
            ("C1", "C", c1),
            ("C2", "C", c2),
            ("CL2", "CL", cl2),
            ("C3", "C", c3),
            ("C4", "C", other(c3)),
            ("C5", "C", other(c2)),
            ("CL5", "CL", other(cl2)),
            ("C6", "C", other(c1)),
        ] {
            structure
                .atoms
                .push(atom(900, "LIG", name, element, turned(position)));
        }
        let definitions = format!(
            "data_LIG\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
             _chem_comp_atom.type_symbol\n_chem_comp_atom.pdbx_aromatic_flag\n\
             LIG C1 C N\nLIG C2 C N\nLIG CL2 CL N\nLIG C3 C N\nLIG C4 C N\n\
             LIG C5 C N\nLIG CL5 CL N\nLIG C6 C N\nloop_\n_chem_comp_bond.comp_id\n\
             _chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
             _chem_comp_bond.value_order\nLIG C1 C2 SING\nLIG C2 CL2 SING\n\
             LIG C2 C3 SING\nLIG C3 C4 {order}\nLIG C4 C5 SING\nLIG C5 CL5 SING\n\
             LIG C5 C6 SING\n"
        );
        structure.components = Dictionary::parse(&definitions);

        structure
    }

    #[test]
    fn a_double_bond_turned_about_leaves_a_meso_ligand_s_centres_alone() {
        for order in ["SING", "DOUB"] {
            let reference = dichlorohexene(false, order);

            let turned = assessment_of(&dichlorohexene(true, order), &reference).unwrap();

            assert_eq!(turned.inverted, [], "C3-C4 {order}");
            let flipped: Vec<[&str; 2]> = turned
                .flipped
                .iter()
                .map(|ends| ends.each_ref().map(|(_, _, name)| name.as_str()))
                .collect();
            let expected: &[[&str; 2]] = if order == "DOUB" {
                &[["C3", "C4"]]
            } else {
                &[]
            };
            assert_eq!(flipped, expected, "C3-C4 {order}");
        }
    }

    #[test]
    fn a_pose_that_inverts_a_stereocentre_is_found_out() {
        let reference = chiral_site();
        assert_eq!(fit_of(&reference, &reference).unwrap().inverted, []);

        // The nitrogen and the oxygen trade places, which turns C1 the other
        // way.
        let mut pose = reference.clone();
        let (n1, o1) = (named(&pose, "N1"), named(&pose, "O1"));
        let (a, b) = (pose.atoms[n1].position, pose.atoms[o1].position);
        pose.atoms[n1].position = b;
        pose.atoms[o1].position = a;

        let c1 = (
            (FixedStr::new("A"), 900, None),
            FixedStr::new("LIG"),
            FixedStr::new("C1"),
        );
        assert_eq!(fit_of(&pose, &reference).unwrap().inverted, [c1]);
    }

    #[test]
    fn a_mirror_image_named_otherwise_is_found_out_through_its_bond_graph() {
        let reference = chiral_site();
        let mut pose = reference.clone();
        let ligand_start = pose.atoms.len() - 5;
        for (i, atom) in pose.atoms[ligand_start..].iter_mut().enumerate() {
            atom.name = FixedStr::new(&format!("X{i}"));
            atom.position.x = -atom.position.x;
        }

        let inverted = fit_of(&pose, &reference).unwrap().inverted;
        let names: Vec<&str> = inverted.iter().map(|(_, _, name)| name.as_str()).collect();
        assert_eq!(names, ["C1"]);
    }

    #[test]
    fn a_meso_ligand_turned_about_inverts_nothing() {
        let reference = meso_site(MESO_CHLORINE);
        let mut turned = reference.clone();
        let ligand_start = turned.atoms.len() - 6;
        for atom in &mut turned.atoms[ligand_start..] {
            atom.position = Vec3::new(-atom.position.x, -atom.position.y, atom.position.z);
        }

        let fit = fit_of(&turned, &reference).unwrap();
        assert!(fit.rmsd > 1.0, "{fit:?}");
        assert_eq!(fit.inverted, []);
    }

    #[test]
    fn one_centre_of_a_meso_ligand_inverted_is_still_found_out() {
        let reference = meso_site(MESO_CHLORINE);
        let pose = meso_site(CHIRAL_CHLORINE);

        let inverted = fit_of(&pose, &reference).unwrap().inverted;
        let names: Vec<&str> = inverted.iter().map(|(_, _, name)| name.as_str()).collect();
        assert_eq!(names, ["C2"]);
    }
}
