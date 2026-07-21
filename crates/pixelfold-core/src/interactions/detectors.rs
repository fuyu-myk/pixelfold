//! Per-type detectors. Each one narrows the atom set to the roles it needs,
//! grids those positions, and only tests candidates inside the cutoff.

use std::collections::HashMap;
use std::hash::Hash;

use glam::Vec3;

use crate::fixed_str::FixedStr;
use crate::spatial::Grid;
use crate::structure::{Atom, Protein};

use super::classify::{
    Charge, charged_group, is_apolar_carbon, is_cysteine_sulfur, is_hbond_acceptor, is_metal,
    is_metal_donor, is_water,
};
use super::params::{
    DISULFIDE_DIST_MAX, HBOND_DIST_MAX, HBOND_DONOR_ANGLE_MIN, HYDROPHOBIC_DIST_MAX,
    METAL_DIST_MAX, MIN_DIST, PICATION_DIST_MAX, PISTACK_ANGLE_DEV, PISTACK_DIST_MAX,
    PISTACK_OFFSET_MAX, SALTBRIDGE_DIST_MAX, WATER_BRIDGE_MAX_DIST, WATER_BRIDGE_MIN_DIST,
    WATER_BRIDGE_OMEGA_MAX, WATER_BRIDGE_OMEGA_MIN, WATER_BRIDGE_THETA_MIN,
};
use super::{Interaction, InteractionKind, aromatic, hydrogens};

type ResidueKey = (FixedStr<4>, u32, Option<char>);

fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

/// Indices of the atoms satisfying `role`, with their positions.
fn collect(protein: &Protein, role: impl Fn(&Atom) -> bool) -> (Vec<usize>, Vec<Vec3>) {
    let indices: Vec<usize> = protein
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| role(atom))
        .map(|(index, _)| index)
        .collect();
    let positions = indices.iter().map(|&i| protein.atoms[i].position).collect();

    (indices, positions)
}

/// Cysteine SG-SG disulfide bridges. Detection is distance-only; the CB-SG-SG-CB
/// dihedral clusters near 90 degrees in real bridges but is not a gate.
pub fn disulfides(protein: &Protein) -> Vec<Interaction> {
    let (sulfurs, positions) = collect(protein, is_cysteine_sulfur);
    if sulfurs.len() < 2 {
        return Vec::new();
    }

    let grid = Grid::build(&positions, DISULFIDE_DIST_MAX);

    let mut out = Vec::new();
    for (a, &atom_a) in sulfurs.iter().enumerate() {
        let center = positions[a];
        grid.for_each_within(center, DISULFIDE_DIST_MAX, |b| {
            if b <= a {
                return; // report each pair once, never a sulfur against itself
            }

            let distance = (center - positions[b]).length();
            if distance <= MIN_DIST {
                return;
            }

            out.push(Interaction {
                kind: InteractionKind::Disulfide,
                atoms_a: vec![atom_a],
                atoms_b: vec![sulfurs[b]],
                distance,
                angle: None,
                bridge: None,
            });
        });
    }

    out
}

/// Metal ions coordinated by an O, N, or S donor.
pub fn metal_coordination(protein: &Protein) -> Vec<Interaction> {
    let (metals, metal_positions) = collect(protein, is_metal);
    if metals.is_empty() {
        return Vec::new();
    }

    let (donors, donor_positions) = collect(protein, is_metal_donor);
    if donors.is_empty() {
        return Vec::new();
    }

    let grid = Grid::build(&donor_positions, METAL_DIST_MAX);

    let mut out = Vec::new();
    for (m, &metal) in metals.iter().enumerate() {
        let center = metal_positions[m];
        grid.for_each_within(center, METAL_DIST_MAX, |d| {
            let distance = (center - donor_positions[d]).length();
            if !(MIN_DIST < distance && distance < METAL_DIST_MAX) {
                return;
            }

            out.push(Interaction {
                kind: InteractionKind::MetalCoordination,
                atoms_a: vec![metal],
                atoms_b: vec![donors[d]],
                distance,
                angle: None,
                bridge: None,
            });
        });
    }

    out
}

/// One hydrogen bond that passed both gates.
struct Bond {
    donor: usize,
    acceptor: usize,
    distance: f32,
    angle: f32,
}

/// An acceptor inside a donor's cutoff, with the angle each of the donor's
/// hydrogens makes to it.
struct Reachable {
    acceptor: usize,
    distance: f32,
    angles: Vec<f32>,
}

/// Give each of a donor's hydrogens its straightest acceptor.
///
/// Taking the best bond, spending the hydrogen that makes it, and repeating is
/// what stops one hydrogen from being counted twice. Capping the number of bonds
/// at the number of hydrogens is not enough on its own: both hydrogens of an
/// amide can be the best route to the same oxygen.
fn assign(donor: &hydrogens::Donor, reachable: Vec<Reachable>) -> Vec<Bond> {
    let mut spent = vec![false; donor.hydrogens.capacity()];
    let mut taken = vec![false; reachable.len()];
    let mut bonds = Vec::new();

    loop {
        let mut best: Option<(usize, usize, f32)> = None;
        for (r, candidate) in reachable.iter().enumerate() {
            if taken[r] {
                continue;
            }

            for (h, &angle) in candidate.angles.iter().enumerate() {
                if spent[h] || angle <= HBOND_DONOR_ANGLE_MIN {
                    continue;
                }

                if best.is_none_or(|(_, _, straightest)| angle > straightest) {
                    best = Some((r, h, angle));
                }
            }
        }

        let Some((r, h, angle)) = best else {
            return bonds;
        };

        spent[h] = true;
        taken[r] = true;
        bonds.push(Bond {
            donor: donor.atom,
            acceptor: reachable[r].acceptor,
            distance: reachable[r].distance,
            angle,
        });
    }
}

/// Hydrogen bonds, on PLIP's geometry: a donor-acceptor heavy-atom distance
/// inside the cutoff and a D-H...A angle at the hydrogen above the minimum.
///
/// The hydrogens come from [`super::hydrogens`], which places a rotatable
/// donor's hydrogen against the acceptor under test rather than at an arbitrary
/// torsion. Each hydrogen then makes at most one bond, straightest first: a
/// donor's hydrogens are *assigned* to acceptors rather than merely counted,
/// because an amide can easily point both of its hydrogens at the same
/// carboxylate and reach one oxygen twice while the other hydrogen reaches
/// nothing.
///
/// Three differences from PLIP are deliberate:
///
/// - PLIP caps a donor at **one** bond however many hydrogens it has, so a
///   lysine ammonium or an arginine reports a single bond at an interface where
///   it makes three. That cap is safe for a small ligand and wrong for a
///   protein, which is what this engine is for. For every single-hydrogen donor,
///   the majority, the two rules agree exactly.
/// - PLIP drops a hydrogen bond whose two ends already form a salt bridge. RING
///   does not, and lets a residue pair carry both edges at once, which says more
///   in a network. Neither does this.
/// - PLIP has a guard its comment calls a backbone-backbone filter. The property
///   it reads is OpenBabel's `SIDECHAIN`, so it actually drops side-chain to
///   side-chain bonds, and only in intra mode. It is not reproduced.
///
/// These are potential hydrogen bonds, in HBPLUS's sense: where a tautomer or a
/// torsion is undetermined by the coordinates, what is reported is what the
/// geometry allows, and some of what it allows is mutually exclusive.
pub fn hydrogen_bonds(protein: &Protein) -> Vec<Interaction> {
    let donors = hydrogens::donors(protein);
    let (acceptors, acceptor_positions) = collect(protein, is_hbond_acceptor);
    if donors.is_empty() || acceptors.is_empty() {
        return Vec::new();
    }

    let grid = Grid::build(&acceptor_positions, HBOND_DIST_MAX);

    let mut bonds = Vec::new();
    for donor in &donors {
        let origin = protein.atoms[donor.atom].position;
        let mut reachable = Vec::new();

        grid.for_each_within(origin, HBOND_DIST_MAX, |i| {
            let acceptor = acceptors[i];
            if residue_key(&protein.atoms[donor.atom]) == residue_key(&protein.atoms[acceptor]) {
                return; // a residue does not bond to itself, nor an atom to itself
            }

            let target = acceptor_positions[i];
            let distance = (target - origin).length();
            if !(MIN_DIST < distance && distance < HBOND_DIST_MAX) {
                return;
            }

            let angles = donor.hydrogens.angles_toward(origin, target);
            if angles.iter().all(|&angle| angle <= HBOND_DONOR_ANGLE_MIN) {
                return; // no hydrogen of this donor can reach it
            }

            reachable.push(Reachable {
                acceptor,
                distance,
                angles,
            });
        });

        bonds.extend(assign(donor, reachable));
    }

    bonds.sort_unstable_by_key(|bond| (bond.donor, bond.acceptor));
    bonds
        .into_iter()
        .map(|bond| Interaction {
            kind: InteractionKind::HydrogenBond,
            atoms_a: vec![bond.donor],
            atoms_b: vec![bond.acceptor],
            distance: bond.distance,
            angle: Some(bond.angle),
            bridge: None,
        })
        .collect()
}

/// A water oxygen, the only atom that bridges two residues in this engine.
fn is_water_oxygen(atom: &Atom) -> bool {
    is_water(atom.residue_name.as_str()) && atom.element.as_str().eq_ignore_ascii_case("O")
}

/// Water bridges: an ordered water linking a donor and an acceptor that do not
/// reach each other directly.
///
/// Both legs must sit in the same distance band, the donor must point a hydrogen
/// at the water, and the two legs must meet at the water within an angular
/// window. That last gate is what separates a bridge from a water that merely
/// happens to sit between two residues: the geometry has to look like two
/// hydrogen bonds through one oxygen, not a coincidence of proximity.
///
/// Two differences from PLIP. This reports one bridge per donor instead of multiple,
/// through whichever hydrogen reaches the water straightest. Additionally, this looks
/// between any two residues over looking only for bridges between its ligand and its
/// binding site.
pub fn water_bridges(protein: &Protein) -> Vec<Interaction> {
    let (waters, water_positions) = collect(protein, is_water_oxygen);
    let (acceptors, acceptor_positions) = collect(protein, is_hbond_acceptor);
    let donors = hydrogens::donors(protein);
    if waters.is_empty() || acceptors.is_empty() || donors.is_empty() {
        return Vec::new();
    }

    let grid = Grid::build(&water_positions, WATER_BRIDGE_MAX_DIST);
    let band = WATER_BRIDGE_MIN_DIST..=WATER_BRIDGE_MAX_DIST;

    // The acceptor legs of each water, so the donor pass can join on them.
    let mut legs: HashMap<usize, Vec<(usize, f32)>> = HashMap::new();
    for (a, &acceptor) in acceptors.iter().enumerate() {
        let position = acceptor_positions[a];
        grid.for_each_within(position, WATER_BRIDGE_MAX_DIST, |w| {
            let distance = (position - water_positions[w]).length();
            if band.contains(&distance) {
                legs.entry(w).or_default().push((acceptor, distance));
            }
        });
    }

    let mut out = Vec::new();
    for donor in &donors {
        let origin = protein.atoms[donor.atom].position;
        let donor_residue = residue_key(&protein.atoms[donor.atom]);

        grid.for_each_within(origin, WATER_BRIDGE_MAX_DIST, |w| {
            let Some(reachable) = legs.get(&w) else {
                return;
            };

            let water = water_positions[w];
            let donor_leg = (origin - water).length();
            if !band.contains(&donor_leg) {
                return;
            }

            // The hydrogen that points most nearly at the water carries the leg.
            let Some((hydrogen, donor_angle)) = donor
                .hydrogens
                .positions_toward(water)
                .into_iter()
                .map(|h| (h, hydrogens::angle_at(h, origin, water)))
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
            else {
                return;
            };
            if donor_angle <= WATER_BRIDGE_THETA_MIN {
                return;
            }

            for &(acceptor, acceptor_leg) in reachable {
                if residue_key(&protein.atoms[acceptor]) == donor_residue {
                    continue; // a residue does not bridge to itself
                }

                let water_angle =
                    hydrogens::angle_at(water, protein.atoms[acceptor].position, hydrogen);
                if !(WATER_BRIDGE_OMEGA_MIN < water_angle && water_angle < WATER_BRIDGE_OMEGA_MAX) {
                    continue;
                }

                out.push(Interaction {
                    kind: InteractionKind::WaterBridge,
                    atoms_a: vec![donor.atom],
                    atoms_b: vec![acceptor],
                    // The longer leg, the one that limits the bridge.
                    distance: donor_leg.max(acceptor_leg),
                    angle: Some(water_angle),
                    bridge: Some(waters[w]),
                });
            }
        });
    }

    out
}

/// One surviving apolar carbon pair.
struct Contact {
    a: usize,
    b: usize,
    distance: f32,
}

/// Hydrophobic contacts between apolar carbons.
///
/// Every apolar carbon of a residue is close to every apolar carbon of the
/// residue packed against it, so the raw pairs overstate one touch many times
/// over. PLIP reduces them (`PLInteraction.refine_hydrophobic`); the reduction
/// used here is the one PLIP itself switches to when its partner is a peptide
/// rather than a small molecule, which is the case that matches a protein.
pub fn hydrophobic_contacts(protein: &Protein) -> Vec<Interaction> {
    let (carbons, positions) = collect(protein, is_apolar_carbon);
    if carbons.len() < 2 {
        return Vec::new();
    }

    let grid = Grid::build(&positions, HYDROPHOBIC_DIST_MAX);

    let mut contacts = Vec::new();
    for (a, &atom_a) in carbons.iter().enumerate() {
        let center = positions[a];
        grid.for_each_within(center, HYDROPHOBIC_DIST_MAX, |b| {
            if b <= a {
                return; // report each pair once, never a carbon against itself
            }

            let atom_b = carbons[b];
            if residue_key(&protein.atoms[atom_a]) == residue_key(&protein.atoms[atom_b]) {
                return; // a residue does not contact itself
            }

            let distance = (center - positions[b]).length();
            if !(MIN_DIST < distance && distance < HYDROPHOBIC_DIST_MAX) {
                return;
            }

            contacts.push(Contact {
                a: atom_a,
                b: atom_b,
                distance,
            });
        });
    }

    reduce_to_closest(protein, contacts)
        .into_iter()
        .map(|contact| Interaction {
            kind: InteractionKind::Hydrophobic,
            atoms_a: vec![contact.a],
            atoms_b: vec![contact.b],
            distance: contact.distance,
            angle: None,
            bridge: None,
        })
        .collect()
}

/// Reduce each side in turn: first the closest carbon of each residue reaching a
/// given atom, then the same with the sides swapped. One pass alone would only
/// thin one direction, leaving the other side's carbons all reported.
fn reduce_to_closest(protein: &Protein, contacts: Vec<Contact>) -> Vec<Contact> {
    let once = keep_shortest(contacts, |c| (c.a, residue_key(&protein.atoms[c.b])));
    let mut twice = keep_shortest(once, |c| (residue_key(&protein.atoms[c.a]), c.b));
    twice.sort_unstable_by_key(|c| (c.a, c.b));

    twice
}

/// Keep the shortest contact per key. Ties keep the contact seen first, as PLIP
/// does by comparing strictly.
fn keep_shortest<K: Eq + Hash>(
    contacts: Vec<Contact>,
    key: impl Fn(&Contact) -> K,
) -> Vec<Contact> {
    let mut best: HashMap<K, Contact> = HashMap::new();
    for contact in contacts {
        match best.entry(key(&contact)) {
            std::collections::hash_map::Entry::Occupied(mut slot) => {
                if slot.get().distance > contact.distance {
                    slot.insert(contact);
                }
            }
            std::collections::hash_map::Entry::Vacant(slot) => {
                slot.insert(contact);
            }
        }
    }

    best.into_values().collect()
}

/// A residue's charged group: the centroid of its charged atoms.
struct ChargedGroup {
    sign: Charge,
    centroid: Vec3,
    atoms: Vec<usize>,
}

/// One charged group per charged residue, built from the atoms actually present
/// (a residue missing its carboxylate oxygens contributes nothing).
fn charged_groups(protein: &Protein) -> Vec<ChargedGroup> {
    let mut groups = Vec::new();
    let mut index = 0;

    while index < protein.atoms.len() {
        let start = index;
        let key = residue_key(&protein.atoms[start]);
        while index < protein.atoms.len() && residue_key(&protein.atoms[index]) == key {
            index += 1;
        }

        let Some((sign, names)) = charged_group(protein.atoms[start].residue_name.as_str()) else {
            continue;
        };
        let atoms: Vec<usize> = (start..index)
            .filter(|&i| names.contains(&protein.atoms[i].name.as_str()))
            .collect();
        if atoms.is_empty() {
            continue;
        }

        let sum = atoms
            .iter()
            .fold(Vec3::ZERO, |acc, &i| acc + protein.atoms[i].position);
        groups.push(ChargedGroup {
            sign,
            centroid: sum / atoms.len() as f32,
            atoms,
        });
    }

    groups
}

/// Salt bridges between oppositely charged groups, measured centroid to
/// centroid.
pub fn salt_bridges(protein: &Protein) -> Vec<Interaction> {
    let groups = charged_groups(protein);
    let (positive, negative): (Vec<&ChargedGroup>, Vec<&ChargedGroup>) = groups
        .iter()
        .partition(|group| group.sign == Charge::Positive);
    if positive.is_empty() || negative.is_empty() {
        return Vec::new();
    }

    let negative_centroids: Vec<Vec3> = negative.iter().map(|group| group.centroid).collect();
    let grid = Grid::build(&negative_centroids, SALTBRIDGE_DIST_MAX);

    let mut out = Vec::new();
    for cation in &positive {
        grid.for_each_within(cation.centroid, SALTBRIDGE_DIST_MAX, |n| {
            let anion = negative[n];
            let distance = (cation.centroid - anion.centroid).length();
            if !(MIN_DIST < distance && distance < SALTBRIDGE_DIST_MAX) {
                return;
            }

            out.push(Interaction {
                kind: InteractionKind::SaltBridge,
                atoms_a: cation.atoms.clone(),
                atoms_b: anion.atoms.clone(),
                distance,
                angle: None,
                bridge: None,
            });
        });
    }

    out
}

/// Pi-stacking between aromatic rings, on PLIP's geometry.
///
/// A pair passes when the ring centroids are within the cutoff, the offset (how
/// far one centroid sits from over the other's centre) is small, and the planes
/// are either near parallel or near perpendicular.
///
/// A tryptophan's two rings are close and near coplanar, so the same-residue
/// guard is what stops a residue from stacking with itself.
pub fn pi_stacking(protein: &Protein) -> Vec<Interaction> {
    let rings = aromatic::rings(protein);
    if rings.len() < 2 {
        return Vec::new();
    }

    let centers: Vec<Vec3> = rings.iter().map(|ring| ring.center).collect();
    let grid = Grid::build(&centers, PISTACK_DIST_MAX);

    let mut out = Vec::new();
    for (a, ring_a) in rings.iter().enumerate() {
        grid.for_each_within(ring_a.center, PISTACK_DIST_MAX, |b| {
            if b <= a {
                return; // report each pair once
            }

            let ring_b = &rings[b];
            if ring_a.residue == ring_b.residue {
                return; // a residue's own rings do not stack
            }

            let distance = (ring_a.center - ring_b.center).length();
            if !(MIN_DIST < distance && distance < PISTACK_DIST_MAX) {
                return;
            }

            let offset = ring_a
                .planar_offset(ring_b.center)
                .min(ring_b.planar_offset(ring_a.center));
            if offset >= PISTACK_OFFSET_MAX {
                return;
            }

            // PLIP guards the parallel case with `0 < angle`, rejecting only
            // bit-identical normals.
            let angle = aromatic::plane_angle(ring_a.normal, ring_b.normal);
            let parallel = angle < PISTACK_ANGLE_DEV;
            let perpendicular =
                90.0 - PISTACK_ANGLE_DEV < angle && angle < 90.0 + PISTACK_ANGLE_DEV;
            if !(parallel || perpendicular) {
                return;
            }

            out.push(Interaction {
                kind: InteractionKind::PiStacking,
                atoms_a: ring_a.atoms.clone(),
                atoms_b: ring_b.atoms.clone(),
                distance,
                angle: Some(angle),
                bridge: None,
            });
        });
    }

    out
}

/// Pi-cation between an aromatic ring and a positive charge centre.
///
/// The cations are the same groups salt bridges use: the arginine guanidinium,
/// the lysine ammonium, and the histidine imidazole. PLIP restricts pi-cation to
/// exactly these three residues and, on the protein path, brings in no metals,
/// so neither does this. A histidine is both a ring and a cation, so the
/// same-residue guard stops its ring from counting its own charge.
pub fn pi_cation(protein: &Protein) -> Vec<Interaction> {
    let rings = aromatic::rings(protein);
    let cations: Vec<ChargedGroup> = charged_groups(protein)
        .into_iter()
        .filter(|group| group.sign == Charge::Positive)
        .collect();
    if rings.is_empty() || cations.is_empty() {
        return Vec::new();
    }

    let centers: Vec<Vec3> = rings.iter().map(|ring| ring.center).collect();
    let grid = Grid::build(&centers, PICATION_DIST_MAX);

    let mut out = Vec::new();
    for cation in &cations {
        let cation_residue = residue_key(&protein.atoms[cation.atoms[0]]);
        grid.for_each_within(cation.centroid, PICATION_DIST_MAX, |r| {
            let ring = &rings[r];
            if ring.residue == cation_residue {
                return; // a histidine ring does not count its own charge
            }

            let distance = (ring.center - cation.centroid).length();
            if !(MIN_DIST < distance && distance < PICATION_DIST_MAX) {
                return;
            }
            if ring.planar_offset(cation.centroid) >= PISTACK_OFFSET_MAX {
                return;
            }

            out.push(Interaction {
                kind: InteractionKind::PiCation,
                atoms_a: ring.atoms.clone(),
                atoms_b: cation.atoms.clone(),
                distance,
                angle: None,
                bridge: None,
            });
        });
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::SecondaryStructure;

    fn atom(residue: &str, seq: u32, name: &str, element: &str, position: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(element),
            residue_name: FixedStr::new(residue),
            residue_seq: seq,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position,
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    fn protein(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
        }
    }

    #[test]
    fn disulfide_found_at_bond_length_and_rejected_when_too_far() {
        // Two Cys gamma sulfurs at the nominal 2.05 A S-S bond length.
        let bonded = protein(vec![
            atom("CYS", 1, "SG", "S", Vec3::ZERO),
            atom("CYS", 2, "SG", "S", Vec3::new(2.05, 0.0, 0.0)),
        ]);
        let found = disulfides(&bonded);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::Disulfide);
        assert!((found[0].distance - 2.05).abs() < 1e-4);

        // Not a bridge beyond the 2.5 A cutoff.
        let apart = protein(vec![
            atom("CYS", 1, "SG", "S", Vec3::ZERO),
            atom("CYS", 2, "SG", "S", Vec3::new(3.0, 0.0, 0.0)),
        ]);
        assert!(disulfides(&apart).is_empty());
    }

    #[test]
    fn disulfide_reports_each_pair_once() {
        let structure = protein(vec![
            atom("CYS", 1, "SG", "S", Vec3::ZERO),
            atom("CYS", 2, "SG", "S", Vec3::new(2.05, 0.0, 0.0)),
        ]);
        assert_eq!(disulfides(&structure).len(), 1); // not two, mirrored
    }

    #[test]
    fn metal_coordinated_by_donor_but_not_by_carbon() {
        // Zinc with an Asp carboxylate oxygen in range, plus an apolar carbon
        // equally close that must not count.
        let structure = protein(vec![
            atom("ZN", 1, "ZN", "Zn", Vec3::ZERO),
            atom("ASP", 2, "OD1", "O", Vec3::new(2.1, 0.0, 0.0)),
            atom("ALA", 3, "CB", "C", Vec3::new(0.0, 2.1, 0.0)),
        ]);
        let found = metal_coordination(&structure);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::MetalCoordination);
        assert_eq!(found[0].atoms_b, vec![1]); // the Asp oxygen

        // Not coordinated beyond 3.0 A.
        let far = protein(vec![
            atom("ZN", 1, "ZN", "Zn", Vec3::ZERO),
            atom("ASP", 2, "OD1", "O", Vec3::new(4.0, 0.0, 0.0)),
        ]);
        assert!(metal_coordination(&far).is_empty());
    }

    #[test]
    fn calcium_ion_coordinates_but_c_alpha_never_does() {
        // Both are named "CA"; only the element distinguishes them.
        let ion = protein(vec![
            atom("CA", 1, "CA", "Ca", Vec3::ZERO),
            atom("GLU", 2, "OE1", "O", Vec3::new(2.4, 0.0, 0.0)),
        ]);
        assert_eq!(metal_coordination(&ion).len(), 1);

        let backbone = protein(vec![
            atom("ALA", 1, "CA", "C", Vec3::ZERO),
            atom("GLU", 2, "OE1", "O", Vec3::new(2.4, 0.0, 0.0)),
        ]);
        assert!(metal_coordination(&backbone).is_empty());
    }

    #[test]
    fn salt_bridge_between_carboxylate_and_guanidinium_centroids() {
        // Asp carboxylate centroid at x=0, Arg guanidinium centroid at x=4.
        let structure = protein(vec![
            atom("ASP", 1, "OD1", "O", Vec3::new(0.0, 1.0, 0.0)),
            atom("ASP", 1, "OD2", "O", Vec3::new(0.0, -1.0, 0.0)),
            atom("ARG", 2, "NE", "N", Vec3::new(4.0, 1.0, 0.0)),
            atom("ARG", 2, "NH1", "N", Vec3::new(4.0, 0.0, 0.0)),
            atom("ARG", 2, "NH2", "N", Vec3::new(4.0, -1.0, 0.0)),
        ]);
        let found = salt_bridges(&structure);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::SaltBridge);
        assert!((found[0].distance - 4.0).abs() < 1e-4); // centroid separation
        assert_eq!(found[0].atoms_a.len(), 3); // whole guanidinium
        assert_eq!(found[0].atoms_b.len(), 2); // whole carboxylate
    }

    /// Give an atom an explicit residue number so contacts can span residues.
    fn at(residue: &str, seq: u32, name: &str, position: Vec3) -> Atom {
        atom(residue, seq, name, "C", position)
    }

    /// Two residues, the first donating its amide hydrogen towards a carbonyl
    /// oxygen placed by the caller. Residue 1's C=O fixes residue 2's amide H.
    fn amide_donor_to(acceptor: Vec3) -> Protein {
        // Residue 1 is a proline: it supplies the carbonyl that fixes residue
        // 2's amide hydrogen without donating a terminal ammonium of its own.
        protein(vec![
            atom("PRO", 1, "N", "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", "C", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("ALA", 2, "N", "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("ALA", 2, "CA", "C", Vec3::new(4.0, 2.8, 0.0)),
            atom("GLY", 9, "O", "O", acceptor),
        ])
    }

    #[test]
    fn hydrogen_bond_found_when_the_amide_points_at_the_acceptor() {
        // Residue 2's amide H sits at N + unit(C_prev - O_prev); put an acceptor
        // out along that direction, at a credible donor-acceptor distance.
        let n = Vec3::new(3.3, 1.5, 0.0);
        let along = (Vec3::new(2.0, 1.4, 0.0) - Vec3::new(1.2, 2.3, 0.0)).normalize();
        let structure = amide_donor_to(n + along * 3.0);

        let found: Vec<_> = hydrogen_bonds(&structure);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::HydrogenBond);
        assert_eq!(found[0].atoms_a, vec![4]); // the amide nitrogen
        assert_eq!(found[0].atoms_b, vec![6]); // the carbonyl oxygen
        assert!((found[0].distance - 3.0).abs() < 1e-4);
        assert!(found[0].angle.unwrap() > 179.0, "expected a straight bond");
    }

    #[test]
    fn hydrogen_bond_rejected_beyond_the_cutoff_and_at_a_bad_angle() {
        let n = Vec3::new(3.3, 1.5, 0.0);
        let along = (Vec3::new(2.0, 1.4, 0.0) - Vec3::new(1.2, 2.3, 0.0)).normalize();

        // Straight on, but too far away.
        assert!(hydrogen_bonds(&amide_donor_to(n + along * 5.0)).is_empty());

        // Close enough, but sitting behind the donor where the hydrogen is not.
        assert!(hydrogen_bonds(&amide_donor_to(n - along * 3.0)).is_empty());
    }

    #[test]
    fn a_donor_makes_no_more_bonds_than_it_has_hydrogens() {
        // One amide hydrogen, three acceptors all within reach: the straightest
        // bond wins and the other two are dropped.
        let n = Vec3::new(3.3, 1.5, 0.0);
        let along = (Vec3::new(2.0, 1.4, 0.0) - Vec3::new(1.2, 2.3, 0.0)).normalize();

        let mut atoms = amide_donor_to(n + along * 3.0).atoms;
        atoms.push(atom(
            "GLY",
            10,
            "O",
            "O",
            n + along * 3.0 + Vec3::new(0.9, 0.0, 0.0),
        ));
        atoms.push(atom(
            "GLY",
            11,
            "O",
            "O",
            n + along * 3.0 + Vec3::new(0.0, 0.9, 0.0),
        ));
        let structure = protein(atoms);

        let found = hydrogen_bonds(&structure);
        assert_eq!(found.len(), 1, "one hydrogen cannot make three bonds");
        assert_eq!(found[0].atoms_b, vec![6], "the straightest acceptor wins");
    }

    /// An amide's two hydrogens are fixed and often point much the same way, so
    /// a carboxylate offering two oxygens can be reachable twice over by one of
    /// them and not at all by the other.
    #[test]
    fn two_hydrogens_cannot_both_be_the_same_hydrogen() {
        // An asparagine amide, with both oxygens of a carboxylate placed out
        // along the direction of a single one of its hydrogens.
        let nd2 = Vec3::new(0.0, 0.0, 0.0);
        let structure = protein(vec![
            atom("ASN", 1, "CB", "C", Vec3::new(-2.5, -1.2, 0.0)),
            atom("ASN", 1, "CG", "C", Vec3::new(-1.5, 0.0, 0.0)),
            atom("ASN", 1, "OD1", "O", Vec3::new(-2.0, 1.2, 0.0)),
            atom("ASN", 1, "ND2", "N", nd2),
            atom("ASP", 2, "OD1", "O", Vec3::new(2.4, -1.9, 0.0)),
            atom("ASP", 2, "OD2", "O", Vec3::new(1.4, -2.6, 0.0)),
        ]);

        let donor = &hydrogens::donors(&structure)[0];
        let hydrogens::Hydrogens::Fixed(hs) = &donor.hydrogens else {
            panic!("an amide is fixed");
        };
        assert_eq!(hs.len(), 2);

        // Both acceptors are reachable, and only ever by the same hydrogen.
        for acceptor in [Vec3::new(2.4, -1.9, 0.0), Vec3::new(1.4, -2.6, 0.0)] {
            let angles = donor.hydrogens.angles_toward(nd2, acceptor);
            let reaching: Vec<usize> = angles
                .iter()
                .enumerate()
                .filter(|&(_, &angle)| angle > HBOND_DONOR_ANGLE_MIN)
                .map(|(h, _)| h)
                .collect();
            assert_eq!(reaching, vec![0], "only one hydrogen reaches: {angles:?}");
        }

        let found = hydrogen_bonds(&structure);
        assert_eq!(found.len(), 1, "one hydrogen made two bonds: {found:?}");
    }

    #[test]
    fn only_typed_acceptors_accept() {
        let n = Vec3::new(3.3, 1.5, 0.0);
        let along = (Vec3::new(2.0, 1.4, 0.0) - Vec3::new(1.2, 2.3, 0.0)).normalize();
        let target = n + along * 3.0;

        // A lysine ammonium is a donor, never an acceptor; a cysteine sulfur is
        // neither, under the typing PLIP inherits from OpenBabel.
        let mut atoms = amide_donor_to(target).atoms;
        atoms.pop();
        atoms.push(atom("LYS", 9, "NZ", "N", target));
        assert!(hydrogen_bonds(&protein(atoms)).is_empty());

        let mut atoms = amide_donor_to(target).atoms;
        atoms.pop();
        atoms.push(atom("CYS", 9, "SG", "S", target));
        assert!(hydrogen_bonds(&protein(atoms)).is_empty());
    }

    #[test]
    fn hydrophobic_contact_found_in_range_and_rejected_beyond_it() {
        let close = protein(vec![
            at("ALA", 1, "CB", Vec3::ZERO),
            at("LEU", 2, "CD1", Vec3::new(3.5, 0.0, 0.0)),
        ]);
        let found = hydrophobic_contacts(&close);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::Hydrophobic);
        assert!((found[0].distance - 3.5).abs() < 1e-4);

        let apart = protein(vec![
            at("ALA", 1, "CB", Vec3::ZERO),
            at("LEU", 2, "CD1", Vec3::new(4.5, 0.0, 0.0)),
        ]);
        assert!(hydrophobic_contacts(&apart).is_empty());
    }

    #[test]
    fn only_apolar_carbons_contact() {
        // A serine CB reaches its hydroxyl oxygen, and a methionine CG its
        // sulfur, so neither is apolar however close it sits.
        let structure = protein(vec![
            at("ALA", 1, "CB", Vec3::ZERO),
            at("SER", 2, "CB", Vec3::new(3.0, 0.0, 0.0)),
            at("MET", 3, "CG", Vec3::new(0.0, 3.0, 0.0)),
            atom("MET", 4, "SD", "S", Vec3::new(0.0, 0.0, 3.0)),
            at("ALA", 5, "CA", Vec3::new(0.0, 0.0, -3.0)),
        ]);
        assert!(hydrophobic_contacts(&structure).is_empty());
    }

    #[test]
    fn atoms_of_one_residue_do_not_contact_each_other() {
        // Leucine's own carbons are all within 4 A of one another.
        let structure = protein(vec![
            at("LEU", 1, "CB", Vec3::ZERO),
            at("LEU", 1, "CG", Vec3::new(1.5, 0.0, 0.0)),
            at("LEU", 1, "CD1", Vec3::new(2.5, 1.0, 0.0)),
        ]);
        assert!(hydrophobic_contacts(&structure).is_empty());
    }

    #[test]
    fn one_residue_pair_reports_one_contact_not_every_carbon_pair() {
        // A leucine packed against an alanine: four apolar carbons all inside
        // the cutoff of the alanine's CB, which is one touch, not four.
        let structure = protein(vec![
            at("ALA", 1, "CB", Vec3::ZERO),
            at("LEU", 2, "CB", Vec3::new(3.0, 0.0, 0.0)),
            at("LEU", 2, "CG", Vec3::new(3.2, 1.0, 0.0)),
            at("LEU", 2, "CD1", Vec3::new(3.4, 2.0, 0.0)),
            at("LEU", 2, "CD2", Vec3::new(3.6, 0.0, 1.0)),
        ]);

        let found = hydrophobic_contacts(&structure);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].atoms_a, vec![0]); // the alanine CB
        assert_eq!(found[0].atoms_b, vec![1]); // the closest leucine carbon
        assert!((found[0].distance - 3.0).abs() < 1e-4);
    }

    #[test]
    fn the_reduction_thins_both_sides() {
        // The mirror of the case above: two carbons of one residue reaching a
        // single carbon of another. Reducing only by the first atom of the pair
        // would leave both, so this is what the second pass is for.
        let structure = protein(vec![
            at("LEU", 1, "CB", Vec3::ZERO),
            at("LEU", 1, "CG", Vec3::new(0.0, 1.5, 0.0)),
            at("ALA", 2, "CB", Vec3::new(3.0, 0.0, 0.0)),
        ]);

        let found = hydrophobic_contacts(&structure);
        assert_eq!(found.len(), 1);
        assert!((found[0].distance - 3.0).abs() < 1e-4); // the shorter leg
    }

    #[test]
    fn distinct_residue_pairs_each_keep_a_contact() {
        // Reduction must not collapse genuinely different partners: one leucine
        // carbon touching three separate alanines is three contacts.
        let structure = protein(vec![
            at("LEU", 1, "CB", Vec3::ZERO),
            at("ALA", 2, "CB", Vec3::new(3.0, 0.0, 0.0)),
            at("ALA", 3, "CB", Vec3::new(0.0, 3.1, 0.0)),
            at("ALA", 4, "CB", Vec3::new(0.0, 0.0, 3.2)),
        ]);
        assert_eq!(hydrophobic_contacts(&structure).len(), 3);
    }

    #[test]
    fn salt_bridge_rejected_beyond_cutoff_and_between_like_charges() {
        let far = protein(vec![
            atom("ASP", 1, "OD1", "O", Vec3::ZERO),
            atom("ASP", 1, "OD2", "O", Vec3::ZERO),
            atom("LYS", 2, "NZ", "N", Vec3::new(7.0, 0.0, 0.0)),
        ]);
        assert!(salt_bridges(&far).is_empty());

        // Two anions close together must not pair.
        let like = protein(vec![
            atom("ASP", 1, "OD1", "O", Vec3::ZERO),
            atom("GLU", 2, "OE1", "O", Vec3::new(3.0, 0.0, 0.0)),
        ]);
        assert!(salt_bridges(&like).is_empty());
    }

    /// A flat six-membered ring, named as a phenylalanine, centred at `center`
    /// and spanned by two in-plane axes so its plane can be aimed.
    fn ring_at(seq: u32, center: Vec3, u: Vec3, v: Vec3) -> Vec<Atom> {
        let names = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"];
        (0..6)
            .map(|k| {
                let t = std::f32::consts::PI / 3.0 * k as f32;
                atom(
                    "PHE",
                    seq,
                    names[k],
                    "C",
                    center + (u * t.cos() + v * t.sin()) * 1.4,
                )
            })
            .collect()
    }

    fn stack(a: Vec<Atom>, b: Vec<Atom>) -> Vec<Interaction> {
        let mut atoms = a;
        atoms.extend(b);
        pi_stacking(&protein(atoms))
    }

    #[test]
    fn face_to_face_rings_stack_and_report_a_parallel_angle() {
        // Two rings in parallel planes, one above the other.
        let found = stack(
            ring_at(1, Vec3::ZERO, Vec3::X, Vec3::Y),
            ring_at(2, Vec3::new(0.0, 0.0, 3.8), Vec3::X, Vec3::Y),
        );
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::PiStacking);
        assert!((found[0].distance - 3.8).abs() < 1e-3);
        assert!(
            found[0].angle.unwrap() < PISTACK_ANGLE_DEV,
            "a parallel stack"
        );
    }

    #[test]
    fn edge_to_face_rings_stack_and_report_a_perpendicular_angle() {
        // The second ring stands on edge, its plane perpendicular to the first.
        let found = stack(
            ring_at(1, Vec3::ZERO, Vec3::X, Vec3::Y),
            ring_at(2, Vec3::new(0.0, 0.0, 4.0), Vec3::X, Vec3::Z),
        );
        assert_eq!(found.len(), 1);
        assert!(
            found[0].angle.unwrap() > 90.0 - PISTACK_ANGLE_DEV,
            "a T-shape"
        );
    }

    #[test]
    fn rings_do_not_stack_beyond_the_distance_or_the_offset() {
        // Parallel but too far apart.
        assert!(
            stack(
                ring_at(1, Vec3::ZERO, Vec3::X, Vec3::Y),
                ring_at(2, Vec3::new(0.0, 0.0, 6.0), Vec3::X, Vec3::Y),
            )
            .is_empty()
        );

        // Close enough, parallel, but slipped sideways past the offset cutoff.
        assert!(
            stack(
                ring_at(1, Vec3::ZERO, Vec3::X, Vec3::Y),
                ring_at(2, Vec3::new(3.0, 0.0, 3.0), Vec3::X, Vec3::Y),
            )
            .is_empty()
        );
    }

    #[test]
    fn rings_tilted_between_parallel_and_perpendicular_do_not_stack() {
        // A 45 degree tilt is neither face-to-face nor edge-to-face.
        let tilt = (Vec3::Y + Vec3::Z).normalize();
        assert!(
            stack(
                ring_at(1, Vec3::ZERO, Vec3::X, Vec3::Y),
                ring_at(2, Vec3::new(0.0, 0.0, 3.8), Vec3::X, tilt),
            )
            .is_empty()
        );
    }

    #[test]
    fn a_tryptophans_two_rings_do_not_stack_with_each_other() {
        let trp = protein(vec![
            atom("TRP", 1, "CG", "C", Vec3::new(0.0, 0.0, 0.0)),
            atom("TRP", 1, "CD1", "C", Vec3::new(1.2, 0.5, 0.0)),
            atom("TRP", 1, "NE1", "N", Vec3::new(2.0, -0.5, 0.0)),
            atom("TRP", 1, "CE2", "C", Vec3::new(1.3, -1.6, 0.0)),
            atom("TRP", 1, "CD2", "C", Vec3::new(0.0, -1.4, 0.0)),
            atom("TRP", 1, "CZ2", "C", Vec3::new(1.6, -2.9, 0.0)),
            atom("TRP", 1, "CH2", "C", Vec3::new(0.6, -3.9, 0.0)),
            atom("TRP", 1, "CZ3", "C", Vec3::new(-0.7, -3.7, 0.0)),
            atom("TRP", 1, "CE3", "C", Vec3::new(-1.0, -2.4, 0.0)),
        ]);
        assert!(pi_stacking(&trp).is_empty());
    }

    fn cation_over(ring: Vec<Atom>, nz: Vec3) -> Vec<Interaction> {
        let mut atoms = ring;
        atoms.push(atom("LYS", 9, "NZ", "N", nz));
        pi_cation(&protein(atoms))
    }

    #[test]
    fn a_cation_over_a_ring_face_is_a_pi_cation() {
        let found = cation_over(
            ring_at(1, Vec3::ZERO, Vec3::X, Vec3::Y),
            Vec3::new(0.0, 0.0, 4.5),
        );
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::PiCation);
        assert!((found[0].distance - 4.5).abs() < 1e-3);
    }

    #[test]
    fn a_cation_is_rejected_beyond_the_distance_or_off_the_face() {
        // Straight over the face but past the six angstrom cutoff.
        assert!(
            cation_over(
                ring_at(1, Vec3::ZERO, Vec3::X, Vec3::Y),
                Vec3::new(0.0, 0.0, 6.5)
            )
            .is_empty()
        );

        // Close, but beside the ring rather than over it.
        assert!(
            cation_over(
                ring_at(1, Vec3::ZERO, Vec3::X, Vec3::Y),
                Vec3::new(3.0, 0.0, 2.0)
            )
            .is_empty()
        );
    }

    #[test]
    fn a_histidine_ring_does_not_count_its_own_charge() {
        // A lone histidine is both a ring and a cation; it must not pi-cation
        // with itself.
        let his = protein(vec![
            atom("HIS", 1, "CG", "C", Vec3::new(0.0, 0.0, 0.0)),
            atom("HIS", 1, "ND1", "N", Vec3::new(1.1, 0.7, 0.0)),
            atom("HIS", 1, "CE1", "C", Vec3::new(2.2, 0.0, 0.0)),
            atom("HIS", 1, "NE2", "N", Vec3::new(1.9, -1.2, 0.0)),
            atom("HIS", 1, "CD2", "C", Vec3::new(0.5, -1.2, 0.0)),
        ]);
        assert!(pi_cation(&his).is_empty());
    }

    /// A donor residue, an acceptor residue, and a water between them. Residue 1
    /// is a proline so it contributes the carbonyl that fixes residue 2's amide
    /// hydrogen without donating a terminal ammonium of its own.
    fn bridged(water: Vec3, acceptor: Vec3) -> Protein {
        protein(vec![
            atom("PRO", 1, "N", "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", "C", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("ALA", 2, "N", "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("ALA", 2, "CA", "C", Vec3::new(4.0, 2.8, 0.0)),
            atom("HOH", 9, "O", "O", water),
            atom("GLY", 8, "O", "O", acceptor),
        ])
    }

    /// Residue 2's amide hydrogen points along this direction from its nitrogen.
    fn amide_direction() -> Vec3 {
        (Vec3::new(2.0, 1.4, 0.0) - Vec3::new(1.2, 2.3, 0.0)).normalize()
    }

    #[test]
    fn a_water_bridges_a_donor_to_an_acceptor() {
        // Water 3 A along the amide hydrogen's direction, with an acceptor off
        // to the side so the angle at the water lands inside the window.
        let nitrogen = Vec3::new(3.3, 1.5, 0.0);
        let water = nitrogen + amide_direction() * 3.0;
        let found = water_bridges(&bridged(water, water + Vec3::new(0.0, -2.9, 0.6)));

        assert_eq!(found.len(), 1);
        assert_eq!(found[0].kind, InteractionKind::WaterBridge);
        assert_eq!(found[0].atoms_a, vec![4]); // the amide nitrogen
        assert_eq!(found[0].atoms_b, vec![7]); // the carbonyl oxygen
        assert_eq!(found[0].bridge, Some(6)); // the water oxygen
        let angle = found[0].angle.unwrap();
        assert!(
            WATER_BRIDGE_OMEGA_MIN < angle && angle < WATER_BRIDGE_OMEGA_MAX,
            "the angle at the water is inside the window: {angle}"
        );
    }

    #[test]
    fn a_leg_outside_the_distance_band_is_not_a_bridge() {
        let nitrogen = Vec3::new(3.3, 1.5, 0.0);

        // The water sits too far from the donor.
        let far = nitrogen + amide_direction() * 5.0;
        assert!(water_bridges(&bridged(far, far + Vec3::new(0.0, -2.9, 0.6))).is_empty());

        // The water is in range of the donor but the acceptor is too far from it.
        let water = nitrogen + amide_direction() * 3.0;
        assert!(water_bridges(&bridged(water, water + Vec3::new(0.0, -6.0, 0.0))).is_empty());
    }

    #[test]
    fn an_acceptor_straight_through_the_water_is_not_a_bridge() {
        // Two hydrogen bonds cannot leave one oxygen in opposite directions.
        let nitrogen = Vec3::new(3.3, 1.5, 0.0);
        let water = nitrogen + amide_direction() * 3.0;
        let beyond = water + amide_direction() * 2.9;

        assert!(water_bridges(&bridged(water, beyond)).is_empty());
    }

    #[test]
    fn a_structure_without_water_has_no_bridges() {
        let dry = protein(vec![
            atom("PRO", 1, "N", "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", "C", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("ALA", 2, "N", "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("GLY", 8, "O", "O", Vec3::new(5.0, 2.0, 0.0)),
        ]);
        assert!(water_bridges(&dry).is_empty());
    }
}
