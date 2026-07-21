//! Polar hydrogen inference.
//!
//! A crystal structure has no hydrogens, and every hydrogen bond criterion worth
//! using needs one. Donors fall into two kinds:
//!
//! - The heavy atoms **fix** the hydrogen (a backbone amide, a guanidinium, an
//!   indole NH). It is placed once, and the resulting angle is real evidence.
//! - The heavy atoms leave a **rotation** free (a hydroxyl, an ammonium). No
//!   single position is right, so the hydrogen is placed against the acceptor
//!   being tested: the point on its circle that comes closest to reaching it.
//!
//! Placing a rotatable hydrogen per acceptor is HBPLUS's rule ("if the hydrogen
//! on the donor is not fixed then the point on the hydrogen's locus that is
//! closest to the acceptor being investigated is selected"), and Reduce does the
//! same thing by scanning torsions for the one that points a hydroxyl hydrogen at
//! a neighbouring acceptor.
//!
//! The rotatable case still gates: the hydrogen is confined to a circle, so an
//! acceptor lying at a bad angle to the donor's bond cannot be reached at any
//! torsion. It is a weaker constraint than the fixed case, not a vacuous one.
//! Both HBPLUS and Chimera call the results "potential" hydrogen bonds for this
//! reason, and so does this.
//!
//! Geometry is HBPLUS's Table III. Hydrogen bond lengths there are 1.00 A for
//! every donor here except the lysine ammonium at 1.01 A; that 0.01 A is far
//! under the coordinate error of the structures this runs on, so one length is
//! used.
//!
//! A standard residue's donors are tabulated. Everything else is read off the
//! component definition the file carries. The heavy atoms come from the bond graph
//! rather than the definition.
//!
//! Sources:
//! - HBPLUS: McDonald and Thornton, J Mol Biol 1994; manual Tables II and III at
//!   `https://www.uoxray.uoregon.edu/local/manuals/hbplus/hbplus.txt`
//! - Reduce: Word et al., J Mol Biol 285:1735 (1999); `RotDonor::finalize` in
//!   `https://github.com/rlabduke/reduce`

use glam::Vec3;

use crate::structure::{Atom, Protein};

use super::classify;
use super::connectivity::Bonds;
use super::topology;

/// Length of every inferred polar hydrogen bond (HBPLUS Table III).
const H_BOND_LENGTH: f32 = 1.0;
/// The DD-D-H angle at an sp2 donor carrying two hydrogens (HBPLUS Table III).
const SP2_H_ANGLE: f32 = 120.0;
/// The DD-D-H angle at an sp3 hydroxyl or ammonium (HBPLUS Table III).
const SP3_H_ANGLE: f32 = 110.0;
/// The H-D-H angle between two hydrogens on one sp3 centre.
const SP3_H_H_ANGLE: f32 = 109.5;
/// A C(prev)-N separation beyond this (A) is a chain break, not a peptide bond.
pub(super) const MAX_PEPTIDE_BOND: f32 = 2.5;

/// A hydrogen bond donor and the polar hydrogens it carries.
pub struct Donor {
    pub atom: usize,
    pub hydrogens: Hydrogens,
}

/// Where a donor's polar hydrogens sit.
pub enum Hydrogens {
    /// Positions the heavy atoms determine.
    Fixed(Vec<Vec3>),
    /// Hydrogens free to rotate about the donor's single bond. They share one
    /// circle, so `count` of them can be satisfied at once but each is placed
    /// against the acceptor under test.
    Rotor { circle: Circle, count: usize },
    /// One hydrogen, and more than one site the heavy atoms allow it: a planar
    /// centre carrying a single hydrogen can take either of the two positions in
    /// its plane, and the coordinates do not say which.
    Undetermined(Vec<Vec3>),
}

impl Hydrogens {
    /// How many hydrogen bonds this donor can make at once: one per hydrogen.
    pub fn capacity(&self) -> usize {
        match self {
            Hydrogens::Fixed(positions) => positions.len(),
            Hydrogens::Rotor { count, .. } => *count,
            Hydrogens::Undetermined(_) => 1,
        }
    }

    /// The D-H...A angle each of this donor's hydrogens makes with `acceptor`,
    /// in degrees, one per hydrogen.
    ///
    /// A caller choosing between acceptors needs the angle per hydrogen rather
    /// than the best of them: a donor carrying two hydrogens can reach two
    /// acceptors only if they are reached by *different* ones, and an amide
    /// pointing both hydrogens the same way often is not.
    ///
    /// A rotor's hydrogens are interchangeable, so each is reported at the best
    /// angle the circle can reach. That understates the coupling between them,
    /// since turning one hydrogen turns them all, and it is the approximation
    /// this makes of a rotor that Reduce models by scanning torsions.
    pub fn angles_toward(&self, donor: Vec3, acceptor: Vec3) -> Vec<f32> {
        self.positions_toward(acceptor)
            .into_iter()
            .map(|hydrogen| angle_at(hydrogen, donor, acceptor))
            .collect()
    }

    /// Where each of this donor's hydrogens sits when reaching toward `target`.
    ///
    /// A rotor's hydrogens share one circle, so each is reported at the point
    /// that best reaches the target and they coincide there. A caller that needs
    /// one hydrogen rather than all of them takes the one making the best angle.
    pub fn positions_toward(&self, target: Vec3) -> Vec<Vec3> {
        match self {
            Hydrogens::Fixed(positions) => positions.clone(),
            Hydrogens::Rotor { circle, count } => vec![circle.toward(target); *count],
            Hydrogens::Undetermined(sites) => sites
                .iter()
                .copied()
                .min_by(|a, b| {
                    (*a - target)
                        .length_squared()
                        .total_cmp(&(*b - target).length_squared())
                })
                .into_iter()
                .collect(),
        }
    }
}

/// The circle a rotatable hydrogen traces about its donor's one bond.
pub struct Circle {
    donor: Vec3,
    /// Unit vector from the antecedent towards the donor.
    axis: Vec3,
    /// The hydrogen's offset from the donor along the axis.
    along: f32,
    /// The circle's radius about the axis.
    radius: f32,
}

impl Circle {
    /// A circle at `angle` (the antecedent-donor-hydrogen angle) about the bond
    /// from `antecedent` to `donor`.
    fn new(donor: Vec3, antecedent: Vec3, angle: f32) -> Option<Self> {
        let axis = (donor - antecedent).try_normalize()?;
        let (sin, cos) = angle.to_radians().sin_cos();

        Some(Circle {
            donor,
            axis,
            along: -cos * H_BOND_LENGTH,
            radius: sin * H_BOND_LENGTH,
        })
    }

    /// The point on the circle closest to the donor-acceptor line, which is the
    /// hydrogen position that best reaches `acceptor`.
    ///
    /// The angle D-H...A shrinks as the hydrogen swings off that line, so the
    /// closest point is also the one giving the largest angle. When the acceptor
    /// lies on the axis every torsion is equivalent and any spoke will do.
    fn toward(&self, acceptor: Vec3) -> Vec3 {
        let to_acceptor = acceptor - self.donor;
        let off_axis = to_acceptor - self.axis * to_acceptor.dot(self.axis);
        let spoke = off_axis
            .try_normalize()
            .unwrap_or_else(|| self.axis.any_orthonormal_vector());

        self.donor + self.axis * self.along + spoke * self.radius
    }
}

/// The angle at `vertex` between the directions to `a` and to `b`, in degrees.
pub fn angle_at(vertex: Vec3, a: Vec3, b: Vec3) -> f32 {
    let (Some(u), Some(v)) = ((a - vertex).try_normalize(), (b - vertex).try_normalize()) else {
        return 0.0;
    };

    u.dot(v).clamp(-1.0, 1.0).acos().to_degrees()
}

/// How a donor's hydrogens are constrained by the heavy atoms around it.
enum Geometry {
    /// The backbone amide: one hydrogen opposite the preceding carbonyl.
    Amide,
    /// One hydrogen on the outward bisector of two heavy neighbours.
    Bisector,
    /// Two hydrogens in the sp2 plane, at `SP2_H_ANGLE` to the single neighbour.
    PlanarPair,
    /// Hydrogens free to rotate about the bond to the single neighbour.
    Rotor { count: usize },
}

/// The donor table.
///
/// Sulfur is absent: PLIP inherits OpenBabel's typing, which admits only N, O
/// and F as donors, so a cysteine thiol never donates. HBPLUS, CHARMM and
/// ChimeraX all disagree and let free cysteine donate. Water is absent too; its
/// hydrogens have no antecedent to place them against, and PLIP routes waters
/// through its water-bridge detector rather than this one.
///
/// Histidine's ND1 and NE2 both appear here, and both also accept. Only one of
/// them carries a hydrogen in a given tautomer, and the coordinates do not say
/// which, so both possibilities are reported.
fn geometry_of(atom: &Atom) -> Option<Geometry> {
    let residue = atom.residue_name.as_str();
    let name = atom.name.as_str();

    if !topology::is_standard_residue(residue) {
        return None;
    }
    if name == "N" {
        // Proline's nitrogen is tertiary and carries no hydrogen.
        return (residue != "PRO").then_some(Geometry::Amide);
    }

    match (residue, name) {
        ("ARG", "NE") | ("TRP", "NE1") | ("HIS", "ND1" | "NE2") => Some(Geometry::Bisector),
        ("ASN", "ND2") | ("GLN", "NE2") | ("ARG", "NH1" | "NH2") => Some(Geometry::PlanarPair),
        ("SER", "OG") | ("THR", "OG1") | ("TYR", "OH") => Some(Geometry::Rotor { count: 1 }),
        ("LYS", "NZ") => Some(Geometry::Rotor { count: 3 }),
        _ => None,
    }
}

/// One hydrogen pointing away from every heavy neighbour at once: the bisector
/// of two of them, or the fourth vertex left by three.
fn outward(donor: Vec3, neighbours: &[Vec3]) -> Option<Vec3> {
    let inward = neighbours.iter().try_fold(Vec3::ZERO, |sum, &neighbour| {
        Some(sum + (neighbour - donor).try_normalize()?)
    })?;

    Some(donor - inward.try_normalize()? * H_BOND_LENGTH)
}

/// Two hydrogens either side of the plane through `donor` and its two
/// neighbours, which is where a protonated secondary amine carries them.
fn straddling(donor: Vec3, first: Vec3, second: Vec3) -> Option<[Vec3; 2]> {
    let away = (outward(donor, &[first, second])? - donor).try_normalize()?;
    let normal = (first - donor).cross(second - donor).try_normalize()?;
    let (sin, cos) = (SP3_H_H_ANGLE / 2.0).to_radians().sin_cos();

    Some([
        donor + (away * cos + normal * sin) * H_BOND_LENGTH,
        donor + (away * cos - normal * sin) * H_BOND_LENGTH,
    ])
}

/// Two hydrogens in the plane of `beyond`-`neighbour`-`donor`, each at
/// `SP2_H_ANGLE` to the bond, which is where an sp2 amine's pair sits.
fn planar_pair(donor: Vec3, neighbour: Vec3, beyond: Vec3) -> Option<[Vec3; 2]> {
    let bond = (donor - neighbour).try_normalize()?;
    let normal = (donor - neighbour)
        .cross(beyond - neighbour)
        .try_normalize()?;
    let in_plane = normal.cross(bond);
    let (sin, cos) = (180.0 - SP2_H_ANGLE).to_radians().sin_cos();

    Some([
        donor + (bond * cos + in_plane * sin) * H_BOND_LENGTH,
        donor + (bond * cos - in_plane * sin) * H_BOND_LENGTH,
    ])
}

/// The preceding residue's carbonyl, which fixes the next amide hydrogen.
struct Carbonyl {
    chain: crate::fixed_str::FixedStr<4>,
    c: Vec3,
    /// The carbonyl oxygen, absent when the model left it out.
    o: Option<Vec3>,
}

/// Every donor in the structure, with its hydrogens placed.
///
/// A donor whose hydrogen cannot be placed is not a donor. That is the residue
/// after a chain break, and the residue whose predecessor was modelled without
/// its carbonyl oxygen: both have an amide hydrogen with no direction to take.
/// It is not the first residue of a chain, which carries an ammonium rather than
/// an amide and needs nothing outside itself.
pub fn donors(protein: &Protein, bonds: &Bonds) -> Vec<Donor> {
    let mut out = Vec::new();
    let mut preceding: Option<Carbonyl> = None;

    let mut index = 0;
    while index < protein.atoms.len() {
        let start = index;
        let key = residue_key(&protein.atoms[start]);
        while index < protein.atoms.len() && residue_key(&protein.atoms[index]) == key {
            index += 1;
        }
        let residue = start..index;

        for i in residue.clone() {
            let placed = match geometry_of(&protein.atoms[i]) {
                Some(geometry) => place(protein, residue.clone(), i, geometry, preceding.as_ref()),
                None => place_from_component(protein, bonds, i),
            };

            if let Some(hydrogens) = placed {
                out.push(Donor { atom: i, hydrogens });
            }
        }

        // Remember the last backbone carbonyl carbon seen, with its oxygen if
        // the model has one.
        if let Some(c) = position_of(protein, residue.clone(), "C", None) {
            preceding = Some(Carbonyl {
                chain: protein.atoms[start].chain_id,
                c,
                o: position_of(protein, residue.clone(), "O", None),
            });
        }
    }

    out
}

/// Place a donor's hydrogens, reading the antecedents it needs out of the
/// residue's own atoms.
fn place(
    protein: &Protein,
    residue: std::ops::Range<usize>,
    donor_index: usize,
    geometry: Geometry,
    preceding: Option<&Carbonyl>,
) -> Option<Hydrogens> {
    let atom = &protein.atoms[donor_index];
    let donor = atom.position;
    let residue_name = atom.residue_name.as_str();
    let name = atom.name.as_str();

    let altloc = atom.altloc;
    let antecedent = |of: &str| {
        let of = topology::heavy_neighbours(residue_name, of).next()?;
        position_of(protein, residue.clone(), of, altloc)
    };

    match geometry {
        Geometry::Amide => {
            match preceding {
                // An interior amide: the hydrogen lies opposite the preceding
                // C=O, as in DSSP.
                Some(previous) if previous.chain == atom.chain_id => {
                    if (previous.c - donor).length() > MAX_PEPTIDE_BOND {
                        return None;
                    }
                    let along = (previous.c - previous.o?).try_normalize()?;

                    Some(Hydrogens::Fixed(vec![donor + along * H_BOND_LENGTH]))
                }
                // The first residue of a chain is an ammonium, not an amide: it
                // carries three hydrogens turning about the N-CA bond, the same
                // class HBPLUS gives a lysine.
                _ => Some(Hydrogens::Rotor {
                    circle: Circle::new(donor, antecedent("N")?, SP3_H_ANGLE)?,
                    count: 3,
                }),
            }
        }
        Geometry::Bisector => {
            let mut neighbours = topology::heavy_neighbours(residue_name, name);
            let first = position_of(protein, residue.clone(), neighbours.next()?, altloc)?;
            let second = position_of(protein, residue.clone(), neighbours.next()?, altloc)?;

            Some(Hydrogens::Fixed(vec![outward(donor, &[first, second])?]))
        }
        Geometry::PlanarPair => {
            let neighbour_name = topology::heavy_neighbours(residue_name, name).next()?;
            let neighbour = position_of(protein, residue.clone(), neighbour_name, altloc)?;

            // Any further neighbour of the sp2 centre defines the same plane.
            let beyond_name =
                topology::heavy_neighbours(residue_name, neighbour_name).find(|&n| n != name)?;
            let beyond = position_of(protein, residue.clone(), beyond_name, altloc)?;

            Some(Hydrogens::Fixed(
                planar_pair(donor, neighbour, beyond)?.to_vec(),
            ))
        }
        Geometry::Rotor { count } => Some(Hydrogens::Rotor {
            circle: Circle::new(donor, antecedent(name)?, SP3_H_ANGLE)?,
            count,
        }),
    }
}

/// Place the hydrogens of a donor the residue tables do not describe, from what
/// the file says its component is.
///
/// The donor set is OpenBabel's, which PLIP inherits (`OBAtom::IsHbondDonor`):
/// nitrogen, oxygen or fluorine bearing at least one hydrogen. Sulfur is absent
/// there and so absent here, as it is for the standard residues.
///
/// The hydrogen count comes from the definition and the heavy neighbours from
/// the bond graph, which is what lets a polymerised component be handled.
fn place_from_component(protein: &Protein, bonds: &Bonds, index: usize) -> Option<Hydrogens> {
    let atom = &protein.atoms[index];
    let residue_name = atom.residue_name.as_str();
    if topology::is_standard_residue(residue_name) || classify::is_water(residue_name) {
        return None;
    }

    let element = atom.element.as_str();
    if !["N", "O", "F"]
        .iter()
        .any(|donor| element.eq_ignore_ascii_case(donor))
    {
        return None;
    }

    let component = protein.components.get(residue_name)?;
    let name = atom.name.as_str();
    if classify::is_acid_oxygen(component, name) {
        return None;
    }

    let neighbours = bonds.of(index);
    let count = hydrogen_count(protein, bonds, component, index)?;

    // A hydroxyl hydrogen turns about its bond whatever the centre is conjugated to.
    let planar = element.eq_ignore_ascii_case("N")
        && classify::hybridisation(component, name) == 2
        && neighbours
            .iter()
            .all(|&n| classify::hybridisation(component, protein.atoms[n].name.as_str()) == 2);

    let positions: Vec<Vec3> = neighbours
        .iter()
        .map(|&n| protein.atoms[n].position)
        .collect();
    // Any second-shell atom fixes the plane of a planar centre.
    let beyond = neighbours
        .first()
        .and_then(|&first| bonds.of(first).iter().find(|&&n| n != index))
        .map(|&n| protein.atoms[n].position);

    place_around(atom.position, &positions, beyond, count, planar)
}

/// How many hydrogens an atom carries in the structure, from the count the
/// definition gives the free component and the bonds the coordinates show.
///
/// The two disagree whenever the component is not free.
///
/// The converse is deliberately not applied. A bond the definition lists and the
/// coordinates lack would, by the same argument, leave a hydrogen behind, and at
/// a chain's first residue it genuinely does.
fn hydrogen_count(
    protein: &Protein,
    bonds: &Bonds,
    component: &crate::components::Component,
    index: usize,
) -> Option<usize> {
    let name = protein.atoms[index].name.as_str();

    let gained = bonds
        .of(index)
        .iter()
        .filter(|&&n| {
            let neighbour = protein.atoms[n].name.as_str();
            !component.neighbours(name).any(|other| other == neighbour)
        })
        .count();

    component
        .hydrogen_count(name)
        .checked_sub(gained)
        .filter(|&count| count > 0)
}

/// Place `count` hydrogens on a donor surrounded by the heavy atoms at
/// `neighbours`.
fn place_around(
    donor: Vec3,
    neighbours: &[Vec3],
    beyond: Option<Vec3>,
    count: usize,
    planar: bool,
) -> Option<Hydrogens> {
    match neighbours {
        [] => None,
        [single] if !planar => Some(Hydrogens::Rotor {
            circle: Circle::new(donor, *single, SP3_H_ANGLE)?,
            count,
        }),
        [single] => {
            let sites = planar_pair(donor, *single, beyond?)?.to_vec();
            match count {
                1 => Some(Hydrogens::Undetermined(sites)),
                2 => Some(Hydrogens::Fixed(sites)),
                _ => None,
            }
        }
        [first, second] if count == 2 => Some(Hydrogens::Fixed(
            straddling(donor, *first, *second)?.to_vec(),
        )),
        _ if count == 1 => Some(Hydrogens::Fixed(vec![outward(donor, neighbours)?])),
        _ => None,
    }
}

/// The position of a named atom in a residue, taken from the same alternate
/// conformer as the atom asking.
///
/// A residue's range spans every conformer, so matching on name alone would read
/// conformer A's carbon while placing conformer B's hydroxyl hydrogen, building
/// the rotor about an axis that is not a bond in either one. Atoms outside any
/// conformer are shared by all of them, so they are the fallback.
fn position_of(
    protein: &Protein,
    residue: std::ops::Range<usize>,
    name: &str,
    altloc: Option<char>,
) -> Option<Vec3> {
    let atoms = &protein.atoms[residue];
    let named = |want: Option<char>| atoms.iter().find(|a| a.name == name && a.altloc == want);

    named(altloc)
        .or_else(|| named(None))
        .or_else(|| atoms.iter().find(|a| a.name == name))
        .map(|atom| atom.position)
}

type ResidueKey = (crate::fixed_str::FixedStr<4>, u32, Option<char>);

fn residue_key(atom: &Atom) -> ResidueKey {
    (atom.chain_id, atom.residue_seq, atom.insertion_code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixed_str::FixedStr;
    use crate::structure::SecondaryStructure;

    fn atom(residue: &str, seq: u32, name: &str, position: Vec3) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new(name),
            element: FixedStr::new(&name[..1]),
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

    /// Every donor of a structure, over its own bond graph.
    fn placed(structure: &Protein) -> Vec<Donor> {
        donors(structure, &super::super::connectivity::bonds(structure))
    }

    fn protein(atoms: Vec<Atom>) -> Protein {
        Protein {
            atoms,
            title: String::new(),
            surface_points: Vec::new(),
            hbonds: Vec::new(),
            assembly: None,
            components: Default::default(),
        }
    }

    /// One component carrying every donor shape the ligand path decides between,
    /// written the way a structure file states it.
    ///
    /// O2 is a hydroxyl, N4 a carboxamide, N6 a secondary amine, N11 an amidine
    /// nitrogen, N20 a protonated secondary amine, and S8 a thiol.
    const DEFINITION: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
LIG C1  C
LIG O2  O
LIG H2  H
LIG C3  C
LIG O9  O
LIG N4  N
LIG H4A H
LIG H4B H
LIG C5  C
LIG N6  N
LIG H6  H
LIG C7  C
LIG S8  S
LIG H8  H
LIG C10 C
LIG C12 C
LIG N11 N
LIG H11 H
LIG C19 C
LIG C21 C
LIG N20 N
LIG H20A H
LIG H20B H
LIG C31 C
LIG O32 O
LIG C33 C
LIG C34 C
LIG N30 N
LIG H30 H
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG C1  O2   sing
LIG O2  H2   sing
LIG C3  O9   doub
LIG C3  N4   sing
LIG N4  H4A  sing
LIG N4  H4B  sing
LIG C5  N6   sing
LIG N6  C7   sing
LIG N6  H6   sing
LIG C7  S8   sing
LIG S8  H8   sing
LIG C10 N11  doub
LIG C10 C12  sing
LIG N11 H11  sing
LIG C19 N20  sing
LIG N20 C21  sing
LIG N20 H20A sing
LIG N20 H20B sing
LIG C31 O32  doub
LIG C31 N30  sing
LIG N30 C33  sing
LIG N30 H30  sing
LIG C33 C34  sing
#
";

    /// An atom of the ligand above, with its element stated rather than guessed.
    fn ligand_atom(name: &str, element: &str, position: Vec3) -> Atom {
        let mut atom = atom("LIG", 1, name, position);
        atom.element = FixedStr::new(element);
        atom.is_hetatm = true;

        atom
    }

    /// The ligand above, with only the atoms listed modelled.
    fn ligand(atoms: Vec<Atom>) -> Protein {
        let mut structure = protein(atoms);
        structure.components = crate::components::Dictionary::parse(DEFINITION);

        structure
    }

    /// The hydrogens placed on the named atom of a structure.
    fn hydrogens_on(structure: &Protein, name: &str) -> Option<Hydrogens> {
        placed(structure)
            .into_iter()
            .find(|donor| structure.atoms[donor.atom].name == name)
            .map(|donor| donor.hydrogens)
    }

    #[test]
    fn a_ligand_hydroxyl_turns_about_its_bond() {
        let structure = ligand(vec![
            ligand_atom("C1", "C", Vec3::new(-1.4, 0.0, 0.0)),
            ligand_atom("O2", "O", Vec3::ZERO),
        ]);

        let Some(Hydrogens::Rotor { count, .. }) = hydrogens_on(&structure, "O2") else {
            panic!("a hydroxyl hydrogen turns about the C-O bond");
        };
        assert_eq!(count, 1);
    }

    #[test]
    fn a_ligand_carboxamide_carries_a_fixed_pair_in_its_plane() {
        let structure = ligand(vec![
            ligand_atom("O9", "O", Vec3::new(-2.5, 1.0, 0.0)),
            ligand_atom("C3", "C", Vec3::new(-1.5, 0.0, 0.0)),
            ligand_atom("N4", "N", Vec3::ZERO),
        ]);

        let Some(Hydrogens::Fixed(positions)) = hydrogens_on(&structure, "N4") else {
            panic!("an sp2 amide nitrogen has both hydrogens fixed");
        };
        assert_eq!(positions.len(), 2);
        for position in positions {
            assert!(position.z.abs() < 1e-5, "the hydrogen left the amide plane");
        }
    }

    #[test]
    fn a_ligand_nitrogen_between_two_carbons_points_its_hydrogen_away_from_both() {
        let structure = ligand(vec![
            ligand_atom("C5", "C", Vec3::new(1.0, 1.0, 0.0)),
            ligand_atom("N6", "N", Vec3::ZERO),
            ligand_atom("C7", "C", Vec3::new(1.0, -1.0, 0.0)),
        ]);

        let Some(Hydrogens::Fixed(positions)) = hydrogens_on(&structure, "N6") else {
            panic!("one hydrogen between two neighbours is fixed");
        };
        assert_eq!(positions.len(), 1);
        assert!((positions[0] - Vec3::new(-1.0, 0.0, 0.0)).length() < 1e-5);
    }

    /// A planar centre carrying one hydrogen has two sites in its plane and the
    /// coordinates do not say which is taken, so the site is chosen against the
    /// acceptor as a rotor's is.
    #[test]
    fn an_amidine_nitrogen_takes_whichever_of_its_two_sites_faces_the_acceptor() {
        let nitrogen = Vec3::ZERO;
        let structure = ligand(vec![
            ligand_atom("C12", "C", Vec3::new(-2.5, 1.0, 0.0)),
            ligand_atom("C10", "C", Vec3::new(-1.5, 0.0, 0.0)),
            ligand_atom("N11", "N", nitrogen),
        ]);

        let Some(hydrogens) = hydrogens_on(&structure, "N11") else {
            panic!("an amidine nitrogen donates");
        };
        assert_eq!(hydrogens.capacity(), 1, "it carries one hydrogen, not two");

        // An acceptor either side of the plane's axis draws the hydrogen to the
        // site on its own side.
        let above = hydrogens.positions_toward(Vec3::new(1.0, 3.0, 0.0))[0];
        let below = hydrogens.positions_toward(Vec3::new(1.0, -3.0, 0.0))[0];
        assert!(above.y > 0.0 && below.y < 0.0, "the site did not follow");
        assert!((above - below).length() > 1.0, "both sites coincide");
    }

    #[test]
    fn a_protonated_secondary_amine_carries_its_pair_either_side_of_the_plane() {
        let structure = ligand(vec![
            ligand_atom("C19", "C", Vec3::new(1.0, 1.0, 0.0)),
            ligand_atom("N20", "N", Vec3::ZERO),
            ligand_atom("C21", "C", Vec3::new(1.0, -1.0, 0.0)),
        ]);

        let Some(Hydrogens::Fixed(positions)) = hydrogens_on(&structure, "N20") else {
            panic!("both hydrogens of an ammonium are fixed");
        };
        assert_eq!(positions.len(), 2);
        // The neighbours lie in z = 0, so the hydrogens straddle it.
        assert!(positions[0].z * positions[1].z < 0.0, "not straddling");
        assert!((positions[0].x + positions[1].x) < 0.0, "not pointing away");
    }

    #[test]
    fn a_ligand_thiol_and_carbon_never_donate() {
        let structure = ligand(vec![
            ligand_atom("C1", "C", Vec3::new(-1.4, 0.0, 0.0)),
            ligand_atom("C7", "C", Vec3::ZERO),
            ligand_atom("S8", "S", Vec3::new(1.8, 0.0, 0.0)),
        ]);

        assert!(hydrogens_on(&structure, "S8").is_none());
        assert!(hydrogens_on(&structure, "C1").is_none());
    }

    #[test]
    fn a_ligand_atom_with_no_neighbour_modelled_has_nothing_to_orient_against() {
        let structure = ligand(vec![ligand_atom("O2", "O", Vec3::ZERO)]);

        assert!(placed(&structure).is_empty());
    }

    /// A definition describes the free component, so a residue in a chain has
    /// one fewer hydrogen than it lists: the peptide bond stands where it was.
    #[test]
    fn a_bond_to_the_next_component_accounts_for_one_of_its_hydrogens() {
        const MODIFIED: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
MSE N   N
MSE H   H
MSE H2  H
MSE CA  C
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
MSE N  H
MSE N  H2
MSE N  CA
#
";
        let mut structure = protein(vec![
            atom("PRO", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("MSE", 2, "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("MSE", 2, "CA", Vec3::new(4.0, 2.8, 0.0)),
        ]);
        structure.components = crate::components::Dictionary::parse(MODIFIED);

        let Some(hydrogens) = hydrogens_on(&structure, "N") else {
            panic!("a residue in a chain still donates");
        };
        assert_eq!(
            hydrogens.capacity(),
            1,
            "the peptide bond replaced the second hydrogen the definition lists"
        );
    }

    /// The same definitions cover the standard residues, and there they describe
    /// the free amino acid rather than the polymer. The residue tables answer for
    /// those, so a definition cannot hand an interior amide a second hydrogen.
    #[test]
    fn a_standard_residue_ignores_the_hydrogen_count_in_the_definitions() {
        const FREE_ALANINE: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
ALA N   N
ALA H   H
ALA H2  H
ALA CA  C
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
ALA N  H
ALA N  H2
ALA N  CA
#
";
        let mut structure = protein(vec![
            atom("PRO", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("ALA", 2, "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("ALA", 2, "CA", Vec3::new(4.0, 2.8, 0.0)),
        ]);
        structure.components = crate::components::Dictionary::parse(FREE_ALANINE);

        let found = placed(&structure);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].hydrogens.capacity(), 1);
    }

    /// A water definition carries two hydrogens, and a water is still routed
    /// through the water-bridge detector rather than donating here.
    #[test]
    fn water_never_donates_however_its_definition_reads() {
        const WATER: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
HOH O   O
HOH H1  H
HOH H2  H
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
HOH O  H1
HOH O  H2
#
";
        let mut structure = protein(vec![atom("HOH", 1, "O", Vec3::ZERO)]);
        structure.components = crate::components::Dictionary::parse(WATER);

        assert!(placed(&structure).is_empty());
    }

    /// A nucleotide, whose definition is the free 5'-monophosphate: the acid
    /// hydrogen the file writes on OP2 is not there at pH 7.4, and the hydroxyl
    /// on O3' is given up to the bond to the next residue.
    const NUCLEOTIDE: &str = "\
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
DG P     P
DG OP1   O
DG OP2   O
DG HOP2  H
DG OP3   O
DG HOP3  H
DG O5'   O
DG C5'   C
DG C3'   C
DG O3'   O
DG HO3'  H
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
DG P   OP1  doub
DG P   OP2  sing
DG OP2 HOP2 sing
DG P   OP3  sing
DG OP3 HOP3 sing
DG P   O5'  sing
DG O5' C5'  sing
DG C3' O3'  sing
DG O3' HO3' sing
#
";

    fn nucleotide(atoms: Vec<Atom>) -> Protein {
        let mut structure = protein(atoms);
        structure.components = crate::components::Dictionary::parse(NUCLEOTIDE);

        structure
    }

    /// A phosphate is ionised at pH 7.4, and its two free oxygens are equivalent
    /// however arbitrarily the definition puts the acid hydrogen on one of them.
    #[test]
    fn a_phosphate_oxygen_does_not_donate() {
        let structure = nucleotide(vec![
            atom("DG", 1, "P", Vec3::ZERO),
            atom("DG", 1, "OP1", Vec3::new(1.5, 0.0, 0.0)),
            atom("DG", 1, "OP2", Vec3::new(-0.5, 1.4, 0.0)),
            atom("DG", 1, "O5'", Vec3::new(-0.5, -1.4, 0.0)),
        ]);

        assert!(hydrogens_on(&structure, "OP2").is_none());
        assert!(hydrogens_on(&structure, "OP1").is_none());
    }

    /// The bridging oxygen of a phosphodiester gave its hydroxyl hydrogen to the
    /// bond. Placing one anyway aims it at the phosphorus it is bonded through
    /// and reports a covalent bond as a hydrogen bond.
    #[test]
    fn a_linkage_oxygen_does_not_keep_the_free_monomers_hydroxyl() {
        let linked = nucleotide(vec![
            atom("DG", 1, "C3'", Vec3::new(-1.4, 0.0, 0.0)),
            atom("DG", 1, "O3'", Vec3::ZERO),
            atom("DG", 2, "P", Vec3::new(1.6, 0.0, 0.0)),
        ]);
        assert!(hydrogens_on(&linked, "O3'").is_none());

        // The 3' terminus keeps its genuine hydroxyl.
        let terminal = nucleotide(vec![
            atom("DG", 1, "C3'", Vec3::new(-1.4, 0.0, 0.0)),
            atom("DG", 1, "O3'", Vec3::ZERO),
        ]);
        assert!(hydrogens_on(&terminal, "O3'").is_some());
    }

    /// A neighbour absent by chemistry cannot be told from one absent through
    /// disorder, so a missing bond never frees a hydrogen.
    #[test]
    fn a_missing_neighbour_does_not_invent_a_hydrogen() {
        let structure = nucleotide(vec![
            atom("DG", 1, "C5'", Vec3::new(-1.4, 0.0, 0.0)),
            atom("DG", 1, "O5'", Vec3::ZERO),
        ]);

        assert!(hydrogens_on(&structure, "O5'").is_none());
    }

    /// A definition can call a centre planar through a partner the coordinates
    /// never modelled.
    #[test]
    fn a_planar_centre_whose_conjugated_partner_is_missing_turns_instead() {
        let structure = ligand(vec![
            // N4 is planar through C3, which carries the carbonyl. Model only
            // the amide nitrogen's other side.
            ligand_atom("C5", "C", Vec3::new(-1.5, 0.0, 0.0)),
            ligand_atom("N6", "N", Vec3::ZERO),
        ]);

        assert!(
            matches!(
                hydrogens_on(&structure, "N6"),
                Some(Hydrogens::Rotor { .. })
            ),
            "a plane taken from an sp3 stand-in is worse than no plane"
        );
    }

    #[test]
    fn bisector_hydrogen_points_away_from_both_neighbours() {
        // Neighbours symmetric about the x axis; the hydrogen must fall on -x.
        let donor = Vec3::ZERO;
        let h = outward(
            donor,
            &[Vec3::new(1.0, 1.0, 0.0), Vec3::new(1.0, -1.0, 0.0)],
        )
        .unwrap();

        assert!((h - Vec3::new(-1.0, 0.0, 0.0)).length() < 1e-5);
        assert!((h - donor).length() - H_BOND_LENGTH < 1e-5);
    }

    #[test]
    fn planar_pair_sits_in_plane_at_the_sp2_angle() {
        // An amide: neighbour on -x, the atom beyond it off in +y.
        let donor = Vec3::ZERO;
        let neighbour = Vec3::new(-1.5, 0.0, 0.0);
        let beyond = Vec3::new(-2.5, 1.0, 0.0);
        let [h1, h2] = planar_pair(donor, neighbour, beyond).unwrap();

        for h in [h1, h2] {
            assert!(((h - donor).length() - H_BOND_LENGTH).abs() < 1e-4);
            assert!((angle_at(donor, neighbour, h) - SP2_H_ANGLE).abs() < 1e-3);
            assert!(h.z.abs() < 1e-5, "hydrogen left the plane: {h:?}");
        }
        assert!((h1 - h2).length() > 1e-3, "the two hydrogens coincide");
    }

    #[test]
    fn a_rotor_aims_its_hydrogen_at_the_acceptor() {
        // A hydroxyl: antecedent on -x, donor at the origin.
        let donor = Vec3::ZERO;
        let circle = Circle::new(donor, Vec3::new(-1.4, 0.0, 0.0), SP3_H_ANGLE).unwrap();

        // An acceptor swung around the axis must drag the hydrogen with it, and
        // the angle at the hydrogen must stay the best the rotor can manage.
        let near = circle.toward(Vec3::new(1.0, 2.8, 0.0));
        assert!(near.z.abs() < 1e-5);
        assert!(near.y > 0.0, "hydrogen did not turn towards the acceptor");

        let far = circle.toward(Vec3::new(1.0, -2.8, 0.0));
        assert!(far.y < 0.0, "hydrogen did not follow the acceptor round");

        // The bond angle is preserved whatever the torsion.
        for h in [near, far] {
            assert!(((h - donor).length() - H_BOND_LENGTH).abs() < 1e-4);
            assert!((angle_at(donor, Vec3::new(-1.4, 0.0, 0.0), h) - SP3_H_ANGLE).abs() < 1e-3);
        }
    }

    /// The rotor is a real gate, not a vacuous one. A hydrogen free to spin is
    /// still pinned to a cone, so how straight a bond it can make depends on
    /// where the acceptor sits relative to that cone: perfect on it, worse off
    /// it, hopeless behind it. Turning the torsion cannot rescue the last case.
    #[test]
    fn what_a_rotor_can_reach_depends_on_the_angle_to_its_bond() {
        let donor = Vec3::ZERO;
        let antecedent = Vec3::new(-1.4, 0.0, 0.0);
        let rotor = || Hydrogens::Rotor {
            circle: Circle::new(donor, antecedent, SP3_H_ANGLE).unwrap(),
            count: 1,
        };

        // On the cone: the hydrogen lands on the donor-acceptor line exactly.
        let cone_direction = Vec3::new(
            (180.0f32 - SP3_H_ANGLE).to_radians().cos(),
            (180.0f32 - SP3_H_ANGLE).to_radians().sin(),
            0.0,
        );
        let angle = rotor().angles_toward(donor, cone_direction * 2.8)[0];
        assert!(angle > 179.0, "expected a straight bond, got {angle}");

        // Off the cone: still acceptable, but measurably bent.
        let angle = rotor().angles_toward(donor, Vec3::new(2.0, 1.5, 0.0))[0];
        assert!(
            (100.0..150.0).contains(&angle),
            "expected a bent but passing bond, got {angle}"
        );

        // Behind the donor: unreachable at any torsion, so the gate rejects it.
        let angle = rotor().angles_toward(donor, Vec3::new(-3.0, 0.2, 0.0))[0];
        assert!(angle < 100.0, "expected a rejected bond, got {angle}");
    }

    #[test]
    fn capacity_is_one_bond_per_hydrogen() {
        let structure = protein(vec![
            atom("LYS", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("LYS", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("LYS", 1, "CB", Vec3::new(2.0, 1.4, 0.0)),
            atom("LYS", 1, "CG", Vec3::new(3.5, 1.4, 0.0)),
            atom("LYS", 1, "CD", Vec3::new(4.0, 2.8, 0.0)),
            atom("LYS", 1, "CE", Vec3::new(5.5, 2.8, 0.0)),
            atom("LYS", 1, "NZ", Vec3::new(6.0, 4.2, 0.0)),
        ]);

        let found = placed(&structure);
        let nz = found
            .iter()
            .find(|d| structure.atoms[d.atom].name == "NZ")
            .expect("the ammonium donates");
        assert_eq!(nz.hydrogens.capacity(), 3);
    }

    /// Proline in the middle of a chain: a preceding carbonyl, but a tertiary
    /// nitrogen with no hydrogen on it to donate.
    #[test]
    fn proline_never_donates() {
        let structure = protein(vec![
            atom("PRO", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("PRO", 2, "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("PRO", 2, "CA", Vec3::new(4.0, 2.8, 0.0)),
        ]);

        assert!(placed(&structure).is_empty());
    }

    /// A chain's first nitrogen is an ammonium, not an amide. It has no
    /// preceding carbonyl and needs none: three hydrogens turn about its N-CA
    /// bond. Dropping it would lose a charged donor from every chain.
    #[test]
    fn a_chain_starts_with_an_ammonium_not_an_amide() {
        let structure = protein(vec![
            atom("ALA", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("ALA", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("ALA", 1, "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("ALA", 1, "O", Vec3::new(1.2, 2.3, 0.0)),
        ]);

        let found = placed(&structure);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].hydrogens.capacity(), 3);
        assert!(
            matches!(found[0].hydrogens, Hydrogens::Rotor { .. }),
            "the terminal ammonium turns about its bond to CA"
        );
    }

    /// A residue whose predecessor is too far for a peptide bond has a genuine
    /// amide hydrogen and nothing to orient it against, so it cannot donate.
    #[test]
    fn a_residue_after_a_chain_break_does_not_donate() {
        let structure = protein(vec![
            atom("PRO", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", Vec3::new(1.2, 2.3, 0.0)),
            // Far past MAX_PEPTIDE_BOND from residue 1's carbonyl carbon.
            atom("ALA", 2, "N", Vec3::new(30.0, 1.5, 0.0)),
            atom("ALA", 2, "CA", Vec3::new(31.0, 2.8, 0.0)),
        ]);

        assert!(placed(&structure).is_empty());
    }

    /// A predecessor bonded but modelled without its carbonyl oxygen leaves the
    /// amide hydrogen with no direction. Treating that as a chain start would
    /// hand a one-hydrogen amide three hydrogens free to swing anywhere, which
    /// reaches acceptors the real hydrogen cannot.
    #[test]
    fn an_amide_whose_preceding_carbonyl_oxygen_is_missing_does_not_donate() {
        let structure = protein(vec![
            atom("PRO", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", Vec3::new(2.0, 1.4, 0.0)),
            // The carbonyl oxygen of residue 1 was not modelled.
            atom("ALA", 2, "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("ALA", 2, "CA", Vec3::new(4.0, 2.8, 0.0)),
        ]);

        assert!(placed(&structure).is_empty());
    }

    /// A water or ligand between two bonded residues must not make the second
    /// look like the start of a chain.
    #[test]
    fn a_residue_without_a_backbone_carbon_does_not_break_the_chain() {
        let structure = protein(vec![
            atom("PRO", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("HOH", 90, "O", Vec3::new(20.0, 20.0, 20.0)),
            atom("ALA", 2, "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("ALA", 2, "CA", Vec3::new(4.0, 2.8, 0.0)),
        ]);

        let found = placed(&structure);
        assert_eq!(found.len(), 1);
        assert!(
            matches!(found[0].hydrogens, Hydrogens::Fixed(_)),
            "the amide is still interior, not a terminal ammonium"
        );
    }

    #[test]
    fn an_interior_amide_donates_opposite_the_preceding_carbonyl() {
        // Residue 1 is a proline so that only residue 2 donates.
        let structure = protein(vec![
            atom("PRO", 1, "N", Vec3::new(0.0, 0.0, 0.0)),
            atom("PRO", 1, "CA", Vec3::new(1.5, 0.0, 0.0)),
            atom("PRO", 1, "C", Vec3::new(2.0, 1.4, 0.0)),
            atom("PRO", 1, "O", Vec3::new(1.2, 2.3, 0.0)),
            atom("ALA", 2, "N", Vec3::new(3.3, 1.5, 0.0)),
            atom("ALA", 2, "CA", Vec3::new(4.0, 2.8, 0.0)),
        ]);

        let found = placed(&structure);
        assert_eq!(found.len(), 1);

        let Hydrogens::Fixed(positions) = &found[0].hydrogens else {
            panic!("an amide hydrogen is fixed by the heavy atoms");
        };
        assert_eq!(positions.len(), 1);

        // H = N + unit(C_prev - O_prev), the rule DSSP uses.
        let expected = Vec3::new(3.3, 1.5, 0.0)
            + (Vec3::new(2.0, 1.4, 0.0) - Vec3::new(1.2, 2.3, 0.0)).normalize();
        assert!((positions[0] - expected).length() < 1e-5);
    }
}
