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
//! Sources:
//! - HBPLUS: McDonald and Thornton, J Mol Biol 1994; manual Tables II and III at
//!   `https://www.uoxray.uoregon.edu/local/manuals/hbplus/hbplus.txt`
//! - Reduce: Word et al., J Mol Biol 285:1735 (1999); `RotDonor::finalize` in
//!   `https://github.com/rlabduke/reduce`

use glam::Vec3;

use crate::structure::{Atom, Protein};

use super::topology;

/// Length of every inferred polar hydrogen bond (HBPLUS Table III).
const H_BOND_LENGTH: f32 = 1.0;
/// The DD-D-H angle at an sp2 donor carrying two hydrogens (HBPLUS Table III).
const SP2_H_ANGLE: f32 = 120.0;
/// The DD-D-H angle at an sp3 hydroxyl or ammonium (HBPLUS Table III).
const SP3_H_ANGLE: f32 = 110.0;
/// A C(prev)-N separation beyond this (A) is a chain break, not a peptide bond.
const MAX_PEPTIDE_BOND: f32 = 2.5;

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
}

impl Hydrogens {
    /// How many hydrogen bonds this donor can make at once: one per hydrogen.
    pub fn capacity(&self) -> usize {
        match self {
            Hydrogens::Fixed(positions) => positions.len(),
            Hydrogens::Rotor { count, .. } => *count,
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

/// One hydrogen on the bisector of `first`-`donor`-`second`, pointing away from
/// both, which is where an sp2 nitrogen's single hydrogen sits.
fn bisector(donor: Vec3, first: Vec3, second: Vec3) -> Option<Vec3> {
    let inward = (first - donor).try_normalize()? + (second - donor).try_normalize()?;

    Some(donor - inward.try_normalize()? * H_BOND_LENGTH)
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
pub fn donors(protein: &Protein) -> Vec<Donor> {
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
            let Some(geometry) = geometry_of(&protein.atoms[i]) else {
                continue;
            };

            if let Some(hydrogens) =
                place(protein, residue.clone(), i, geometry, preceding.as_ref())
            {
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

            Some(Hydrogens::Fixed(vec![bisector(donor, first, second)?]))
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
    fn bisector_hydrogen_points_away_from_both_neighbours() {
        // Neighbours symmetric about the x axis; the hydrogen must fall on -x.
        let donor = Vec3::ZERO;
        let h = bisector(donor, Vec3::new(1.0, 1.0, 0.0), Vec3::new(1.0, -1.0, 0.0)).unwrap();

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

        let found = donors(&structure);
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

        assert!(donors(&structure).is_empty());
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

        let found = donors(&structure);
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

        assert!(donors(&structure).is_empty());
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

        assert!(donors(&structure).is_empty());
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

        let found = donors(&structure);
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

        let found = donors(&structure);
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
