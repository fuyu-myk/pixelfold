//! Framing a pose: the camera turned so its ligand is seen past the protein
//! around it.

use glam::{Quat, Vec3};
use pixelfold_core::Protein;
use pixelfold_render::renderer::{self, Camera};

/// How many directions a pose's ligand is weighed for viewing from.
const VIEWS: usize = 256;

/// How far beyond the ligand an atom still hides it: a Cα sphere's radius and
/// a little.
const SIGHT_MARGIN: f32 = 2.0;

/// Frame `protein` in `width` by `height` pixels, turned so the `ligand`
/// atoms are seen along the line of sight the trace blocks least, which a
/// pocket otherwise hides them from.
pub(crate) fn frame_pose(
    camera: &mut Camera,
    protein: &Protein,
    ligand: &[usize],
    width: f32,
    height: f32,
) {
    renderer::auto_frame_protein(protein, camera, width, height);
    let atoms: Vec<Vec3> = ligand
        .iter()
        .filter_map(|&atom| protein.atoms.get(atom))
        .map(|atom| atom.position)
        .collect();
    if atoms.is_empty() {
        return;
    }

    let centre = atoms.iter().sum::<Vec3>() / atoms.len() as f32;
    let reach = atoms
        .iter()
        .map(|&atom| atom.distance(centre))
        .fold(0.0, f32::max)
        + SIGHT_MARGIN;
    let trace: Vec<Vec3> = protein
        .atoms
        .iter()
        .enumerate()
        .filter(|&(index, atom)| atom.name == "CA" && ligand.binary_search(&index).is_err())
        .map(|(_, atom)| atom.position)
        .collect();
    let view = clearest_view(centre, reach, &trace, centre - camera.center);
    camera.orientation = Quat::from_rotation_arc(view, Vec3::Z);
}

/// Of [`VIEWS`] directions spread over a sphere, the one to view a ligand at
/// `centre` along whose line of sight the fewest `blockers` sit within `reach`
/// of it; of those as clear, the one fewest sit within twice that of, which
/// looks down the middle of a pocket's mouth rather than past its rim; then
/// the nearest `outward`.
fn clearest_view(centre: Vec3, reach: f32, blockers: &[Vec3], outward: Vec3) -> Vec3 {
    let blocked = |view: Vec3| {
        let near = |within: f32| {
            blockers
                .iter()
                .filter(|&&atom| {
                    let offset = atom - centre;
                    let along = offset.dot(view);
                    along > 0.0 && (offset - view * along).length_squared() < within * within
                })
                .count()
        };
        (near(reach), near(2.0 * reach))
    };
    let outward = outward.normalize_or_zero();

    directions()
        .map(|view| (view, blocked(view)))
        .min_by(|(a, a_blocked), (b, b_blocked)| {
            a_blocked
                .cmp(b_blocked)
                .then(b.dot(outward).total_cmp(&a.dot(outward)))
        })
        .map(|(view, _)| view)
        .expect("there are views to weigh")
}

/// [`VIEWS`] unit vectors spread evenly over a sphere (a Fibonacci lattice).
fn directions() -> impl Iterator<Item = Vec3> {
    let turn = std::f32::consts::PI * (3.0 - 5.0_f32.sqrt());
    (0..VIEWS).map(move |i| {
        let z = 1.0 - 2.0 * (i as f32 + 0.5) / VIEWS as f32;
        let radius = (1.0 - z * z).sqrt();
        let angle = turn * i as f32;
        Vec3::new(radius * angle.cos(), radius * angle.sin(), z)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_clearest_view_of_a_ligand_is_through_its_pocket_s_mouth() {
        // A bowl around the origin open towards +x, and far off along -y the
        // rest of the protein, which pulls its centre that way.
        let bowl = directions()
            .filter(|view| view.x < 0.5)
            .map(|view| view * 6.0);
        let bulk = (0..200).map(|i| Vec3::new(0.0, -30.0 - i as f32 * 0.1, 0.0));
        let atoms: Vec<Vec3> = bowl.chain(bulk).collect();
        let centre = atoms.iter().sum::<Vec3>() / atoms.len() as f32;

        let view = clearest_view(Vec3::ZERO, 2.0, &atoms, Vec3::ZERO - centre);

        assert!(view.x > 0.7, "{view}");
    }
}
