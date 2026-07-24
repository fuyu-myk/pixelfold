use glam::{Quat, Vec2, Vec3};

use pixelfold_core::{FixedStr, Protein, SecondaryStructure};

pub struct Camera {
    /// Accumulated orientation as a unit quaternion.
    pub orientation: Quat,
    /// Rotation pivot: the structure center, set by `auto_frame_protein`.
    pub center: Vec3,
    /// Screen-space pan offset in pixels, independent of orientation.
    pub pan: Vec2,
    /// Current zoom (screen pixels per Angstrom).
    pub zoom: f32,
    /// Zoom that fits the structure, set by `auto_frame_protein`. The manual
    /// zoom clamp is relative to this, so it works at any structure scale.
    pub fit_zoom: f32,
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            orientation: Quat::IDENTITY,
            center: Vec3::ZERO,
            pan: Vec2::ZERO,
            zoom: 1.0,
            fit_zoom: 1.0,
        }
    }
}

impl Camera {
    pub fn new() -> Self {
        Self::default()
    }

    /// Rotate about the screen axes (trackball feel). The incremental rotation
    /// is pre-multiplied, so it composes in camera space regardless of the
    /// current orientation, and the quaternion is renormalised to avoid drift.
    pub fn rotate(&mut self, pitch: f32, yaw: f32, roll: f32) {
        let delta =
            Quat::from_rotation_x(pitch) * Quat::from_rotation_y(yaw) * Quat::from_rotation_z(roll);
        self.orientation = (delta * self.orientation).normalize();
    }

    /// Zoom in (positive delta) or out (negative), scaling the current zoom so
    /// the step stays proportional at any scale, clamped around the fit zoom.
    pub fn adjust_zoom(&mut self, delta: f32) {
        self.zoom = (self.zoom * (1.0 + delta)).clamp(self.fit_zoom * 0.1, self.fit_zoom * 10.0);
    }

    /// Pan the view in screen space (pixels), independent of orientation.
    pub fn pan_camera(&mut self, dx: f32, dy: f32) {
        self.pan.x += dx;
        self.pan.y += dy;
    }

    /// Project a 3D point to screen coordinates plus a depth value for sorting.
    ///
    /// The point is centered on `center`, rotated by `orientation`, scaled by
    /// `zoom`, then offset by the screen-space `pan`. Because `pan` is added
    /// after projection, panning always moves the view along the screen axes.
    /// Screen y increases upward (the ratatui Canvas has a lower-left origin).
    pub fn project_point_with_depth(
        &self,
        point: Vec3,
        width: f32,
        height: f32,
    ) -> (f32, f32, f32) {
        let view = self.orientation * (point - self.center);
        let x = view.x * self.zoom + width / 2.0 + self.pan.x;
        let y = view.y * self.zoom + height / 2.0 + self.pan.y;

        (x, y, view.z)
    }
}

#[derive(Clone)]
/// Represents a projected atom for rendering
pub struct ProjectedAtom {
    pub x: f32,
    pub y: f32,
    pub depth: f32,
    pub name: FixedStr<4>,
    pub residue_name: FixedStr<6>,
    pub secondary_structure: SecondaryStructure,
}

/// Represents a projected surface point for rendering
pub struct ProjectedSurfacePoint {
    pub x: f32,
    pub y: f32,
    pub depth: f32,
    pub hydrophobicity: f32,
}

/// Project all atoms in a protein to 2D screen space
pub fn project_protein(
    protein: &Protein,
    camera: &Camera,
    width: f32,
    height: f32,
) -> Vec<ProjectedAtom> {
    protein
        .atoms
        .iter()
        .map(|atom| {
            let (x, y, depth) = camera.project_point_with_depth(atom.position, width, height);
            ProjectedAtom {
                x,
                y,
                depth,
                name: atom.name,
                residue_name: atom.residue_name,
                secondary_structure: atom.secondary_structure,
            }
        })
        .collect()
}

/// Project all surface points in a protein to 2D screen space
pub fn project_surface(
    protein: &Protein,
    camera: &Camera,
    width: f32,
    height: f32,
) -> Vec<ProjectedSurfacePoint> {
    protein
        .surface_points
        .iter()
        .map(|surface_point| {
            let (x, y, depth) =
                camera.project_point_with_depth(surface_point.position, width, height);
            ProjectedSurfacePoint {
                x,
                y,
                depth,
                hydrophobicity: surface_point.hydrophobicity,
            }
        })
        .collect()
}

/// Calculate the center of mass of a protein
pub fn calculate_center_of_mass(protein: &Protein) -> Vec3 {
    if protein.atoms.is_empty() {
        return Vec3::ZERO;
    }

    let sum = protein
        .atoms
        .iter()
        .fold(Vec3::ZERO, |acc, atom| acc + atom.position);

    sum / protein.atoms.len() as f32
}

/// Calculate bounding box of protein
pub fn calculate_bounds(protein: &Protein) -> (Vec3, Vec3) {
    if protein.atoms.is_empty() {
        return (Vec3::ZERO, Vec3::ZERO);
    }

    let mut min = Vec3::splat(f32::INFINITY);
    let mut max = Vec3::splat(f32::NEG_INFINITY);

    for atom in &protein.atoms {
        min = min.min(atom.position);
        max = max.max(atom.position);
    }

    (min, max)
}

/// Auto-frame the camera to fit the protein: reset the orientation and pan,
/// pivot on the center of mass, and zoom to fit with padding.
pub fn auto_frame_protein(protein: &Protein, camera: &mut Camera, width: f32, height: f32) {
    let positions: Vec<Vec3> = protein.atoms.iter().map(|a| a.position).collect();
    auto_frame_points(&positions, camera, width, height);
}

/// Auto-frame the camera to fit an arbitrary set of points: reset orientation
/// and pan, pivot on their mean, and zoom to fit with padding. This is what
/// lets `render` frame a selected subset instead of the whole structure.
pub fn auto_frame_points(positions: &[Vec3], camera: &mut Camera, width: f32, height: f32) {
    if positions.is_empty() {
        return;
    }

    let center = positions.iter().copied().sum::<Vec3>() / positions.len() as f32;
    let mut min = Vec3::splat(f32::INFINITY);
    let mut max = Vec3::splat(f32::NEG_INFINITY);
    for &p in positions {
        min = min.min(p);
        max = max.max(p);
    }

    let extent = max - min;
    let max_extent = extent.x.max(extent.y).max(extent.z);

    camera.orientation = Quat::IDENTITY;
    camera.center = center;
    camera.pan = Vec2::ZERO;

    // 1.4 leaves ~20% space around the structure on each side.
    let padding_factor = 1.4;
    let viewport_size = width.min(height);
    if max_extent > 0.0 {
        camera.fit_zoom = viewport_size / (max_extent * padding_factor);
        camera.zoom = camera.fit_zoom;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pixelfold_core::{Atom, SecondaryStructure};
    use std::f32::consts::PI;

    fn atom_at(x: f32) -> Atom {
        Atom {
            serial: 0,
            name: FixedStr::new("CA"),
            element: FixedStr::new("C"),
            residue_name: FixedStr::new("ALA"),
            residue_seq: 1,
            chain_id: FixedStr::new("A"),
            is_hetatm: false,
            altloc: None,
            occupancy: 1.0,
            insertion_code: None,
            model_number: 1,
            position: Vec3::new(x, 0.0, 0.0),
            b_factor: 0.0,
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    #[test]
    fn center_projects_to_screen_center_regardless_of_orientation() {
        let mut cam = Camera::new();
        cam.center = Vec3::new(3.0, -2.0, 5.0);
        cam.zoom = 2.5;
        cam.rotate(0.7, 1.3, 0.2);

        let (x, y, _) = cam.project_point_with_depth(cam.center, 100.0, 80.0);
        assert!((x - 50.0).abs() < 1e-3, "x = {x}");
        assert!((y - 40.0).abs() < 1e-3, "y = {y}");
    }

    #[test]
    fn pan_is_screen_space_and_orientation_independent() {
        let mut cam = Camera::new();
        cam.center = Vec3::new(1.0, 2.0, 3.0);
        cam.zoom = 1.7;
        cam.rotate(0.9, -0.4, 1.1);

        let p = Vec3::new(10.0, -5.0, 4.0);
        let (x0, y0, _) = cam.project_point_with_depth(p, 200.0, 150.0);
        cam.pan_camera(12.0, -7.0);
        let (x1, y1, _) = cam.project_point_with_depth(p, 200.0, 150.0);

        assert!((x1 - x0 - 12.0).abs() < 1e-3, "dx = {}", x1 - x0);
        assert!((y1 - y0 + 7.0).abs() < 1e-3, "dy = {}", y1 - y0);
    }

    #[test]
    fn rotation_stays_normalized() {
        let mut cam = Camera::new();
        for _ in 0..100 {
            cam.rotate(0.3, 0.5, 0.2);
        }
        assert!((cam.orientation.length() - 1.0).abs() < 1e-4);
    }

    #[test]
    fn full_turn_returns_to_start() {
        let mut cam = Camera::new();
        let steps = 1000;
        let step = 2.0 * PI / steps as f32;
        for _ in 0..steps {
            cam.rotate(0.0, step, 0.0);
        }
        let dot = cam.orientation.dot(Quat::IDENTITY).abs();
        assert!(dot > 0.999, "dot = {dot}");
    }

    #[test]
    fn zoom_scales_offset_from_center() {
        let mut cam = Camera::new();
        let p = Vec3::new(4.0, 0.0, 0.0);
        cam.zoom = 1.0;
        let (x1, _, _) = cam.project_point_with_depth(p, 100.0, 100.0);
        cam.zoom = 2.0;
        let (x2, _, _) = cam.project_point_with_depth(p, 100.0, 100.0);
        assert!(((x2 - 50.0) - 2.0 * (x1 - 50.0)).abs() < 1e-3);
    }

    #[test]
    fn auto_frame_resets_orientation_and_pan() {
        let protein = Protein {
            atoms: vec![atom_at(-5.0), atom_at(5.0)],
            title: String::new(),
            surface_points: vec![],
            hbonds: vec![],
            assembly: None,
            components: Default::default(),
        };
        let mut cam = Camera::new();
        cam.rotate(1.0, 2.0, 3.0);
        cam.pan_camera(50.0, 30.0);

        auto_frame_protein(&protein, &mut cam, 100.0, 100.0);

        assert!(cam.orientation.dot(Quat::IDENTITY).abs() > 0.999);
        assert_eq!(cam.pan, Vec2::ZERO);
        assert_eq!(cam.center, Vec3::ZERO);
    }

    #[test]
    fn manual_zoom_is_relative_to_fit() {
        let protein = Protein {
            atoms: vec![atom_at(-4.0), atom_at(4.0)],
            title: String::new(),
            surface_points: vec![],
            hbonds: vec![],
            assembly: None,
            components: Default::default(),
        };
        let mut cam = Camera::new();
        auto_frame_protein(&protein, &mut cam, 152.0, 152.0);

        let fit = cam.zoom;
        assert!(fit > 10.0, "fit = {fit}"); // exceeds the old fixed clamp band
        cam.adjust_zoom(0.1);
        assert!(cam.zoom > fit, "zoom in should increase from the fit");
    }
}
