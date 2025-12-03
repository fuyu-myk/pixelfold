use glam::{Mat4, Vec3, Vec4};

use crate::{Protein, SecondaryStructure};


pub struct Camera {
    pub position: Vec3,
    pub rotation: Vec3, // Euler angles (pitch, yaw, roll)
    pub zoom: f32,
    pub pan: Vec3,
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            position: Vec3::new(0.0, 0.0, 100.0),
            rotation: Vec3::ZERO,
            zoom: 1.0,
            pan: Vec3::ZERO,
        }
    }
}

impl Camera {
    pub fn new() -> Self {
        Self::default()
    }

    /// Rotate camera by delta angles (in radians)
    pub fn rotate(&mut self, pitch: f32, yaw: f32, roll: f32) {
        self.rotation.x += pitch;
        self.rotation.y += yaw;
        self.rotation.z += roll;
    }

    /// Zoom in (positive) or out (negative)
    pub fn adjust_zoom(&mut self, delta: f32) {
        self.zoom = (self.zoom + delta).max(0.1).min(10.0);
    }

    /// Pan camera in screen space
    pub fn pan_camera(&mut self, dx: f32, dy: f32) {
        self.pan.x += dx;
        self.pan.y += dy;
    }

    /// Get the view-projection matrix for transforming 3D points
    pub fn get_view_matrix(&self) -> Mat4 {
        // Create rotation matrix from Euler angles
        let rotation_x = Mat4::from_rotation_x(self.rotation.x);
        let rotation_y = Mat4::from_rotation_y(self.rotation.y);
        let rotation_z = Mat4::from_rotation_z(self.rotation.z);
        let rotation = rotation_z * rotation_y * rotation_x;

        let center_translation = Mat4::from_translation(self.pan);
        let camera_translation = Mat4::from_translation(self.position);

        camera_translation * rotation * center_translation
    }

    /// Project a 3D point to 2D screen coordinates
    pub fn project_point(&self, point: Vec3, width: f32, height: f32) -> (f32, f32) {
        let view_matrix = self.get_view_matrix();
        
        // Transform point by view matrix
        let transformed = view_matrix * Vec4::new(point.x, point.y, point.z, 1.0);
        
        // Apply orthographic projection with zoom
        let x = (transformed.x * self.zoom) + width / 2.0;
        let y = (transformed.y * self.zoom) + height / 2.0;
        
        (x, y)
    }

    /// Project a 3D point with depth information (for depth sorting)
    pub fn project_point_with_depth(&self, point: Vec3, width: f32, height: f32) -> (f32, f32, f32) {
        let view_matrix = self.get_view_matrix();
        
        let transformed = view_matrix * Vec4::new(point.x, point.y, point.z, 1.0);
        
        let x = (transformed.x * self.zoom) + width / 2.0;
        let y = (transformed.y * self.zoom) + height / 2.0;
        let z = transformed.z; // Keep depth for sorting
        
        (x, y, z)
    }
}

/// Represents a projected atom for rendering
pub struct ProjectedAtom {
    pub x: f32,
    pub y: f32,
    pub depth: f32,
    pub name: String,
    pub residue_name: String,
    pub secondary_structure: SecondaryStructure,
}

/// Project all atoms in a protein to 2D screen space
pub fn project_protein(protein: &Protein, camera: &Camera, width: f32, height: f32) -> Vec<ProjectedAtom> {
    protein
        .atoms
        .iter()
        .map(|atom| {
            let (x, y, depth) = camera.project_point_with_depth(atom.position, width, height);
            ProjectedAtom {
                x,
                y,
                depth,
                name: atom.name.clone(),
                residue_name: atom.residue_name.clone(),
                secondary_structure: atom.secondary_structure,
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

/// Auto-frame the camera to fit the protein in view
pub fn auto_frame_protein(protein: &Protein, camera: &mut Camera, width: f32, height: f32) {
    let center = calculate_center_of_mass(protein);
    let (min, max) = calculate_bounds(protein);
    
    // Calculate the extent of the protein
    let extent = max - min;
    let max_extent = extent.x.max(extent.y).max(extent.z);
    
    camera.pan = -center;
    
    // Adjust zoom to fit protein in view with padding
    // Use 1.4 as padding factor to leave ~20% space around the protein on each side
    let padding_factor = 1.4;
    let viewport_size = width.min(height);
    if max_extent > 0.0 {
        camera.zoom = viewport_size / (max_extent * padding_factor);
    }
}
