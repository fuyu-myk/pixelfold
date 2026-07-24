//! The single RGBA8 framebuffer every render path writes into.
//!
//! Colour, depth, and atom id are parallel buffers in row-major order with a
//! top-left origin (row 0 is the top). The rasteriser owns the flip from the
//! camera's lower-left screen space; this type only stores pixels.

/// Atom id written for pixels no fragment has covered.
pub const NO_ID: u32 = u32::MAX;

/// Depth of an uncovered pixel. Larger depth is nearer the viewer, so the far
/// plane is the most negative value and any real fragment passes the z-test.
const FAR_DEPTH: f32 = f32::NEG_INFINITY;

pub struct Framebuffer {
    width: u32,
    height: u32,
    background: [u8; 4],
    color: Vec<[u8; 4]>,
    depth: Vec<f32>,
    ids: Vec<u32>,
}

impl Framebuffer {
    /// A framebuffer cleared to `background`, the far depth, and [`NO_ID`].
    pub fn new(width: u32, height: u32, background: [u8; 4]) -> Self {
        let len = width as usize * height as usize;
        Self {
            width,
            height,
            background,
            color: vec![background; len],
            depth: vec![FAR_DEPTH; len],
            ids: vec![NO_ID; len],
        }
    }

    pub fn width(&self) -> u32 {
        self.width
    }

    pub fn height(&self) -> u32 {
        self.height
    }

    pub fn background(&self) -> [u8; 4] {
        self.background
    }

    /// RGBA8 pixels, row-major, top-left origin. The natural PNG layout.
    pub fn color(&self) -> &[[u8; 4]] {
        &self.color
    }

    pub fn ids(&self) -> &[u32] {
        &self.ids
    }

    pub fn depth(&self) -> &[f32] {
        &self.depth
    }

    fn index(&self, x: u32, y: u32) -> usize {
        y as usize * self.width as usize + x as usize
    }

    pub fn pixel(&self, x: u32, y: u32) -> [u8; 4] {
        self.color[self.index(x, y)]
    }

    /// The atom id at a pixel, or [`NO_ID`] where the background shows through.
    /// Exact pixel picking is one call to this against the click position.
    pub fn id_at(&self, x: u32, y: u32) -> u32 {
        self.ids[self.index(x, y)]
    }

    pub fn depth_at(&self, x: u32, y: u32) -> f32 {
        self.depth[self.index(x, y)]
    }

    /// Depth-test a fragment and, if it is nearer than what is stored, write its
    /// colour, depth, and atom id. Returns whether the fragment was kept.
    pub fn test_and_set(&mut self, x: u32, y: u32, depth: f32, color: [u8; 4], id: u32) -> bool {
        let i = self.index(x, y);
        if depth > self.depth[i] {
            self.depth[i] = depth;
            self.color[i] = color;
            self.ids[i] = id;
            true
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_clears_to_background_far_and_no_id() {
        let bg = [10, 20, 30, 255];
        let fb = Framebuffer::new(4, 3, bg);

        assert_eq!(fb.width(), 4);
        assert_eq!(fb.height(), 3);
        assert_eq!(fb.color().len(), 12);
        assert!(fb.color().iter().all(|&p| p == bg));
        assert!(fb.ids().iter().all(|&id| id == NO_ID));
        assert!(fb.depth().iter().all(|&d| d == FAR_DEPTH));
    }

    #[test]
    fn nearer_fragment_wins_and_farther_is_rejected() {
        let mut fb = Framebuffer::new(2, 2, [0, 0, 0, 255]);

        assert!(fb.test_and_set(1, 1, 0.0, [255, 0, 0, 255], 7));
        assert_eq!(fb.pixel(1, 1), [255, 0, 0, 255]);
        assert_eq!(fb.id_at(1, 1), 7);

        // Farther (smaller depth) loses.
        assert!(!fb.test_and_set(1, 1, -1.0, [0, 255, 0, 255], 9));
        assert_eq!(fb.pixel(1, 1), [255, 0, 0, 255]);
        assert_eq!(fb.id_at(1, 1), 7);

        // Nearer (larger depth) wins.
        assert!(fb.test_and_set(1, 1, 5.0, [0, 0, 255, 255], 3));
        assert_eq!(fb.pixel(1, 1), [0, 0, 255, 255]);
        assert_eq!(fb.id_at(1, 1), 3);
    }

    #[test]
    fn pixels_are_addressed_row_major() {
        let mut fb = Framebuffer::new(3, 2, [0, 0, 0, 255]);
        fb.test_and_set(2, 1, 0.0, [1, 2, 3, 255], 1);

        // The written pixel is the last one in row-major order.
        assert_eq!(fb.color()[5], [1, 2, 3, 255]);
        assert_eq!(fb.pixel(2, 1), [1, 2, 3, 255]);
    }
}
