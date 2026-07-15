//! Uniform grid spatial index.
//!
//! Points are bucketed into cubic cells of a fixed size. When the cell size is
//! at least the largest query radius, a neighbour query only inspects a small
//! fixed block of cells, giving O(1) queries after an O(N) build. This powers
//! hydrogen-bond and interaction detection, `within` selections, and clash
//! checks.

use std::collections::HashMap;

use glam::Vec3;

type Cell = (i32, i32, i32);

/// A uniform grid over 3D points. The grid keeps its own copy of the positions
/// so queries are self-contained; index `i` in a query result refers to point
/// `i` in the slice the grid was built from.
pub struct Grid {
    inv_cell_size: f32,
    positions: Vec<Vec3>,
    cells: HashMap<Cell, Vec<usize>>,
}

impl Grid {
    /// Build a grid, bucketing every point index into its cell. Choose
    /// `cell_size` no smaller than the largest radius that will be queried.
    pub fn build(positions: &[Vec3], cell_size: f32) -> Self {
        assert!(cell_size > 0.0, "grid cell size must be positive");
        let inv_cell_size = 1.0 / cell_size;

        let mut cells: HashMap<Cell, Vec<usize>> = HashMap::new();
        for (i, &p) in positions.iter().enumerate() {
            cells.entry(cell_of(p, inv_cell_size)).or_default().push(i);
        }

        Self {
            inv_cell_size,
            positions: positions.to_vec(),
            cells,
        }
    }

    /// Collect the indices of points within `radius` of `center`.
    pub fn within(&self, center: Vec3, radius: f32) -> Vec<usize> {
        let mut out = Vec::new();
        self.for_each_within(center, radius, |i| out.push(i));
        out
    }

    /// Invoke `f` with the index of each point within `radius` of `center`, in
    /// no particular order.
    pub fn for_each_within(&self, center: Vec3, radius: f32, mut f: impl FnMut(usize)) {
        let r2 = radius * radius;
        let base = cell_of(center, self.inv_cell_size);
        // Enough cells so no point within `radius` is missed for any cell size.
        let reach = (radius * self.inv_cell_size).ceil() as i32;

        for dx in -reach..=reach {
            for dy in -reach..=reach {
                for dz in -reach..=reach {
                    let cell = (base.0 + dx, base.1 + dy, base.2 + dz);
                    let Some(indices) = self.cells.get(&cell) else {
                        continue;
                    };

                    for &i in indices {
                        if (self.positions[i] - center).length_squared() <= r2 {
                            f(i);
                        }
                    }
                }
            }
        }
    }
}

fn cell_of(p: Vec3, inv_cell_size: f32) -> Cell {
    (
        (p.x * inv_cell_size).floor() as i32,
        (p.y * inv_cell_size).floor() as i32,
        (p.z * inv_cell_size).floor() as i32,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sorted(mut v: Vec<usize>) -> Vec<usize> {
        v.sort_unstable();
        v
    }

    #[test]
    fn finds_points_inside_and_excludes_points_outside() {
        let points = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),   // just outside r=4
            Vec3::new(0.0, 3.0, 0.0),   // inside r=4
            Vec3::new(20.0, 20.0, 0.0), // far away
        ];
        let grid = Grid::build(&points, 4.0);
        assert_eq!(sorted(grid.within(Vec3::ZERO, 4.0)), vec![0, 1, 3]);
    }

    #[test]
    fn boundary_distance_is_inclusive() {
        let points = vec![Vec3::ZERO, Vec3::new(4.0, 0.0, 0.0)];
        let grid = Grid::build(&points, 4.0);
        assert_eq!(sorted(grid.within(Vec3::ZERO, 4.0)), vec![0, 1]);
    }

    #[test]
    fn handles_negative_coordinates() {
        let points = vec![
            Vec3::new(-10.0, -10.0, -10.0),
            Vec3::new(-11.0, -10.0, -10.0),
            Vec3::new(10.0, 10.0, 10.0),
        ];
        let grid = Grid::build(&points, 6.0);
        assert_eq!(
            sorted(grid.within(Vec3::new(-10.0, -10.0, -10.0), 3.0)),
            vec![0, 1]
        );
    }

    #[test]
    fn radius_larger_than_cell_size_scans_enough_cells() {
        // Points 9 A away must be found even with a 3 A cell size.
        let points = vec![Vec3::ZERO, Vec3::new(9.0, 0.0, 0.0)];
        let grid = Grid::build(&points, 3.0);
        assert_eq!(sorted(grid.within(Vec3::ZERO, 10.0)), vec![0, 1]);
        assert_eq!(grid.within(Vec3::ZERO, 8.0), vec![0]);
    }

    #[test]
    fn empty_grid_returns_nothing() {
        let grid = Grid::build(&[], 5.0);
        assert!(grid.within(Vec3::ZERO, 100.0).is_empty());
    }
}
