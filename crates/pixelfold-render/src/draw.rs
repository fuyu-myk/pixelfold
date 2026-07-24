/// Map B-factor value to RGB color using blue –> cyan –> yellow –> red gradient
/// Normalizes B-factors using percentile-based scaling
pub fn bfactor_to_color(b_factor: f32, b_min: f32, b_max: f32) -> (u8, u8, u8) {
    if b_max <= b_min {
        return (128, 128, 128); // Gray for invalid range
    }

    let t = ((b_factor - b_min) / (b_max - b_min)).clamp(0.0, 1.0);

    // Blue –> Cyan –> Yellow –> Red gradient
    // Blue (0, 100, 255) at t = 0
    // Cyan (0, 255, 255) at t = 0.33
    // Yellow (255, 255, 0) at t = 0.67
    // Red (255, 0, 0) at t = 1.0

    let (r, g, b) = if t < 0.33 {
        // Blue to Cyan
        let local_t = t / 0.33;
        let r = 0.0;
        let g = 100.0 + (255.0 - 100.0) * local_t;
        let b = 255.0;
        (r, g, b)
    } else if t < 0.67 {
        // Cyan to Yellow
        let local_t = (t - 0.33) / 0.34;
        let r = 255.0 * local_t;
        let g = 255.0;
        let b = 255.0 * (1.0 - local_t);
        (r, g, b)
    } else {
        // Yellow to Red
        let local_t = (t - 0.67) / 0.33;
        let r = 255.0;
        let g = 255.0 * (1.0 - local_t);
        let b = 0.0;
        (r, g, b)
    };

    (r as u8, g as u8, b as u8)
}

/// CPK/Jmol colour for an element symbol, keyed on the element column (so a
/// calcium ion is green, not carbon grey). Values are the Jmol defaults; an
/// unrecognised element falls back to the Jmol "unknown" pink.
pub fn element_to_color(element: &str) -> (u8, u8, u8) {
    match element.to_uppercase().as_str() {
        "H" => (255, 255, 255),
        "C" => (144, 144, 144),
        "N" => (48, 80, 248),
        "O" => (255, 13, 13),
        "F" => (144, 224, 80),
        "P" => (255, 128, 0),
        "S" => (255, 255, 48),
        "CL" => (31, 240, 31),
        "BR" => (166, 41, 41),
        "I" => (148, 0, 148),
        "FE" => (224, 102, 51),
        "ZN" => (125, 128, 176),
        "MG" => (138, 255, 0),
        "CA" => (61, 255, 0),
        "NA" => (171, 92, 242),
        "K" => (143, 64, 212),
        "MN" => (156, 122, 199),
        _ => (255, 192, 203),
    }
}

/// Draws a line between two points using Bresenham-like algorithm
pub fn draw_line(x0: f32, y0: f32, x1: f32, y1: f32) -> Vec<(f64, f64)> {
    let mut points = Vec::new();

    let dx = (x1 - x0).abs();
    let dy = (y1 - y0).abs();

    let sx = if x0 < x1 { 1.0 } else { -1.0 };
    let sy = if y0 < y1 { 1.0 } else { -1.0 };

    let mut err = dx - dy;
    let mut x = x0;
    let mut y = y0;

    // Starting point
    points.push((x as f64, y as f64));

    // Draw line using Bresenham's algorithm
    let max_steps = (dx.max(dy) as usize).max(1) + 2; // Safety limit

    for _ in 0..max_steps {
        if (x - x1).abs() < 0.5 && (y - y1).abs() < 0.5 {
            // Add endpoint if not already close to last point
            if let Some(&last) = points.last() {
                let dist = ((last.0 - x1 as f64).powi(2) + (last.1 - y1 as f64).powi(2)).sqrt();
                if dist > 1.0 {
                    points.push((x1 as f64, y1 as f64));
                }
            }
            break;
        }

        let e2 = 2.0 * err;

        // Move in x direction
        if e2 > -dy {
            err -= dy;
            x += sx;
        }

        // Move in y direction
        if e2 < dx {
            err += dx;
            y += sy;
        }

        points.push((x as f64, y as f64));
    }

    points
}

/// Map H-bond energy to color (cyan -> yellow -> orange)
///
/// Weak bonds (~ -0.5 kcal/mol): cyan;
/// Medium bonds (~ -2.0 kcal/mol): yellow;
/// Strong bonds (~ -5.0 kcal/mol): orange
pub fn hbond_energy_to_color(energy: f32) -> (u8, u8, u8) {
    // Normalization
    let t = ((energy.abs() - 0.5) / 4.5).clamp(0.0, 1.0);

    if t < 0.5 {
        // Cyan (0, 255, 255) to Yellow (255, 255, 0)
        let local_t = t / 0.5;
        let r = (255.0 * local_t) as u8;
        let g = 255;
        let b = (255.0 * (1.0 - local_t)) as u8;
        (r, g, b)
    } else {
        // Yellow (255, 255, 0) to Orange (255, 165, 0)
        let local_t = (t - 0.5) / 0.5;
        let r = 255;
        let g = (255.0 - 90.0 * local_t) as u8;
        let b = 0;
        (r, g, b)
    }
}

/// Draw dashed line by sampling every nth point
pub fn draw_dashed_line(
    x0: f32,
    y0: f32,
    x1: f32,
    y1: f32,
    dash_spacing: usize,
) -> Vec<(f64, f64)> {
    let full_line = draw_line(x0, y0, x1, y1);
    full_line
        .into_iter()
        .enumerate()
        .filter(|(i, _)| i % dash_spacing == 0)
        .map(|(_, p)| p)
        .collect()
}
