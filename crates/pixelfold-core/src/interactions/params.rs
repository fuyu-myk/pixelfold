//! Geometric parameters for interaction detection.
//!
//! Every value is taken from PLIP's `plip/basic/config.py`, with the exact PLIP
//! variable name in each doc comment. Where RING 4.0 uses a different value,
//! the difference is recorded: it is the expected reason pixelfold and RING
//! disagree on an edge, not a bug. Distances are Angstroms, angles are degrees.
//!
//! RING's values are quoted at its strict defaults. They are not in its paper:
//! the paper describes the method and leaves every threshold to the web
//! documentation, and RING's engine is closed source. Note also that RING 2.0
//! and RING 4.0 differ on both distance and angle, so a criterion attributed to
//! "RING" without a version is ambiguous; the version is named here each time.
//!
//! Sources:
//! - PLIP config.py: `https://raw.githubusercontent.com/pharmai/plip/master/plip/basic/config.py`
//! - RING 4.0 criteria and strict/relaxed defaults: `https://ring.biocomputingup.it/about`
//! - RING 4.0: Clementel et al., Nucleic Acids Research 2024 (gkae337)
//! - RING 2.0 (superseded, for the version differences): Piovesan et al., NAR 2016 (gkw315)

/// Shared lower distance bound applied to every distance gate (PLIP `MIN_DIST`).
pub const MIN_DIST: f32 = 0.5;

/// Max donor-acceptor heavy-atom distance D...A (PLIP `HBOND_DIST_MAX`, which is
/// Hubbard and Haider's 3.5 A plus 0.6). RING 4.0 gates the same distance at
/// 3.9 A and adds two constraints PLIP has no equivalent of: an H...A leg under
/// 2.5 A and a D-H...A angle over 90 degrees, so RING reports fewer.
pub const HBOND_DIST_MAX: f32 = 4.1;
/// Min D-H...A angle (PLIP `HBOND_DON_ANGLE_MIN`, Hubbard and Haider's 90
/// degrees plus 10).
///
/// Despite the name, and despite PLIP's own comment calling it the angle "at the
/// hydrogen bond donor", `detection.py` measures it at the **hydrogen**: it takes
/// the angle between the vectors H->D and H->A. RING measures the same angle at
/// the same vertex, at 90 degrees.
pub const HBOND_DONOR_ANGLE_MIN: f32 = 100.0;

/// Max charge-centre to charge-centre distance for a salt bridge (PLIP
/// `SALTBRIDGE_DIST_MAX`). RING uses a stricter ~4.0 A, atom-to-atom.
pub const SALTBRIDGE_DIST_MAX: f32 = 5.5;

/// Max apolar carbon-carbon distance for a hydrophobic contact (PLIP
/// `HYDROPH_DIST_MAX`), applied between atom centres.
///
/// RING has no comparable centre-to-centre threshold: it folds these into van
/// der Waals contacts measured surface to surface.
pub const HYDROPHOBIC_DIST_MAX: f32 = 4.0;

/// Max aromatic ring centroid distance for pi-stacking (PLIP `PISTACK_DIST_MAX`).
/// RING is more permissive at ~6.5 A.
pub const PISTACK_DIST_MAX: f32 = 5.5;
/// Max deviation from parallel or perpendicular ring orientation (PLIP
/// `PISTACK_ANG_DEV`).
pub const PISTACK_ANGLE_DEV: f32 = 30.0;
/// Max lateral centroid offset (PLIP `PISTACK_OFFSET_MAX`, benzene radius + 0.5).
pub const PISTACK_OFFSET_MAX: f32 = 2.0;

/// Max cation to aromatic-ring-centroid distance (PLIP `PICATION_DIST_MAX`). RING
/// is stricter at ~5.0 A.
pub const PICATION_DIST_MAX: f32 = 6.0;

/// Max halogen X...A distance (PLIP `HALOGEN_DIST_MAX`, Auffinger's van der
/// Waals sums widened by 0.5).
pub const HALOGEN_DIST_MAX: f32 = 4.0;

/// Elements that donate a halogen bond, each through the single carbon it hangs
/// from.
///
/// Fluorine is deliberately absent. PLIP admits F, Cl, Br and I alike
/// (`Mol.is_functional_group`, atomic numbers 9/17/35/53), but the paper PLIP
/// cites for these very thresholds excludes fluorine on physical grounds:
/// "Fluorine atoms remain entirely electronegative, whereas each of the other
/// three halogen atoms shows the emergence of an electropositive crown ...
/// fluorines are more likely to serve as hydrogen bond acceptors in F...H-O-type
/// interactions" (Auffinger et al., PNAS 2004).
pub const HALOGEN_ELEMENTS: &[&str] = &["CL", "BR", "I"];

/// Elements that accept a halogen bond, from PLIP `BindingSite.find_hal`
/// (atomic numbers 8, 7, 16).
pub const HALOGEN_ACCEPTOR_ELEMENTS: &[&str] = &["O", "N", "S"];

/// Elements that may stand as the acceptor's antecedent Y, from the same
/// function (atomic numbers 6, 7, 15, 16).
pub const HALOGEN_ANTECEDENT_ELEMENTS: &[&str] = &["C", "N", "P", "S"];
/// Optimal C-X...A donor angle (PLIP `HALOGEN_DON_ANGLE`).
pub const HALOGEN_DONOR_ANGLE: f32 = 165.0;
/// Optimal Y-A...X acceptor angle (PLIP `HALOGEN_ACC_ANGLE`).
pub const HALOGEN_ACCEPTOR_ANGLE: f32 = 120.0;
/// Tolerance applied to both halogen angles (PLIP `HALOGEN_ANGLE_DEV`).
pub const HALOGEN_ANGLE_DEV: f32 = 30.0;

/// Min leg distance through the bridging water (PLIP `WATER_BRIDGE_MINDIST`).
pub const WATER_BRIDGE_MIN_DIST: f32 = 2.5;
/// Max leg distance through the bridging water (PLIP `WATER_BRIDGE_MAXDIST`).
pub const WATER_BRIDGE_MAX_DIST: f32 = 4.1;
/// Water-bridge angle bounds (PLIP `WATER_BRIDGE_OMEGA_MIN`).
pub const WATER_BRIDGE_OMEGA_MIN: f32 = 71.0;
/// Water-bridge angle bounds (PLIP `WATER_BRIDGE_OMEGA_MAX`).
pub const WATER_BRIDGE_OMEGA_MAX: f32 = 140.0;
/// Water-bridge angle bounds (PLIP `WATER_BRIDGE_THETA_MIN`).
pub const WATER_BRIDGE_THETA_MIN: f32 = 100.0;

/// Max metal-ion to coordinating O/N/S atom distance (PLIP `METAL_DIST_MAX`,
/// Harding 2001).
pub const METAL_DIST_MAX: f32 = 3.0;

/// Max S(gamma)-S(gamma) distance for a disulfide. PLIP does not detect
/// disulfides; this is RING's SSBOND cutoff (nominal S-S bond length ~2.05 A).
pub const DISULFIDE_DIST_MAX: f32 = 2.5;

/// The cell size for the interaction spatial grid; the largest interaction cutoff.
pub const GRID_CELL_SIZE: f32 = PICATION_DIST_MAX;

/// Element symbols treated as coordinating metal ions, from PLIP `METAL_IONS`
/// (reduced to element symbols; matched case-insensitively against the element
/// column, which is why a calcium ion "CA" is distinct from a C-alpha "C").
pub const METAL_ELEMENTS: &[&str] = &[
    "CA", "CO", "MG", "MN", "FE", "CU", "ZN", "LI", "NA", "K", "RB", "SR", "CS", "BA", "CR", "NI",
    "RU", "RH", "PD", "AG", "CD", "LA", "W", "OS", "IR", "PT", "AU", "HG", "CE", "PR", "SM", "EU",
    "GD", "TB", "YB", "LU", "AL", "GA", "IN", "SB", "TL", "PB",
];
