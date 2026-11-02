use glam::Vec3;
use pixelfold_core::interactions::InteractionKind;
use pixelfold_core::{Atom, FixedStr, Interaction, Protein, SecondaryStructure, rin};
use ratatui::style::Color;

use super::draw::normalize_positions;
use super::palette::{burial_color, edge_intensity, scale_color};
use super::*;

fn atom(chain: &str, seq: i32, name: &str, residue: &str) -> Atom {
    Atom {
        serial: 0,
        name: FixedStr::new(name),
        element: FixedStr::new(&name[..1]),
        residue_name: FixedStr::new(residue),
        residue_seq: seq,
        chain_id: FixedStr::new(chain),
        is_hetatm: false,
        molecule: pixelfold_core::Molecule::Unknown,
        altloc: None,
        occupancy: 1.0,
        insertion_code: None,
        model_number: 1,
        position: Vec3::ZERO,
        b_factor: 0.0,
        secondary_structure: SecondaryStructure::Coil,
    }
}

fn hbond(a: usize, b: usize) -> Interaction {
    Interaction {
        kind: InteractionKind::HydrogenBond,
        atoms_a: vec![a],
        atoms_b: vec![b],
        distance: 3.0,
        angle: None,
        bridge: None,
    }
}

#[test]
fn normalize_fits_positions_into_the_unit_square() {
    assert!(normalize_positions(vec![]).is_empty());
    assert_eq!(normalize_positions(vec![(5.0, 5.0)]), vec![(0.5, 0.5)]);

    // Aspect-preserving: the larger (y) axis spans the frame, both centred.
    let out = normalize_positions(vec![(2.0, 2.0), (4.0, 6.0)]);
    for &(x, y) in &out {
        assert!(
            (0.0..=1.0).contains(&x) && (0.0..=1.0).contains(&y),
            "{x},{y}"
        );
    }
    assert!((out[0].1 - 0.0).abs() < 1e-6 && (out[1].1 - 1.0).abs() < 1e-6);
}

#[test]
fn the_residue_atom_mapping_is_consistent_both_ways() {
    // Residue A/1 (atoms 0,1) hydrogen-bonds A/2 (atoms 2,3); atom 4 is an
    // isolated residue that never interacts.
    let protein = Protein {
        atoms: vec![
            atom("A", 1, "N", "LYS"),
            atom("A", 1, "NZ", "LYS"),
            atom("A", 2, "O", "GLU"),
            atom("A", 2, "OE1", "GLU"),
            atom("B", 5, "CB", "LEU"),
        ],
        title: String::new(),
        hbonds: Vec::new(),
        assembly: None,
        components: Default::default(),
    };
    let network = rin::build(&protein, &[hbond(0, 2)]);
    let view = NetworkView::from_network(&protein, network);

    // Both atoms of A/1 map to A/1's node, and that node's representative
    // atom belongs to A/1.
    let node = view.node_of_atom(0).expect("atom 0 is in the network");
    assert_eq!(
        view.node_of_atom(1),
        Some(node),
        "both A/1 atoms share a node"
    );
    assert_ne!(view.node_of_atom(2), Some(node), "A/2 is a different node");
    let repr = view.atom_of_node(node);
    assert!(
        repr == Some(0) || repr == Some(1),
        "the node's atom is one of A/1's"
    );

    // The non-interacting residue is not in the network.
    assert_eq!(view.node_of_atom(4), None);

    // One laid-out position per node.
    assert_eq!(view.positions.len(), view.network.nodes.len());
    // And one burial value per node.
    assert_eq!(view.rsa.len(), view.network.nodes.len());
}

#[test]
fn toggle_color_switches_the_node_scheme() {
    let protein = Protein {
        atoms: vec![atom("A", 1, "N", "LYS"), atom("A", 2, "O", "GLU")],
        title: String::new(),
        hbonds: Vec::new(),
        assembly: None,
        components: Default::default(),
    };
    let network = rin::build(&protein, &[hbond(0, 1)]);
    let mut view = NetworkView::from_network(&protein, network);

    assert_eq!(view.color_mode, NodeColorMode::Structure);
    view.toggle_color();
    assert_eq!(view.color_mode, NodeColorMode::Burial);
    view.toggle_color();
    assert_eq!(view.color_mode, NodeColorMode::Structure);
}

#[test]
fn burial_color_separates_buried_from_exposed() {
    // Buried reads warm (red over blue), exposed reads cool (blue over red).
    let Color::Rgb(br, _, bb) = burial_color(Some(0.0)) else {
        panic!("expected an rgb colour");
    };
    assert!(br > bb, "a buried residue should read warm");

    let Color::Rgb(er, _, eb) = burial_color(Some(1.0)) else {
        panic!("expected an rgb colour");
    };
    assert!(eb > er, "an exposed residue should read cool");

    // Above the reference maximum still reads fully exposed, not past it.
    assert_eq!(burial_color(Some(2.0)), burial_color(Some(1.0)));
    // A residue with no reference maximum is neutral.
    assert_eq!(burial_color(None), Color::DarkGray);
}

#[test]
fn edge_intensity_rises_with_energy_but_compresses_the_disulfide() {
    // Disulfide (251) sets the top; a hydrogen bond (20) and hydrophobic
    // contact (0.3) sit below it, brighter with energy but never dark.
    let max = 251.0;
    let weak = edge_intensity(0.3, max);
    let hbond = edge_intensity(20.0, max);
    let disulfide = edge_intensity(251.0, max);

    assert!(weak < hbond && hbond < disulfide, "intensity tracks energy");
    assert!(weak >= 0.4, "the weakest edge still draws above the floor");
    assert!((disulfide - 1.0).abs() < 1e-6, "the strongest edge is full");
    // Log compression: the 12x energy gap is nowhere near a 12x brightness
    // gap, so one covalent bond does not blow out the scale.
    assert!(
        disulfide / hbond < 2.0,
        "disulfide does not dominate visually"
    );
}

#[test]
fn scale_color_dims_toward_black() {
    assert_eq!(
        scale_color(Color::Rgb(200, 100, 50), 0.5),
        Color::Rgb(100, 50, 25)
    );
    assert_eq!(
        scale_color(Color::Rgb(200, 100, 50), 1.0),
        Color::Rgb(200, 100, 50)
    );
    // A non-rgb colour is left untouched.
    assert_eq!(scale_color(Color::Gray, 0.5), Color::Gray);
}

#[test]
fn edge_passes_honours_kind_and_chain_filters() {
    let intra_hbond = EdgeVis {
        u: 0,
        v: 1,
        kind: InteractionKind::HydrogenBond,
        energy: 20.0,
        interchain: false,
        frequency: 1.0,
    };
    let inter_salt = EdgeVis {
        u: 0,
        v: 2,
        kind: InteractionKind::SaltBridge,
        energy: 17.0,
        interchain: true,
        frequency: 1.0,
    };

    // No filter: everything passes.
    assert!(edge_passes(&intra_hbond, None, ChainFilter::All, 0.0));
    assert!(edge_passes(&inter_salt, None, ChainFilter::All, 0.0));

    // Kind filter keeps only the chosen kind.
    let only_hbond = Some(InteractionKind::HydrogenBond);
    assert!(edge_passes(&intra_hbond, only_hbond, ChainFilter::All, 0.0));
    assert!(!edge_passes(&inter_salt, only_hbond, ChainFilter::All, 0.0));

    // Chain filter splits intra from inter.
    assert!(edge_passes(&intra_hbond, None, ChainFilter::IntraOnly, 0.0));
    assert!(!edge_passes(
        &intra_hbond,
        None,
        ChainFilter::InterOnly,
        0.0
    ));
    assert!(edge_passes(&inter_salt, None, ChainFilter::InterOnly, 0.0));
    assert!(!edge_passes(&inter_salt, None, ChainFilter::IntraOnly, 0.0));
}

#[test]
fn cycle_chain_filter_rotates_all_intra_inter() {
    let protein = Protein {
        atoms: vec![atom("A", 1, "N", "LYS"), atom("A", 2, "O", "GLU")],
        title: String::new(),
        hbonds: Vec::new(),
        assembly: None,
        components: Default::default(),
    };
    let network = rin::build(&protein, &[hbond(0, 1)]);
    let mut view = NetworkView::from_network(&protein, network);

    assert_eq!(view.chain_filter, ChainFilter::All);
    view.cycle_chain_filter();
    assert_eq!(view.chain_filter, ChainFilter::IntraOnly);
    view.cycle_chain_filter();
    assert_eq!(view.chain_filter, ChainFilter::InterOnly);
    view.cycle_chain_filter();
    assert_eq!(view.chain_filter, ChainFilter::All);
}

#[test]
fn cycle_kind_filter_walks_present_kinds_then_clears() {
    // One residue pair joined by both a hydrogen bond and a salt bridge, so
    // two kinds are present.
    let protein = Protein {
        atoms: vec![
            atom("A", 1, "N", "LYS"),
            atom("A", 1, "NZ", "LYS"),
            atom("A", 2, "O", "GLU"),
            atom("A", 2, "OE1", "GLU"),
        ],
        title: String::new(),
        hbonds: Vec::new(),
        assembly: None,
        components: Default::default(),
    };
    let salt = Interaction {
        kind: InteractionKind::SaltBridge,
        atoms_a: vec![1],
        atoms_b: vec![3],
        distance: 3.5,
        angle: None,
        bridge: None,
    };
    let network = rin::build(&protein, &[hbond(0, 2), salt]);
    let mut view = NetworkView::from_network(&protein, network);

    assert_eq!(
        view.present_kinds,
        vec![InteractionKind::HydrogenBond, InteractionKind::SaltBridge]
    );
    assert_eq!(view.kind_filter, None);
    view.cycle_kind_filter();
    assert_eq!(view.kind_filter, Some(InteractionKind::HydrogenBond));
    view.cycle_kind_filter();
    assert_eq!(view.kind_filter, Some(InteractionKind::SaltBridge));
    view.cycle_kind_filter();
    assert_eq!(
        view.kind_filter, None,
        "past the last kind clears the filter"
    );
}

/// Three residues, A/1 hydrogen-bonded to A/2 and in hydrophobic contact
/// with B/5, at `x` along the x axis.
fn residues(x: f32) -> Protein {
    let at = |mut atom: Atom| {
        atom.position.x += x;
        atom
    };
    Protein {
        atoms: vec![
            at(atom("A", 1, "N", "LYS")),
            at(atom("A", 2, "O", "GLU")),
            at(atom("B", 5, "CB", "LEU")),
        ],
        title: String::new(),
        hbonds: Vec::new(),
        assembly: None,
        components: Default::default(),
    }
}

fn contact(a: usize, b: usize) -> Interaction {
    Interaction {
        kind: InteractionKind::Hydrophobic,
        ..hbond(a, b)
    }
}

/// Four models: the hydrogen bond in every one, the contact in one, and in
/// the last A/1 and A/2 in hydrophobic contact with no hydrogen bond.
fn ensemble() -> rin::Ensemble {
    let mut ensemble = rin::Ensemble::default();
    ensemble.add(&residues(0.0), &[hbond(0, 1), contact(0, 2)]);
    for _ in 0..2 {
        ensemble.add(&residues(0.0), &[hbond(0, 1)]);
    }
    ensemble.add(&residues(0.0), &[contact(0, 1)]);
    ensemble
}

#[test]
fn an_ensemble_brightens_frequent_edges_and_a_floor_hides_rare_ones() {
    let mut view = NetworkView::from_ensemble(&residues(0.0), &ensemble());
    let edge = |kind, v| {
        view.edges
            .iter()
            .find(|e| e.kind == kind && e.v == v)
            .expect("an edge of the kind")
    };
    let (hbond, contact) = (
        edge(InteractionKind::HydrogenBond, 1),
        edge(InteractionKind::Hydrophobic, 2),
    );
    assert_eq!((hbond.frequency, contact.frequency), (0.75, 0.25));
    assert!(view.brightness(hbond, 0.0) > view.brightness(contact, 0.0));
    let shown = |view: &NetworkView| {
        view.edges
            .iter()
            .filter(|e| edge_passes(e, None, ChainFilter::All, view.min_frequency))
            .count()
    };
    assert_eq!(shown(&view), 3);

    view.set_min_frequency(0.5);
    assert_eq!(shown(&view), 1);
}

#[test]
fn the_scope_steps_through_the_floors_and_back_to_the_model() {
    let mut scope = NetworkScope::default();
    let mut floors = Vec::new();
    for _ in 0..4 {
        scope = scope.next();
        floors.push(scope.is_ensemble().then(|| scope.min_frequency()));
    }

    assert_eq!(floors, [Some(0.0), Some(0.5), Some(0.9), None]);
}

// A/1 and A/2 hydrogen-bond in three models and touch only hydrophobically
// in the fourth: they interact in every model, though no edge says so.
#[test]
fn a_residue_s_partners_count_the_models_they_touch_it_in_by_any_kind() {
    let view = NetworkView::from_ensemble(&residues(0.0), &ensemble());
    let a1 = view.node_of_atom(0).expect("A/1 is in the network");

    assert_eq!(view.partners(a1), [("A/2", 4), ("B/5", 1)]);
}

#[test]
fn a_network_laid_over_another_model_keeps_its_layout_and_settings() {
    let mut view = NetworkView::from_ensemble(&residues(0.0), &ensemble());
    view.toggle_mode();
    view.set_min_frequency(0.5);
    let positions = view.positions.clone();

    let moved = view.remapped(&residues(10.0));

    assert_eq!(moved.positions, positions);
    assert_eq!(moved.mode, LayoutMode::Spatial);
    assert_eq!(moved.min_frequency, 0.5);
    assert_eq!(moved.partners(0), [("A/2", 4), ("B/5", 1)]);
    assert_eq!(moved.node_pos3[0].x, 10.0, "the new model's coordinates");
    assert_eq!(moved.node_of_atom(2), Some(2));
}

#[test]
fn a_residue_the_shown_model_lacks_maps_to_no_atom() {
    let view = NetworkView::from_ensemble(&residues(0.0), &ensemble());
    let without_b5 = Protein {
        atoms: residues(0.0).atoms[..2].to_vec(),
        ..residues(0.0)
    };

    let moved = view.remapped(&without_b5);

    assert_eq!(moved.atom_of_node(2), None);
    assert_eq!(moved.rsa[2], None);
    assert_eq!(moved.atom_of_node(1), Some(1));
}
