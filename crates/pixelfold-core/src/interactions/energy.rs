//! Per-interaction energy weights for the residue interaction network.
//!
//! Each interaction kind carries a representative interaction energy, a positive
//! magnitude in kilojoules per mole, for one detected contact of that kind. An
//! edge sums the weights of the contacts folded into it (weight times count), so
//! a residue pair joined by three hydrogen bonds scores three times one.
//!
//! The values are a literature-backed scoring table on one consistent basis: the
//! intrinsic strength of a single contact in a protein / condensed-phase setting,
//! not gas-phase extremes and not the much smaller net free energy of folding
//! (which cancels most of the raw interaction against desolvation). That is the
//! same basis RING's own weights use, and the four kinds RING scores land inside
//! these independently-sourced ranges (RING: H-bond 17, ionic 20, van der Waals
//! 6 per residue pair, pi-stacking 9.4 kJ/mol), so the two tables corroborate.
//!
//! Two values carry caveats worth stating up front. The disulfide is a covalent
//! bond (~251 kJ/mol), an order of magnitude above every non-covalent contact,
//! so it dominates any additive score it joins; the pane compresses the scale
//! perceptually rather than distort the number. The hydrophobic weight is per
//! atom-atom contact, the scale this engine detects them at, which is why it is
//! small and why an edge's many contacts sum to the per-residue-pair magnitude
//! RING reports.

use super::InteractionKind;

impl InteractionKind {
    /// The representative energy of one detected contact of this kind, a positive
    /// magnitude in kJ/mol. See the module documentation for the basis and the
    /// per-value sources.
    ///
    /// The match is exhaustive on purpose: a new interaction kind cannot compile
    /// without choosing a weight.
    pub fn contact_energy_kj(self) -> f32 {
        match self {
            // Intrinsic backbone N-H...O=C bond strength, moderate and neutral:
            // 4.79 kcal/mol (beta-sheet) / 5.57 (alpha-helix), i.e. ~20 kJ/mol.
            // Sheu, Yang, Selzle & Schlag, "Energetics of hydrogen bonds in
            // peptides", PNAS 100(22):12683-12687 (2003). Agrees with RING's 17.
            InteractionKind::HydrogenBond => 20.0,

            // A single well-formed protein salt bridge contributes 3-5 kcal/mol
            // to folding stability; 4 kcal/mol ~ 17 kJ/mol. Anderson, Becktel &
            // Dahlquist, Biochemistry 29(9):2403-2408 (1990). Agrees with RING's
            // 20.
            InteractionKind::SaltBridge => 17.0,

            // One apolar atom-atom van der Waals contact at the Lennard-Jones
            // minimum: OPLS-AA carbon well depth epsilon ~ 0.07 kcal/mol ~ 0.3
            // kJ/mol. Jorgensen, Maxwell & Tirado-Rives, JACS 118:11225 (1996).
            // Per atom-atom, not per residue pair: ~15-20 such contacts sum to
            // RING's per-pair 6 kJ/mol, which is how the two reconcile.
            InteractionKind::Hydrophobic => 0.3,

            // A stacked aromatic pair, between the protein estimate (Burley &
            // Petsko, Science 229:23 (1985), ~4-8 kJ/mol) and the benzene-dimer
            // minima (Sinnokrot, Valeev & Sherrill, JACS 124:10887 (2002),
            // ~11.5 kJ/mol). Agrees with RING's 9.4.
            InteractionKind::PiStacking => 10.0,

            // A cation-pi contact (Lys/Arg over an aromatic face) in an aqueous /
            // protein context, ~3 kcal/mol. Dougherty, Acc. Chem. Res.
            // 46(4):885 (2013); Gallivan & Dougherty, PNAS 96:9459 (1999).
            InteractionKind::PiCation => 12.5,

            // A C-X...O=C halogen bond to a backbone carbonyl, mid-range across
            // Cl/Br/I (per-halogen ~8/10/16 kJ/mol). Auffinger et al., PNAS
            // 101:16789 (2004); Nishikawa et al., J. Chem. Inf. Model. 64(15)
            // (2024).
            InteractionKind::HalogenBond => 12.0,

            // A water-mediated bridge: geometrically two hydrogen bonds, but the
            // net per-bridge stabilisation is about one H-bond once the trapped
            // water's entropy is paid. Ben-Naim, J. Chem. Phys. 90:7412 (1989);
            // Petukhov et al., Protein Sci. 8:1982 (1999). The softest of these
            // constants; the ~3 kcal/mol estimate was later contested.
            InteractionKind::WaterBridge => 12.5,

            // One metal-ligand coordination bond (His to Zn2+, aqueous, screened),
            // a distinctly strong non-covalent contact. Liao et al., BMC Chem.
            // 7:44 (2013); Rulisek & Havlas, JACS 122:10428 (2000). Varies
            // strongly with metal and donor (defensible band ~50-150 kJ/mol).
            InteractionKind::MetalCoordination => 70.0,

            // The covalent S-S bond of a cystine, ~60 kcal/mol ~ 251 kJ/mol.
            // Benson, Chem. Rev. 78:23 (1978); Luo, Comprehensive Handbook of
            // Chemical Bond Energies (2007). Covalent, so ~12x a hydrogen bond:
            // it dominates any additive score and is compressed apart in the pane.
            InteractionKind::Disulfide => 251.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL: [InteractionKind; 9] = [
        InteractionKind::HydrogenBond,
        InteractionKind::SaltBridge,
        InteractionKind::Hydrophobic,
        InteractionKind::PiStacking,
        InteractionKind::PiCation,
        InteractionKind::HalogenBond,
        InteractionKind::WaterBridge,
        InteractionKind::MetalCoordination,
        InteractionKind::Disulfide,
    ];

    #[test]
    fn every_kind_has_a_positive_weight() {
        for kind in ALL {
            assert!(kind.contact_energy_kj() > 0.0, "{kind:?} has no weight");
        }
    }

    #[test]
    fn the_table_matches_the_signed_off_values() {
        assert_eq!(InteractionKind::HydrogenBond.contact_energy_kj(), 20.0);
        assert_eq!(InteractionKind::SaltBridge.contact_energy_kj(), 17.0);
        assert_eq!(InteractionKind::Hydrophobic.contact_energy_kj(), 0.3);
        assert_eq!(InteractionKind::PiStacking.contact_energy_kj(), 10.0);
        assert_eq!(InteractionKind::PiCation.contact_energy_kj(), 12.5);
        assert_eq!(InteractionKind::HalogenBond.contact_energy_kj(), 12.0);
        assert_eq!(InteractionKind::WaterBridge.contact_energy_kj(), 12.5);
        assert_eq!(InteractionKind::MetalCoordination.contact_energy_kj(), 70.0);
        assert_eq!(InteractionKind::Disulfide.contact_energy_kj(), 251.0);
    }

    #[test]
    fn the_covalent_disulfide_dominates_every_non_covalent_contact() {
        // The disulfide must stay far above the strongest non-covalent weight, so
        // the display layer knows it needs perceptual compression.
        let non_covalent_max = ALL
            .iter()
            .filter(|k| **k != InteractionKind::Disulfide)
            .map(|k| k.contact_energy_kj())
            .fold(0.0f32, f32::max);
        assert!(InteractionKind::Disulfide.contact_energy_kj() > 3.0 * non_covalent_max);
    }
}
