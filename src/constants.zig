// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

/// Universal gas constant in J/(mol·K)
pub const R: f64 = 8.31446261815324;

/// Standard pressure
pub const P_std: f64 = 101325.0;

/// Avogadro's number, 1/mol
pub const N_A: f64 = 6.02214076e23;

/// Planck constant, J·s
pub const h: f64 = 6.62607015e-34;

/// Boltzmann constant, J/K
pub const k_B: f64 = 1.380649e-23;

/// Water autodissociation constant at 25°C
pub const Kw: f64 = 1e-14;

/// Temperature for standard conditions (K)
pub const T_std: f64 = 298.15;

/// Minimum volume threshold below which a substance is considered absent.
pub const MIN_VOLUME: f64 = 1e-12;

/// Minimum mole threshold below which a substance is considered absent.
///
/// 1e-12 mol ≈ 1 picomole ≈ ~6*10¹¹ molecules.
///
/// Rationale:
/// - Large enough to avoid floating-point noise (f64 epsilon ≈ 1e-16)
/// - Small enough to not lose pharmacologically relevant amounts
///   (typical drug dose ~1e-3 mol, receptor binding ~1e-9 mol)
/// - Matches analytical chemistry detection limits (LC-MS ~picomolar)
pub const MIN_MOLES: f64 = 1e-12;

pub const MIN_AREA: f64 = 1e-15;

pub const MIN_CONCENTRATION: f64 = 1e-15;

/// Returns true if the amount is too small to be chemically meaningful.
pub inline fn isNegligible(moles: f64) bool {
    return moles < MIN_MOLES;
}

/// Minimum odor intensity required for an odor to be considered.
pub const MIN_ODOR_INTENSITY: f64 = 1e-6;

/// Minimum color intensity required for a color to be considered.
pub const MIN_COLOR_INTENSITY: f64 = 1e-6;

pub inline fn kToC(k: f64) f64 {
    return k - 273.15;
}

pub inline fn cToK(c: f64) f64 {
    return c + 273.15;
}

/// Maximum relative difference in particle diameter for solid phase merging.
/// Phases with |d1 - d2| / max(d1, d2) < this threshold are considered
/// similar enough to combine into a single phase.
pub const PARTICLE_SIZE_MERGE_RATIO: f64 = 0.2;

/// Default particle diameter (meters) for newly precipitated solids.
/// 100 μm = 0.1 mm, typical for moderate precipitation rates.
pub const DEFAULT_PRECIPITATION_DIAMETER: f64 = 1e-4;

/// Default particle diameter (meters) for solids formed by freezing.
/// 1 cm, representing bulk solidification rather than fine crystals.
pub const DEFAULT_FREEZING_DIAMETER: f64 = 0.01;
