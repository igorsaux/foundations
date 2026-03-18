// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const std = @import("std");

const constants = @import("constants.zig");
const Flavor = @import("flavor.zig").Flavor;
const LiquidPhase = @import("liquid_phase.zig").LiquidPhase;
const MoleculeId = @import("molecule_id.zig").MoleculeId;
const ScoredAttribute = @import("scored_attribute.zig").ScoredAttribute;

pub const SolidPhase = struct {
    /// The molecular species of this solid.
    molecule: MoleculeId,
    /// Amount of substance in moles.
    moles: f64,
    /// Characteristic particle diameter in meters.
    /// Typical range: 1e-6 m (fine powder) to 0.1 m (large crystal).
    /// Used for surface area calculations, dissolution rates, and phase merging.
    particle_diameter: f64,
    /// Solution trapped within the crystal lattice during crystallization.
    /// Contains dissolved impurities that were incorporated during growth.
    occluded_solution: ?LiquidPhase,
    /// Whether this solid has a regular crystal lattice (true) or is
    /// amorphous (false). Crystalline solids have sharper melting points
    /// and more predictable solubility.
    crystalline: bool,

    pub inline fn deinit(this: *SolidPhase, allocator: std.mem.Allocator) void {
        if (this.occluded_solution) |*s| {
            s.deinit(allocator);
            this.occluded_solution = null;
        }
    }

    /// Returns the total volume of this solid phase in liters.
    ///
    /// V = n * MW / ρ [cm³] → / 1000 [L]
    ///
    /// This is the actual material volume, not the bulk volume
    /// (which would include void space between particles).
    pub inline fn getVolume(this: *const SolidPhase) f64 {
        const mw = this.molecule.getWeight();
        const density = this.molecule.getDensity();

        // moles * (g/mol) / (g/cm³) = cm³ = mL, then / 1000 for L
        const volume = this.moles * mw / density / 1000.0;

        if (volume < constants.MIN_VOLUME) {
            return 0.0;
        }

        return volume;
    }

    /// Returns the approximate number of particles in this solid phase.
    ///
    /// Assumes roughly spherical particles:
    /// V_particle = (π/6) * d³
    /// N = V_total / V_particle
    ///
    /// Returns 0 if particle_diameter is zero or negative.
    pub inline fn getParticleCount(this: *const SolidPhase) f64 {
        if (this.particle_diameter <= 0.0) {
            return 0.0;
        }

        // Volume in m³ (getVolume returns liters, 1 L = 0.001 m³)
        const total_volume_m3 = this.getVolume() * 0.001;
        const d = this.particle_diameter;
        const particle_volume = (std.math.pi / 6.0) * d * d * d;

        return total_volume_m3 / particle_volume;
    }

    /// Returns the total surface area of all particles in m².
    ///
    /// For N spherical particles of diameter d:
    /// A_total = N * π * d²
    ///
    /// Surface area affects dissolution rate, heat transfer, and reactivity.
    pub inline fn getSurfaceArea(this: *const SolidPhase) f64 {
        const count = this.getParticleCount();
        const d = this.particle_diameter;

        return count * std.math.pi * d * d;
    }

    /// Returns the specific surface area in m²/g.
    ///
    /// SSA = A_total / mass
    ///
    /// Higher values indicate finer particles with greater reactivity.
    /// Typical values: 0.1 m²/g (coarse) to 100+ m²/g (nanoparticles).
    pub inline fn getSpecificSurfaceArea(this: *const SolidPhase) f64 {
        const mass = this.moles * this.molecule.getWeight(); // grams

        if (mass <= 0.0) {
            return 0.0;
        }

        return this.getSurfaceArea() / mass;
    }

    /// Returns the color of this solid phase.
    ///
    /// Solid color depends on:
    /// 1. Intrinsic color of the molecular species
    /// 2. Amount of substance (larger crystals appear more intensely colored)
    /// 3. Crystallinity (amorphous solids may appear different)
    ///
    /// For solids with occluded solution, the trapped liquid may
    /// contribute to color if the solid itself is transparent.
    pub fn getColors(this: *const SolidPhase) Flavor {
        var flavor = Flavor{};
        const meta = this.molecule.molecule();

        // Primary color from the solid itself
        if (meta.color != .transparent and meta.color_intensity > constants.MIN_COLOR_INTENSITY) {
            flavor.colors[0] = meta.color;

            // Intensity scales with mass^(1/3) for visual size effect
            // Normalize to 1 gram reference
            const mass = this.moles * meta.weight;
            const size_factor = std.math.pow(f64, @max(mass, 0.001), 1.0 / 3.0);

            flavor.color_intensities[0] = @floatCast(@min(@as(f64, meta.color_intensity) * size_factor, 1.0));
        }

        // Secondary color from occluded solution (if solid is transparent)
        if (this.occluded_solution) |*occ| {
            const occ_flavor = occ.getColors();

            if (meta.color == .transparent) {
                // Occluded solution shows through
                flavor.colors[0] = occ_flavor.colors[0];
                flavor.color_intensities[0] = occ_flavor.color_intensities[0] * 0.5; // Dimmed
                flavor.colors[1] = occ_flavor.colors[1];
                flavor.color_intensities[1] = occ_flavor.color_intensities[1] * 0.5;
            } else if (occ_flavor.colors[0] != .transparent) {
                // Solid is colored; occluded adds secondary tint
                flavor.colors[1] = occ_flavor.colors[0];
                flavor.color_intensities[1] = occ_flavor.color_intensities[0] * 0.3;
            }
        }

        return flavor;
    }

    /// Returns the odors emanating from this solid phase.
    ///
    /// Solids smell via:
    /// 1. Sublimation (direct solid -> gas, especially low-melting solids)
    /// 2. Surface desorption of volatile impurities
    /// 3. Occluded solution evaporation
    ///
    /// Sublimation rate approximated from melting point:
    /// lower melting point = higher sublimation = stronger odor.
    pub fn getOdors(this: *const SolidPhase, temperature: f64) Flavor {
        var flavor = Flavor{};
        const meta = this.molecule.molecule();

        // Direct sublimation odor from the solid itself
        if (meta.odor != .none and meta.odor_intensity > constants.MIN_ODOR_INTENSITY) {
            // Sublimation factor: lower melting point = more sublimation
            var sublimation_factor: f64 = 0.1; // Base: most solids barely sublime

            if (meta.melting_point) |mp| {
                // Sublimation enthalpy ≈ ΔH_fus + ΔH_vap, roughly 1.2x ΔH_vap
                const k: f64 = if (meta.delta_hvap) |hvap|
                    hvap * 1.2 * 1000.0 / 8.314
                else
                    6000.0; // Higher barrier for sublimation

                // Reference: at T = 0.7 * T_melt, very little sublimation
                // At T = 0.95 * T_melt, noticeable sublimation
                sublimation_factor = @exp(k * (1.0 / mp - 1.0 / temperature));
                sublimation_factor = @min(@max(sublimation_factor, 0.001), 1.0);
            }

            // Surface area factor: more moles = more surface = more smell
            const surface_factor = std.math.pow(f64, this.moles * 100.0, 2.0 / 3.0);

            const score = @as(f64, meta.odor_intensity) * sublimation_factor *
                @min(surface_factor, 10.0);

            if (score > 0.01) {
                flavor.odors[0] = meta.odor;
                flavor.odor_intensities[0] = @floatCast(@min(score, 1.0));
            }
        }

        // Odor from occluded solution (wet crystals smell)
        if (this.occluded_solution) |*occ| {
            const occ_flavor = occ.getOdors(temperature);

            if (flavor.odors[0] == .none) {
                flavor.odors[0] = occ_flavor.odors[0];
                flavor.odor_intensities[0] = occ_flavor.odor_intensities[0] * 0.7;
                flavor.odors[1] = occ_flavor.odors[1];
                flavor.odor_intensities[1] = occ_flavor.odor_intensities[1] * 0.7;
            } else if (occ_flavor.odors[0] != .none and occ_flavor.odors[0] != flavor.odors[0]) {
                flavor.odors[1] = occ_flavor.odors[0];
                flavor.odor_intensities[1] = occ_flavor.odor_intensities[0] * 0.5;
            }
        }

        return flavor;
    }

    /// Returns combined color and odor profile for this solid phase.
    pub inline fn getFlavor(this: *const SolidPhase, temperature: f64) Flavor {
        var flavor = this.getColors();
        const odor_flavor = this.getOdors(temperature);

        flavor.odors = odor_flavor.odors;
        flavor.odor_intensities = odor_flavor.odor_intensities;

        return flavor;
    }
};
