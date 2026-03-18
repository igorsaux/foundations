// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const std = @import("std");

const Color = @import("color.zig").Color;
const constants = @import("constants.zig");
const Flavor = @import("flavor.zig").Flavor;
const MoleculeId = @import("molecule_id.zig").MoleculeId;
const Odor = @import("odor.zig").Odor;
const ScoredAttribute = @import("scored_attribute.zig").ScoredAttribute;

pub const LiquidPhase = struct {
    /// Array of mole amounts indexed by MoleculeId enum ordinal.
    solutes: []f64,

    pub inline fn init(allocator: std.mem.Allocator) error{OutOfMemory}!LiquidPhase {
        const solutes = try allocator.alloc(f64, MoleculeId.count());
        @memset(solutes, 0.0);

        return .{
            .solutes = solutes,
        };
    }

    pub inline fn deinit(this: *LiquidPhase, allocator: std.mem.Allocator) void {
        allocator.free(this.solutes);
        this.solutes = &.{};
    }

    pub inline fn addMolecule(this: *LiquidPhase, molecule: MoleculeId, moles: f64) void {
        this.solutes[@intFromEnum(molecule)] += moles;
    }

    /// Add moles to the phase (element-wise accumulation).
    pub inline fn addMoles(this: *LiquidPhase, moles: []const f64) void {
        std.debug.assert(this.solutes.len == moles.len);

        for (0..this.solutes.len) |i| {
            this.solutes[i] += moles[i];
        }
    }

    pub inline fn getMoles(this: *const LiquidPhase) f64 {
        var total: f64 = 0.0;

        for (this.solutes) |m| {
            total += m;
        }

        return total;
    }

    /// Computes total liquid volume in liters.
    ///
    /// For each component:
    /// mass_i = moles_i * MW_i [g]
    /// volume_i = mass_i / density_i [cm³ = mL]
    /// volume_i_L = volume_i / 1000 [L]
    pub inline fn getVolume(this: *const LiquidPhase) f64 {
        var total: f64 = 0.0;

        for (0..this.solutes.len) |i| {
            const mid: MoleculeId = @enumFromInt(i);
            const moles = this.solutes[i];

            // For gases and solids use partial molar volume
            if (mid.getPartialMolarVolume()) |pmv| {
                total += moles * pmv / 1000.0; // cm³ -> L
            } else {
                // V_i = n_i * MW_i / ρ_i [g / (g/cm³) = cm³ = mL]
                // Convert mL to L by dividing by 1000.
                total += moles * mid.getWeight() / mid.getDensity() / 1000.0;
            }
        }

        if (total < constants.MIN_VOLUME) {
            return 0.0;
        }

        return total;
    }

    /// Returns molar concentration (molarity, mol/L) of the given molecule
    /// in this liquid phase.
    ///  C = n / V_total
    pub inline fn getConcentration(this: *const LiquidPhase, molecule: MoleculeId) f64 {
        const volume = this.getVolume();

        if (volume == 0.0) {
            return 0.0;
        }

        return this.solutes[@intFromEnum(molecule)] / volume;
    }

    /// Returns the mole fraction of the given molecule:
    ///  x_i = n_i / Σ n_j
    /// Returns 0 if the phase is empty.
    pub inline fn getMoleFraction(this: *const LiquidPhase, molecule: MoleculeId) f64 {
        const total = this.getMoles();

        if (total == 0.0) {
            return 0.0;
        }

        return this.solutes[@intFromEnum(molecule)] / total;
    }

    /// Returns the average mass density of the liquid phase in g/cm³.
    ///
    /// Computed as total mass / total volume, assuming ideal mixing:
    /// ρ = (Σ n_i * MW_i) / (Σ n_i * MW_i / ρ_i)
    ///
    /// This is the volume-weighted harmonic mean of component densities.
    pub fn getDensity(this: *const LiquidPhase) f64 {
        var mass: f64 = 0.0;
        var volume_cm3: f64 = 0.0;

        for (0..this.solutes.len) |i| {
            const molecule: MoleculeId = @enumFromInt(i);
            const moles = this.solutes[i];
            const m = moles * molecule.getWeight();

            mass += m;
            volume_cm3 += m / molecule.getDensity();
        }

        if (volume_cm3 == 0.0) {
            return 0.0;
        }

        // density in g/cm³
        return mass / volume_cm3;
    }

    /// Returns total mass of all components in grams.
    ///  m = Σ n_i * MW_i
    pub inline fn getMass(this: *const LiquidPhase) f64 {
        var total: f64 = 0.0;

        for (0..this.solutes.len) |i| {
            const molecule: MoleculeId = @enumFromInt(i);

            total += this.solutes[i] * molecule.getWeight();
        }

        return total;
    }

    /// Returns the molecule that occupies the largest volume in the phase.
    pub inline fn getPrimaryComponent(this: *const LiquidPhase) ?MoleculeId {
        var max_volume: f64 = 0.0;
        var primary: ?MoleculeId = null;

        for (0..this.solutes.len) |i| {
            const molecule: MoleculeId = @enumFromInt(i);
            const volume = this.solutes[i] * molecule.getWeight() / molecule.getDensity() / 1000.0;

            if (volume > max_volume) {
                primary = molecule;
                max_volume = volume;
            }
        }

        return primary;
    }

    /// Returns the volume-weighted average Hildebrand solubility parameter
    /// of the phase in MPa^½.
    pub inline fn getHildebrand(this: *const LiquidPhase) f64 {
        var total_volume: f64 = 0.0;
        var weighted_delta: f64 = 0.0;

        for (0..this.solutes.len) |i| {
            const molecule: MoleculeId = @enumFromInt(i);
            const moles = this.solutes[i];

            if (constants.isNegligible(moles)) {
                continue;
            }

            const molar_volume = molecule.getWeight() / molecule.getDensity();
            const volume = moles * molar_volume;

            total_volume += volume;
            weighted_delta += volume * molecule.getHildebrand();
        }

        if (total_volume == 0.0) {
            return 0.0;
        }

        return weighted_delta / total_volume;
    }

    /// Estimates the saturation concentration (mol/L) of a molecule in this
    /// liquid phase - the maximum dissolved concentration before precipitation.
    ///
    /// Uses aqueous solubility as a reference and adjusts via regular solution
    /// theory (Scatchard-Hildebrand):
    ///
    /// ln(γ_i) = V_m_i * (δ_i - δ_solvent)² / RT
    ///
    /// Since solubility ∝ 1/γ:
    ///
    /// S_phase = S_water * exp[V_m * ((δ_mol - δ_water)² - (δ_mol - δ_phase)²) / RT]
    ///
    /// Returns +inf if water_solubility is null (fully miscible, no solid form).
    pub inline fn getSaturationLimit(this: *const LiquidPhase, molecule: MoleculeId, temperature: f64) f64 {
        const ws = molecule.getWaterSolubility() orelse {
            return std.math.inf(f64);
        };

        const delta_phase: f64 = this.getHildebrand();
        const delta_mol: f64 = molecule.getHildebrand();
        const delta_water: f64 = 47.8;
        const vm: f64 = molecule.getWeight() / molecule.getDensity();
        const RT: f64 = 8.314 * temperature;

        // Regular solution theory: activity coefficient ratio.
        // V_m [cm³/mol] * δ² [J/cm³] / RT [J/mol] -> dimensionless.
        const water_dd = delta_mol - delta_water;
        const phase_dd = delta_mol - delta_phase;
        const exponent = vm * (water_dd * water_dd - phase_dd * phase_dd) / RT;

        return ws * @exp(exponent);
    }

    /// Estimates the boiling point of the liquid phase at a given pressure
    /// in Pascals.
    ///
    /// Strategy:
    /// 1. Find the component with the lowest normal boiling point (at 1 atm)
    ///    among components present above trace amounts (mole fraction > 0.001).
    ///
    /// 2. Apply the Clausius-Clapeyron equation to correct for pressure:
    ///    1/T₂ = 1/T₁ - R * ln(P₂/P₁) / ΔH_vap
    ///
    /// 3. Apply Raoult's law correction for mixtures:
    ///    The effective pressure is P / x_i, meaning the pure component must
    ///    reach a higher vapor pressure to boil in solution.
    ///
    /// Limitations:
    /// - Assumes ideal solution behavior (Raoult's law)
    /// - Assumes ΔH_vap is constant with temperature
    /// - Ignores azeotrope formation
    pub inline fn getBoilingPoint(this: *const LiquidPhase, pressure: f64) f64 {
        var lowest_bp = std.math.inf(f64);
        var best_molecule: ?MoleculeId = null;

        // Step 1: Find the component with the lowest normal boiling point
        // that is present in non-trace amounts.
        for (0..this.solutes.len) |i| {
            const molecule: MoleculeId = @enumFromInt(i);

            if (this.getMoleFraction(molecule) < 0.001) {
                continue;
            }

            if (molecule.getBoilingPoint()) |bp| {
                if (bp < lowest_bp) {
                    lowest_bp = bp;
                    best_molecule = molecule;
                }
            }
        }

        // If no component has a known boiling point, return +inf
        const mol = best_molecule orelse {
            return std.math.inf(f64);
        };

        // Step 2: Apply Clausius-Clapeyron correction for pressure
        if (mol.getDeltaHvap()) |delta_hvap_kj| {
            // Convert kJ/mol to J/mol
            const delta_hvap: f64 = delta_hvap_kj * 1000.0;

            // Step 3: Apply Raoult's law correction.
            // In an ideal mixture, the vapor pressure of component i is:
            // P_i = x_i * P_i*
            // The solution boils when the total vapor pressure equals external pressure.
            // For the dominant volatile component, the effective pressure it must
            // overcome is P_actual / x_i (it needs higher pure-component vapor pressure).
            const mole_frac = this.getMoleFraction(mol);
            const effective_pressure = if (mole_frac > 0.001)
                pressure / mole_frac
            else
                pressure;

            // Clausius-Clapeyron: 1/T₂ = 1/T₁ - (R/ΔH_vap) * ln(P₂/P₁)
            const ln_ratio = @log(effective_pressure / constants.P_std);
            const inv_T2 = (1.0 / lowest_bp) - (constants.R * ln_ratio / delta_hvap);

            if (inv_T2 > 0.0) {
                return 1.0 / inv_T2;
            }
        }

        // Fallback: no ΔH_vap data, return the normal boiling point uncorrected.
        return lowest_bp;
    }

    /// Estimates the freezing point of the liquid phase.
    ///
    /// Strategy:
    /// 1. Find the component with the highest normal melting point (the solvent).
    ///
    /// 2. Apply the cryoscopic equation (freezing point depression):
    ///   ΔT_f = K_f * m * i
    ///
    ///   Where:
    ///   - K_f = R * T_f² * M_solvent / ΔH_fus (cryoscopic constant)
    ///   - m = molality of solutes (mol solute / kg solvent)
    ///   - i = van't Hoff factor (1 for non-electrolytes)
    ///
    /// Simplified approach using Blagden's law:
    /// ΔT_f ≈ K_f * Σ(n_solute) / m_solvent
    ///
    /// For ideal solutions without ΔH_fus data, uses Raoult's law approximation:
    /// T_f = T_f° * x_solvent
    /// (valid for dilute solutions)
    ///
    /// Limitations:
    /// - Assumes ideal solution behavior
    /// - Ignores eutectic formation
    /// - Uses simplified cryoscopic model
    pub inline fn getFreezingPoint(this: *const LiquidPhase) f64 {
        var highest_mp = -std.math.inf(f64);
        var best_molecule: ?MoleculeId = null;
        var best_mole_fraction: f64 = 0.0;

        // Step 1: Find the solvent (component with highest melting point
        // that is present in significant amounts).
        for (0..this.solutes.len) |i| {
            const molecule: MoleculeId = @enumFromInt(i);
            const mole_fraction = this.getMoleFraction(molecule);

            if (mole_fraction < 0.001) {
                continue;
            }

            if (molecule.getMeltingPoint()) |mp| {
                if (mp > highest_mp) {
                    highest_mp = mp;
                    best_molecule = molecule;
                    best_mole_fraction = mole_fraction;
                }
            }
        }

        // If no component has a known melting point, return -inf (never freezes).
        const solvent = best_molecule orelse {
            return -std.math.inf(f64);
        };

        // Step 2: Calculate freezing point depression.
        //
        // Full cryoscopic equation:
        // ΔT_f = (R * T_f² * M_solvent / ΔH_fus) * m
        //
        // Raoult's law approximation (for ideal solutions):
        // ln(x_solvent) ≈ -ΔH_fus / R * (1/T_f - 1/T_f°)
        //
        // For simplicity, using the approximation:
        // T_f_mixture ≈ T_f° + (R * T_f°² / ΔH_fus) * ln(x_solvent)

        if (solvent.getDeltaHfus()) |delta_hfus_kj| {
            // Convert kJ/mol to J/mol
            const delta_hfus: f64 = delta_hfus_kj * 1000.0;

            // Calculate mole fraction of solvent
            const x_solvent: f64 = best_mole_fraction;

            if (x_solvent > 0.0 and x_solvent < 1.0) {
                // Using thermodynamically rigorous formula:
                // ln(x_solvent) = -ΔH_fus/R * (1/T_f - 1/T_f°)
                // Solving for T_f:
                // 1/T_f = 1/T_f° - R * ln(x_solvent) / ΔH_fus

                const ln_x = @log(x_solvent);
                const inv_Tf_new = (1.0 / highest_mp) - (constants.R * ln_x / delta_hfus);

                if (inv_Tf_new > 0.0) {
                    return 1.0 / inv_Tf_new;
                }
            }
        }

        // Fallback: Simple Raoult's law approximation
        // T_f_mixture ≈ T_f° * x_solvent^0.5 (empirical adjustment for better fit)
        // Or just return pure solvent melting point if it dominates
        if (best_mole_fraction > 0.9) {
            return highest_mp;
        }

        // For significant mixtures without ΔH_fus data, estimate depression
        // using simplified linear approximation.
        // Typical K_f values are 1.5-5 K·kg/mol, using 3 as average.
        const total_moles = this.getMoles();
        const solvent_moles = this.solutes[@intFromEnum(solvent)];
        const solute_moles = total_moles - solvent_moles;

        if (solvent_moles > 0.0) {
            const solvent_mass_kg = solvent_moles * solvent.getWeight() / 1000.0;
            const molality = solute_moles / solvent_mass_kg;

            // Estimate K_f ≈ R * T_f² * M / ΔH_fus
            // Assume ΔH_fus ≈ 10 kJ/mol (typical for organic compounds)
            const assumed_delta_hfus: f64 = 10000.0; // J/mol
            const M_solvent: f64 = solvent.getWeight() / 1000.0; // kg/mol
            const K_f = constants.R * highest_mp * highest_mp * M_solvent / assumed_delta_hfus;

            const delta_T = K_f * molality;

            return highest_mp - delta_T;
        }

        return highest_mp;
    }

    /// Returns the dominant colors of this liquid phase.
    ///
    /// Color perception depends on:
    /// 1. Molar concentration of each colored species
    /// 2. Color intensity coefficient (molar absorptivity analog)
    /// 3. Path length approximation (phase volume^(1/3))
    ///
    /// Uses Beer-Lambert-like scoring: score = C * ε * l
    /// where C = concentration, ε = color_intensity, l = path length.
    ///
    /// Transparent molecules are ignored. Returns up to 2 colors
    /// sorted by visual dominance.
    pub fn getColors(this: *const LiquidPhase) Flavor {
        var flavor = Flavor{};
        const volume = this.getVolume();

        if (volume < constants.MIN_VOLUME) {
            return flavor;
        }

        // Path length approximation: cube root of volume in cm
        // (1 L = 1000 cm³, so 1 L -> ~10 cm path)
        const path_length: f64 = std.math.pow(f64, volume * 1000.0, 1.0 / 3.0);

        // Accumulate scores for each color type
        var color_scores: [@typeInfo(Color).@"enum".fields.len]f64 =
            .{0.0} ** @typeInfo(Color).@"enum".fields.len;

        for (0..this.solutes.len) |i| {
            const moles = this.solutes[i];

            if (constants.isNegligible(moles)) {
                continue;
            }

            const molecule: MoleculeId = @enumFromInt(i);
            const meta = molecule.molecule();

            // Skip transparent molecules
            if (meta.color == .transparent) {
                continue;
            }

            if (meta.color_intensity < constants.MIN_COLOR_INTENSITY) {
                continue;
            }

            // Beer-Lambert-like score: concentration * absorptivity * path
            const concentration = moles / volume;
            const score = concentration * @as(f64, meta.color_intensity) * path_length;

            color_scores[@intFromEnum(meta.color)] += score;
        }

        // Find top 2 colors (excluding transparent at index 0)
        var top: [2]ScoredAttribute = .{
            .{ .index = 0, .score = 0.0 },
            .{ .index = 0, .score = 0.0 },
        };

        for (1..color_scores.len) |i| {
            const score = color_scores[i];

            if (score > top[0].score) {
                top[1] = top[0];
                top[0] = .{ .index = @intCast(i), .score = score };
            } else if (score > top[1].score) {
                top[1] = .{ .index = @intCast(i), .score = score };
            }
        }

        // Normalize intensities relative to strongest
        const max_score = @max(top[0].score, 1e-9);

        if (top[0].score > 1e-9) {
            flavor.colors[0] = @enumFromInt(top[0].index);
            flavor.color_intensities[0] = @floatCast(@min(top[0].score / max_score, 1.0));
        }

        if (top[1].score > 1e-9) {
            flavor.colors[1] = @enumFromInt(top[1].index);
            flavor.color_intensities[1] = @floatCast(@min(top[1].score / max_score, 1.0));
        }

        return flavor;
    }

    /// Returns the dominant odors of this liquid phase.
    ///
    /// Odor perception from liquids depends on:
    /// 1. Vapor pressure (volatile molecules smell more)
    /// 2. Mole fraction in solution (Raoult's law)
    /// 3. Odor intensity coefficient (odor threshold analog)
    /// 4. Surface area approximation (volume^(2/3))
    ///
    /// Lower boiling point = higher vapor pressure = stronger smell.
    /// Score = x_i * P_vap_rel * odor_intensity * surface_area
    pub fn getOdors(this: *const LiquidPhase, temperature: f64) Flavor {
        var flavor = Flavor{};
        const volume = this.getVolume();

        if (volume < constants.MIN_VOLUME) {
            return flavor;
        }

        // Total moles for mole fraction calculation
        var total_moles: f64 = 0.0;
        for (this.solutes) |moles| {
            total_moles += moles;
        }

        if (constants.isNegligible(total_moles)) {
            return flavor;
        }

        // Surface area approximation: volume^(2/3) in cm²
        const surface_area: f64 = std.math.pow(f64, volume * 1000.0, 2.0 / 3.0);

        // Accumulate scores for each odor type
        var odor_scores: [@typeInfo(Odor).@"enum".fields.len]f64 =
            .{0.0} ** @typeInfo(Odor).@"enum".fields.len;

        for (0..this.solutes.len) |i| {
            const moles = this.solutes[i];
            if (moles < 1e-9) continue;

            const molecule: MoleculeId = @enumFromInt(i);
            const meta = molecule.molecule();

            // Skip odorless molecules
            if (meta.odor == .none) continue;
            if (meta.odor_intensity < 1e-6) continue;

            // Mole fraction (Raoult's law contribution)
            const mole_fraction = moles / total_moles;

            // Relative vapor pressure from boiling point using Clausius-Clapeyron.
            // Lower BP = higher vapor pressure at room temperature.
            //
            // P_vap(T) ∝ exp(-ΔH_vap / R * (1/T - 1/T_bp))
            //
            // Simplified: P_rel = exp(k * (1/T_bp - 1/T))
            // where k ≈ ΔH_vap / R ≈ 4000-5000 K for typical organics.
            //
            // At T = 298 K:
            //   - BP = 298 K → P_rel = 1.0 (boiling, max vapor)
            //   - BP = 373 K → P_rel ≈ 0.07 (water-like, low vapor)
            //   - BP = 313 K → P_rel ≈ 0.6 (DCM-like, high vapor)
            var vapor_pressure_rel: f64 = 1.0;

            if (meta.boiling_point) |bp| {
                // k ≈ 4500 K (average ΔH_vap ≈ 37 kJ/mol)
                const k: f64 = 4500.0;
                // Note: (1/T_bp - 1/T) is negative when T < T_bp,
                // giving P_rel < 1, which is correct.
                vapor_pressure_rel = @exp(k * (1.0 / bp - 1.0 / temperature));
                // Clamp to reasonable range [0.001, 10]
                vapor_pressure_rel = @min(@max(vapor_pressure_rel, 0.001), 10.0);
            }

            const score = mole_fraction * vapor_pressure_rel *
                @as(f64, meta.odor_intensity) * surface_area;

            odor_scores[@intFromEnum(meta.odor)] += score;
        }

        // Find top 2 odors (excluding .none at index 0)
        var top: [2]ScoredAttribute = .{
            .{ .index = 0, .score = 0.0 },
            .{ .index = 0, .score = 0.0 },
        };

        for (1..odor_scores.len) |i| {
            const score = odor_scores[i];

            if (score > top[0].score) {
                top[1] = top[0];
                top[0] = .{ .index = @intCast(i), .score = score };
            } else if (score > top[1].score) {
                top[1] = .{ .index = @intCast(i), .score = score };
            }
        }

        // Normalize intensities
        const max_score = @max(top[0].score, 1e-9);

        if (top[0].score > 1e-9) {
            flavor.odors[0] = @enumFromInt(top[0].index);
            flavor.odor_intensities[0] = @floatCast(@min(top[0].score / max_score, 1.0));
        }
        if (top[1].score > 1e-9) {
            flavor.odors[1] = @enumFromInt(top[1].index);
            flavor.odor_intensities[1] = @floatCast(@min(top[1].score / max_score, 1.0));
        }

        return flavor;
    }

    /// Returns combined color and odor profile for this liquid phase.
    pub inline fn getFlavor(this: *const LiquidPhase, temperature: f64) Flavor {
        var flavor = this.getColors();
        const odor_flavor = this.getOdors(temperature);

        flavor.odors = odor_flavor.odors;
        flavor.odor_intensities = odor_flavor.odor_intensities;

        return flavor;
    }
};
