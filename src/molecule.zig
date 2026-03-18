// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const std = @import("std");

const Color = @import("color.zig").Color;
const constants = @import("constants.zig");
const Odor = @import("odor.zig").Odor;

pub const Molecule = struct {
    pub const State = enum(u8) {
        liquid,
        gas,
        solid,
    };

    /// Molecular weight in g/mol (Daltons).
    weight: f64,
    /// Mass density in g/cm³ at 25°C and 1 atm.
    density: f64,
    /// Acid dissociation constants (pKa values), sorted in ascending order.
    /// Used for Henderson-Hasselbalch charge calculations.
    pka: []const f64 = &.{},
    /// Net charge of the fully protonated form of the molecule.
    /// For each pKa, one proton can be lost, decreasing charge by 1.
    fully_protonated_charge: i8 = 0,
    /// Octanol-water partition coefficient (log₁₀ P).
    /// Measures hydrophobicity: positive = lipophilic, negative = hydrophilic.
    log_p: f64 = -7.0,
    /// Number of hydrogen bond donors (N-H and O-H bonds).
    hbd: u8 = 0,
    /// Number of hydrogen bond acceptors (N and O atoms).
    hba: u8 = 0,
    /// Topological polar surface area in Ų.
    /// If 0, estimated from hbd/hba using Ertl's approximation.
    psa: f64 = 0,
    /// Physical state at 25°C and 1 atm.
    state: State = .liquid,
    /// Normal boiling point in Kelvin at 1 atm.
    /// Used for Clausius-Clapeyron pressure corrections.
    boiling_point: ?f64 = null,
    /// Melting point in Kelvin.
    melting_point: ?f64 = null,
    /// Aqueous solubility in mol/L (M) at 25°C.
    ///
    /// water_solubility = null means "no saturation limit"
    ///  -> The molecule is either:
    /// (a) a liquid miscible with water in all proportions, or
    /// (b) a gas (solubility handled differently), or
    /// (c) data unknown (safe default: won't precipitate)
    ///
    /// water_solubility = X means "precipitates above X mol/L"
    ///  -> The molecule has a SOLID FORM at 25°C that can crash out
    /// RULE OF THUMB:
    ///  state == .solid  -> MUST set water_solubility (it will precipitate)
    ///  state == .liquid -> usually null (miscible or forms separate layer)
    ///  state == .gas    -> optionally set for dissolved gas limit (Henry's law)
    water_solubility: ?f64 = null,
    /// Standard enthalpy of formation in kJ/mol at 25°C.
    delta_hf: f64 = 0.0,
    /// Hildebrand solubility parameter δ in MPa^½ (≡ (J/cm³)^½).
    /// Measures cohesive energy density; "like dissolves like" = similar δ.
    /// Typical range: ~12 (fluorocarbons) – ~48 (water).
    /// If null, estimated from delta_hvap, density, and weight.
    hildebrand: ?f64 = null,
    /// Plasma protein binding fraction (0.0 = unbound, 1.0 = fully bound).
    /// Affects free drug concentration in pharmacokinetic models.
    protein_binding: f64 = 0.0,
    /// Gas-lipid partition coefficient (dimensionless).
    /// Ratio of equilibrium concentration in lipid to gas phase.
    /// Relevant for volatile anesthetics and inhaled substances.
    gas_lipid_partition: f64 = 1.0,
    /// Enthalpy of vaporization in kJ/mol at the normal boiling point.
    /// Required for Clausius-Clapeyron boiling point correction with pressure.
    delta_hvap: ?f64 = null,
    /// Enthalpy of fusion (melting) in kJ/mol at the melting point.
    /// Required for freezing point depression calculations in mixtures.
    /// Used in the cryoscopic equation: ΔT_f = K_f * m
    /// where K_f = R * T_f² * M / ΔH_fus
    delta_hfus: ?f64 = null,
    /// Molar heat capacity at constant pressure (C_p) in J/(mol·K).
    /// Determines how much thermal energy is required to raise the temperature.
    /// Default is approx 50 J/(mol·K) (typical for organic substances).
    heat_capacity: f64 = 50.0,
    /// Color of the molecule.
    color: Color = .transparent,
    /// Color intensity coefficient (0.0 = transparent, 1.0 = deeply colored).
    /// Reflects molar absorptivity (extinction coefficient).
    /// Dyes ~1.0, water ~0.0, dilute solutions ~0.1-0.3.
    color_intensity: f32 = 0.0,
    /// Odor of the molecule.
    odor: Odor = .none,
    /// Odor intensity coefficient (0.0 = odorless, 1.0 = very strong).
    /// Reflects vapor pressure, odor threshold, and perception strength.
    /// Thiols ~1.0, water ~0.0, alcohols ~0.3-0.5.
    odor_intensity: f32 = 0.0,
    /// If the molecule is an ion.
    is_ion: bool = false,
    /// Bulk Modulus (Volume elasticity) in Pascals (Pa).
    /// Resistance to compression.
    /// Used to calculate pressure spikes when a rigid container is overfilled.
    /// Only relevant for liquids and solids.
    bulk_modulus: f64 = 2.0e9, // Default to water-like
    /// Partial molar volume in cm³/mol at infinite dilution in water at 25°C.
    /// Used to calculate the volume occupied by dissolved gases/solids.
    partial_molar_volume: ?f64 = null,

    /// Calculates the net molecular charge at a given pH using the
    /// Henderson-Hasselbalch equation.
    ///
    /// For each ionizable group with a given pKa, the fraction deprotonated
    /// (alpha) is: α = 1 / (1 + 10^(pKa - pH))
    ///
    /// Starting from the fully protonated charge, each deprotonation event
    /// decreases the charge by 1.
    pub inline fn getCharge(this: Molecule, ph: f64) f64 {
        var charge: f64 = @floatFromInt(this.fully_protonated_charge);

        for (this.pka) |pka| {
            // Fraction of molecules that have lost this proton at the given pH.
            // When pH >> pKa, alpha -> 1 (fully deprotonated).
            // When pH << pKa, alpha -> 0 (fully protonated).
            const pka_f64: f64 = @floatCast(pka);
            const ph_f64: f64 = @floatCast(ph);
            const alpha = 1.0 / (1.0 + std.math.pow(f64, 10.0, pka_f64 - ph_f64));

            charge -= alpha;
        }

        return charge;
    }

    /// Returns the topological polar surface area in Ų.
    /// If not explicitly set, estimates using Ertl's approximation:
    /// PSA ≈ 23.0 * HBD + 9.2 * HBA
    pub inline fn getPSA(this: Molecule) f64 {
        const psa = this.psa;

        if (psa > 0) {
            return psa;
        }

        return @as(f64, @floatFromInt(this.hbd)) * 23.0 +
            @as(f64, @floatFromInt(this.hba)) * 9.2;
    }

    pub inline fn getPressureCorrectedBoilingPoint(this: Molecule, pressure: f64) ?f64 {
        const bp_normal = this.boiling_point orelse {
            return null;
        };

        var bp_corr: f64 = @floatCast(bp_normal);

        if (this.delta_hvap) |hvap_kj| {
            const hvap_j = hvap_kj * 1000.0;
            const ln_ratio = @log(@as(f64, @floatCast(pressure)) / constants.P_std);
            const inv_T2 = (1.0 / bp_corr) - (constants.R * ln_ratio / hvap_j);

            if (inv_T2 > 0.0) {
                bp_corr = 1.0 / inv_T2;
            }
        }

        return bp_corr;
    }

    /// Returns the enthalpy of vaporization in J/mol, or 30000.0 as a fallback.
    pub inline fn getDeltaHvapJ(this: Molecule) f64 {
        return if (this.delta_hvap) |h|
            h * 1000.0
        else
            30000.0;
    }

    /// Returns the enthalpy of fusion in J/mol, or 10000.0 as a fallback.
    pub inline fn getDeltaHfusJ(this: Molecule) f64 {
        return if (this.delta_hfus) |h|
            h * 1000.0
        else
            10000.0;
    }

    /// Returns the Hildebrand solubility parameter in MPa^½.
    /// If not set explicitly, estimates from ΔH_vap:
    /// δ = √((ΔH_vap - RT) / V_m)
    /// where V_m = M / ρ (molar volume in cm³/mol).
    pub inline fn getHildebrand(this: Molecule) f64 {
        if (this.hildebrand) |h| {
            return h;
        }

        if (this.delta_hvap) |hvap| {
            const rt = 2.478; // kJ/mol at 298 K
            const vm = this.weight / this.density; // cm³/mol
            // (kJ/mol * 1000) / (cm³/mol) = J/cm³ = MPa
            const delta_sq = (hvap - rt) * 1000.0 / vm;

            if (delta_sq > 0) {
                return @sqrt(delta_sq);
            }
        }

        return 20.0; // fallback: typical organic molecule
    }

    /// Returns the equilibrium phase of a pure substance at the given T, P.
    ///
    /// The boiling point is corrected for pressure via Clausius-Clapeyron:
    ///
    /// 1/T_b(P) = 1/T_b° - (R / ΔH_vap) * ln(P / P°)
    ///
    /// When T_b(P) drops below T_melt (very low pressure) the liquid
    /// window vanishes - the substance sublimes. Falls back to the
    /// molecule's catalogued standard state when no transition data exists.
    pub inline fn targetState(this: Molecule, temperature: f64, pressure: f64) Molecule.State {
        const T: f64 = @floatCast(temperature);

        const mp: ?f64 = if (this.melting_point) |m|
            @as(f64, @floatCast(m))
        else
            null;

        const bp: ?f64 = this.getPressureCorrectedBoilingPoint(pressure);

        // No thermodynamic data - fall back to the catalogued state at 25 °C.
        if (mp == null and bp == null) {
            return this.state;
        }

        // Sublimation regime: liquid window has collapsed.
        if (mp != null and bp != null and bp.? <= mp.?) {
            return if (T < bp.?)
                .solid
            else
                .gas;
        }

        // Normal sequence: solid -> liquid -> gas.
        if (mp) |melting| {
            if (T < melting) {
                return .solid;
            }
        }

        if (bp) |boiling| {
            if (T >= boiling) {
                return .gas;
            }
        }

        return .liquid;
    }

    /// Determines whether two liquids form a single phase (miscible) or
    /// separate into immiscible layers.
    ///
    /// Uses a modified Hildebrand criterion with hydrogen-bonding correction:
    ///
    /// 1. Base criterion: |δ_a - δ_b| < threshold (7.0 MPa^½)
    ///    Works for non-polar/weakly-polar pairs (alkanes, ethers, DCM, etc.)
    ///
    /// 2. H-bond correction: if both molecules can donate AND accept H-bonds,
    ///    they form a hydrogen-bonded network that stabilizes mixing beyond
    ///    what δ alone predicts. The effective δ distance is reduced.
    ///
    /// 3. Full miscibility override: if both molecules have no known solid
    ///    form limit (water_solubility = null) AND share H-bonding capability,
    ///    they are considered miscible (e.g., water + ethanol, water + methanol).
    pub inline fn areMiscible(a: Molecule, b: Molecule) bool {
        const delta_a: f64 = @floatCast(a.getHildebrand());
        const delta_b: f64 = @floatCast(b.getHildebrand());
        const delta_diff = @abs(delta_a - delta_b);

        // Both can form hydrogen bonds (have both donors and acceptors)?
        const a_hbond = a.hbd > 0 and a.hba > 0;
        const b_hbond = b.hbd > 0 and b.hba > 0;
        const mutual_hbond = a_hbond and b_hbond;

        // If both are fully water-miscible liquids with mutual H-bonding,
        // they are miscible regardless of δ difference.
        // (water_solubility == null for liquids means "miscible in all proportions")
        if (mutual_hbond) {
            const a_miscible = a.state == .liquid and a.water_solubility == null;
            const b_miscible = b.state == .liquid and b.water_solubility == null;

            if (a_miscible and b_miscible) {
                return true;
            }
        }

        // Base Hildebrand criterion.
        const threshold: f64 = if (mutual_hbond)
            25.0
        else
            7.0;

        if (delta_diff < threshold) {
            return true;
        }

        return false;
    }

    /// Calculates the saturation vapor pressure of the pure substance at a given temperature
    /// using the Clausius-Clapeyron equation.
    /// Returns 0.0 if thermodynamic data (boiling_point, delta_hvap) is missing.
    pub inline fn getVaporPressure(this: Molecule, temperature: f64) f64 {
        const bp_normal = this.boiling_point orelse {
            return 0.0;
        };

        const dH_vap = this.delta_hvap orelse {
            return 0.0;
        };

        const T_b: f64 = @floatCast(bp_normal);

        if (T_b <= 0.0 or temperature <= 0.0) {
            return 0.0;
        }

        const hvap_j = dH_vap * 1000.0;
        // ln(P / P_std) = - (dH_vap / R) * (1/T - 1/T_b)
        const ln_ratio = -hvap_j / constants.R * (1.0 / temperature - 1.0 / T_b);

        return constants.P_std * @exp(ln_ratio);
    }
};
