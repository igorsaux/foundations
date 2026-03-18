// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const std = @import("std");

const Color = @import("color.zig").Color;
const constants = @import("constants.zig");
const Molecule = @import("molecule.zig").Molecule;
const Odor = @import("odor.zig").Odor;

pub const MoleculeId = enum(u16) {
    /// Returns the total number of defined molecules (compile-time).
    pub inline fn count() comptime_int {
        return @typeInfo(MoleculeId).@"enum".fields.len;
    }

    /// Returns molecular weight in g/mol.
    pub inline fn getWeight(this: MoleculeId) f64 {
        return this.molecule().weight;
    }

    /// Returns mass density in g/cm³ at 25°C.
    pub inline fn getDensity(this: MoleculeId) f64 {
        return this.molecule().density;
    }

    /// Returns the list of pKa values (ascending order).
    pub inline fn getPka(this: MoleculeId) []const f64 {
        return this.molecule().pka;
    }

    /// Returns the net charge of the fully protonated species.
    pub inline fn getFullyProtonatedCharge(this: MoleculeId) i8 {
        return this.molecule().fully_protonated_charge;
    }

    /// Calculates the net molecular charge at a given pH using the
    /// Henderson-Hasselbalch equation.
    ///
    /// For each ionizable group with a given pKa, the fraction deprotonated
    /// (alpha) is: α = 1 / (1 + 10^(pKa - pH))
    ///
    /// Starting from the fully protonated charge, each deprotonation event
    /// decreases the charge by 1.
    pub inline fn getCharge(this: MoleculeId, ph: f64) f64 {
        return this.molecule().getCharge(ph);
    }

    /// Returns the octanol-water partition coefficient (log10 P).
    pub inline fn getLogP(this: MoleculeId) f64 {
        return this.molecule().log_p;
    }

    /// Returns the number of hydrogen bond donors.
    pub inline fn getHBD(this: MoleculeId) u8 {
        return this.molecule().hbd;
    }

    /// Returns the number of hydrogen bond acceptors.
    pub inline fn getHBA(this: MoleculeId) u8 {
        return this.molecule().hba;
    }

    /// Returns the topological polar surface area in Ų.
    /// If not explicitly set, estimates using Ertl's approximation:
    /// PSA ≈ 23.0 * HBD + 9.2 * HBA
    pub inline fn getPSA(this: MoleculeId) f64 {
        return this.molecule().getPSA();
    }

    /// Returns the normal boiling point in Kelvin (at 1 atm), or null if unknown.
    pub inline fn getBoilingPoint(this: MoleculeId) ?f64 {
        return this.molecule().boiling_point;
    }

    pub inline fn getPressureCorrectedBoilingPoint(this: MoleculeId, pressure: f64) ?f64 {
        return this.molecule().getPressureCorrectedBoilingPoint(pressure);
    }

    /// Returns the melting point in Kelvin, or null if unknown.
    pub inline fn getMeltingPoint(this: MoleculeId) ?f64 {
        return this.molecule().melting_point;
    }

    /// Returns aqueous solubility in mol/L, or null if unknown / miscible.
    pub inline fn getWaterSolubility(this: MoleculeId) ?f64 {
        return this.molecule().water_solubility;
    }

    /// Returns the standard enthalpy of formation in kJ/mol.
    pub inline fn getDeltaHf(this: MoleculeId) f64 {
        return this.molecule().delta_hf;
    }

    /// Returns the enthalpy of vaporization in kJ/mol, or null if unknown.
    pub inline fn getDeltaHvap(this: MoleculeId) ?f64 {
        return this.molecule().delta_hvap;
    }

    /// Returns the enthalpy of vaporization in J/mol, or 30000.0 as a fallback.
    pub inline fn getDeltaHvapJ(this: MoleculeId) f64 {
        return this.molecule().getDeltaHvapJ();
    }

    /// Returns the enthalpy of fusion (melting) in kJ/mol, or null if unknown.
    pub inline fn getDeltaHfus(this: MoleculeId) ?f64 {
        return this.molecule().delta_hfus;
    }

    /// Returns the enthalpy of fusion in J/mol, or 10000.0 as a fallback.
    pub inline fn getDeltaHfusJ(this: MoleculeId) f64 {
        return this.molecule().getDeltaHfusJ();
    }

    /// Returns the molar heat capacity in J/(mol·K).
    pub inline fn getHeatCapacity(this: MoleculeId) f64 {
        return @as(f64, @floatCast(this.molecule().heat_capacity));
    }

    /// Returns the Hildebrand solubility parameter in MPa^½.
    /// If not set explicitly, estimates from ΔH_vap:
    /// δ = √((ΔH_vap - RT) / V_m)
    /// where V_m = M / ρ (molar volume in cm³/mol).
    pub inline fn getHildebrand(this: MoleculeId) f64 {
        return this.molecule().getHildebrand();
    }

    /// Returns the physical state at 25°C and 1 atm.
    pub inline fn getState(this: MoleculeId) Molecule.State {
        return this.molecule().state;
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
    pub inline fn targetState(this: MoleculeId, temperature: f64, pressure: f64) Molecule.State {
        return this.molecule().targetState(temperature, pressure);
    }

    pub inline fn isLiquid(this: MoleculeId) bool {
        return this.getState() == .liquid;
    }

    pub inline fn isGas(this: MoleculeId) bool {
        return this.getState() == .gas;
    }

    pub inline fn isSolid(this: MoleculeId) bool {
        return this.getState() == .solid;
    }

    /// Returns plasma protein binding fraction (0.0–1.0).
    pub inline fn getProteinBinding(this: MoleculeId) f32 {
        return this.molecule().protein_binding;
    }

    /// Returns gas-lipid partition coefficient (dimensionless).
    pub inline fn getGasLipidPartition(this: MoleculeId) f64 {
        return this.molecule().gas_lipid_partition;
    }

    pub inline fn isIon(this: MoleculeId) bool {
        return this.molecule().is_ion;
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
    pub inline fn areMiscible(a: MoleculeId, b: MoleculeId) bool {
        return a.molecule().areMiscible(b.molecule());
    }

    /// Calculates the saturation vapor pressure of the pure substance at a given temperature
    /// using the Clausius-Clapeyron equation.
    /// Returns 0.0 if thermodynamic data (boiling_point, delta_hvap) is missing.
    pub inline fn getVaporPressure(this: MoleculeId, temperature: f64) f64 {
        return this.molecule().getVaporPressure(temperature);
    }

    pub inline fn getPartialMolarVolume(this: MoleculeId) ?f64 {
        return this.molecule().partial_molar_volume;
    }

    /// Returns a human-readable display name for the molecule.
    pub inline fn getDisplayName(this: MoleculeId) []const u8 {
        return switch (this) {
            .acetate_ion => "Acetate ion",
            .acetic_acid => "Acetic Acid",
            .argon => "Argon",
            .azane => "Ammonia",
            .calcium_carbonate => "Calcium Carbonate",
            .calcium_ion => "Ca ion",
            .carbon_dioxide => "Carbon Dioxide",
            .carbonate_ion => "CO ion",
            .chloride_ion => "Cl ion",
            .dichloromethane => "Dichloromethane",
            .ethanol => "Ethanol",
            .graphite => "Graphite",
            .hydrogen_ion => "Hydrogen ion",
            .nitrogen => "Nitrogen",
            .oxidane => "Water",
            .oxygen => "Oxygen",
            .sodium_chloride => "Sodium Chloride",
            .sodium_ion => "Na ion",
        };
    }

    pub inline fn molesPerLiter(this: MoleculeId) f64 {
        return (1000.0 * this.getDensity()) / this.getWeight();
    }

    pub inline fn molesPerG(this: MoleculeId) f64 {
        return 1.0 / this.getWeight();
    }

    pub inline fn molecule(this: MoleculeId) Molecule {
        return switch (this) {
            .acetate_ion => .{ .weight = 59.044, .density = 1.0, .fully_protonated_charge = -1, .pka = &.{4.76}, .log_p = -3.72, .hbd = 0, .hba = 2, .state = .solid, .water_solubility = 5.7, .delta_hf = -486.0, .heat_capacity = 80.0, .is_ion = true },
            .acetic_acid => .{ .weight = 60.052, .density = 1.049, .pka = &.{4.76}, .log_p = -0.17, .hbd = 1, .hba = 2, .state = .liquid, .boiling_point = 391.1, .melting_point = 289.8, .water_solubility = null, .delta_hf = -484.5, .delta_hvap = 23.7, .delta_hfus = 11.73, .hildebrand = 20.7, .heat_capacity = 123.1, .odor = .acrid, .odor_intensity = 0.7 },
            .argon => .{ .weight = 39.948, .density = 0.001633, .log_p = 0.0, .hbd = 0, .hba = 0, .state = .gas, .boiling_point = 87.30, .melting_point = 83.81, .water_solubility = 0.0014, .delta_hf = 0.0, .delta_hvap = 6.45, .delta_hfus = 1.18, .heat_capacity = 20.79, .partial_molar_volume = 32.0 },
            .azane => .{ .weight = 17.031, .density = 0.73, .log_p = -1.38, .hbd = 3, .hba = 1, .state = .gas, .boiling_point = 239.8, .melting_point = 195.4, .water_solubility = null, .delta_hf = -45.9, .delta_hvap = 23.35, .delta_hfus = 5.66, .heat_capacity = 35.06, .odor = .ammoniacal, .odor_intensity = 1.0, .partial_molar_volume = 24.6 },
            .calcium_carbonate => .{ .weight = 100.086, .density = 2.71, .log_p = -2.12, .hbd = 0, .hba = 3, .psa = 63.19, .state = .solid, .melting_point = 1100.0, .water_solubility = 0.00013, .delta_hf = -1207.6, .delta_hfus = 36.0, .hildebrand = 25.0, .heat_capacity = 81.9, .color = .white, .color_intensity = 0.1, .odor = .none, .odor_intensity = 0.0, .bulk_modulus = 15.0e9 },
            .calcium_ion => .{ .weight = 40.078, .density = 1.55, .fully_protonated_charge = 2, .log_p = -2.5, .hbd = 0, .hba = 0, .state = .solid, .melting_point = 1115.0, .boiling_point = 1757.0, .water_solubility = null, .delta_hf = -542.8, .heat_capacity = 25.0, .is_ion = true },
            .carbon_dioxide => .{ .weight = 44.01, .density = 0.001799, .log_p = 0.0, .hbd = 0, .hba = 2, .state = .gas, .boiling_point = 194.7, .melting_point = 216.55, .water_solubility = 0.034, .delta_hf = -393.5, .delta_hvap = 16.7, .delta_hfus = 8.33, .heat_capacity = 37.11, .partial_molar_volume = 33.0 },
            .carbonate_ion => .{ .weight = 60.008, .density = 1.0, .fully_protonated_charge = -2, .pka = &.{ 6.35, 10.33 }, .log_p = -4.0, .hbd = 0, .hba = 3, .state = .solid, .water_solubility = 2.0, .delta_hf = -677.1, .heat_capacity = 50.0, .is_ion = true },
            .chloride_ion => .{ .weight = 35.453, .density = 1.56, .fully_protonated_charge = -1, .log_p = -3.5, .hbd = 0, .hba = 0, .psa = 0, .state = .liquid, .boiling_point = null, .melting_point = null, .water_solubility = null, .delta_hf = -167.08, .hildebrand = null, .protein_binding = 0.0, .gas_lipid_partition = 0.0, .delta_hvap = null, .heat_capacity = 25.0, .color = .transparent, .color_intensity = 0.0, .odor = .none, .odor_intensity = 0.0, .is_ion = true },
            .dichloromethane => .{ .weight = 84.93, .density = 1.322, .log_p = 1.2, .hbd = 0, .hba = 0, .boiling_point = 312.8, .melting_point = 178.0, .delta_hf = 124.3, .delta_hvap = 28.06, .delta_hfus = 4.60, .heat_capacity = 100.0, .odor = .sweet, .odor_intensity = 0.7 },
            .ethanol => .{ .weight = 46.07, .density = 0.79, .pka = &.{15.9}, .log_p = -0.31, .hbd = 1, .hba = 1, .boiling_point = 351.5, .melting_point = 159.0, .delta_hf = -277.6, .delta_hvap = 38.56, .delta_hfus = 4.93, .hildebrand = 26.5, .heat_capacity = 112.4, .odor = .ethereal, .odor_intensity = 0.4, .bulk_modulus = 1.1e9 },
            .graphite => .{ .weight = 12.011, .density = 2.26, .log_p = 0.0, .hbd = 0, .hba = 0, .psa = 0.0, .state = .solid, .melting_point = 3800.0, .water_solubility = 0.0, .delta_hf = 0.0, .delta_hfus = 117.0, .hildebrand = 0.0, .heat_capacity = 8.52, .color = .black, .color_intensity = 1.0, .odor = .none, .odor_intensity = 0.0, .bulk_modulus = 15.0e9 },
            .hydrogen_ion => .{ .weight = 1.008, .density = 1.0, .fully_protonated_charge = 1, .log_p = -6.0, .hbd = 1, .hba = 0, .state = .gas, .water_solubility = null, .delta_hf = 0.0, .heat_capacity = 10.0, .is_ion = true },
            .nitrogen => .{ .weight = 28.014, .density = 0.001145, .log_p = 0.0, .hbd = 0, .hba = 0, .state = .gas, .boiling_point = 77.36, .melting_point = 63.15, .water_solubility = 0.0006, .delta_hf = 0.0, .delta_hvap = 5.56, .delta_hfus = 0.72, .heat_capacity = 29.12, .partial_molar_volume = 34.0 },
            .oxidane => .{ .weight = 18.015, .density = 0.9970, .log_p = -1.38, .hbd = 2, .hba = 1, .boiling_point = 373.13, .melting_point = 273.15, .delta_hf = -285.8, .delta_hvap = 40.65, .delta_hfus = 6.01, .hildebrand = 47.8, .heat_capacity = 75.38, .bulk_modulus = 2.2e9 },
            .oxygen => .{ .weight = 31.998, .density = 0.001308, .log_p = 0.0, .hbd = 0, .hba = 0, .state = .gas, .boiling_point = 90.19, .melting_point = 54.36, .water_solubility = 0.0012, .delta_hf = 0.0, .delta_hvap = 6.82, .delta_hfus = 0.44, .heat_capacity = 29.38, .partial_molar_volume = 31.0 },
            .sodium_chloride => .{ .weight = 58.44, .density = 2.165, .log_p = -3.0, .hbd = 0, .hba = 0, .psa = 0, .state = .solid, .boiling_point = 1686.0, .melting_point = 1074.0, .water_solubility = 6.15, .delta_hf = -411.12, .delta_hfus = 28.16, .hildebrand = 25.0, .protein_binding = 0.0, .gas_lipid_partition = 0.0, .delta_hvap = 207.0, .heat_capacity = 50.5, .color = .white, .color_intensity = 0.8, .odor = .none, .odor_intensity = 0.0, .bulk_modulus = 25.0e9 },
            .sodium_ion => .{ .weight = 22.990, .density = 0.968, .fully_protonated_charge = 1, .log_p = -3.5, .hbd = 0, .hba = 0, .psa = 0, .state = .liquid, .boiling_point = null, .melting_point = null, .water_solubility = null, .delta_hf = -240.34, .hildebrand = null, .protein_binding = 0.0, .gas_lipid_partition = 0.0, .delta_hvap = null, .heat_capacity = 25.0, .color = .transparent, .color_intensity = 0.0, .odor = .none, .odor_intensity = 0.0, .is_ion = true },
        };
    }

    acetate_ion,
    acetic_acid,
    argon,
    azane,
    calcium_carbonate,
    calcium_ion,
    carbon_dioxide,
    carbonate_ion,
    chloride_ion,
    dichloromethane,
    ethanol,
    graphite,
    hydrogen_ion,
    nitrogen,
    oxidane,
    oxygen,
    sodium_chloride,
    sodium_ion,
};
