// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const std = @import("std");

const Color = @import("color.zig").Color;
const constants = @import("constants.zig");
const Flavor = @import("flavor.zig").Flavor;
const MoleculeId = @import("molecule_id.zig").MoleculeId;
const Odor = @import("odor.zig").Odor;
const ScoredAttribute = @import("scored_attribute.zig").ScoredAttribute;

pub const GasPhase = struct {
    /// Array of mole amounts indexed by Molecule enum ordinal.
    molecules: []f64,

    pub inline fn init(allocator: std.mem.Allocator) error{OutOfMemory}!GasPhase {
        const molecules = try allocator.alloc(f64, MoleculeId.count());
        @memset(molecules, 0.0);

        return .{
            .molecules = molecules,
        };
    }

    pub inline fn deinit(this: *GasPhase, allocator: std.mem.Allocator) void {
        allocator.free(this.molecules);
        this.molecules = &.{};
    }

    pub inline fn addMolecule(this: *GasPhase, molecule: MoleculeId, moles: f64) void {
        this.molecules[@intFromEnum(molecule)] += moles;
    }

    pub inline fn addMoles(this: *GasPhase, moles: []const f64) void {
        std.debug.assert(this.molecules.len == moles.len);

        for (0..this.molecules.len) |i| {
            this.molecules[i] += moles[i];
        }
    }

    /// Returns the color of this gas phase.
    ///
    /// Most gases are colorless. Colored gases (Cl2, Br2, NO2, I2 vapor)
    /// are relatively rare. Color visibility depends on:
    /// 1. Partial pressure of the colored species
    /// 2. Path length through the gas (volume^(1/3))
    /// 3. Intrinsic color intensity
    ///
    /// Uses Beer-Lambert law for gases: A = ε * c * l
    /// where c is derived from partial pressure via ideal gas law.
    pub fn getColors(this: *const GasPhase, temperature: f64, pressure: f64) Flavor {
        var flavor = Flavor{};

        // Total moles for partial pressure calculation
        var total_moles: f64 = 0.0;

        for (this.molecules) |moles| {
            total_moles += moles;
        }

        if (constants.isNegligible(total_moles)) {
            return flavor;
        }

        // Gas volume via ideal gas law: V = nRT/P [L]
        const R: f64 = 8.314; // J/(mol·K)
        const volume_m3 = total_moles * R * temperature / pressure;
        const volume_L = volume_m3 * 1000.0;

        // Path length approximation: cube root of volume in cm
        const path_length: f64 = std.math.pow(f64, volume_L * 1000.0, 1.0 / 3.0);

        var color_scores: [@typeInfo(Color).@"enum".fields.len]f64 =
            .{0.0} ** @typeInfo(Color).@"enum".fields.len;

        for (0..this.molecules.len) |i| {
            const moles = this.molecules[i];

            if (constants.isNegligible(moles)) {
                continue;
            }

            const mid: MoleculeId = @enumFromInt(i);
            const meta = mid.molecule();

            if (meta.color == .transparent) {
                continue;
            }

            if (meta.color_intensity < constants.MIN_COLOR_INTENSITY) {
                continue;
            }

            // Concentration from partial pressure: c = n/V
            const concentration = moles / volume_L;

            // Gases are generally much less intensely colored than liquids
            // Scale down by factor of 10
            const score = concentration * @as(f64, meta.color_intensity) * path_length * 0.1;

            color_scores[@intFromEnum(meta.color)] += score;
        }

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

    /// Returns the dominant odors of this gas phase.
    ///
    /// Gas-phase odor is the most direct form of smell perception.
    /// Odor intensity depends on:
    /// 1. Partial pressure of the odorant (mole fraction * total pressure)
    /// 2. Odor intensity coefficient (inverse of odor threshold)
    /// 3. Absolute amount (more gas = more molecules reaching nose)
    ///
    /// Score = (n_i / n_total) * P_total * odor_intensity * volume_factor
    pub fn getOdors(this: *const GasPhase, temperature: f64, pressure: f64) Flavor {
        var flavor = Flavor{};

        var total_moles: f64 = 0.0;

        for (this.molecules) |moles| {
            total_moles += moles;
        }

        if (constants.isNegligible(total_moles)) {
            return flavor;
        }

        // Volume for scaling (larger volumes = more dilute but more total molecules)
        const R: f64 = 8.314;
        const volume_L = total_moles * R * temperature / pressure * 1000.0;

        // Volume factor: cube root for realistic scaling
        const volume_factor = std.math.pow(f64, @max(volume_L, 0.001), 1.0 / 3.0);

        var odor_scores: [@typeInfo(Odor).@"enum".fields.len]f64 =
            .{0.0} ** @typeInfo(Odor).@"enum".fields.len;

        for (0..this.molecules.len) |i| {
            const moles = this.molecules[i];

            if (constants.isNegligible(moles)) {
                continue;
            }

            const mid: MoleculeId = @enumFromInt(i);
            const meta = mid.molecule();

            if (meta.odor == .none) {
                continue;
            }

            if (meta.odor_intensity < constants.MIN_ODOR_INTENSITY) {
                continue;
            }

            // Partial pressure contribution (in atm for intuitive scaling)
            const mole_fraction = moles / total_moles;
            const partial_pressure_atm = mole_fraction * pressure / 101325.0;

            // Gas-phase odors are perceived directly, no evaporation barrier
            // Scale by 10x compared to liquid evaporation
            const score = partial_pressure_atm * @as(f64, meta.odor_intensity) *
                volume_factor * 10.0;

            odor_scores[@intFromEnum(meta.odor)] += score;
        }

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

    /// Returns combined color and odor profile for this gas phase.
    pub inline fn getFlavor(this: *const GasPhase, temperature: f64, pressure: f64) Flavor {
        var flavor = this.getColors(temperature, pressure);
        const odor_flavor = this.getOdors(temperature, pressure);

        flavor.odors = odor_flavor.odors;
        flavor.odor_intensities = odor_flavor.odor_intensities;

        return flavor;
    }
};
