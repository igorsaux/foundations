// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const std = @import("std");

const constants = @import("constants.zig");
const Contents = @import("contents.zig").Contents;
const Molecule = @import("molecule_id.zig").MoleculeId;

/// Compute pH of a liquid phase from its composition.
pub inline fn computePhasePH(moles: []const f64, volume: f64) f64 {
    std.debug.assert(moles.len == Molecule.count());

    if (volume < constants.MIN_VOLUME) {
        return 7.0;
    }

    const Kw: f64 = 1e-14;
    const max_iterations: usize = 50;
    const tolerance: f64 = 1e-10;

    var static_charge: f64 = 0.0;

    for (0..moles.len) |i| {
        const mol: Molecule = @enumFromInt(i);
        const conc = moles[i] / volume;

        if (mol == .oxidane or mol == .hydrogen_ion or mol.getPka().len > 0) {
            continue;
        }

        static_charge += @as(f64, @floatFromInt(mol.getFullyProtonatedCharge())) * conc;
    }

    var ph: f64 = 7.0;

    for (0..max_iterations) |_| {
        const h_conc = std.math.pow(f64, 10.0, -ph);
        const oh_conc = Kw / h_conc;

        var total_charge: f64 = static_charge + h_conc - oh_conc;
        var d_charge_d_ph: f64 = -@log(10.0) * h_conc - @log(10.0) * oh_conc;

        for (0..moles.len) |i| {
            const molecule: Molecule = @enumFromInt(i);
            const n = moles[i];

            if (constants.isNegligible(n)) {
                continue;
            }

            const pkas = molecule.getPka();

            if (pkas.len == 0) {
                continue;
            }

            const conc = n / volume;
            var charge: f64 = @floatFromInt(molecule.getFullyProtonatedCharge());
            var d_charge: f64 = 0.0;

            for (pkas) |pka| {
                const exp_term = std.math.pow(f64, 10.0, @as(f64, @floatCast(pka)) - ph);
                const alpha = 1.0 / (1.0 + exp_term);
                const d_alpha = @log(10.0) * alpha * (1.0 - alpha);

                charge -= alpha;
                d_charge -= d_alpha;
            }

            total_charge += charge * conc;
            d_charge_d_ph += d_charge * conc;
        }

        if (@abs(d_charge_d_ph) < 1e-30) {
            break;
        }

        var delta = -total_charge / d_charge_d_ph;
        delta = std.math.clamp(delta, -2.0, 2.0);

        ph += delta;
        ph = std.math.clamp(ph, -1.0, 15.0);

        if (@abs(delta) < tolerance) {
            break;
        }
    }

    return @floatCast(ph);
}

pub inline fn resetStdAir(dst: *Contents, max_volume: f64) void {
    const headspace = dst.getGasVolume(max_volume);

    @memset(dst.gas.molecules, 0.0);

    const total_moles = (constants.P_std * (headspace / 1000.0)) / (8.314 * constants.T_std);

    if (!constants.isNegligible(total_moles)) {
        dst.gas.addMolecule(.nitrogen, total_moles * 0.78084);
        dst.gas.addMolecule(.oxygen, total_moles * 0.20946);
        dst.gas.addMolecule(.argon, total_moles * 0.00934);
        dst.gas.addMolecule(.carbon_dioxide, total_moles * 0.0004);
        dst.gas.addMolecule(.oxidane, total_moles * 0.015);
    }
}
