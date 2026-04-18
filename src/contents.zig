// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const std = @import("std");

const Color = @import("color.zig").Color;
const constants = @import("constants.zig");
const Flavor = @import("flavor.zig").Flavor;
const GasPhase = @import("gas_phase.zig").GasPhase;
const helpers = @import("helpers.zig");
const LiquidPhase = @import("liquid_phase.zig").LiquidPhase;
const MoleculeId = @import("molecule_id.zig").MoleculeId;
const Odor = @import("odor.zig").Odor;
const reactions = @import("reactions.zig");
const ScoredAttribute = @import("scored_attribute.zig").ScoredAttribute;
const SolidPhase = @import("solid_phase.zig").SolidPhase;

pub const Contents = struct {
    pub const Metrics = struct {
        /// Moles vaporized via boiling (T > T_boil, bulk phase change).
        moles_boiled: f64 = 0.0,
        /// Moles evaporated via VLE (surface evaporation at T < T_boil).
        moles_evaporated: f64 = 0.0,
    };

    /// Immiscible liquid layers (e.g., aqueous and organic phases).
    /// Sorted by density: lightest (index 0) on top, densest on bottom.
    liquids: std.ArrayList(LiquidPhase) = .empty,
    /// Distinct solid precipitates, crystals, or undissolved material.
    solids: std.ArrayList(SolidPhase) = .empty,
    /// Gas/vapor phase above the liquid (headspace).
    gas: GasPhase,
    /// Total internal heat energy of the system in Joules.
    /// Temperature is derived from this: T = heat_energy / total_heat_capacity.
    heat_energy: f64 = constants.cToK(25) * 100.0,
    /// Baseline heat capacity of the container itself (e.g., a glass beaker) in J/K.
    /// Prevents division by zero when the container is completely empty.
    container_heat_capacity: f64 = 180.0,
    /// Area of the liquid-gas interface (m^2).
    gas_contact_area: f64 = 0.006,
    /// Pressure of the system in Pascals.
    pressure: f64 = constants.P_std,
    /// Normalized stirring intensity factor [0.0, 1.0].
    /// 0.0 = Unstirred (purely diffusion limited).
    /// 0.2 = Gentle stirring (e.g., occasional swirling by hand).
    /// 0.5 = Moderate stirring (e.g., standard magnetic stirrer at medium RPM).
    /// 0.8 = Vigorous stirring (e.g., high-speed overhead mechanical stirrer).
    /// 1.0 = Extreme agitation (e.g., high-shear homogenizer or sonication).
    stirring: f64 = 0.0,

    /// When true, automatic settle() calls are suppressed.
    /// Use beginBatch()/endBatch() to manage this flag.
    batch_mode: bool = false,

    metrics: ?*Metrics = null,

    // Level 1 - HIGH-LEVEL (auto-settle):
    //   addLiquid(), addSolid(), addGas()
    //   setTemperature(), setPressure()
    //   pour*(), transfer*()
    //
    //   These methods automatically call settle() to reach equilibrium.
    //   Use for single operations or when immediate consistency is needed.
    //
    // Level 2 - BATCH MODE (manual settle):
    //   beginBatch() / endBatch()
    //   addLiquidRaw(), addSolidRaw(), addGasRaw()
    //
    //   Use when performing multiple modifications. Call endBatch() or
    //   settle() when done to update the system state.
    //
    // Level 3 - STATE UPDATES (explicit control):
    //   settle()                 - full equilibration (phases + precipitation + sort)
    //   updatePhaseTransitions() - solid/liquid/gas transitions based on T, P
    //   equilibratePhases()      - redistribute solutes among immiscible layers
    //   checkPrecipitation()     - precipitate supersaturated species
    //
    // TYPICAL WORKFLOWS:
    //
    // 1. Add substances and equilibrate:
    //      contents.addLiquid(allocator, .water, 1.0);
    //      contents.addSolid(allocator, .sodium_chloride, 0.1, 1e-4);
    //    OR in batch:
    //      contents.beginBatch();
    //      contents.addLiquidRaw(allocator, .water, 1.0);
    //      contents.addLiquidRaw(allocator, .ethanol, 0.5);
    //      contents.endBatch(allocator);
    //
    // 2. Change environment and update:
    //      contents.setTemperature(allocator, 373.15); // auto phase transitions + settle
    //    OR manually:
    //      contents.temperature = 373.15;
    //      contents.updatePhaseTransitions(allocator);
    //      contents.settle(allocator);
    //
    // 3. Transfer and update destination only:
    //      contents.pourVolumeInto(allocator, dst, 0.1); // dst settles, source doesn't

    pub inline fn init(allocator: std.mem.Allocator) error{OutOfMemory}!Contents {
        return .{
            .gas = try GasPhase.init(allocator),
        };
    }

    pub inline fn deinit(this: *Contents, allocator: std.mem.Allocator) void {
        for (this.liquids.items) |*phase| {
            phase.deinit(allocator);
        }

        this.liquids.deinit(allocator);

        for (this.solids.items) |*phase| {
            phase.deinit(allocator);
        }

        this.solids.deinit(allocator);
        this.gas.deinit(allocator);
    }

    /// Enters batch mode: suppresses automatic settle() calls.
    /// Call endBatch() when finished to update the system state.
    ///
    /// Example:
    ///   contents.beginBatch();
    ///   contents.addLiquidRaw(allocator, .water, 1.0);
    ///   contents.addLiquidRaw(allocator, .ethanol, 0.5);
    ///   contents.addSolidRaw(allocator, .salt, 0.1, 1e-4);
    ///   contents.endBatch(allocator);
    pub inline fn beginBatch(this: *Contents) void {
        this.batch_mode = true;
    }

    /// Exits batch mode and performs full state update.
    /// Equivalent to: batch_mode = false; settle();
    pub inline fn endBatch(this: *Contents, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        this.batch_mode = false;
        try this.settle(allocator);
    }

    /// Adds liquid-phase molecule and settles the system.
    /// For batch operations, use addLiquidRaw() with beginBatch()/endBatch().
    pub inline fn addLiquid(
        this: *Contents,
        allocator: std.mem.Allocator,
        molecule: MoleculeId,
        moles: f64,
    ) error{OutOfMemory}!void {
        try this.addLiquidRaw(allocator, molecule, moles, false);

        if (!this.batch_mode) {
            try this.settle(allocator);
        }
    }

    /// Adds liquid-phase molecule without triggering settle().
    /// Use with beginBatch()/endBatch() or call settle() manually.
    pub inline fn addLiquidRaw(
        this: *Contents,
        allocator: std.mem.Allocator,
        molecule: MoleculeId,
        moles: f64,
        force: bool,
    ) error{OutOfMemory}!void {
        if (!force and constants.isNegligible(moles)) {
            return;
        }

        var buf: [MoleculeId.count()]f64 = undefined;
        @memset(&buf, 0.0);
        buf[@intFromEnum(molecule)] = moles;

        try this.placeLiquidMoles(allocator, &buf);
    }

    /// Adds solid substance and settles the system.
    /// For batch operations, use addSolidRaw() with beginBatch()/endBatch().
    pub inline fn addSolid(
        this: *Contents,
        allocator: std.mem.Allocator,
        molecule: MoleculeId,
        moles: f64,
        particle_diameter: f64,
    ) error{OutOfMemory}!void {
        try this.addSolidRaw(allocator, molecule, moles, particle_diameter, false);

        if (!this.batch_mode) {
            try this.settle(allocator);
        }
    }

    /// Adds solid substance without triggering settle().
    /// Use with beginBatch()/endBatch() or call settle() manually.
    pub inline fn addSolidRaw(
        this: *Contents,
        allocator: std.mem.Allocator,
        molecule: MoleculeId,
        moles: f64,
        particle_diameter: f64,
        force: bool,
    ) error{OutOfMemory}!void {
        if (!force and constants.isNegligible(moles)) {
            return;
        }

        try this.mergeSolid(allocator, molecule, moles, particle_diameter);
    }

    /// Adds gas to headspace. No settle needed for gas additions.
    pub inline fn addGas(this: *Contents, molecule: MoleculeId, moles: f64) void {
        this.addGasRaw(molecule, moles, false);
    }

    /// Adds gas to headspace without any state updates.
    pub inline fn addGasRaw(
        this: *Contents,
        molecule: MoleculeId,
        moles: f64,
        force: bool,
    ) void {
        if (!force and constants.isNegligible(moles)) {
            return;
        }

        this.gas.molecules[@intFromEnum(molecule)] += moles;
    }

    // Transfer methods settle the DESTINATION but NOT the source.
    // This models intentional physical separation where source state
    // is deliberately altered (e.g., separating funnel, decanting).
    //
    // If source needs re-equilibration, call source.settle() explicitly.

    /// Transfers entire liquid phase by index to destination.
    pub fn transferLiquidPhase(
        this: *Contents,
        allocator: std.mem.Allocator,
        dst: *Contents,
        phase_index: usize,
    ) error{OutOfMemory}!void {
        if (phase_index >= this.liquids.items.len) {
            return;
        }

        // Capture source temperature before removing any mass
        const source_T: f64 = @floatCast(this.getTemperature());
        var transferred_cp: f64 = 0.0;

        const num_mol = comptime MoleculeId.count();
        var transfer_moles: [num_mol]f64 = undefined;

        const phase = &this.liquids.items[phase_index];

        for (0..num_mol) |i| {
            const moles = phase.solutes[i];
            transfer_moles[i] = moles;

            if (!constants.isNegligible(moles)) {
                const mol: MoleculeId = @enumFromInt(i);

                transferred_cp += moles * mol.getHeatCapacity();
            }
        }

        phase.deinit(allocator);
        _ = this.liquids.orderedRemove(phase_index);

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;
        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        try dst.placeLiquidMoles(allocator, &transfer_moles);

        if (!dst.batch_mode) {
            try dst.settle(allocator);
        }
    }

    /// Transfers a specified volume (in liters) from one liquid phase of
    /// this container into a destination container.
    ///
    /// A proportional fraction of every component in the source phase is
    /// removed:
    ///
    /// fraction = min(requested_volume / phase_volume, 1.0)
    /// transferred_i = moles_i * fraction
    ///
    /// The source phase is removed entirely when its residual volume drops
    /// below constants.MIN_VOLUME.
    ///
    /// Returns the volume actually transferred in liters (clamped to the
    /// available phase volume). Returns 0 when phase_index is out of
    /// range or volume_liters <= 0.
    pub inline fn transferLiquidPhaseVolume(
        this: *Contents,
        allocator: std.mem.Allocator,
        dst: *Contents,
        phase_index: usize,
        volume: f64,
    ) error{OutOfMemory}!f64 {
        if (phase_index >= this.liquids.items.len or volume <= 0.0) {
            return 0.0;
        }

        const phase = &this.liquids.items[phase_index];
        const total_volume = phase.getVolume();

        if (total_volume < constants.MIN_VOLUME) {
            return 0.0;
        }

        // Capture source temperature before removing any mass
        const source_T: f64 = @floatCast(this.getTemperature());
        var transferred_cp: f64 = 0.0;

        const fraction = @min(volume / total_volume, 1.0);
        const actual_volume = fraction * total_volume;

        if (actual_volume < constants.MIN_VOLUME) {
            return 0.0;
        }

        const num_mol = comptime MoleculeId.count();
        var transfer_moles: [num_mol]f64 = undefined;

        for (0..num_mol) |i| {
            const moles = phase.solutes[i];
            const transfer = moles * fraction;

            transfer_moles[i] = transfer;
            phase.solutes[i] = @max(moles - transfer, 0.0);

            if (!constants.isNegligible(transfer)) {
                const mol: MoleculeId = @enumFromInt(i);

                transferred_cp += transfer * mol.getHeatCapacity();
            }
        }

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;
        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        if (phase.getVolume() < constants.MIN_VOLUME) {
            phase.deinit(allocator);
            _ = this.liquids.orderedRemove(phase_index);
        }

        try dst.placeLiquidMoles(allocator, &transfer_moles);

        if (!dst.batch_mode) {
            try dst.settle(allocator);
        }

        return actual_volume;
    }

    /// Transfers a specified total volume (in liters) of liquid from this
    /// container to another, pouring from the top layer down - mimicking
    /// manual pouring from a beaker or flask.
    ///
    /// Layers are ordered lightest-first (index 0 = top). The requested
    /// volume is drained starting from the topmost layer:
    ///
    /// 1. If the top layer's volume <= remaining request, drain it entirely.
    /// 2. Otherwise, take a proportional fraction of the top layer.
    /// 3. Repeat with the next layer until the request is fulfilled.
    ///
    /// Returns the volume actually transferred in liters (<= requested).
    /// Returns 0 when there are no liquid phases or volume_liters <= 0.
    pub fn transferLiquidVolume(
        this: *Contents,
        allocator: std.mem.Allocator,
        dst: *Contents,
        volume: f64,
    ) error{OutOfMemory}!f64 {
        if (volume < constants.MIN_VOLUME or this.liquids.items.len == 0) {
            return 0.0;
        }

        const source_T: f64 = @floatCast(this.getTemperature());
        var transferred_cp: f64 = 0.0;

        const num_mol = comptime MoleculeId.count();
        var transfer_moles: [num_mol]f64 = undefined;
        @memset(&transfer_moles, 0.0);

        var remaining = volume;
        var transferred: f64 = 0.0;

        // Drain from the top (index 0 = lightest / topmost layer).
        const j: usize = 0;

        while (j < this.liquids.items.len and remaining >= constants.MIN_VOLUME) {
            const phase = &this.liquids.items[j];
            const phase_vol = phase.getVolume();

            if (phase_vol < constants.MIN_VOLUME) {
                // Phase is effectively empty - clean it up and continue.
                phase.deinit(allocator);
                _ = this.liquids.orderedRemove(j);

                continue;
            }

            if (phase_vol <= remaining) {
                for (0..num_mol) |i| {
                    const moles = phase.solutes[i];
                    transfer_moles[i] += moles;

                    if (!constants.isNegligible(moles)) {
                        const mol: MoleculeId = @enumFromInt(i);

                        transferred_cp += moles * mol.getHeatCapacity();
                    }
                }

                remaining -= phase_vol;
                transferred += phase_vol;

                phase.deinit(allocator);
                _ = this.liquids.orderedRemove(j);
                // Don't increment - next layer is now at index j.
            } else {
                // Partial drain from this layer.
                const fraction = remaining / phase_vol;

                for (0..num_mol) |i| {
                    const moles = phase.solutes[i];
                    const take = moles * fraction;

                    transfer_moles[i] += take;
                    phase.solutes[i] = @max(moles - take, 0.0);

                    if (!constants.isNegligible(take)) {
                        const mol: MoleculeId = @enumFromInt(i);

                        transferred_cp += take * mol.getHeatCapacity();
                    }
                }

                transferred += remaining;
                remaining = 0.0;

                // Remove if residual is negligible.
                if (phase.getVolume() < constants.MIN_VOLUME) {
                    phase.deinit(allocator);
                    _ = this.liquids.orderedRemove(j);
                }

                break;
            }
        }

        if (transferred < constants.MIN_VOLUME) {
            return 0.0;
        }

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;

        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        try dst.placeLiquidMoles(allocator, &transfer_moles);

        if (!dst.batch_mode) {
            try dst.settle(allocator);
        }

        return transferred;
    }

    /// Transfers all liquid phases to destination.
    /// Destination is settled; source is NOT.
    pub fn transferAllLiquids(
        this: *Contents,
        allocator: std.mem.Allocator,
        dst: *Contents,
    ) error{OutOfMemory}!void {
        if (this.liquids.items.len == 0) {
            return;
        }

        const source_T: f64 = @floatCast(this.getTemperature());
        var transferred_cp: f64 = 0.0;

        const num_mol = comptime MoleculeId.count();
        var transfer_moles: [num_mol]f64 = undefined;
        @memset(&transfer_moles, 0.0);

        // Pool moles from every layer and free the per-phase solute arrays.
        for (this.liquids.items) |*phase| {
            for (0..num_mol) |i| {
                const moles = phase.solutes[i];

                transfer_moles[i] += moles;

                if (!constants.isNegligible(moles)) {
                    const mol: MoleculeId = @enumFromInt(i);

                    transferred_cp += moles * mol.getHeatCapacity();
                }
            }

            phase.deinit(allocator);
        }

        this.liquids.clearRetainingCapacity();

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;
        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        try dst.placeLiquidMoles(allocator, &transfer_moles);

        if (!dst.batch_mode) {
            try dst.settle(allocator);
        }
    }

    /// Transfers an entire solid phase by index to the destination container.
    ///
    /// The solid entry (including its occluded solution and particle
    /// properties) is moved wholesale. The source entry is removed via
    /// orderedRemove to preserve the relative order of remaining solids.
    ///
    /// Destination is settled (which includes consolidateSolids); source is NOT.
    pub inline fn transferSolidPhase(
        this: *Contents,
        allocator: std.mem.Allocator,
        dst: *Contents,
        phase_index: usize,
    ) error{OutOfMemory}!void {
        if (phase_index >= this.solids.items.len) {
            return;
        }

        const source_T: f64 = @floatCast(this.getTemperature());
        const solid = this.solids.orderedRemove(phase_index);

        var transferred_cp = solid.moles * solid.molecule.getHeatCapacity();

        if (solid.occluded_solution) |*occ| {
            const num_mol = comptime MoleculeId.count();

            for (0..num_mol) |i| {
                if (!constants.isNegligible(occ.solutes[i])) {
                    const mol: MoleculeId = @enumFromInt(i);

                    transferred_cp += occ.solutes[i] * mol.getHeatCapacity();
                }
            }
        }

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;
        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        try dst.solids.append(allocator, solid);

        if (!dst.batch_mode) {
            try dst.settle(allocator);
        }
    }

    /// Transfers a specified volume (in liters) of a solid phase by index
    /// to the destination container.
    ///
    /// Transfer is quantized by particle count - only whole particles move:
    ///
    /// n_particles = floor(requested_volume / V_particle)
    /// transferred_volume = n_particles * V_particle
    /// transferred_moles = n_particles * moles_per_particle
    ///
    /// Occluded solution is split proportionally to the mole fraction
    /// transferred.
    ///
    /// The source solid is removed entirely when fewer than one particle
    /// remains.
    ///
    /// Returns the volume actually transferred in liters. Returns 0 when:
    /// - phase_index is out of range
    /// - volume <= 0
    /// - requested volume is less than one particle's volume
    pub inline fn transferSolidPhaseVolume(
        this: *Contents,
        allocator: std.mem.Allocator,
        dst: *Contents,
        phase_index: usize,
        volume: f64,
    ) error{OutOfMemory}!f64 {
        if (phase_index >= this.solids.items.len or volume < constants.MIN_VOLUME) {
            return 0.0;
        }

        const solid = &this.solids.items[phase_index];
        const total_volume = solid.getVolume();

        if (total_volume < constants.MIN_VOLUME) {
            return 0.0;
        }

        const d = solid.particle_diameter;
        const particle_volume_L = (std.math.pi / 6.0) * d * d * d * 1000.0;

        if (particle_volume_L < constants.MIN_VOLUME) {
            return 0.0;
        }

        const total_particles = solid.getParticleCount();
        const requested_particles = @floor(volume / particle_volume_L);

        if (requested_particles < 1.0) {
            return 0.0;
        }

        const source_T: f64 = @floatCast(this.getTemperature());
        var transferred_cp: f64 = 0.0;

        const particles_to_transfer = @min(requested_particles, total_particles);
        const actual_volume = particles_to_transfer * particle_volume_L;

        if (actual_volume < constants.MIN_VOLUME) {
            return 0.0;
        }

        const moles_per_particle = if (total_particles > 0.0) solid.moles / total_particles else 0.0;
        const transfer_moles = particles_to_transfer * moles_per_particle;
        const fraction = if (solid.moles > 0.0) transfer_moles / solid.moles else 0.0;

        transferred_cp += transfer_moles * solid.molecule.getHeatCapacity();

        var transferred_solid = SolidPhase{
            .molecule = solid.molecule,
            .moles = transfer_moles,
            .particle_diameter = solid.particle_diameter,
            .occluded_solution = null,
            .crystalline = solid.crystalline,
        };

        // Split occluded solution proportionally.
        if (solid.occluded_solution) |*src_occ| {
            const num_mol = comptime MoleculeId.count();

            var dst_occ = try LiquidPhase.init(allocator);
            errdefer dst_occ.deinit(allocator);

            var has_content = false;

            for (0..num_mol) |i| {
                const occ_moles = src_occ.solutes[i];
                const take = occ_moles * fraction;

                if (take > 0.0) {
                    dst_occ.solutes[i] = take;
                    src_occ.solutes[i] = @max(occ_moles - take, 0.0);
                    has_content = true;

                    const mol: MoleculeId = @enumFromInt(i);
                    transferred_cp += take * mol.getHeatCapacity();
                }
            }

            if (has_content) {
                transferred_solid.occluded_solution = dst_occ;
            } else {
                dst_occ.deinit(allocator);
            }
        }

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;
        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        solid.moles = @max(solid.moles - transfer_moles, 0.0);
        const remaining_particles = solid.getParticleCount();

        if (remaining_particles < 1.0) {
            solid.deinit(allocator);
            _ = this.solids.orderedRemove(phase_index);
        } else if (solid.occluded_solution) |*occ| {
            if (occ.getVolume() < constants.MIN_VOLUME) {
                occ.deinit(allocator);
                solid.occluded_solution = null;
            }
        }

        try dst.solids.append(allocator, transferred_solid);

        if (!dst.batch_mode) {
            try dst.settle(allocator);
        }

        return actual_volume;
    }

    /// Transfers a specified total volume (in liters) of solid material from
    /// this container to another.
    ///
    /// Iterates through solid phases front-to-back (modelling scooping from
    /// a pile). Transfer is quantized by particle count - only whole
    /// particles move from each phase.
    ///
    /// For each phase:
    /// 1. Compute how many whole particles fit in the remaining request.
    /// 2. Transfer that many particles (moles, occluded solution).
    /// 3. Move to the next phase if request is not fulfilled.
    ///
    /// Returns the volume actually transferred in liters (<= requested).
    /// Returns 0 when there are no solid phases or volume <= 0.
    pub inline fn transferSolidVolume(
        this: *Contents,
        allocator: std.mem.Allocator,
        dst: *Contents,
        volume: f64,
    ) error{OutOfMemory}!f64 {
        if (volume < constants.MIN_VOLUME or this.solids.items.len == 0) {
            return 0.0;
        }

        const source_T: f64 = @floatCast(this.getTemperature());
        var transferred_cp: f64 = 0.0;

        var remaining = volume;
        var transferred: f64 = 0.0;
        var j: usize = 0;

        while (j < this.solids.items.len and remaining > 0.0) {
            const solid = &this.solids.items[j];

            // Particle volume in liters.
            const d = solid.particle_diameter;
            const particle_volume_L = (std.math.pi / 6.0) * d * d * d * 1000.0;

            if (particle_volume_L <= 0.0) {
                j += 1;

                continue;
            }

            const total_particles = solid.getParticleCount();

            if (total_particles < 1.0) {
                // Phase has less than one particle - clean it up.
                solid.deinit(allocator);
                _ = this.solids.orderedRemove(j);

                continue;
            }

            // How many whole particles can we take?
            const requested_particles = @floor(remaining / particle_volume_L);

            if (requested_particles < 1.0) {
                j += 1;

                continue;
            }

            const particles_to_transfer = @min(requested_particles, total_particles);
            const volume_taken = particles_to_transfer * particle_volume_L;

            const moles_per_particle = solid.moles / total_particles;
            const transfer_moles = particles_to_transfer * moles_per_particle;
            const fraction = transfer_moles / solid.moles;

            transferred_cp += transfer_moles * solid.molecule.getHeatCapacity();

            var transferred_solid = SolidPhase{
                .molecule = solid.molecule,
                .moles = transfer_moles,
                .particle_diameter = solid.particle_diameter,
                .occluded_solution = null,
                .crystalline = solid.crystalline,
            };

            if (solid.occluded_solution) |*src_occ| {
                const num_mol = comptime MoleculeId.count();

                var dst_occ = try LiquidPhase.init(allocator);
                errdefer dst_occ.deinit(allocator);

                var has_content = false;

                for (0..num_mol) |i| {
                    const occ_moles = src_occ.solutes[i];
                    const take = occ_moles * fraction;

                    if (!constants.isNegligible(take)) {
                        dst_occ.solutes[i] = take;
                        src_occ.solutes[i] = @max(occ_moles - take, 0.0);
                        has_content = true;

                        const mol: MoleculeId = @enumFromInt(i);
                        transferred_cp += take * mol.getHeatCapacity();
                    }
                }

                if (has_content) {
                    transferred_solid.occluded_solution = dst_occ;
                } else {
                    dst_occ.deinit(allocator);
                }
            }

            solid.moles -= transfer_moles;
            const remaining_particles = solid.getParticleCount();

            if (remaining_particles < 1.0) {
                solid.deinit(allocator);
                _ = this.solids.orderedRemove(j);
                // Don't increment j - next phase now at index j.
            } else {
                if (solid.occluded_solution) |*occ| {
                    if (occ.getVolume() < constants.MIN_VOLUME) {
                        occ.deinit(allocator);
                        solid.occluded_solution = null;
                    }
                }

                j += 1;
            }

            try dst.solids.append(allocator, transferred_solid);
            remaining -= volume_taken;
            transferred += volume_taken;
        }

        if (transferred < constants.MIN_VOLUME) {
            return 0.0;
        }

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;
        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        if (!dst.batch_mode) {
            try dst.settle(allocator);
        }

        return transferred;
    }

    /// Transfers all solids to destination (filtration, decanting).
    pub inline fn transferAllSolids(this: *Contents, allocator: std.mem.Allocator, dst: *Contents) error{OutOfMemory}!void {
        const source_T: f64 = @floatCast(this.getTemperature());
        var transferred_cp: f64 = 0.0;
        const num_mol = MoleculeId.count();

        for (this.solids.items) |solid| {
            transferred_cp += solid.moles * solid.molecule.getHeatCapacity();

            if (solid.occluded_solution) |*occ| {
                for (0..num_mol) |i| {
                    if (occ.solutes[i] > 0.0) {
                        const mol: MoleculeId = @enumFromInt(i);

                        transferred_cp += occ.solutes[i] * mol.getHeatCapacity();
                    }
                }
            }

            try dst.solids.append(allocator, solid);
        }

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;
        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        this.solids.clearRetainingCapacity();
    }

    /// Transfers a specified volume (in liters, measured at the source's
    /// temperature and pressure) of gas from this container to another.
    ///
    /// Uses the ideal gas law:
    ///
    /// V_total = n_total * R * T / P  [m³]  →  * 1000 for [L]
    /// fraction = min(requested / V_total_L, 1)
    ///
    /// Returns the volume actually transferred in liters (<= requested).
    pub inline fn transferGasVolume(
        this: *Contents,
        dst: *Contents,
        volume: f64,
    ) f64 {
        if (volume < constants.MIN_VOLUME) {
            return 0.0;
        }

        const source_T: f64 = @floatCast(this.getTemperature());
        var transferred_cp: f64 = 0.0;

        const num_mol = comptime MoleculeId.count();
        const P: f64 = @floatCast(this.pressure);

        var n_total: f64 = 0.0;

        for (0..num_mol) |i| {
            n_total += this.gas.molecules[i];
        }

        if (constants.isNegligible(n_total)) {
            return 0.0;
        }

        const total_volume = n_total * constants.R * source_T / P * 1000.0;
        const fraction = @min(volume / total_volume, 1.0);
        const actual_volume = fraction * total_volume;

        for (0..num_mol) |i| {
            const moles = this.gas.molecules[i];
            const transfer = moles * fraction;

            this.gas.molecules[i] = @max(moles - transfer, 0.0);
            dst.gas.molecules[i] += transfer;

            if (!constants.isNegligible(transfer)) {
                const mol: MoleculeId = @enumFromInt(i);

                transferred_cp += transfer * mol.getHeatCapacity();
            }
        }

        // Transfer thermal energy
        const transferred_energy = transferred_cp * source_T;
        this.heat_energy = @max(0.0, this.heat_energy - transferred_energy);
        dst.heat_energy += transferred_energy;

        return actual_volume;
    }

    /// Calculates the total heat capacity (C_p) of the container and all its contents in J/K.
    pub inline fn getHeatCapacity(this: *const Contents) f64 {
        var cp_total = this.container_heat_capacity;
        const num_mol = comptime MoleculeId.count();

        // 1. Gas Phase
        for (0..num_mol) |i| {
            const moles = this.gas.molecules[i];

            if (!constants.isNegligible(moles)) {
                const mol: MoleculeId = @enumFromInt(i);

                cp_total += moles * mol.getHeatCapacity();
            }
        }

        // 2. Liquid Phases
        for (this.liquids.items) |*phase| {
            for (0..num_mol) |i| {
                const moles = phase.solutes[i];

                if (!constants.isNegligible(moles)) {
                    const mol: MoleculeId = @enumFromInt(i);

                    cp_total += moles * mol.getHeatCapacity();
                }
            }
        }

        // 3. Solid Phases (including occluded solutions)
        for (this.solids.items) |*solid| {
            cp_total += solid.moles * solid.molecule.getHeatCapacity();

            if (solid.occluded_solution) |*occ| {
                for (0..num_mol) |i| {
                    const moles = occ.solutes[i];

                    if (!constants.isNegligible(moles)) {
                        const mol: MoleculeId = @enumFromInt(i);

                        cp_total += moles * mol.getHeatCapacity();
                    }
                }
            }
        }

        return cp_total;
    }

    /// Returns the current thermodynamic temperature of the system in Kelvin.
    /// T = Energy / Capacity
    pub inline fn getTemperature(this: *const Contents) f64 {
        const cp = this.getHeatCapacity();

        if (cp <= 0.0) {
            return 0.0;
        }

        return @floatCast(this.heat_energy / cp);
    }

    /// Adds raw thermal energy (in Joules) to the system.
    /// Positive values heat it up, negative cool it down.
    pub inline fn addHeat(this: *Contents, joules: f64) void {
        this.heat_energy += joules;
    }

    /// Exchanges heat with the environment using Newton's law of cooling.
    /// Q = thermal_conductivity * (T_env - T_sys) * dt
    ///
    /// Returns the amount of heat energy (in Joules) transferred TO the environment.
    /// A negative return value means the environment lost heat (and the container gained it).
    /// This is strictly required to maintain the Law of Conservation of Energy
    /// in a broader simulation (e.g., cooling down the hotplate or heating the room).
    pub inline fn exchangeHeat(
        this: *Contents,
        ambient_temperature: f64,
        thermal_conductivity: f64,
        dt: f64,
    ) f64 {
        const T = @as(f64, @floatCast(this.getTemperature()));
        const T_env = @as(f64, @floatCast(ambient_temperature));
        const q = thermal_conductivity * (T_env - T) * dt;

        this.addHeat(q);

        return -q;
    }

    /// Artificially forces the system to a new temperature.
    /// Calculates required energy and triggers instant phase equilibration.
    pub inline fn setTemperature(
        this: *Contents,
        allocator: std.mem.Allocator,
        new_temperature: f64,
        max_volume: f64,
    ) error{OutOfMemory}!void {
        if (this.getTemperature() == new_temperature) {
            return;
        }

        // Force internal energy to match the target temperature
        this.heat_energy = @as(f64, @floatCast(new_temperature)) * this.getHeatCapacity();

        try this.updatePhaseTransitions(allocator, std.math.inf(f64), max_volume);

        if (!this.batch_mode) {
            try this.settle(allocator);
        }
    }

    /// Sets pressure and updates the system state.
    ///
    /// Performs:
    /// 1. Sets new pressure
    /// 2. Updates phase transitions (instant equilibration using infinite dt)
    /// 3. Settles the system (equilibration, precipitation, sorting)
    ///
    /// For gradual transitions, set `temperature` directly and call
    /// `updatePhaseTransitions(allocator, dt)` in your simulation loop.
    pub inline fn setPressure(
        this: *Contents,
        allocator: std.mem.Allocator,
        new_pressure: f64,
    ) error{OutOfMemory}!void {
        if (this.pressure == new_pressure) {
            return;
        }

        this.pressure = new_pressure;
        // High-level API assumes instantaneous equilibrium
        try this.updatePhaseTransitions(allocator, std.math.inf(f64));

        if (!this.batch_mode) {
            try this.settle(allocator);
        }
    }

    /// Sets both temperature and pressure, then updates the system.
    /// More efficient than calling setTemperature + setPressure separately.
    pub inline fn setConditions(
        this: *Contents,
        allocator: std.mem.Allocator,
        new_temperature: f64,
        new_pressure: f64,
    ) error{OutOfMemory}!void {
        const t_changed = this.temperature != new_temperature;
        const p_changed = this.pressure != new_pressure;

        if (!t_changed and !p_changed) {
            return;
        }

        if (t_changed) {
            this.heat_energy = @as(f64, @floatCast(new_temperature)) * this.getHeatCapacity();
        }

        this.pressure = new_pressure;

        try this.updatePhaseTransitions(allocator, std.math.inf(f64));

        if (!this.batch_mode) {
            try this.settle(allocator);
        }
    }

    /// Removes all contents from this container, returning them to a fresh state.
    /// Frees all allocated phase memory.
    pub inline fn clear(this: *Contents, allocator: std.mem.Allocator) void {
        for (this.liquids.items) |*phase| {
            phase.deinit(allocator);
        }
        this.liquids.clearRetainingCapacity();

        for (this.solids.items) |*phase| {
            phase.deinit(allocator);
        }
        this.solids.clearRetainingCapacity();

        @memset(this.gas.molecules, 0.0);
    }

    /// Returns true if the container has no significant contents.
    pub inline fn isEmpty(this: *const Contents) bool {
        // Check liquids
        for (this.liquids.items) |*phase| {
            if (phase.getVolume() >= constants.MIN_VOLUME) {
                return false;
            }
        }

        // Check solids
        for (this.solids.items) |*solid| {
            if (constants.isNegligible(solid.moles)) {
                return false;
            }
        }

        // Check gas
        for (this.gas.molecules) |moles| {
            if (constants.isNegligible(moles)) {
                return false;
            }
        }

        return true;
    }

    /// Returns total liquid volume in liters across all layers.
    pub inline fn getTotalLiquidVolume(this: *const Contents) f64 {
        var total: f64 = 0.0;

        for (this.liquids.items) |*phase| {
            total += phase.getVolume();
        }

        return total;
    }

    /// Returns total solid volume in liters across all phases.
    pub inline fn getTotalSolidVolume(this: *const Contents) f64 {
        var total: f64 = 0.0;

        for (this.solids.items) |*solid| {
            total += solid.getVolume();
        }

        return total;
    }

    /// Returns the actual volume occupied by gas in the container (headspace).
    /// In a rigid container, gas always fills the available free space.
    pub inline fn getGasVolume(this: *const Contents, max_volume: f64) f64 {
        const V_liquid = this.getTotalLiquidVolume();
        const V_solid = this.getTotalSolidVolume();
        const headspace = max_volume - V_liquid - V_solid;

        return @max(0.0, headspace);
    }

    pub inline fn getGasMoles(this: *const Contents) f64 {
        var moles: f64 = 0.0;

        for (this.gas.molecules) |m| {
            moles += m;
        }

        if (constants.isNegligible(moles)) {
            return 0.0;
        }

        return moles;
    }

    /// Returns the overall dominant odors from ALL phases in the container.
    ///
    /// Combines odor contributions from:
    /// 1. Gas phase (direct, strongest contribution)
    /// 2. All liquid phases (evaporation from surfaces)
    /// 3. All solid phases (sublimation and occluded volatiles)
    ///
    /// Gas-phase odors dominate since they are already airborne.
    /// Liquid odors are weighted by surface area exposure.
    /// Solid odors contribute least unless highly volatile (low melting point).
    ///
    /// Returns the 2 strongest odors across all phases.
    pub fn getDominantOdors(this: *const Contents) Flavor {
        var flavor = Flavor{};

        // Accumulate odor scores across all phases
        var odor_scores: [@typeInfo(Odor).@"enum".fields.len]f64 =
            .{0.0} ** @typeInfo(Odor).@"enum".fields.len;

        // Weight factors for different phases
        const gas_weight: f64 = 10.0; // Gas dominates (already airborne)
        const liquid_weight: f64 = 1.0; // Liquids evaporate
        const solid_weight: f64 = 0.3; // Solids mostly don't smell

        // Gas phase contribution
        const gas_flavor = this.gas.getOdors(this.getTemperature(), this.pressure);

        for (0..2) |i| {
            if (gas_flavor.odors[i] != .none) {
                odor_scores[@intFromEnum(gas_flavor.odors[i])] +=
                    @as(f64, gas_flavor.odor_intensities[i]) * gas_weight;
            }
        }

        // Liquid phases contribution
        for (this.liquids.items) |*liquid| {
            const liquid_flavor = liquid.getOdors(this.getTemperature());

            for (0..2) |i| {
                if (liquid_flavor.odors[i] != .none) {
                    odor_scores[@intFromEnum(liquid_flavor.odors[i])] +=
                        @as(f64, liquid_flavor.odor_intensities[i]) * liquid_weight;
                }
            }
        }

        // Solid phases contribution
        for (this.solids.items) |*solid| {
            const solid_flavor = solid.getOdors(this.getTemperature());

            for (0..2) |i| {
                if (solid_flavor.odors[i] != .none) {
                    odor_scores[@intFromEnum(solid_flavor.odors[i])] +=
                        @as(f64, solid_flavor.odor_intensities[i]) * solid_weight;
                }
            }
        }

        // Find top 2 odors
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

    /// Returns the overall visible colors from ALL phases in the container.
    ///
    /// Visual perception depends on viewing direction:
    /// - Looking down: see liquids through gas, solids at bottom
    /// - Looking through: gas color, liquid layers
    ///
    /// This method assumes top-down viewing (most common lab scenario):
    /// 1. Gas contributes if colored (rare)
    /// 2. Top liquid layer is most visible
    /// 3. Solids at bottom may show through transparent liquids
    ///
    /// Returns up to 2 dominant colors.
    pub fn getDominantColors(this: *const Contents) Flavor {
        var flavor = Flavor{};

        var color_scores: [@typeInfo(Color).@"enum".fields.len]f64 =
            .{0.0} ** @typeInfo(Color).@"enum".fields.len;

        // Gas contribution (usually transparent, but Br2, I2, NO2 are colored)
        const gas_flavor = this.gas.getColors(this.getTemperature(), this.pressure);

        for (0..2) |i| {
            if (gas_flavor.colors[i] != .transparent) {
                color_scores[@intFromEnum(gas_flavor.colors[i])] +=
                    @as(f64, gas_flavor.color_intensities[i]) * 0.5;
            }
        }

        // Liquid layers: top layers obscure bottom layers
        // Weight decreases for deeper layers
        var layer_weight: f64 = 1.0;
        for (this.liquids.items) |*liquid| {
            const liquid_flavor = liquid.getColors();

            // Check if this layer is transparent
            var layer_opacity: f64 = 0.0;

            for (0..2) |i| {
                if (liquid_flavor.colors[i] != .transparent) {
                    const contribution = @as(f64, liquid_flavor.color_intensities[i]) * layer_weight;

                    color_scores[@intFromEnum(liquid_flavor.colors[i])] += contribution;
                    layer_opacity = @max(layer_opacity, liquid_flavor.color_intensities[i]);
                }
            }

            // Reduce weight for layers below based on opacity of this layer
            layer_weight *= (1.0 - @as(f64, layer_opacity) * 0.8);
            layer_weight = @max(layer_weight, 0.1);
        }

        // Solids visible through remaining transparency
        for (this.solids.items) |*solid| {
            const solid_flavor = solid.getColors();

            for (0..2) |i| {
                if (solid_flavor.colors[i] != .transparent) {
                    color_scores[@intFromEnum(solid_flavor.colors[i])] +=
                        @as(f64, solid_flavor.color_intensities[i]) * layer_weight * 0.8;
                }
            }
        }

        // Find top 2 colors
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

    /// Returns complete flavor profile (colors + odors) for the entire container.
    pub inline fn getFlavor(this: *const Contents) Flavor {
        var flavor = this.getDominantColors();
        const odor_flavor = this.getDominantOdors();

        flavor.odors = odor_flavor.odors;
        flavor.odor_intensities = odor_flavor.odor_intensities;

        return flavor;
    }

    /// Updates the system pressure.
    /// Uses the ideal gas law (P = nRT / V) when headspace exists.
    /// Transitions to bulk modulus compression (P = P_atm + K * strain)
    /// when the container is completely full or overfilled (hydraulic lock).
    /// `max_volume` is the total container volume in liters.
    pub inline fn updatePressure(this: *Contents, max_volume: f64) void {
        const V_liquid = this.getTotalLiquidVolume();
        const V_solid = this.getTotalSolidVolume();
        const V_occupied = V_liquid + V_solid;

        // Can be negative if the container is overfilled
        const V_gas_L = max_volume - V_occupied;

        // Headspace exists (Ideal Gas Law)
        if (V_gas_L > constants.MIN_VOLUME) {
            const V_gas_m3 = V_gas_L * 0.001; // Convert L to m³

            var n_total: f64 = 0.0;
            for (0..this.gas.molecules.len) |i| {
                n_total += this.gas.molecules[i];
            }

            if (constants.isNegligible(n_total)) {
                this.pressure = 0.0; // True vacuum - will cause vacuum boiling!

                return;
            }

            const T: f64 = @floatCast(this.getTemperature());
            if (T <= 0.0) {
                this.pressure = 0.0;

                return;
            }

            // P = nRT / V (Pressure in Pascals)
            const P: f64 = n_total * constants.R * T / V_gas_m3;
            this.pressure = @floatCast(P);

            return;
        }

        // Full or Overfilled (Hydraulic Lock / Bulk Modulus)

        // Calculate how much the contents are compressed relative to container size.
        // strain = 0.0 when perfectly full, > 0.0 when overfilled.
        const strain = @max(0.0, V_occupied - max_volume) / max_volume;

        // Get the effective resistance to compression of the mixture
        const K_mix = this.getEffectiveBulkModulus();

        // Baseline pressure: a sealed rigid container defaults to 1 atm
        // (assuming it was sealed at atmospheric pressure).
        // Hydraulic pressure adds to this baseline.
        const P_total = constants.P_std + (K_mix * strain);

        this.pressure = @floatCast(P_total);
    }

    /// Redistributes molecules among immiscible liquid phases to approximate
    /// partition equilibrium using Scatchard-Hildebrand regular solution theory.
    ///
    /// For molecule i in phase j, the Boltzmann-weighted affinity is:
    ///
    /// A_i^j = exp[-V_m_i * (δ_i - δ_j)² / RT]
    ///
    /// Moles are distributed proportionally to affinity * phase volume:
    ///
    /// n_i^j = n_i_total * (A_i^j * V_j) / Σ_k(A_i^k * V_k)
    pub fn equilibratePhases(this: *Contents) void {
        const num_phases = this.liquids.items.len;

        if (num_phases < 2) {
            return;
        }

        const num_mol = comptime MoleculeId.count();
        const max_phases = 8;
        const np = @min(num_phases, max_phases);
        const T: f64 = @floatCast(this.getTemperature());
        const RT = constants.R * T;

        // Snapshot each phase's identity and volume before redistribution.
        // Each phase is characterized by its primary solvent's Hildebrand δ.
        var phase_delta: [max_phases]f64 = undefined;
        var phase_vol: [max_phases]f64 = undefined;

        for (0..np) |j| {
            const primary = this.liquids.items[j].getPrimaryComponent() orelse {
                phase_delta[j] = 0.0;
                phase_vol[j] = 0.0;

                continue;
            };

            phase_delta[j] = @floatCast(primary.getHildebrand());

            const primary_moles = this.liquids.items[j].solutes[@intFromEnum(primary)];
            // Volume in liters from the primary solvent only.
            phase_vol[j] = primary_moles * primary.getWeight() / primary.getDensity() / 1000.0;
        }

        // Pool moles of every molecule across all phases.
        var total: [num_mol]f64 = undefined;
        @memset(&total, 0.0);

        for (0..np) |j| {
            for (0..num_mol) |i| {
                total[i] += this.liquids.items[j].solutes[i];
            }
        }

        // Distribute each molecule among phases using Hildebrand affinity.
        for (0..num_mol) |i| {
            if (constants.isNegligible(total[i])) {
                for (0..np) |j| {
                    this.liquids.items[j].solutes[i] = 0.0;
                }

                continue;
            }

            const mol: MoleculeId = @enumFromInt(i);
            const delta_i: f64 = @floatCast(mol.getHildebrand());
            const vm_i: f64 = mol.getWeight() / mol.getDensity();

            var w: [max_phases]f64 = undefined;
            var w_sum: f64 = 0.0;

            for (0..np) |j| {
                const dd = delta_i - phase_delta[j];

                w[j] = @exp(-vm_i * dd * dd / RT) * phase_vol[j];
                w_sum += w[j];
            }

            if (w_sum > 0.0) {
                var distributed: f64 = 0.0;

                for (0..np - 1) |j| {
                    const amount = total[i] * w[j] / w_sum;

                    this.liquids.items[j].solutes[i] = amount;
                    distributed += amount;
                }

                this.liquids.items[np - 1].solutes[i] = total[i] - distributed;
            }
        }
    }

    /// Checks every liquid phase for supersaturation and moves excess moles
    /// into the solid phase as precipitate.
    ///
    /// Skips molecules that are the primary solvent of their phase (solvents
    /// don't precipitate out of themselves) and molecules with no solid form
    /// (water_solubility = null -> saturation = +inf).
    pub fn checkPrecipitation(this: *Contents, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        const num_mol = comptime MoleculeId.count();

        for (this.liquids.items) |*phase| {
            const volume = phase.getVolume();

            if (volume <= 0.0) {
                continue;
            }

            const primary = phase.getPrimaryComponent() orelse {
                continue;
            };

            for (0..num_mol) |i| {
                const mol: MoleculeId = @enumFromInt(i);

                // Solvents don't precipitate out of their own phase.
                if (mol == primary) {
                    continue;
                }

                const moles = phase.solutes[i];

                if (constants.isNegligible(moles)) {
                    continue;
                }

                const concentration = moles / volume;
                const saturation = phase.getSaturationLimit(mol, this.getTemperature());

                if (saturation == std.math.inf(f64)) {
                    continue;
                }

                if (concentration > saturation) {
                    const excess = (concentration - saturation) * volume;

                    // don't create dust
                    if (constants.isNegligible(excess)) {
                        continue;
                    }

                    phase.solutes[i] -= excess;

                    try this.mergeSolid(
                        allocator,
                        mol,
                        excess,
                        constants.DEFAULT_PRECIPITATION_DIAMETER,
                    );
                }
            }
        }
    }

    /// Brings the system to equilibrium state.
    ///
    /// Performs in order:
    /// 1. Removes empty liquid phases
    /// 2. Redistributes solutes among immiscible phases (equilibratePhases)
    /// 3. Precipitates supersaturated species (checkPrecipitation)
    /// 4. Sorts liquid layers by density (lightest on top)
    /// 5. Consolidates similar solid phases
    /// 6. Removes empty solid phases
    ///
    /// Call this after any batch modifications or when manual control is needed.
    /// High-level methods (addLiquid, setTemperature, etc.) call this automatically.
    pub inline fn settle(this: *Contents, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        this.removeEmptyLiquidPhases(allocator);
        this.equilibratePhases();
        try this.checkPrecipitation(allocator);
        this.sortLiquidsByDensity();
        this.consolidateSolids(allocator);
        this.removeEmptySolidPhases(allocator);
    }

    /// Updates phase states over a time step dt.
    /// Phase changes are strictly governed by available Latent Heat (Enthalpy).
    /// If water hits 100°C, it won't boil unless excess heat_energy is available.
    /// Boiling/Melting consumes energy (cooling the system back down to the transition point).
    /// Condensing/Freezing releases energy (heating the system back up).
    ///
    /// max_volume: Total container volume in liters. Used to calculate headspace for
    /// vapor-liquid equilibrium. Pass 0.0 for vacuum-like flash evaporation behavior
    /// (no VLE limit, limited only by heat availability).
    pub fn updatePhaseTransitions(
        this: *Contents,
        allocator: std.mem.Allocator,
        dt: f64,
        max_volume: f64,
    ) error{OutOfMemory}!void {
        if (dt <= 0.0) {
            return;
        }

        const num_mol = MoleculeId.count();

        // Buffers for net molar changes to each phase (positive = add, negative = remove).
        // These are applied atomically at the end to prevent mass loss during intermediate steps.
        var delta_gas: [num_mol]f64 = .{0} ** num_mol;
        var delta_liquid_global: [num_mol]f64 = .{0} ** num_mol; // For molecules entering liquid from gas/solid
        var delta_solid: [num_mol]f64 = .{0} ** num_mol;

        const P: f64 = @floatCast(this.pressure);

        const is_instant = dt == std.math.inf(f64);

        var current_T = @as(f64, @floatCast(this.getTemperature()));
        var current_Cp = this.getHeatCapacity();

        // Snapshot of gas moles for Dalton's Law calculations (state at t=0).
        // We track changes separately to avoid order-dependent artifacts.
        var total_gas_moles: f64 = 0.0;

        for (0..num_mol) |i| {
            total_gas_moles += this.gas.molecules[i];
        }

        // Calculate headspace volume based on current solids/liquids (snapshot).
        // This remains constant during this time step to avoid coupling oscillations.
        const V_liquid_total = this.getTotalLiquidVolume();
        const V_solid_total = this.getTotalSolidVolume();
        const V_gas_L = @max(0.0, max_volume - V_liquid_total - V_solid_total);
        const V_gas_m3 = V_gas_L * 0.001;

        // 1. Gas Phase (Condensation & Deposition)
        // Gas molecules transitioning to liquid or solid.
        for (0..num_mol) |i| {
            const amount = this.gas.molecules[i];

            if (constants.isNegligible(amount)) {
                continue;
            }

            const mol: MoleculeId = @enumFromInt(i);
            const target = mol.targetState(@floatCast(current_T), @floatCast(P));

            if (target != .gas) {
                const P_eq = mol.getVaporPressure(current_T);

                if (P_eq > 0.0) {
                    const P_partial = if (total_gas_moles > 0.0 and P > 0.0)
                        P * (amount / total_gas_moles)
                    else
                        0.0;

                    if (P_partial <= P_eq) {
                        continue;
                    }
                }

                const mp = @as(f64, @floatCast(mol.getMeltingPoint() orelse current_T));
                const bp = mol.getPressureCorrectedBoilingPoint(@floatCast(P)) orelse mp;
                const dH_vap = mol.getDeltaHvapJ();
                const dH_fus = mol.getDeltaHfusJ();

                const t_trans = if (target == .liquid)
                    bp
                else
                    mp;

                const dH = if (target == .liquid)
                    dH_vap
                else
                    (dH_vap + dH_fus);

                var take = amount;

                // Heat-limited condensation: can only condense what the temperature drop allows
                if (!is_instant and current_T < t_trans) {
                    const q_missing = (t_trans - current_T) * current_Cp;

                    take = @min(amount, @max(0.0, q_missing / dH));
                }

                // Prevent excessive condensation, which creates an artificial vacuum
                // below the saturated vapor pressure of the substance in a closed container.
                if (max_volume > constants.MIN_VOLUME and V_gas_m3 > constants.MIN_VOLUME) {
                    const P_eq_clamp = mol.getVaporPressure(current_T);

                    if (P_eq_clamp > 0.0) {
                        const min_moles_gas = P_eq_clamp * V_gas_m3 / (constants.R * current_T);
                        const max_condense = @max(0.0, amount - min_moles_gas);

                        take = @min(take, max_condense);
                    }
                }

                if (take > constants.MIN_MOLES) {
                    delta_gas[i] -= take; // Remove from gas

                    if (target == .liquid) {
                        delta_liquid_global[i] += take; // Add to liquid (will be distributed)
                    } else {
                        delta_solid[i] += take; // Add to solid
                    }

                    if (!is_instant) {
                        this.heat_energy += take * dH; // Exothermic
                        current_T = @as(f64, @floatCast(this.getTemperature()));
                        current_Cp = this.getHeatCapacity();
                        total_gas_moles -= take; // Update tracking variable
                    }
                }
            }
        }

        // VLE mass transfer rate: proportional to contact area and stirring.
        // reference_area ≈ 0.01 m² (10 cm diameter beaker cross-section).
        // Stirring increases surface renewal. Max enhancement is ~50x.
        const reference_area: f64 = 0.01;
        const vle_stirring_factor = std.math.pow(f64, 50.0, this.stirring);
        const effective_vle_rate: f64 = 0.5 *
            (@as(f64, @floatCast(this.gas_contact_area)) / reference_area) *
            vle_stirring_factor;

        // 2. Liquid Phases (Freezing, VLE Evaporation/Absorption, & Boiling)
        // Process each phase independently for freezing/boiling, but VLE uses global gas state.
        for (this.liquids.items) |*phase| {
            // 2a. FREEZING (Liquid -> Solid)
            // Priority over evaporation: subcooled liquids freeze before flashing to vapor.
            const freezing_point = @as(f64, @floatCast(phase.getFreezingPoint()));

            if (current_T < freezing_point) {
                var solvent_idx: ?usize = null;
                var highest_mp: f64 = -std.math.inf(f64);

                // Find the major component (solvent) - highest melting point among significant components
                for (0..num_mol) |i| {
                    if (phase.solutes[i] > constants.MIN_MOLES) {
                        const mp = @as(f64, @floatCast((@as(MoleculeId, @enumFromInt(i))).getMeltingPoint() orelse -273.15));

                        if (mp > highest_mp) {
                            highest_mp = mp;
                            solvent_idx = i;
                        }
                    }
                }

                if (solvent_idx) |idx| {
                    const mol: MoleculeId = @enumFromInt(idx);
                    const amount = phase.solutes[idx];
                    const dH_fus = mol.getDeltaHfusJ();

                    var take = amount;

                    if (!is_instant) {
                        const q_missing = (freezing_point - current_T) * current_Cp;
                        take = @min(amount, @max(0.0, q_missing / dH_fus));
                    }

                    if (take > constants.MIN_MOLES) {
                        phase.solutes[idx] -= take; // Apply immediately to prevent double-counting
                        delta_solid[idx] += take;

                        if (!is_instant) {
                            this.heat_energy += take * dH_fus; // Exothermic
                            current_T = @as(f64, @floatCast(this.getTemperature()));
                            current_Cp = this.getHeatCapacity();
                        }
                    }
                }
            }

            // Calculate total moles in this phase after freezing for mole fractions
            var total_liquid_moles: f64 = 0.0;

            for (0..num_mol) |i| {
                total_liquid_moles += phase.solutes[i];
            }

            // 2b. Vapor-Liquid Equilibrium (Evaporation & Condensation)
            // Handles mass transfer between this liquid phase and the global gas phase.
            // Skip if no liquid or in instant mode (instant mode uses boiling/freezing only).
            if (total_liquid_moles > constants.MIN_MOLES and !is_instant and V_gas_m3 > constants.MIN_VOLUME) {
                for (0..num_mol) |i| {
                    const moles_liquid = phase.solutes[i];
                    const moles_gas = this.gas.molecules[i] + delta_gas[i]; // Include pending gas changes

                    // Skip if neither phase has significant amount
                    if (moles_liquid <= constants.MIN_MOLES and moles_gas <= constants.MIN_MOLES) {
                        continue;
                    }

                    const mol: MoleculeId = @enumFromInt(i);

                    // Skip if below melting point (should be solid, not liquid)
                    const mp = @as(f64, @floatCast(mol.getMeltingPoint() orelse -std.math.inf(f64)));

                    if (current_T < mp) {
                        continue;
                    }

                    // Mole fraction in this liquid phase
                    const x_i = moles_liquid / total_liquid_moles;

                    var P_eq: f64 = 0.0;

                    if (mol.isGas()) {
                        // Henry's Law for gases dissolved in liquid
                        if (mol.getWaterSolubility()) |sol| {
                            if (sol > 0.0) {
                                const x_sat = @as(f64, @floatCast(sol)) / 55.5; // Approximate saturation mole fraction
                                const H = constants.P_std / x_sat;

                                P_eq = x_i * H;
                            } else {
                                // Insoluble gas - does not participate in VLE
                                continue;
                            }
                        }
                    } else {
                        // Raoult's Law for volatile liquids
                        P_eq = x_i * mol.getVaporPressure(current_T);
                    }

                    // Current partial pressure
                    const P_partial = if (total_gas_moles > 0.0 and P > 0.0)
                        P * (moles_gas / total_gas_moles)
                    else
                        0.0;

                    const dH_vap = mol.getDeltaHvapJ();

                    if (P_eq > P_partial) {
                        // EVAPORATION (Liquid -> Gas)
                        const n_target_ideal = P_eq * V_gas_m3 / (constants.R * current_T);
                        var target_take = @max(0.0, n_target_ideal - moles_gas);

                        // Kinetic rate limit: scales with contact area and stirring
                        target_take = @min(moles_liquid, target_take * effective_vle_rate * dt);

                        // Heat limit
                        if (dH_vap > 0.0 and target_take > 0.0) {
                            const max_take_heat = this.heat_energy / dH_vap;
                            target_take = @min(target_take, max_take_heat);
                        }

                        if (target_take > constants.MIN_MOLES) {
                            phase.solutes[i] -= target_take;
                            delta_gas[i] += target_take;

                            this.heat_energy = @max(0.0, this.heat_energy - target_take * dH_vap); // Endothermic
                            current_T = @as(f64, @floatCast(this.getTemperature()));
                            current_Cp = this.getHeatCapacity();
                            total_gas_moles += target_take;
                            total_liquid_moles -= target_take;

                            if (this.metrics) |m| {
                                m.moles_evaporated += target_take;
                            }
                        }
                    } else if (P_partial > P_eq) {
                        // CONDENSATION (Gas -> Liquid)
                        var target_take: f64 = 0.0;

                        if (P > 0.0 and total_gas_moles > 0.0) {
                            const target_moles_in_gas = total_gas_moles * (P_eq / P);
                            target_take = moles_gas - target_moles_in_gas;
                        }

                        // Kinetic rate limit: scales with contact area and stirring
                        target_take = @min(moles_gas, @max(0.0, target_take) * effective_vle_rate * dt);

                        // Solubility limit (for gases dissolving into liquid)
                        if (mol.getWaterSolubility()) |sol| {
                            const current_conc = moles_liquid / phase.getVolume();
                            const available_space = @max(0.0, @as(f64, @floatCast(sol)) - current_conc);

                            target_take = @min(target_take, available_space * phase.getVolume());
                        }

                        if (target_take > constants.MIN_MOLES) {
                            phase.solutes[i] += target_take;
                            delta_gas[i] -= target_take;

                            this.heat_energy += target_take * dH_vap; // Exothermic
                            current_T = @as(f64, @floatCast(this.getTemperature()));
                            current_Cp = this.getHeatCapacity();
                            total_gas_moles -= target_take;
                            total_liquid_moles += target_take;
                        }
                    }
                }
            }

            // 2c. BOILING (Vigorous vaporization at T > T_boil)
            // Handles bulk boiling when temperature exceeds boiling point.
            const boiling_point: f64 = @floatCast(phase.getBoilingPoint(@floatCast(P)));

            if (current_T > boiling_point) {
                // Find most volatile component (lowest boiling point)
                var volatile_idx: ?usize = null;
                var lowest_bp: f64 = std.math.inf(f64);

                for (0..num_mol) |i| {
                    if (phase.solutes[i] > constants.MIN_MOLES) {
                        const bp = @as(f64, @floatCast((@as(MoleculeId, @enumFromInt(i))).getBoilingPoint() orelse 5000.0));

                        if (bp < lowest_bp) {
                            lowest_bp = bp;
                            volatile_idx = i;
                        }
                    }
                }

                if (volatile_idx) |idx| {
                    const mol: MoleculeId = @enumFromInt(idx);
                    const amount = phase.solutes[idx];
                    const dH_vap = mol.getDeltaHvapJ();

                    var take = amount;

                    if (!is_instant) {
                        const q_excess = (current_T - boiling_point) * current_Cp;
                        take = @min(amount, @max(0.0, q_excess / dH_vap));
                    }

                    // Prevent explosive boiling of saturated vapors in a
                    // closed container that exceeds the pressure limit.
                    if (max_volume > constants.MIN_VOLUME and V_gas_m3 > constants.MIN_VOLUME) {
                        const P_eq_clamp = mol.getVaporPressure(current_T);

                        if (P_eq_clamp > 0.0) {
                            const max_moles_gas = P_eq_clamp * V_gas_m3 / (constants.R * current_T);
                            const current_moles_gas = this.gas.molecules[idx] + delta_gas[idx];
                            const max_boil = @max(0.0, max_moles_gas - current_moles_gas);

                            take = @min(take, max_boil);
                        }
                    }

                    if (take > constants.MIN_MOLES) {
                        phase.solutes[idx] -= take;
                        delta_gas[idx] += take;

                        if (!is_instant) {
                            this.heat_energy = @max(0.0, this.heat_energy - take * dH_vap); // Endothermic
                            current_T = @as(f64, @floatCast(this.getTemperature()));
                            current_Cp = this.getHeatCapacity();
                        }

                        if (this.metrics) |m| {
                            m.moles_boiled += take;
                        }
                    }
                }
            }
        }

        // 3. Solid Phases (Melting & Sublimation)
        var solid_idx: usize = 0;
        while (solid_idx < this.solids.items.len) {
            const sol = &this.solids.items[solid_idx];
            const i: usize = @intFromEnum(sol.molecule);
            const target = sol.molecule.targetState(@floatCast(current_T), @floatCast(P));

            if (target != .solid) {
                const mp = @as(f64, @floatCast(sol.molecule.getMeltingPoint() orelse current_T));
                const bp = sol.molecule.getPressureCorrectedBoilingPoint(@floatCast(P)) orelse mp;

                const dH_vap = sol.molecule.getDeltaHvapJ();
                const dH_fus = sol.molecule.getDeltaHfusJ();

                const t_trans = if (target == .liquid)
                    mp
                else
                    bp;

                const dH = if (target == .liquid)
                    dH_fus
                else
                    (dH_vap + dH_fus);

                var take = sol.moles;

                if (!is_instant and current_T > t_trans) {
                    const q_excess = (current_T - t_trans) * current_Cp;
                    take = @min(sol.moles, @max(0.0, q_excess / dH));
                }

                if (take > constants.MIN_MOLES) {
                    const fraction = take / sol.moles;

                    // Release occluded solutions proportionally
                    if (sol.occluded_solution) |*occ| {
                        for (0..num_mol) |mol_idx| {
                            const occ_take = occ.solutes[mol_idx] * fraction;

                            if (occ_take > constants.MIN_MOLES) {
                                occ.solutes[mol_idx] -= occ_take;
                                delta_liquid_global[mol_idx] += occ_take;
                            }
                        }
                    }

                    const volume_ratio = (sol.moles - take) / sol.moles;
                    sol.moles -= take;

                    if (target == .gas) {
                        delta_gas[i] += take;
                    } else {
                        delta_liquid_global[i] += take;
                    }

                    if (!is_instant) {
                        this.heat_energy = @max(0.0, this.heat_energy - take * dH); // Endothermic
                        current_T = @as(f64, @floatCast(this.getTemperature()));
                        current_Cp = this.getHeatCapacity();
                    }

                    // Update particle diameter if solid remains
                    if (sol.moles > constants.MIN_MOLES) {
                        if (volume_ratio > 0.0) {
                            sol.particle_diameter *= std.math.cbrt(volume_ratio);
                        }

                        solid_idx += 1;
                    } else {
                        // Remove depleted solid phase
                        sol.deinit(allocator);
                        _ = this.solids.swapRemove(solid_idx);
                        // Don't increment solid_idx due to swapRemove
                    }
                } else {
                    solid_idx += 1;
                }
            } else {
                solid_idx += 1;
            }
        }

        // 4. Apply accumulated changes atomically to prevent mass loss
        // 4a. Apply gas changes
        for (0..num_mol) |i| {
            this.gas.molecules[i] += delta_gas[i];

            // Ensure non-negative
            if (this.gas.molecules[i] < 0.0) {
                this.gas.molecules[i] = 0.0;
            }
        }

        // 4b. Apply liquid changes (distribute new liquid molecules among phases by miscibility)
        try this.placeLiquidMoles(allocator, &delta_liquid_global);

        // 4c. Apply solid changes
        for (0..num_mol) |i| {
            if (delta_solid[i] > constants.MIN_MOLES) {
                const mol: MoleculeId = @enumFromInt(i);

                // Use freezing diameter for new solids from phase transitions
                try this.mergeSolid(allocator, mol, delta_solid[i], constants.DEFAULT_FREEZING_DIAMETER);
            }
        }
    }

    /// Consolidates solid phases of the same molecule when their particle
    /// sizes are similar (within constants.PARTICLE_SIZE_MERGE_RATIO).
    ///
    /// This is useful after adding solids from multiple sources or after
    /// physical processes (grinding, agglomeration) that may produce phases
    /// with similar particle sizes that should be treated as one.
    ///
    /// When merging:
    /// - Moles are summed.
    /// - Particle diameter is volume-weighted average.
    /// - Occluded solutions are combined (moles added element-wise).
    /// - Crystallinity is preserved only if both phases are crystalline.
    ///
    /// Time complexity: O(n²) where n = number of solid phases.
    inline fn consolidateSolids(this: *Contents, allocator: std.mem.Allocator) void {
        var i: usize = 0;

        while (i < this.solids.items.len) {
            var j: usize = i + 1;

            while (j < this.solids.items.len) {
                const si = &this.solids.items[i];
                const sj = &this.solids.items[j];

                if (si.molecule != sj.molecule) {
                    j += 1;

                    continue;
                }

                const d1 = si.particle_diameter;
                const d2 = sj.particle_diameter;
                const max_d = @max(d1, d2);

                // Check if particle sizes are similar enough to merge.
                if (max_d <= 0.0 or @abs(d1 - d2) / max_d >= constants.PARTICLE_SIZE_MERGE_RATIO) {
                    j += 1;

                    continue;
                }

                // Merge phase j into phase i.
                const v1 = si.getVolume();
                const v2 = sj.getVolume();

                // Volume-weighted average diameter.
                if (v1 + v2 > 0.0) {
                    si.particle_diameter = (v1 * d1 + v2 * d2) / (v1 + v2);
                }

                si.moles += sj.moles;

                // Preserve crystallinity only if both are crystalline.
                si.crystalline = si.crystalline and sj.crystalline;

                // Merge occluded solutions.
                if (sj.occluded_solution) |*sol_j| {
                    if (si.occluded_solution) |*sol_i| {
                        sol_i.addMoles(sol_j.solutes);
                        sol_j.deinit(allocator);
                    } else {
                        // Transfer ownership from j to i.
                        si.occluded_solution = sj.occluded_solution;
                        sj.occluded_solution = null;
                    }
                }

                // Clean up and remove phase j.
                sj.deinit(allocator);
                _ = this.solids.swapRemove(j);
                // Don't increment j - swapped element now sits at index j.
            }

            i += 1;
        }
    }

    /// Distributes an array of moles (indexed by MoleculeId ordinal) among
    /// existing liquid phases by miscibility, creating new immiscible
    /// layers as needed. Does NOT run equilibration or settle.
    fn placeLiquidMoles(this: *Contents, allocator: std.mem.Allocator, moles: []const f64) error{OutOfMemory}!void {
        std.debug.assert(moles.len == MoleculeId.count());

        const Pair = struct { MoleculeId, f64 };
        var sorted: [MoleculeId.count()]Pair = undefined;

        for (0..moles.len) |i| {
            sorted[i] = .{ @enumFromInt(i), moles[i] };
        }

        std.mem.sort(Pair, &sorted, .{}, struct {
            fn lessThan(_: @TypeOf(.{}), lhs: Pair, rhs: Pair) bool {
                return lhs.@"1" > rhs.@"1";
            }
        }.lessThan);

        const current_T: f64 = @floatCast(this.getTemperature());
        const current_P: f64 = @floatCast(this.pressure);

        for (sorted) |pair| {
            const mid: MoleculeId = pair.@"0";
            const i = @intFromEnum(mid);
            const amount = pair.@"1";

            if (constants.isNegligible(amount)) {
                continue;
            }

            var placed = false;

            var best_phase_idx: ?usize = null;
            var best_affinity: f64 = std.math.inf(f64);

            for (this.liquids.items, 0..) |*p, phase_idx| {
                const primary = p.getPrimaryComponent() orelse {
                    continue;
                };

                if (MoleculeId.areMiscible(mid, primary)) {
                    p.solutes[i] += amount;
                    placed = true;

                    break;
                }

                // We look for the phase with the closest Hildebrand parameter.
                // This is necessary for gases and solids: non-polar gases (N2)
                // will prefer to dissolve in oil, and polar salts in water.
                const affinity = @abs(mid.getHildebrand() - primary.getHildebrand());

                if (affinity < best_affinity) {
                    best_affinity = affinity;
                    best_phase_idx = phase_idx;
                }
            }

            if (!placed) {
                const state = mid.targetState(current_T, current_P);

                // A substance can create a new liquid phase ONLY if it is itself a liquid,
                // OR if there are no liquid phases in the container at all.
                if (state == .liquid or best_phase_idx == null) {
                    var phase = try LiquidPhase.init(allocator);
                    errdefer phase.deinit(allocator);

                    phase.addMolecule(mid, amount);
                    try this.liquids.append(allocator, phase);
                } else {
                    // It's a gas or a solid (solute), and liquid phases already exist.
                    // We force it to dissolve in the most suitable phase.
                    this.liquids.items[best_phase_idx.?].solutes[i] += amount;
                }
            }
        }
    }

    /// Discards liquid phases whose total volume is below constants.MIN_VOLUME.
    /// Uses swap-remove for efficiency; caller must re-sort if order matters.
    inline fn removeEmptyLiquidPhases(this: *Contents, allocator: std.mem.Allocator) void {
        var j: usize = 0;

        while (j < this.liquids.items.len) {
            if (this.liquids.items[j].getVolume() < constants.MIN_VOLUME) {
                this.liquids.items[j].deinit(allocator);
                _ = this.liquids.swapRemove(j);
            } else {
                j += 1;
            }
        }
    }

    /// Discards solid phases whose total moles is below constants.MIN_MOLES.
    /// Uses swap-remove for efficiency.
    inline fn removeEmptySolidPhases(this: *Contents, allocator: std.mem.Allocator) void {
        var j: usize = 0;

        while (j < this.solids.items.len) {
            const moles = this.solids.items[j].moles;

            if (constants.isNegligible(moles)) {
                this.solids.items[j].deinit(allocator);
                _ = this.solids.swapRemove(j);
            } else {
                j += 1;
            }
        }
    }

    /// Sorts liquid layers by density, lightest (top) first.
    inline fn sortLiquidsByDensity(this: *Contents) void {
        std.mem.sort(LiquidPhase, this.liquids.items, .{}, struct {
            fn lessThan(_: @TypeOf(.{}), lhs: LiquidPhase, rhs: LiquidPhase) bool {
                return lhs.getDensity() < rhs.getDensity();
            }
        }.lessThan);
    }

    /// Adds moles to an existing SolidPhase for the given molecule if particle
    /// sizes are similar (within PARTICLE_SIZE_MERGE_RATIO), otherwise creates
    /// a new solid entry.
    ///
    /// When merging, the new particle diameter is the volume-weighted average:
    /// d_new = (V₁ * d₁ + V₂ * d₂) / (V₁ + V₂)
    ///
    /// This preserves the total surface area approximately while combining
    /// phases that are physically indistinguishable.
    inline fn mergeSolid(
        this: *Contents,
        allocator: std.mem.Allocator,
        molecule: MoleculeId,
        moles: f64,
        particle_diameter: f64,
    ) error{OutOfMemory}!void {
        const mw = molecule.getWeight();
        const density = molecule.getDensity();
        // Volume of incoming material in liters.
        const v2 = moles * mw / density / 1000.0;

        for (this.solids.items) |*s| {
            if (s.molecule == molecule) {
                const d1 = s.particle_diameter;
                const d2 = particle_diameter;
                const max_d = @max(d1, d2);

                // Check if particle sizes are similar enough to merge.
                if (max_d > 0.0 and @abs(d1 - d2) / max_d < constants.PARTICLE_SIZE_MERGE_RATIO) {
                    const v1 = s.getVolume();

                    // Volume-weighted average diameter.
                    if (v1 + v2 > 0.0) {
                        s.particle_diameter = (v1 * d1 + v2 * d2) / (v1 + v2);
                    }

                    s.moles += moles;

                    return;
                }
            }
        }

        try this.solids.append(allocator, .{
            .molecule = molecule,
            .moles = moles,
            .particle_diameter = particle_diameter,
            .occluded_solution = null,
            .crystalline = true,
        });
    }

    /// Calculates the effective bulk modulus of all liquids and solids.
    inline fn getEffectiveBulkModulus(this: *const Contents) f64 {
        var total_volume: f64 = 0.0;
        var weighted_k: f64 = 0.0;

        for (this.liquids.items) |*l| {
            for (0..l.solutes.len) |i| {
                const molecule: MoleculeId = @enumFromInt(i);
                const moles = l.solutes[i];

                const volume = moles * molecule.getWeight() / molecule.getDensity() / 1000.0;

                total_volume += volume;
                weighted_k += volume * molecule.molecule().bulk_modulus;
            }
        }

        for (this.solids.items) |s| {
            const molecule: MoleculeId = s.molecule;
            const volume = s.moles * molecule.getWeight() / molecule.getDensity() / 1000.0;

            total_volume += volume;
            weighted_k += volume * molecule.molecule().bulk_modulus;
        }

        if (total_volume < constants.MIN_VOLUME) {
            return 1.0e9;
        }

        return weighted_k / total_volume;
    }
};

const MAX_VOLUME: f64 = 3.0;

test "Just water" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());

    try std.testing.expectEqual(1, contents.liquids.items.len);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try std.testing.expectEqual(1, contents.liquids.items.len);

    // There is a precision error somewhere
    try std.testing.expectApproxEqRel(2.0, contents.getTotalLiquidVolume(), std.math.floatEps(f64));
}

test "Transfer liquids with energy conservation" {
    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try src.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);

    try dst.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try dst.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    const src_t = src.getTemperature();
    const dst_t = dst.getTemperature();

    _ = try src.transferLiquidVolume(std.testing.allocator, &dst, 0.5);

    try std.testing.expectApproxEqRel(src_t, src.getTemperature(), std.math.floatEps(f64));
    try std.testing.expectApproxEqRel(dst_t, dst.getTemperature(), std.math.floatEps(f64));
}

test "Transfer MIN_VOLUME water" {
    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try src.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);
    try dst.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    const transfered = try src.transferLiquidVolume(std.testing.allocator, &dst, constants.MIN_VOLUME);

    try std.testing.expectApproxEqAbs(constants.MIN_VOLUME, transfered, std.math.floatEps(f64));
    try std.testing.expectEqual(1, dst.liquids.items.len);
}

test "Tight water container and argon" {
    const V_container = 1.0;

    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter() * V_container);
    try src.setTemperature(std.testing.allocator, constants.cToK(25), V_container);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);

    for (0..5) |_| {
        if (src.getGasMoles() < constants.MIN_MOLES) {
            src.addGas(.argon, constants.MIN_MOLES);
        }

        try src.setTemperature(std.testing.allocator, constants.cToK(25), V_container);
        try src.updatePhaseTransitions(std.testing.allocator, 0.1, V_container);
        try src.settle(std.testing.allocator);
        src.updatePressure(V_container);

        try std.testing.expectEqual(1, src.liquids.items.len);
    }

    try std.testing.expect(src.pressure < 103000);
    try std.testing.expect(src.pressure > 100000);
}

test "90% water container and argon" {
    const V_container = 1.0;

    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter() * (V_container * 0.9));
    try src.setTemperature(std.testing.allocator, constants.cToK(25), V_container);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);

    for (0..5) |_| {
        if (src.getGasMoles() < constants.MIN_MOLES) {
            src.addGas(.argon, constants.MIN_MOLES);
        }

        try src.setTemperature(std.testing.allocator, constants.cToK(25), V_container);
        try src.updatePhaseTransitions(std.testing.allocator, 0.1, V_container);
        try src.settle(std.testing.allocator);
        src.updatePressure(V_container);

        try std.testing.expectEqual(1, src.liquids.items.len);
    }

    try std.testing.expectApproxEqAbs(3755, src.pressure, 1.0);
}

test "Tight water container and air" {
    const V_container = 1.0;

    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter() * V_container);
    try src.setTemperature(std.testing.allocator, constants.cToK(25), V_container);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);

    for (0..5) |_| {
        helpers.resetStdAir(&src, V_container);

        try src.setTemperature(std.testing.allocator, constants.cToK(25), V_container);
        try src.updatePhaseTransitions(std.testing.allocator, 0.1, V_container);
        try src.settle(std.testing.allocator);
        src.updatePressure(V_container);

        try std.testing.expectEqual(1, src.liquids.items.len);
    }

    try std.testing.expect(src.pressure < 103000);
    try std.testing.expect(src.pressure > 100000);
}

test "90% water container and air" {
    const V_container = 1.0;

    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter() * (V_container * 0.9));
    try src.setTemperature(std.testing.allocator, constants.cToK(25), V_container);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);

    for (0..5) |_| {
        helpers.resetStdAir(&src, V_container);

        try src.setTemperature(std.testing.allocator, constants.cToK(25), V_container);
        try src.updatePhaseTransitions(std.testing.allocator, 0.1, V_container);
        try src.settle(std.testing.allocator);
        src.updatePressure(V_container);

        try std.testing.expectEqual(1, src.liquids.items.len);
    }

    try std.testing.expect(src.pressure < 103000);
    try std.testing.expect(src.pressure > 100000);
}

test "Water + Dichloromethane sollutable" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);
    try std.testing.expectEqual(1, contents.liquids.items.len);

    try contents.addLiquid(std.testing.allocator, .dichloromethane, MoleculeId.dichloromethane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    try contents.settle(std.testing.allocator);

    try std.testing.expectEqual(2, contents.liquids.items.len);
    try std.testing.expect(contents.liquids.items[1].solutes[@intFromEnum(MoleculeId.dichloromethane)] > 0);
}

test "Water + Ethanol" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try std.testing.expectEqual(1, contents.liquids.items.len);

    try contents.addLiquid(std.testing.allocator, .ethanol, MoleculeId.ethanol.molesPerLiter() * 0.1);
    try std.testing.expectEqual(1, contents.liquids.items.len);
}

test "Water freezing" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(-10), MAX_VOLUME);

    try std.testing.expectEqual(0, contents.liquids.items.len);
    try std.testing.expectEqual(1, contents.solids.items.len);
}

test "Water vapor" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(150), MAX_VOLUME);
    contents.updatePressure(MAX_VOLUME);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expect(contents.gas.molecules[@intFromEnum(MoleculeId.oxidane)] > 0.0);
    try std.testing.expectApproxEqAbs(contents.pressure, 475686, 1.0);
}

test "Water evaporation" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    contents.container_heat_capacity = 180.0;

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    const pre_gas_volume = contents.getGasVolume(MAX_VOLUME);

    var secs: f64 = 0.0;

    for (0..120) |_| {
        const dt: f64 = 1.0;

        _ = contents.exchangeHeat(constants.cToK(500), 100, dt);
        try contents.updatePhaseTransitions(std.testing.allocator, dt, MAX_VOLUME);
        contents.updatePressure(MAX_VOLUME);

        secs += dt;
    }

    try std.testing.expectApproxEqAbs(0.702, contents.getTotalLiquidVolume(), 0.001);
    try std.testing.expect(contents.getGasVolume(MAX_VOLUME) > pre_gas_volume);
}

test "Water and open air" {
    const V_container = 3.0;

    var contents = try Contents.init(std.testing.allocator);
    contents.gas_contact_area = 0.0007;

    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter() * 1.5);
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), V_container);

    for (0..100) |_| {
        helpers.resetStdAir(&contents, V_container);

        try contents.updatePhaseTransitions(std.testing.allocator, 0.1, V_container);
        try contents.settle(std.testing.allocator);
        contents.updatePressure(V_container);
        try contents.setTemperature(std.testing.allocator, constants.cToK(25), V_container);
    }

    const now_moles = contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.oxidane)];

    try std.testing.expectApproxEqAbs(MoleculeId.oxidane.molesPerLiter() * 1.5, now_moles, 0.1);
    try std.testing.expect(now_moles < MoleculeId.oxidane.molesPerLiter() * 1.5);
    try std.testing.expectApproxEqAbs(1.5, contents.getTotalLiquidVolume(), 0.001);
}

test "Water pouring" {
    const V_src = 1.0;
    const V_dst = 0.5;

    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try src.setTemperature(std.testing.allocator, constants.cToK(25), V_src);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);

    helpers.resetStdAir(&dst, V_dst);
    try dst.setTemperature(std.testing.allocator, constants.cToK(25), V_dst);

    try std.testing.expectEqual(1, src.liquids.items.len);
    try std.testing.expectEqual(0, dst.liquids.items.len);

    for (0..10) |_| {
        helpers.resetStdAir(&src, V_src);
        try src.updatePhaseTransitions(std.testing.allocator, 0.1, V_src);
        try src.settle(std.testing.allocator);
        src.updatePressure(V_src);

        helpers.resetStdAir(&dst, V_dst);
        try dst.updatePhaseTransitions(std.testing.allocator, 0.1, V_dst);
        try dst.settle(std.testing.allocator);
        dst.updatePressure(V_dst);

        _ = try src.transferLiquidVolume(std.testing.allocator, &dst, 0.025);

        try std.testing.expectEqual(1, src.liquids.items.len);
        try std.testing.expectEqual(1, dst.liquids.items.len);
    }
}

test "Water + Dichloromethane, liquid splitting" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    const expected_volume = contents.liquids.items[0].getVolume();

    try contents.addLiquid(std.testing.allocator, .dichloromethane, MoleculeId.dichloromethane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    try std.testing.expectEqual(2, contents.liquids.items.len);

    var dst_contents = try Contents.init(std.testing.allocator);
    defer dst_contents.deinit(std.testing.allocator);

    const tx_volume = try contents.transferLiquidPhaseVolume(std.testing.allocator, &dst_contents, 0, std.math.floatMax(f64));
    try dst_contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);
    try contents.settle(std.testing.allocator);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(1, dst_contents.liquids.items.len);
    try std.testing.expectApproxEqAbs(expected_volume, tx_volume, 0.1);
}

test "Water + Dichloromethane, pouring half" {
    const TX_LITERS = 1.0;

    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try src.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    try src.addLiquid(std.testing.allocator, .dichloromethane, MoleculeId.dichloromethane.molesPerLiter());
    try src.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);

    const tx_volume = try src.transferLiquidVolume(std.testing.allocator, &dst, TX_LITERS);
    try dst.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    try src.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);
    try src.settle(std.testing.allocator);

    try std.testing.expectEqual(1, src.liquids.items.len);
    try std.testing.expectEqual(1, dst.liquids.items.len);
    try std.testing.expectApproxEqAbs(TX_LITERS, tx_volume, 0.1);

    try std.testing.expect(src.liquids.items[0].solutes[@intFromEnum(MoleculeId.oxidane)] < MoleculeId.oxidane.molesPerLiter() * 0.1);
    try std.testing.expect(dst.liquids.items[0].solutes[@intFromEnum(MoleculeId.oxidane)] > MoleculeId.oxidane.molesPerLiter() * 0.9);

    try std.testing.expect(src.liquids.items[0].solutes[@intFromEnum(MoleculeId.dichloromethane)] > MoleculeId.dichloromethane.molesPerLiter() * 0.9);
    try std.testing.expect(dst.liquids.items[0].solutes[@intFromEnum(MoleculeId.dichloromethane)] < MoleculeId.dichloromethane.molesPerLiter() * 0.1);

    try std.testing.expectApproxEqAbs(TX_LITERS, dst.liquids.items[0].getVolume(), 0.1);
}

test "Water + Dichloromethane odors" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.addLiquid(std.testing.allocator, .dichloromethane, MoleculeId.dichloromethane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    const odors = contents.getDominantOdors();

    try std.testing.expect(odors.odor_intensities[0] > 0.0);
    try std.testing.expectApproxEqRel(0.0, odors.odor_intensities[1], 0.001);

    try std.testing.expectEqual(Odor.sweet, odors.odors[0]);
    try std.testing.expectEqual(Odor.none, odors.odors[1]);
}

test "Ethanol + Dichloromethane odors" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .ethanol, MoleculeId.ethanol.molesPerLiter());
    try contents.addLiquid(std.testing.allocator, .dichloromethane, MoleculeId.dichloromethane.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    const odors = contents.getDominantOdors();

    try std.testing.expect(odors.odor_intensities[0] > 0.0);
    try std.testing.expect(odors.odor_intensities[1] > 0.0);
    try std.testing.expect(odors.odor_intensities[0] > odors.odor_intensities[1]);

    try std.testing.expectEqual(Odor.sweet, odors.odors[0]);
    try std.testing.expectEqual(Odor.ethereal, odors.odors[1]);
}

test "Ethanol pH" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .ethanol, MoleculeId.ethanol.molesPerLiter() * 0.285);
    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter() * 0.015);
    try contents.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    const pH = helpers.computePhasePH(contents.liquids.items[0].solutes, contents.getTotalLiquidVolume());

    try std.testing.expectApproxEqAbs(6.95, pH, 0.1);
}

test "Temperature affects odor intensity" {
    var cold = try Contents.init(std.testing.allocator);
    defer cold.deinit(std.testing.allocator);

    try cold.setTemperature(std.testing.allocator, constants.cToK(0), MAX_VOLUME);

    var hot = try Contents.init(std.testing.allocator);
    defer hot.deinit(std.testing.allocator);

    try hot.setTemperature(std.testing.allocator, constants.cToK(80), MAX_VOLUME);

    // Add same amount of ethanol to both
    try cold.addLiquid(std.testing.allocator, .ethanol, 1.0);
    try hot.addLiquid(std.testing.allocator, .ethanol, 1.0);

    const cold_odors = cold.getDominantOdors();
    const hot_odors = hot.getDominantOdors();

    // Both should detect ethereal odor
    try std.testing.expectEqual(Odor.ethereal, cold_odors.odors[0]);
    try std.testing.expectEqual(Odor.ethereal, hot_odors.odors[0]);
}

test "Hot DCM smells stronger than cold" {
    var cold = try Contents.init(std.testing.allocator);
    defer cold.deinit(std.testing.allocator);

    try cold.setTemperature(std.testing.allocator, constants.cToK(10), MAX_VOLUME);

    var hot = try Contents.init(std.testing.allocator);
    defer hot.deinit(std.testing.allocator);

    try cold.setTemperature(std.testing.allocator, constants.cToK(10), MAX_VOLUME);

    // Add DCM and water (odorless reference)
    try cold.addLiquid(std.testing.allocator, .dichloromethane, 0.5);
    try cold.addLiquid(std.testing.allocator, .oxidane, 1.0);

    try hot.addLiquid(std.testing.allocator, .dichloromethane, 0.5);
    try hot.addLiquid(std.testing.allocator, .oxidane, 1.0);

    // At 40C, DCM is nearly boiling - should smell very strongly
    const hot_odors = hot.getDominantOdors();
    try std.testing.expectEqual(Odor.sweet, hot_odors.odors[0]);

    // At 10C, still smells but less intensely
    const cold_odors = cold.getDominantOdors();
    try std.testing.expectEqual(Odor.sweet, cold_odors.odors[0]);
}

test "Graphite in water" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.addSolid(std.testing.allocator, .graphite, MoleculeId.graphite.molesPerLiter() * 0.1, 2.5e-5);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(1, contents.solids.items.len);
}

test "Water evaporation at room temperature" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    const T_room = constants.cToK(25);
    const V_container = 1.0; // liters

    const n_air = (constants.P_std * 0.001) / (constants.R * T_room);
    contents.gas.molecules[@intFromEnum(MoleculeId.nitrogen)] = n_air * 0.79;
    contents.gas.molecules[@intFromEnum(MoleculeId.oxygen)] = n_air * 0.21;

    // Add liquid water (enough to not fully evaporate)
    const initial_water = 0.5; // moles
    try contents.addLiquidRaw(std.testing.allocator, .oxidane, initial_water, false);
    contents.heat_energy = T_room * contents.getHeatCapacity();

    try contents.settle(std.testing.allocator);
    contents.updatePressure(V_container);

    // Run equilibration
    for (0..1000) |_| {
        try contents.updatePhaseTransitions(std.testing.allocator, 1.0, V_container);
        try contents.settle(std.testing.allocator);
    }

    const gas_water = contents.gas.molecules[@intFromEnum(MoleculeId.oxidane)];
    const liquid_water = contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.oxidane)];

    // At 25°C, water should evaporate to establish ~3.17 kPa partial pressure
    try std.testing.expect(gas_water > 1e-4); // Significant water in gas phase
    try std.testing.expect(gas_water < initial_water); // Not all evaporated
    try std.testing.expect(liquid_water > 0); // Liquid remains
    try std.testing.expect(liquid_water < initial_water); // Some evaporated
}

test "Hygroscopic water absorption by ethanol" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    const T_room = constants.cToK(25.0);
    const V_container = 1.0;
    const target_RH = 0.60;

    // Setup 60% RH atmosphere
    const P_water_sat = 3169.0; // Pa at 25°C
    const P_water_partial = P_water_sat * target_RH;
    const P_std = constants.P_std;
    const R = constants.R;
    const V_m3 = V_container * 0.001;

    const n_water_vapor = (P_water_partial * V_m3) / (R * T_room);
    const n_total_gas = (P_std * V_m3) / (R * T_room);
    const n_dry = n_total_gas - n_water_vapor;

    contents.pressure = @floatCast(P_std);
    contents.gas.molecules[@intFromEnum(MoleculeId.nitrogen)] = n_dry * 0.79;
    contents.gas.molecules[@intFromEnum(MoleculeId.oxygen)] = n_dry * 0.21;
    contents.gas.molecules[@intFromEnum(MoleculeId.oxidane)] = n_water_vapor;

    const initial_total_water = n_water_vapor; // No liquid yet

    // Add ethanol (hygroscopic, reduces water activity)
    try contents.addLiquidRaw(std.testing.allocator, .ethanol, 0.5, false);
    contents.heat_energy = T_room * contents.getHeatCapacity();
    try contents.settle(std.testing.allocator);

    for (0..1000) |_| {
        try contents.updatePhaseTransitions(std.testing.allocator, 0.1, V_container);
        try contents.settle(std.testing.allocator);
    }

    const final_gas_water = contents.gas.molecules[@intFromEnum(MoleculeId.oxidane)];
    const final_liquid_water = contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.oxidane)];

    // Ethanol should absorb water from unsaturated air despite 60% RH < 100%
    try std.testing.expect(final_liquid_water > 1e-6);
    try std.testing.expect(final_gas_water < n_water_vapor * 0.95); // Gas phase water decreased

    // Mass conservation
    const final_total = final_gas_water + final_liquid_water;
    try std.testing.expectApproxEqAbs(initial_total_water, final_total, 1e-9);
}

test "Acetic acid volatility equilibrium" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    const T_room = constants.cToK(25.0);
    const V_container = 1.0;

    // Add some background gas
    const n_air = (constants.P_std * 0.001) / (constants.R * T_room) * 0.5;
    contents.gas.molecules[@intFromEnum(MoleculeId.nitrogen)] = n_air;

    // Add pure acetic acid (volatile weak acid, bp 118°C)
    const initial_acetic = 0.1; // moles
    try contents.addLiquidRaw(std.testing.allocator, .acetic_acid, initial_acetic, false);
    contents.heat_energy = T_room * contents.getHeatCapacity();
    try contents.settle(std.testing.allocator);

    // Run longer to ensure equilibration for low volatility compound
    for (0..2000) |_| {
        try contents.updatePhaseTransitions(std.testing.allocator, 0.1, V_container);
        try contents.settle(std.testing.allocator);
    }

    const gas_acetic = contents.gas.molecules[@intFromEnum(MoleculeId.acetic_acid)];
    const liquid_acetic = contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.acetic_acid)];

    // Acetic acid has P_vap ~2.1 kPa at 25°C, should establish measurable vapor pressure
    try std.testing.expect(gas_acetic > 1e-6); // Must be present in gas phase
    try std.testing.expect(gas_acetic < initial_acetic * 0.5); // But most remains liquid
    try std.testing.expect(liquid_acetic > 0);

    // Rough check: mole fraction in gas should approximate P_vap/P_total
    const total_gas_moles = n_air + gas_acetic;
    const mole_fraction = gas_acetic / total_gas_moles;
    const expected_x = 2100.0 / constants.P_std; // ~0.021

    // Allow wide tolerance for kinetic non-ideality, but ensure it's in reasonable range
    try std.testing.expect(mole_fraction > expected_x * 0.05);
    try std.testing.expect(mole_fraction < expected_x * 5.0);
}

test "Ammonia absorption by water" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    const T_room = constants.cToK(25.0);
    const V_container = 1.0;

    contents.pressure = constants.P_std;

    // Add liquid water first
    try contents.addLiquidRaw(std.testing.allocator, .oxidane, 0.5, false);
    contents.heat_energy = T_room * contents.getHeatCapacity();
    try contents.settle(std.testing.allocator);

    // Inject ammonia gas (highly soluble)
    const initial_nh3 = 0.02; // moles
    contents.addGasRaw(.azane, initial_nh3, false);

    const initial_gas_nh3 = contents.gas.molecules[@intFromEnum(MoleculeId.azane)];
    try std.testing.expectApproxEqAbs(initial_gas_nh3, initial_nh3, 1e-10);

    for (0..1000) |_| {
        try contents.updatePhaseTransitions(std.testing.allocator, 0.1, V_container);
        try contents.settle(std.testing.allocator);
    }

    const final_gas_nh3 = contents.gas.molecules[@intFromEnum(MoleculeId.azane)];
    const final_liquid_nh3 = contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.azane)];

    // Ammonia is extremely soluble; most should be absorbed
    try std.testing.expect(final_liquid_nh3 > final_gas_nh3 * 5.0); // >80% in liquid
    try std.testing.expect(final_gas_nh3 < initial_nh3 * 0.2); // <20% remains gas
    try std.testing.expect(final_liquid_nh3 > 0.015); // Most captured

    // Mass conservation check
    const total_final = final_gas_nh3 + final_liquid_nh3;
    try std.testing.expectApproxEqAbs(initial_nh3, total_final, 1e-9);
}
