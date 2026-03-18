// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const Color = @import("color.zig").Color;
const MoleculeId = @import("molecule_id.zig").MoleculeId;
const Odor = @import("odor.zig").Odor;

/// Flavor profile containing up to 2 dominant colors and 2 dominant odors.
pub const Flavor = struct {
    /// Up to 2 dominant colors, sorted by intensity (strongest first).
    /// Unused slots are set to .transparent.
    colors: [2]Color = .{ .transparent, .transparent },
    /// Up to 2 dominant odors, sorted by intensity (strongest first).
    /// Unused slots are set to .none.
    odors: [2]Odor = .{ .none, .none },
    /// Relative intensity of each color (0.0 to 1.0).
    color_intensities: [2]f32 = .{ 0.0, 0.0 },
    /// Relative intensity of each odor (0.0 to 1.0).
    odor_intensities: [2]f32 = .{ 0.0, 0.0 },
};
