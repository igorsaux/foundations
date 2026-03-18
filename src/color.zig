// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

pub const Color = enum(u8) {
    transparent,
    white,
    black,
    red,
    green,
    yellow,
    blue,
    brown,
    orange,
    pink,
    purple,
    gray,

    pub inline fn getHex(this: Color) []const u8 {
        return switch (this) {
            .transparent => "#FFFFFF00",
            .white => "#FFFFFF",
            .black => "#1A1A1A",
            .red => "#CC3333",
            .green => "#33CC33",
            .yellow => "#CCCC33",
            .blue => "#3333CC",
            .brown => "#8B4513",
            .orange => "#FF8C00",
            .pink => "#FF69B4",
            .purple => "#9932CC",
            .gray => "#808080",
        };
    }
};
