// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

/// Internal struct for scoring colors/odors during accumulation.
pub const ScoredAttribute = struct {
    index: u8,
    score: f64,
};
