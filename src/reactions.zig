// Copyright (C) 2026 Igor Spichkin
// SPDX-License-Identifier: MPL-2.0

const std = @import("std");

const constants = @import("constants.zig");
const Contents = @import("contents.zig").Contents;
const helpers = @import("helpers.zig");
const MoleculeId = @import("molecule_id.zig").MoleculeId;

/// Stoichiometric coefficient for a reaction participant
pub const Stoich = struct {
    molecule: MoleculeId,
    /// Stoichiometric coefficient (positive number)
    coefficient: f64 = 1.0,
};

/// Kinetic rate law types
pub const KineticsType = union(enum) {
    /// Mass action: v = k_forward * ∏[S]^n - k_reverse * ∏[P]^n
    /// Use for: fast equilibria, simple binding, non-enzymatic reactions
    mass_action: MassAction,

    /// Michaelis-Menten: v = kcat * [E] * [S] / (Km + [S])
    /// Use for: most single-substrate enzymatic reactions
    michaelis_menten: MichaelisMenten,

    /// Hill equation: v = kcat * [E] * [S]^n / (K^n + [S]^n)
    /// Use for: cooperative binding, allosteric enzymes, hemoglobin
    hill: Hill,

    /// Ordered Bi-Bi: two substrates bind in order
    /// Use for: NAD-dependent dehydrogenases
    ordered_bi_bi: OrderedBiBi,

    /// Ping-Pong Bi-Bi: substrates bind/release alternately
    /// Use for: transaminases, some kinases
    ping_pong: PingPong,

    /// Random Bi-Bi: substrates can bind in any order
    /// Use for: creatine kinase, hexokinase4
    random_bi_bi: RandomBiBi,

    /// First-order decay: v = k * [S]
    /// Use for: protein degradation, mRNA decay, radioactive decay
    first_order: FirstOrder,

    /// Arrhenius temperature-dependent kinetics
    /// k(T) = A * exp(-Ea / RT)
    /// Use for: any reaction where temperature dependence matters
    arrhenius: Arrhenius,

    /// Eyring (Transition State Theory): more accurate than Arrhenius
    /// k(T) = (kB*T/h) * exp(ΔS‡/R) * exp(-ΔH‡/RT)
    /// Use for: precise temperature modeling, when activation entropy matters
    eyring: Eyring,

    /// Surface reaction (Langmuir-Hinshelwood)
    /// v = k * θ_A * θ_B  where θ = K*P / (1 + K*P)
    /// Use for: heterogeneous catalysis, adsorption-limited reactions
    langmuir_hinshelwood: LangmuirHinshelwood,

    /// Eley-Rideal mechanism
    /// v = k * θ_A * P_B  (one adsorbed, one from gas/solution)
    /// Use for: gas-surface reactions, some catalytic processes
    eley_rideal: EleyRideal,

    /// Diffusion-limited reaction (Smoluchowski)
    /// k = 4π * D * R * NA  where D = diffusion coefficient
    /// Use for: very fast reactions limited by molecular encounter rate
    diffusion_limited: DiffusionLimited,

    /// Autocatalytic: A + B -> 2B (product catalyzes own formation)
    /// v = k * [A] * [B]
    /// Use for: combustion, polymerization, oscillating reactions
    autocatalytic: Autocatalytic,

    /// nth-order kinetics: v = k * [A]^n
    /// Use for: empirical fits, complex mechanisms simplified to single step
    nth_order: NthOrder,

    /// Equilibrium-controlled (uses Keq instead of k_reverse)
    /// v = k_forward * (∏[S]^n - ∏[P]^n / Keq)
    /// Use for: when equilibrium constant is known but not individual rates
    equilibrium: Equilibrium,

    /// Ionic dissociation: AB ⇌ A⁺ + B⁻
    /// For salts, acids, and bases in solution.
    /// Uses Debye-Hückel activity corrections for ionic strength effects.
    dissociation: Dissociation,

    /// Strong electrolyte (complete dissociation): AB → A⁺ + B⁻
    /// Instantaneous, diffusion-limited dissolution.
    /// Use for: NaCl, HCl, NaOH, KNO₃, etc.
    strong_electrolyte: StrongElectrolyte,

    /// Weak electrolyte equilibrium: HA ⇌ H⁺ + A⁻
    /// pH-dependent dissociation with Ka/Kb.
    /// Use for: acetic acid, ammonia, carbonic acid, amino acids.
    weak_electrolyte: WeakElectrolyte,

    /// Salt dissolution from solid: AB(s) ⇌ A⁺(aq) + B⁻(aq)
    /// Surface-area dependent with Ksp.
    /// Use for: CaCO₃, BaSO₄, AgCl precipitation/dissolution.
    dissolution: Dissolution,

    pub const MassAction = struct {
        /// Forward rate constant (units depend on reaction order)
        /// First order: 1/s
        /// Second order: 1/(M*s)
        k_forward: f64,
        /// Reverse rate constant (same units as forward)
        k_reverse: f64 = 0.0,
        /// Substrates consumed
        substrates: []const Stoich,
        /// Products formed
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const MassAction) []const Stoich {
            return this.substrates;
        }

        pub inline fn getProducts(this: *const MassAction) []const Stoich {
            return this.products;
        }
    };

    pub const MichaelisMenten = struct {
        /// Catalytic constant (turnover number) in 1/s
        /// kcat = Vmax / [E]total
        /// Typical range: 1 - 10^6 s^-1
        kcat: f64,
        /// Michaelis constant in M (mol/L)
        /// = (k-1 + kcat) / k1
        /// Approximates substrate concentration at half-maximal velocity
        /// Typical range: 10^-6 to 10^-2 M
        km: f64,
        /// The enzyme catalyzing this reaction
        enzyme: MoleculeId,
        /// Substrate (single substrate reaction)
        substrate: Stoich,
        /// Product(s) formed
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const MichaelisMenten) []const Stoich {
            return &.{this.substrate};
        }

        pub inline fn getProducts(this: *const MichaelisMenten) []const Stoich {
            return this.products;
        }
    };

    pub const Hill = struct {
        /// Catalytic constant in 1/s
        kcat: f64,
        /// Half-saturation constant (K0.5) in M
        k: f64,
        /// Hill coefficient
        /// n = 1: no cooperativity (reduces to MM)
        /// n > 1: positive cooperativity (sigmoidal)
        /// n < 1: negative cooperativity
        /// Typical: 1.5 - 4.0 for allosteric enzymes
        n: f64,
        /// The enzyme
        enzyme: MoleculeId,
        /// Substrate
        substrate: Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const Hill) []const Stoich {
            return &.{this.substrate};
        }

        pub inline fn getProducts(this: *const Hill) []const Stoich {
            return this.products;
        }
    };

    pub const OrderedBiBi = struct {
        /// Catalytic constant in 1/s
        kcat: f64,
        /// Km for substrate A (binds first)
        km_a: f64,
        /// Km for substrate B (binds second)
        km_b: f64,
        /// The enzyme
        enzyme: MoleculeId,
        /// Two substrates: [A, B] where A binds first
        substrates: [2]Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const OrderedBiBi) []const Stoich {
            return &this.substrates;
        }

        pub inline fn getProducts(this: *const OrderedBiBi) []const Stoich {
            return this.products;
        }
    };

    pub const PingPong = struct {
        kcat: f64,
        km_a: f64,
        km_b: f64,
        enzyme: MoleculeId,
        /// Two substrates: [A, B]
        substrates: [2]Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const PingPong) []const Stoich {
            return &this.substrates;
        }

        pub inline fn getProducts(this: *const PingPong) []const Stoich {
            return this.products;
        }
    };

    pub const RandomBiBi = struct {
        kcat: f64,
        km_a: f64,
        km_b: f64,
        /// Interaction factor (alpha)
        /// < 1: substrates enhance each other's binding
        /// > 1: substrates inhibit each other's binding
        /// = 1: independent binding
        alpha: f64 = 1.0,
        enzyme: MoleculeId,
        /// Two substrates: [A, B]
        substrates: [2]Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const RandomBiBi) []const Stoich {
            return &this.substrates;
        }

        pub inline fn getProducts(this: *const RandomBiBi) []const Stoich {
            return this.products;
        }
    };

    pub const FirstOrder = struct {
        /// Rate constant in 1/s
        /// Half-life = ln(2) / k ≈ 0.693 / k
        k: f64,
        /// Single substrate undergoing decay
        substrate: Stoich,
        /// Products (can be empty for pure decay)
        products: []const Stoich = &.{},

        pub inline fn getSubstrates(this: *const FirstOrder) []const Stoich {
            return &.{this.substrate};
        }

        pub inline fn getProducts(this: *const FirstOrder) []const Stoich {
            return this.products;
        }
    };

    pub const Arrhenius = struct {
        /// Pre-exponential factor A (same units as rate constant)
        /// Also called frequency factor
        /// Typical: 10^10 - 10^14 s^-1 for unimolecular
        A: f64,
        /// Activation energy in J/mol (NOT kJ/mol)
        /// Typical: 40,000 - 200,000 J/mol for most reactions
        Ea: f64,
        /// Reverse reaction parameters (null = irreversible)
        A_reverse: ?f64 = null,
        Ea_reverse: ?f64 = null,
        /// Substrates
        substrates: []const Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const Arrhenius) []const Stoich {
            return this.substrates;
        }

        pub inline fn getProducts(this: *const Arrhenius) []const Stoich {
            return this.products;
        }
    };

    pub const Eyring = struct {
        /// Activation enthalpy ΔH‡ in J/mol
        delta_H: f64,
        /// Activation entropy ΔS‡ in J/(mol·K)
        /// Negative = associative (two molecules come together)
        /// Positive = dissociative (molecule falls apart)
        delta_S: f64,
        /// Transmission coefficient κ (usually 0.5 - 1.0)
        kappa: f64 = 1.0,
        /// Reverse parameters (null = irreversible)
        delta_H_reverse: ?f64 = null,
        delta_S_reverse: ?f64 = null,
        /// Substrates
        substrates: []const Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const Eyring) []const Stoich {
            return this.substrates;
        }

        pub inline fn getProducts(this: *const Eyring) []const Stoich {
            return this.products;
        }
    };

    pub const LangmuirHinshelwood = struct {
        /// Surface reaction rate constant (1/s at full coverage)
        k_surface: f64,
        /// Adsorption equilibrium constant for species A (1/Pa or L/mol)
        K_ads_a: f64,
        /// Adsorption equilibrium constant for species B
        K_ads_b: f64,
        /// Two substrates that adsorb: [A, B]
        substrates: [2]Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const LangmuirHinshelwood) []const Stoich {
            return &this.substrates;
        }

        pub inline fn getProducts(this: *const LangmuirHinshelwood) []const Stoich {
            return this.products;
        }
    };

    pub const EleyRideal = struct {
        /// Reaction rate constant
        k: f64,
        /// Adsorption equilibrium constant for adsorbed species
        K_ads: f64,
        /// Species that adsorbs to surface (substrate[0])
        /// Species that reacts from gas/solution phase (substrate[1])
        substrates: [2]Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const EleyRideal) []const Stoich {
            return &this.substrates;
        }

        pub inline fn getProducts(this: *const EleyRideal) []const Stoich {
            return this.products;
        }
    };

    pub const DiffusionLimited = struct {
        /// Combined diffusion coefficient D_A + D_B in m²/s
        /// Water at 25°C: ~1-2 * 10^-9 m²/s for small molecules
        diffusion_coeff: f64,
        /// Reaction radius (sum of molecular radii) in meters
        /// Typical: 0.2 - 1.0 nm = 2e-10 to 1e-9 m
        reaction_radius: f64,
        /// Two substrates that collide: [A, B]
        substrates: [2]Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const DiffusionLimited) []const Stoich {
            return &this.substrates;
        }

        pub inline fn getProducts(this: *const DiffusionLimited) []const Stoich {
            return this.products;
        }
    };

    pub const Autocatalytic = struct {
        /// Rate constant for A + B -> 2B
        k: f64,
        /// Optional saturation (logistic growth limit)
        /// If set, rate = k * [A] * [B] * (1 - [B]/capacity)
        capacity: ?f64 = null,
        /// The consumed substrate (A)
        substrate: Stoich,
        /// The autocatalytic product (B) - appears both as catalyst and product
        catalyst: Stoich,

        pub inline fn getSubstrates(this: *const Autocatalytic) []const Stoich {
            // Both substrate and catalyst are consumed in A + B -> 2B
            // But catalyst is regenerated (net: A -> B)
            return &.{this.substrate};
        }

        pub inline fn getProducts(this: *const Autocatalytic) []const Stoich {
            return &.{this.catalyst};
        }
    };

    pub const NthOrder = struct {
        /// Rate constant (units: M^(1-n) / s)
        k: f64,
        /// Reaction order (can be fractional)
        /// 0 = zero order, 0.5 = half order, 2 = second order, etc.
        n: f64,
        /// Substrate the order applies to
        substrate: Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const NthOrder) []const Stoich {
            return &.{this.substrate};
        }

        pub inline fn getProducts(this: *const NthOrder) []const Stoich {
            return this.products;
        }
    };

    pub const Equilibrium = struct {
        /// Forward rate constant
        k_forward: f64,
        /// Equilibrium constant Keq = [products] / [reactants] at 298.15 K
        Keq: f64,
        /// Substrates
        substrates: []const Stoich,
        /// Products
        products: []const Stoich,

        pub inline fn getSubstrates(this: *const Equilibrium) []const Stoich {
            return this.substrates;
        }

        pub inline fn getProducts(this: *const Equilibrium) []const Stoich {
            return this.products;
        }
    };

    pub const Dissociation = struct {
        /// Dissociation constant Ka (dimensionless or M for 1:1)
        /// Ka = [A⁺][B⁻] / [AB]
        Ka: f64,
        /// The undissociated species (substrate)
        parent: Stoich,
        /// The dissociated ions (products)
        ions: []const Stoich,

        pub inline fn getSubstrates(this: *const Dissociation) []const Stoich {
            return &.{this.parent};
        }

        pub inline fn getProducts(this: *const Dissociation) []const Stoich {
            return this.ions;
        }
    };

    pub const StrongElectrolyte = struct {
        /// Dissolution rate constant (mol/(m²·s))
        k_dissolve: f64 = 1e6,
        /// Saturation concentration in mol/L
        saturation: f64,
        /// The solid salt (substrate)
        solid: Stoich,
        /// The ions produced: [cation, anion]
        ions: [2]Stoich,

        pub inline fn getSubstrates(this: *const StrongElectrolyte) []const Stoich {
            return &.{this.solid};
        }

        pub inline fn getProducts(this: *const StrongElectrolyte) []const Stoich {
            return &this.ions;
        }
    };

    pub const WeakElectrolyte = struct {
        /// Acid dissociation constant pKa (or pKb for bases)
        pKa: f64,
        /// Forward rate constant for proton transfer (1/s)
        k_forward: f64 = 1e5,
        /// Is this a base (uses Kb instead of Ka)?
        is_base: bool = false,
        /// The weak acid/base (substrate)
        parent: Stoich,
        /// Products: [conjugate, proton] for acid, [conjugate, hydroxide] for base
        products: [2]Stoich,

        pub inline fn getSubstrates(this: *const WeakElectrolyte) []const Stoich {
            return &.{this.parent};
        }

        pub inline fn getProducts(this: *const WeakElectrolyte) []const Stoich {
            return &this.products;
        }
    };

    pub const Dissolution = struct {
        /// Solubility product Ksp
        Ksp: f64,
        /// Dissolution rate constant (mol/(m²·s))
        k_dissolve: f64,
        /// Reference temperature for Ksp (K)
        T_ref: f64 = 298.15,
        /// The solid salt (substrate)
        solid: Stoich,
        /// The ions produced: [cation, anion]
        ions: [2]Stoich,

        pub inline fn getSubstrates(this: *const Dissolution) []const Stoich {
            return &.{this.solid};
        }

        pub inline fn getProducts(this: *const Dissolution) []const Stoich {
            return &this.ions;
        }
    };

    /// Get substrates for any kinetics type
    pub inline fn getSubstrates(this: *const KineticsType) []const Stoich {
        return switch (this.*) {
            .mass_action => |*k| k.getSubstrates(),
            .michaelis_menten => |*k| k.getSubstrates(),
            .hill => |*k| k.getSubstrates(),
            .ordered_bi_bi => |*k| k.getSubstrates(),
            .ping_pong => |*k| k.getSubstrates(),
            .random_bi_bi => |*k| k.getSubstrates(),
            .first_order => |*k| k.getSubstrates(),
            .arrhenius => |*k| k.getSubstrates(),
            .eyring => |*k| k.getSubstrates(),
            .langmuir_hinshelwood => |*k| k.getSubstrates(),
            .eley_rideal => |*k| k.getSubstrates(),
            .diffusion_limited => |*k| k.getSubstrates(),
            .autocatalytic => |*k| k.getSubstrates(),
            .nth_order => |*k| k.getSubstrates(),
            .equilibrium => |*k| k.getSubstrates(),
            .dissociation => |*k| k.getSubstrates(),
            .strong_electrolyte => |*k| k.getSubstrates(),
            .weak_electrolyte => |*k| k.getSubstrates(),
            .dissolution => |*k| k.getSubstrates(),
        };
    }

    /// Get products for any kinetics type
    pub inline fn getProducts(this: *const KineticsType) []const Stoich {
        return switch (this.*) {
            .mass_action => |*k| k.getProducts(),
            .michaelis_menten => |*k| k.getProducts(),
            .hill => |*k| k.getProducts(),
            .ordered_bi_bi => |*k| k.getProducts(),
            .ping_pong => |*k| k.getProducts(),
            .random_bi_bi => |*k| k.getProducts(),
            .first_order => |*k| k.getProducts(),
            .arrhenius => |*k| k.getProducts(),
            .eyring => |*k| k.getProducts(),
            .langmuir_hinshelwood => |*k| k.getProducts(),
            .eley_rideal => |*k| k.getProducts(),
            .diffusion_limited => |*k| k.getProducts(),
            .autocatalytic => |*k| k.getProducts(),
            .nth_order => |*k| k.getProducts(),
            .equilibrium => |*k| k.getProducts(),
            .dissociation => |*k| k.getProducts(),
            .strong_electrolyte => |*k| k.getProducts(),
            .weak_electrolyte => |*k| k.getProducts(),
            .dissolution => |*k| k.getProducts(),
        };
    }
};

/// Modulator that affects reaction rate
pub const Modulator = struct {
    molecule: MoleculeId,
    effect: Effect,
    mechanism: Mechanism,
    /// Inhibition/activation constant in M
    ki: f64,
    /// For competitive inhibition: factor by which Km increases at [I] = Ki
    /// Default = 1.0 (Km doubles at [I] = Ki)
    /// For allosteric: max fold-change in rate
    max_effect: f64 = 1.0,

    pub const Effect = enum {
        inhibitor,
        activator,
    };

    pub const Mechanism = enum {
        /// Competitive: increases apparent Km, no effect on Vmax.
        /// WARNING: post-hoc modulation approximates this as non-competitive.
        /// For accurate modeling, incorporate Ki into the kinetics struct directly.
        competitive,

        /// Uncompetitive: reduces both apparent Km and Vmax.
        /// Approximated as pure Vmax reduction (accurate at saturating [S]).
        uncompetitive,

        /// Non-competitive (pure): reduces Vmax, no effect on Km.
        /// Post-hoc scaling is exact for this mechanism.
        non_competitive,

        /// Mixed: affects both Km and Vmax differently.
        /// Uses max_effect as potency scaling factor.
        mixed,

        /// Allosteric: sigmoidal dose-response.
        /// max_effect = maximum fractional change (use <= 1.0 for inhibitors).
        allosteric,

        /// Product inhibition (common in metabolism).
        /// Approximated as non-competitive.
        product,

        /// Irreversible (suicide/mechanism-based).
        irreversible,
    };
};

pub const Reaction = struct {
    /// Human-readable name for debugging
    name: []const u8 = "",

    /// Rate law and kinetic parameters
    kinetics: KineticsType,

    /// Molecules that modulate reaction rate
    modulators: []const Modulator = &.{},

    /// True if reaction runs in both directions significantly
    /// For MM kinetics: adds product inhibition term
    reversible: bool = false,

    /// pH optimum (center of activity bell curve)
    /// null = pH-independent
    ph_optimum: ?f64 = null,

    /// pH sensitivity (width of bell curve)
    /// Lower = more sensitive to pH changes
    /// Typical: 0.5 - 2.0
    ph_width: f64 = 1.0,

    /// Temperature coefficient (Q10)
    /// Rate multiplier per 10°C increase
    /// Most enzymes: 2.0 - 3.0
    q10: f64 = 2.0,

    /// Get substrates from kinetics
    pub inline fn getSubstrates(this: *const Reaction) []const Stoich {
        return this.kinetics.getSubstrates();
    }

    /// Get products from kinetics
    pub inline fn getProducts(this: *const Reaction) []const Stoich {
        return this.kinetics.getProducts();
    }
};

/// Where a reaction takes place
pub const ReactionLocus = union(enum) {
    /// Homogeneous reaction within a single liquid phase
    /// The phase is identified by its primary solvent
    pub const LiquidBulk = struct {
        /// Primary solvent that identifies the phase (e.g., .oxidane for aqueous)
        solvent: MoleculeId,
    };

    /// Reaction at liquid-liquid interface (e.g., extraction, phase-transfer catalysis)
    pub const LiquidInterface = struct {
        phase_a_solvent: MoleculeId,
        phase_b_solvent: MoleculeId,
    };

    /// Gas-liquid interface reaction (e.g., CO2 absorption, evaporation)
    pub const GasLiquidInterface = struct {
        liquid_solvent: MoleculeId,
    };

    /// Heterogeneous reaction at solid-liquid interface
    /// (e.g., dissolution, precipitation, surface catalysis)
    pub const SolidLiquidInterface = struct {
        solid: MoleculeId,
        liquid_solvent: MoleculeId,
    };

    liquid_bulk: LiquidBulk,

    liquid_interface: LiquidInterface,

    gas_liquid_interface: GasLiquidInterface,

    solid_liquid_interface: SolidLiquidInterface,

    /// Gas phase reaction
    gas_bulk,
};

/// Extended reaction definition for multi-phase systems
pub const MultiPhaseReaction = struct {
    /// Base reaction parameters (kinetics, modulators, pH effects, etc.)
    base: Reaction,

    /// Where this reaction occurs
    locus: ReactionLocus,

    /// For interfacial reactions: rate depends on interfacial area
    /// Units: mol/(m²·s) instead of mol/(L·s)
    /// If null, uses volumetric rate from base kinetics
    area_rate_constant: ?f64 = null,

    /// Mass transfer coefficient for interfacial reactions (m/s)
    /// Limits rate when diffusion to interface is slow
    mass_transfer_coeff: ?f64 = null,

    /// Get substrates from base reaction
    pub inline fn getSubstrates(this: *const MultiPhaseReaction) []const Stoich {
        return this.base.getSubstrates();
    }

    /// Get products from base reaction
    pub inline fn getProducts(this: *const MultiPhaseReaction) []const Stoich {
        return this.base.getProducts();
    }
};

/// State vector for ODE integration across all phases.
pub const PhaseState = struct {
    /// Moles in each liquid phase [phase_index][molecule_index]
    liquid_moles: [][]f64 = &.{},
    /// Moles for ALL possible solid phases [molecule_index]
    solid_moles: []f64 = &.{},
    /// Moles in gas phase [molecule_index]
    gas_moles: []f64 = &.{},

    /// Cached phase volumes for rate calculations (liters)
    liquid_volumes: []f64 = &.{},
    /// Cached solid surface areas for ALL possible solids (m²)
    solid_surface_areas: []f64 = &.{},
    /// Cached gas volume (liters)
    gas_volume: f64 = 0.0,

    /// Cached pH values per liquid phase (avoids recalculation)
    liquid_ph: []f64 = &.{},

    /// Total internal heat energy (Joules)
    heat_energy: f64 = 0.0,
    /// Container baseline heat capacity (J/K)
    container_heat_capacity: f64 = 100.0,

    /// Cached heat capacity (J/K) updated each RK stage
    heat_capacity: f64 = 0.0,
    /// Area of the liquid-gas interface (m^2).
    gas_contact_area: f64 = 0.0,
    /// Cached thermodynamic temperature (K) updated each RK stage (T = E/Cp)
    temperature: f64 = 298.15,
    /// Pressure (Pa)
    pressure: f64 = 101325.0,
    /// Normalized stirring intensity factor [0.0, 1.0].
    /// 0.0 = Unstirred (purely diffusion limited).
    /// 0.2 = Gentle stirring (e.g., occasional swirling by hand).
    /// 0.5 = Moderate stirring (e.g., standard magnetic stirrer at medium RPM).
    /// 0.8 = Vigorous stirring (e.g., high-speed overhead mechanical stirrer).
    /// 1.0 = Extreme agitation (e.g., high-shear homogenizer or sonication).
    stirring: f64 = 1.0,

    pub inline fn init(allocator: std.mem.Allocator, contents: *const Contents, max_volume: f64) !PhaseState {
        const n_liquids = contents.liquids.items.len;
        const n_molecules = MoleculeId.count();

        var state: PhaseState = .{};

        // Allocate liquid phase arrays
        if (n_liquids > 0) {
            state.liquid_moles = try allocator.alloc([]f64, n_liquids);
            state.liquid_volumes = try allocator.alloc(f64, n_liquids);
            state.liquid_ph = try allocator.alloc(f64, n_liquids);

            for (0..n_liquids) |i| {
                state.liquid_moles[i] = try allocator.alloc(f64, n_molecules);
                @memcpy(state.liquid_moles[i], contents.liquids.items[i].solutes);

                state.liquid_volumes[i] = contents.liquids.items[i].getVolume();
                state.liquid_ph[i] = helpers.computePhasePH(contents.liquids.items[i].solutes, state.liquid_volumes[i]);
            }
        }

        state.solid_moles = try allocator.alloc(f64, n_molecules);
        @memset(state.solid_moles, 0.0);
        state.solid_surface_areas = try allocator.alloc(f64, n_molecules);
        @memset(state.solid_surface_areas, 0.0);

        for (contents.solids.items) |solid| {
            const idx = @intFromEnum(solid.molecule);

            state.solid_moles[idx] = solid.moles;
            state.solid_surface_areas[idx] = solid.getSurfaceArea();
        }

        // Copy gas phase
        state.gas_moles = try allocator.alloc(f64, n_molecules);
        @memcpy(state.gas_moles, contents.gas.molecules);

        // Initialize thermodynamic values
        state.heat_energy = contents.heat_energy;
        state.container_heat_capacity = contents.container_heat_capacity;
        state.heat_capacity = contents.getHeatCapacity();
        state.gas_contact_area = contents.gas_contact_area;
        state.temperature = @as(f64, @floatCast(contents.getTemperature()));

        state.gas_volume = contents.getGasVolume(max_volume);
        state.pressure = contents.pressure;
        state.stirring = contents.stirring;

        return state;
    }

    pub inline fn deinit(this: *PhaseState, allocator: std.mem.Allocator) void {
        for (this.liquid_moles) |phase| {
            allocator.free(phase);
        }

        allocator.free(this.liquid_moles);
        allocator.free(this.liquid_volumes);
        allocator.free(this.liquid_ph);

        allocator.free(this.solid_moles);
        allocator.free(this.solid_surface_areas);

        allocator.free(this.gas_moles);
    }

    /// Apply state back to Contents after integration.
    pub inline fn applyTo(this: *const PhaseState, allocator: std.mem.Allocator, contents: *Contents) !void {
        for (contents.liquids.items, 0..) |*phase, i| {
            @memcpy(phase.solutes, this.liquid_moles[i]);
        }

        // Update existing solids, or dynamically create newly precipitated ones
        for (0..MoleculeId.count()) |i| {
            const new_moles = this.solid_moles[i];
            const mol: MoleculeId = @enumFromInt(i);

            var found = false;

            for (contents.solids.items) |*solid| {
                if (solid.molecule == mol) {
                    solid.moles = new_moles;
                    found = true;

                    break;
                }
            }

            if (!found and new_moles > constants.MIN_MOLES) {
                try contents.addSolidRaw(allocator, mol, new_moles, 1e-6, true);
            }
        }

        @memcpy(contents.gas.molecules, this.gas_moles);
        // Push the final heat energy back to the physical simulation
        contents.heat_energy = this.heat_energy;
    }

    /// Get concentration (mol/L) of molecule in a liquid phase
    pub inline fn getLiquidConcentration(this: *const PhaseState, phase_idx: usize, mol: MoleculeId) f64 {
        const vol = this.liquid_volumes[phase_idx];

        if (vol < constants.MIN_VOLUME) {
            return 0.0;
        }

        const conc = this.liquid_moles[phase_idx][@intFromEnum(mol)] / vol;

        return @min(conc, 100.0);
    }

    /// Get partial pressure (Pa) of molecule in gas phase
    pub inline fn getPartialPressure(this: *const PhaseState, mol: MoleculeId) f64 {
        var total_moles: f64 = 0;

        for (this.gas_moles) |n| {
            total_moles += n;
        }

        if (constants.isNegligible(total_moles)) {
            return 0.0;
        }

        const mole_fraction = this.gas_moles[@intFromEnum(mol)] / total_moles;
        return mole_fraction * this.pressure;
    }

    /// Get total moles of molecule across all liquid phases
    pub inline fn getTotalLiquidMoles(this: *const PhaseState, mol: MoleculeId) f64 {
        var total: f64 = 0;
        const idx = @intFromEnum(mol);

        for (this.liquid_moles) |phase| {
            total += phase[idx];
        }

        return total;
    }

    /// Get gas concentration (mol/L) using ideal gas law
    pub inline fn getGasConcentration(this: *const PhaseState, mol: MoleculeId) f64 {
        if (this.gas_volume < constants.MIN_VOLUME) {
            return 0.0;
        }

        return this.gas_moles[@intFromEnum(mol)] / this.gas_volume;
    }

    /// Get cached pH for a liquid phase
    pub inline fn getPH(this: *const PhaseState, phase_idx: usize) f64 {
        if (phase_idx >= this.liquid_ph.len) {
            return 7.0;
        }

        return this.liquid_ph[phase_idx];
    }

    /// Update cached pH values after mole changes
    pub inline fn updatePH(this: *PhaseState) void {
        for (0..this.liquid_moles.len) |i| {
            this.liquid_ph[i] = helpers.computePhasePH(this.liquid_moles[i], this.liquid_volumes[i]);
        }
    }

    /// Get surface coverage θ for adsorbed species (0 to 1)
    /// Uses Langmuir isotherm: θ = K*C / (1 + K*C)
    pub inline fn getSurfaceCoverage(
        this: *const PhaseState,
        mol: MoleculeId,
        K_ads: f64,
        phase_idx: usize,
    ) f64 {
        const c = this.getLiquidConcentration(phase_idx, mol);

        return K_ads * c / (1.0 + K_ads * c);
    }
};

/// ODE solver for multi-phase reaction systems.
/// Manages its own memory arena to minimize per-step allocations.
pub const MultiPhaseODESolver = struct {
    /// Minimum timestep (seconds)
    min_dt: f64 = 1e-9,
    /// Maximum timestep (seconds)
    max_dt: f64 = 0.1,
    /// Error tolerance for adaptive stepping
    tolerance: f64 = 1e-6,

    /// Internal arena allocator for RK workspaces
    arena: std.heap.ArenaAllocator,

    // Working arrays for RK stages, sized per integrate() call
    k_liquid: [6][][]f64 = undefined,
    k_solid: [6][]f64 = undefined,
    k_gas: [6][]f64 = undefined,
    k_heat: [6]f64 = .{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    temp_state: PhaseState = undefined,
    trial_state: PhaseState = undefined,
    contents: *const Contents = undefined,

    const StepResult = struct {
        accepted: bool,
        next_dt: f64,
    };

    pub inline fn init(allocator: std.mem.Allocator) MultiPhaseODESolver {
        return .{
            .arena = std.heap.ArenaAllocator.init(allocator),
        };
    }

    pub inline fn deinit(this: *MultiPhaseODESolver) void {
        this.arena.deinit();
    }

    /// Prepare workspace arrays for the current phase topology
    inline fn prepareWorkspace(this: *MultiPhaseODESolver, contents: *const Contents, max_volume: f64) !void {
        _ = this.arena.reset(.retain_capacity);
        const alloc = this.arena.allocator();

        const n_liquids = contents.liquids.items.len;
        const n_mol = MoleculeId.count();

        this.contents = contents;

        inline for (0..6) |stage| {
            if (n_liquids > 0) {
                this.k_liquid[stage] = try alloc.alloc([]f64, n_liquids);

                for (0..n_liquids) |i| {
                    this.k_liquid[stage][i] = try alloc.alloc(f64, n_mol);
                }
            } else {
                this.k_liquid[stage] = &.{};
            }

            // k_solid is now sized for all molecules, identical to k_gas
            this.k_solid[stage] = try alloc.alloc(f64, n_mol);
            this.k_gas[stage] = try alloc.alloc(f64, n_mol);
        }

        this.temp_state = try PhaseState.init(alloc, contents, max_volume);
        this.trial_state = try PhaseState.init(alloc, contents, max_volume);
    }

    /// Integrate for a total time duration
    pub inline fn integrate(
        this: *MultiPhaseODESolver,
        allocator: std.mem.Allocator,
        reactions: []const MultiPhaseReaction,
        contents: *Contents,
        total_time: f64,
        max_volume: f64,
    ) error{OutOfMemory}!void {
        @setEvalBranchQuota(std.math.maxInt(u32));

        // Setup memory for exactly the phases present at THIS tick
        try this.prepareWorkspace(contents, max_volume);
        const arena_alloc = this.arena.allocator();

        var state = try PhaseState.init(arena_alloc, contents, max_volume);

        this.enforceEquilibria(reactions, &state);
        this.updateCachedQuantities(&state);
        state.updatePH();

        var t: f64 = 0;
        var dt: f64 = @min(this.min_dt * 100, total_time / 10);

        while (t < total_time) {
            const remaining = total_time - t;
            var step_dt = @min(dt, remaining);

            // Retry rejected adaptive steps until one is accepted.
            while (true) {
                const result = this.step(reactions, &state, step_dt);
                dt = result.next_dt;

                if (result.accepted) {
                    t += step_dt;

                    break;
                }

                step_dt = @min(dt, remaining);
            }

            // After each accepted ODE step, enforce fast equilibria algebraically.
            this.enforceEquilibria(reactions, &state);

            // Recompute cached thermodynamic quantities after algebraic projection.
            this.updateCachedQuantities(&state);
            state.updatePH();
        }

        // Apply state and dynamically create new solids if needed
        try state.applyTo(allocator, contents);
    }

    /// Computes the standard enthalpy of reaction (Joules/mol) using Hess's Law.
    /// This guarantees strict thermodynamic consistency across the simulation.
    inline fn getReactionEnthalpy(
        rxn: *const MultiPhaseReaction,
    ) f64 {
        var delta_H: f64 = 0.0;

        for (rxn.getProducts()) |p| {
            delta_H += p.coefficient * p.molecule.getDeltaHf();
        }

        for (rxn.getSubstrates()) |s| {
            delta_H -= s.coefficient * s.molecule.getDeltaHf();
        }

        return delta_H;
    }

    inline fn accumulateReactionHeat(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        extent_rate_mol_s: f64,
        k_heat: *f64,
    ) void {
        _ = this;
        const delta_H = getReactionEnthalpy(rxn);

        // Constant-pressure approximation: dQ/dt = -ΔH_rxn * dξ/dt
        k_heat.* += -delta_H * extent_rate_mol_s;
    }

    inline fn applyEquilibriumHeat(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        delta_extent_mol: f64,
        state: *PhaseState,
    ) void {
        _ = this;
        const delta_H = getReactionEnthalpy(rxn);

        // Algebraic equilibrium projection changes composition,
        // updating thermal energy consistently.
        state.heat_energy += -delta_H * delta_extent_mol;
    }

    /// Enforce all fast equilibria algebraically (no ODE, no stiffness).
    /// This projects the current state onto the equilibrium manifold.
    /// Called after each ODE step.
    inline fn enforceEquilibria(
        this: *MultiPhaseODESolver,
        reactions: []const MultiPhaseReaction,
        state: *PhaseState,
    ) void {
        // Iterate a few times to handle coupled equilibria
        for (0..5) |_| {
            var max_change: f64 = 0;

            for (reactions) |rxn| {
                const delta_extent = this.enforceOneEquilibrium(&rxn, state);

                max_change = @max(max_change, @abs(delta_extent));
            }

            // Converged
            if (max_change < 1e-12) {
                break;
            }
        }
    }

    /// Enforce a single equilibrium reaction algebraically.
    /// Returns the magnitude of the change (for convergence check).
    inline fn enforceOneEquilibrium(
        this: *MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        state: *PhaseState,
    ) f64 {
        var max_delta: f64 = 0.0;

        switch (rxn.locus) {
            .liquid_bulk => |loc| {
                for (0..state.liquid_moles.len) |phase_idx| {
                    if (this.isEffectivelyDryLiquidPhase(state, phase_idx, loc.solvent)) {
                        continue;
                    }

                    const vol = state.liquid_volumes[phase_idx];

                    if (vol < constants.MIN_VOLUME) {
                        continue;
                    }

                    switch (rxn.base.kinetics) {
                        .weak_electrolyte => |we| {
                            const d = this.enforceWeakElectrolyteEquilibrium(we, state, phase_idx);
                            this.applyEquilibriumHeat(rxn, d, state);

                            if (@abs(d) > max_delta) {
                                max_delta = @abs(d);
                            }
                        },
                        .dissociation => |dis| {
                            const d = this.enforceDissociationEquilibrium(dis, state, phase_idx);
                            this.applyEquilibriumHeat(rxn, d, state);

                            if (@abs(d) > max_delta) {
                                max_delta = @abs(d);
                            }
                        },
                        .strong_electrolyte => {
                            const d = this.enforceCompleteDissociation(state, phase_idx, rxn.getSubstrates(), rxn.getProducts());
                            this.applyEquilibriumHeat(rxn, d, state);

                            if (@abs(d) > max_delta) {
                                max_delta = @abs(d);
                            }
                        },
                        else => {},
                    }
                }
            },
            .solid_liquid_interface => |loc| {
                for (0..state.liquid_moles.len) |liquid_idx| {
                    switch (rxn.base.kinetics) {
                        .strong_electrolyte => |se| {
                            const solid_idx = @intFromEnum(se.solid.molecule);
                            const d = this.enforceStrongElectrolyteSolubility(se, state, liquid_idx, solid_idx, loc.liquid_solvent);

                            this.applyEquilibriumHeat(rxn, d, state);

                            if (@abs(d) > max_delta) {
                                max_delta = @abs(d);
                            }
                        },
                        .dissolution => |dsl| {
                            const solid_idx = @intFromEnum(dsl.solid.molecule);
                            const d = this.enforceDissolutionSolubility(dsl, state, liquid_idx, solid_idx, loc.liquid_solvent, rxn);

                            this.applyEquilibriumHeat(rxn, d, state);

                            if (@abs(d) > max_delta) {
                                max_delta = @abs(d);
                            }
                        },
                        else => {},
                    }
                }
            },
            else => {},
        }

        return max_delta;
    }

    /// Enforce weak electrolyte equilibrium: HA ⇌ H⁺ + A⁻
    /// Solves the quadratic equation for extent of dissociation.
    inline fn enforceWeakElectrolyteEquilibrium(
        this: *const MultiPhaseODESolver,
        we: KineticsType.WeakElectrolyte,
        state: *PhaseState,
        phase_idx: usize,
    ) f64 {
        _ = this;

        const vol = state.liquid_volumes[phase_idx];

        if (vol < constants.MIN_VOLUME) {
            return 0;
        }

        // Ka from pKa (no temperature correction for simplicity)
        const Ka = std.math.pow(f64, 10.0, -we.pKa);

        const parent_idx = @intFromEnum(we.parent.molecule);
        const conjugate_idx = @intFromEnum(we.products[0].molecule);
        const proton_idx = @intFromEnum(we.products[1].molecule);

        // Current moles
        const n_parent = state.liquid_moles[phase_idx][parent_idx];
        const n_conjugate = state.liquid_moles[phase_idx][conjugate_idx];

        // Total analytical concentration of the acid/base (conserved quantity)
        const n_total = n_parent + n_conjugate;

        if (n_total < constants.MIN_MOLES) {
            return 0;
        }

        const C_total = n_total / vol;

        if (we.is_base) {
            // B + H₂O ⇌ BH⁺ + OH⁻
            const Kw: f64 = 1e-14;
            const Kb = Kw / Ka;

            // Kb = x² / (C_total - x) where x = [BH⁺] = [OH⁻] from this base
            // x² + Kb*x - Kb*C_total = 0
            const disc = Kb * Kb + 4.0 * Kb * C_total;
            const x = (-Kb + @sqrt(@max(0, disc))) / 2.0;
            const x_clamped = std.math.clamp(x, 0, C_total);

            // conjugate = BH⁺, parent = B
            const new_n_conjugate = x_clamped * vol;
            const new_n_parent = (C_total - x_clamped) * vol;
            const delta_extent = (new_n_conjugate - n_conjugate) / we.products[0].coefficient;

            state.liquid_moles[phase_idx][parent_idx] = new_n_parent;
            state.liquid_moles[phase_idx][conjugate_idx] = new_n_conjugate;

            // OH⁻ is handled through pH recalculation elsewhere.
            return delta_extent;
        } else {
            // HA ⇌ H⁺ + A⁻
            // Include the current proton pool from other sources as a background term.
            const current_h = state.liquid_moles[phase_idx][proton_idx];
            const h_from_others = @max(0, current_h - n_conjugate);
            const h_other_conc = h_from_others / vol;

            // Ka = (h_other_conc + x) * x / (C_total - x)
            // where x = equilibrium [A⁻] = H⁺ contributed by this acid
            // x² + (h_other_conc + Ka) * x - Ka * C_total = 0
            const b = h_other_conc + Ka;
            const disc = b * b + 4.0 * Ka * C_total;
            const x = (-b + @sqrt(@max(0, disc))) / 2.0;
            const x_clamped = std.math.clamp(x, 0, C_total);

            const new_n_conjugate = x_clamped * vol;
            const new_n_parent = (C_total - x_clamped) * vol;
            const delta_conjugate = new_n_conjugate - n_conjugate;
            const delta_extent = delta_conjugate / we.products[0].coefficient;

            state.liquid_moles[phase_idx][parent_idx] = new_n_parent;
            state.liquid_moles[phase_idx][conjugate_idx] = new_n_conjugate;
            // H⁺ changes by the same amount as conjugate
            state.liquid_moles[phase_idx][proton_idx] = @max(0, current_h + delta_conjugate);

            return delta_extent;
        }
    }

    /// Enforce dissociation equilibrium: AB ⇌ A⁺ + B⁻
    inline fn enforceDissociationEquilibrium(
        this: *const MultiPhaseODESolver,
        dis: KineticsType.Dissociation,
        state: *PhaseState,
        phase_idx: usize,
    ) f64 {
        _ = this;

        const vol = state.liquid_volumes[phase_idx];

        if (vol < constants.MIN_VOLUME) {
            return 0;
        }

        const parent_idx = @intFromEnum(dis.parent.molecule);
        const n_parent = state.liquid_moles[phase_idx][parent_idx];

        // For simple 1:1 dissociation AB ⇌ A⁺ + B⁻
        // Ka = [A⁺][B⁻] / [AB]
        // This algebraic projection assumes the parent is the dominant source
        // of the two ions in the current phase.

        if (dis.ions.len != 2) {
            // Complex dissociation - skip algebraic enforcement, let ODE handle it.
            return 0;
        }

        const ion0_idx = @intFromEnum(dis.ions[0].molecule);
        const ion1_idx = @intFromEnum(dis.ions[1].molecule);

        const n_ion0 = state.liquid_moles[phase_idx][ion0_idx];
        const n_ion1 = state.liquid_moles[phase_idx][ion1_idx];

        const n_dissociated = @min(n_ion0, n_ion1);
        const n_total = n_parent + n_dissociated;

        if (n_total < constants.MIN_MOLES) {
            return 0;
        }

        const C_total = n_total / vol;

        // Ka = x² / (C_total - x)
        // x² + Ka*x - Ka*C_total = 0
        const disc = dis.Ka * dis.Ka + 4.0 * dis.Ka * C_total;
        const x = (-dis.Ka + @sqrt(@max(0, disc))) / 2.0;
        const x_clamped = std.math.clamp(x, 0, C_total);

        const new_n_parent = (C_total - x_clamped) * vol;
        const delta_extent = x_clamped * vol - n_dissociated;

        state.liquid_moles[phase_idx][parent_idx] = new_n_parent;
        state.liquid_moles[phase_idx][ion0_idx] = @max(0, n_ion0 + delta_extent);
        state.liquid_moles[phase_idx][ion1_idx] = @max(0, n_ion1 + delta_extent);

        return delta_extent;
    }

    /// Enforce complete dissociation for strong electrolyte in liquid bulk
    inline fn enforceCompleteDissociation(
        this: *const MultiPhaseODESolver,
        state: *PhaseState,
        phase_idx: usize,
        substrates: []const Stoich,
        products: []const Stoich,
    ) f64 {
        _ = this;

        // Move all undissociated parent moles into products
        var total_consumed: f64 = 0;

        for (substrates) |s| {
            const idx = @intFromEnum(s.molecule);
            const n = state.liquid_moles[phase_idx][idx];

            if (n > constants.MIN_MOLES) {
                const consumed = n; // Complete dissociation

                state.liquid_moles[phase_idx][idx] = 0;
                total_consumed += consumed / s.coefficient;
            }
        }

        if (total_consumed < constants.MIN_MOLES) {
            return 0;
        }

        for (products) |p| {
            const idx = @intFromEnum(p.molecule);

            state.liquid_moles[phase_idx][idx] += total_consumed * p.coefficient;
        }

        return total_consumed;
    }

    inline fn getPureSolventVolume(
        this: *const MultiPhaseODESolver,
        state: *const PhaseState,
        phase_idx: usize,
        solvent: MoleculeId,
    ) f64 {
        _ = this;

        const n_solvent = state.liquid_moles[phase_idx][@intFromEnum(solvent)];

        if (n_solvent < constants.MIN_MOLES) {
            return 0.0;
        }

        return n_solvent * solvent.getWeight() / solvent.getDensity() / 1000.0;
    }

    inline fn isEffectivelyDryLiquidPhase(
        this: *const MultiPhaseODESolver,
        state: *const PhaseState,
        phase_idx: usize,
        solvent: MoleculeId,
    ) bool {
        const n_solvent = state.liquid_moles[phase_idx][@intFromEnum(solvent)];

        if (constants.isNegligible(n_solvent)) {
            return true;
        }

        const solvent_volume = this.getPureSolventVolume(state, phase_idx, solvent);

        return solvent_volume < constants.MIN_VOLUME;
    }

    inline fn enforceDissolutionSolubility(
        this: *const MultiPhaseODESolver,
        dsl: KineticsType.Dissolution,
        state: *PhaseState,
        liquid_idx: usize,
        solid_idx: usize,
        solvent: MoleculeId,
        rxn: *const MultiPhaseReaction,
    ) f64 {
        if (this.isEffectivelyDryLiquidPhase(state, liquid_idx, solvent)) {
            return this.enforceDrySaltPrecipitation(state, liquid_idx, solid_idx, &dsl.ions, dsl.solid.coefficient);
        }

        const vol = state.liquid_volumes[liquid_idx];

        if (vol < constants.MIN_VOLUME) {
            return 0.0;
        }

        const delta_H = getReactionEnthalpy(rxn);
        const safe_temp = @max(1.0, state.temperature);
        const Ksp_T = @max(1e-300, dsl.Ksp * @exp(-delta_H / constants.R * (1.0 / safe_temp - 1.0 / dsl.T_ref)));

        const cat_idx = @intFromEnum(dsl.ions[0].molecule);
        const an_idx = @intFromEnum(dsl.ions[1].molecule);

        const n_cat = state.liquid_moles[liquid_idx][cat_idx];
        const n_an = state.liquid_moles[liquid_idx][an_idx];

        const c_cat = n_cat / vol;
        const c_an = n_an / vol;

        const Q = c_cat * c_an;

        // 1% tolerance
        if (Q <= Ksp_T * 1.01) {
            return 0.0;
        }

        const b = -(c_cat + c_an);
        const c = Q - Ksp_T;
        const disc = b * b - 4.0 * c;

        if (disc < 0) {
            return 0.0;
        }

        const x_conc = (-b - @sqrt(disc)) / 2.0;
        const precip_moles = x_conc * vol;

        if (precip_moles <= constants.MIN_MOLES) {
            return 0.0;
        }

        state.liquid_moles[liquid_idx][cat_idx] = @max(0.0, n_cat - precip_moles * dsl.ions[0].coefficient);
        state.liquid_moles[liquid_idx][an_idx] = @max(0.0, n_an - precip_moles * dsl.ions[1].coefficient);
        state.solid_moles[solid_idx] += precip_moles * dsl.solid.coefficient;

        return -precip_moles;
    }

    inline fn enforceDrySaltPrecipitation(
        this: *const MultiPhaseODESolver,
        state: *PhaseState,
        liquid_idx: usize,
        solid_idx: usize,
        ions: []const Stoich,
        solid_coefficient: f64,
    ) f64 {
        _ = this;

        if (ions.len == 0) {
            return 0;
        }

        var precip_extent = std.math.inf(f64);

        for (ions) |ion| {
            const n_ion = state.liquid_moles[liquid_idx][@intFromEnum(ion.molecule)];
            precip_extent = @min(precip_extent, n_ion / ion.coefficient);
        }

        if (!std.math.isFinite(precip_extent) or precip_extent < constants.MIN_MOLES) {
            return 0;
        }

        for (ions) |ion| {
            const idx = @intFromEnum(ion.molecule);
            state.liquid_moles[liquid_idx][idx] = @max(
                0.0,
                state.liquid_moles[liquid_idx][idx] - precip_extent * ion.coefficient,
            );
        }

        state.solid_moles[solid_idx] += precip_extent * solid_coefficient;

        // Negative extent means reverse direction: ions -> solid.
        return -precip_extent;
    }

    inline fn enforceStrongElectrolyteSolubility(
        this: *const MultiPhaseODESolver,
        se: KineticsType.StrongElectrolyte,
        state: *PhaseState,
        liquid_idx: usize,
        solid_idx: usize,
        solvent: MoleculeId,
    ) f64 {
        // In a dry phase, collapse dissolved ions fully into the solid.
        if (this.isEffectivelyDryLiquidPhase(state, liquid_idx, solvent)) {
            return this.enforceDrySaltPrecipitation(
                state,
                liquid_idx,
                solid_idx,
                &se.ions,
                se.solid.coefficient,
            );
        }

        const solvent_volume = this.getPureSolventVolume(state, liquid_idx, solvent);

        if (solvent_volume <= 0.0) {
            return this.enforceDrySaltPrecipitation(
                state,
                liquid_idx,
                solid_idx,
                &se.ions,
                se.solid.coefficient,
            );
        }

        const cation_idx = @intFromEnum(se.ions[0].molecule);
        const anion_idx = @intFromEnum(se.ions[1].molecule);

        const nu_cat = se.ions[0].coefficient;
        const nu_an = se.ions[1].coefficient;

        const n_cat = state.liquid_moles[liquid_idx][cation_idx];
        const n_an = state.liquid_moles[liquid_idx][anion_idx];

        const dissolved_formula_units = @min(n_cat / nu_cat, n_an / nu_an);
        const max_dissolved_formula_units = se.saturation * solvent_volume;

        if (dissolved_formula_units <= max_dissolved_formula_units + constants.MIN_MOLES) {
            return 0;
        }

        // Project oversaturated solution back to the solubility ceiling.
        const precip_extent = dissolved_formula_units - max_dissolved_formula_units;

        state.liquid_moles[liquid_idx][cation_idx] = @max(
            0.0,
            n_cat - precip_extent * nu_cat,
        );
        state.liquid_moles[liquid_idx][anion_idx] = @max(
            0.0,
            n_an - precip_extent * nu_an,
        );
        state.solid_moles[solid_idx] += precip_extent * se.solid.coefficient;

        return -precip_extent;
    }

    /// Perform one adaptive RK4(5) step
    pub inline fn step(
        this: *MultiPhaseODESolver,
        reactions: []const MultiPhaseReaction,
        state: *PhaseState,
        dt: f64,
    ) StepResult {
        const n_liquids = state.liquid_moles.len;
        const n_solids = state.solid_moles.len;
        const n_mol = MoleculeId.count();

        // Stage 1: k1 = f(t, y)
        this.computeRates(reactions, state, this.k_liquid[0], this.k_solid[0], this.k_gas[0], &this.k_heat[0]);

        // Stage 2: k2 = f(t + dt/5, y + dt/5 * k1)
        this.advanceTemp(state, &.{0}, &.{1.0 / 5.0}, dt);
        this.computeRates(reactions, &this.temp_state, this.k_liquid[1], this.k_solid[1], this.k_gas[1], &this.k_heat[1]);

        // Stage 3
        this.advanceTemp(state, &.{ 0, 1 }, &.{ 3.0 / 40.0, 9.0 / 40.0 }, dt);
        this.computeRates(reactions, &this.temp_state, this.k_liquid[2], this.k_solid[2], this.k_gas[2], &this.k_heat[2]);

        // Stage 4
        this.advanceTemp(state, &.{ 0, 1, 2 }, &.{ 3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0 }, dt);
        this.computeRates(reactions, &this.temp_state, this.k_liquid[3], this.k_solid[3], this.k_gas[3], &this.k_heat[3]);

        // Stage 5
        this.advanceTemp(state, &.{ 0, 1, 2, 3 }, &.{ -11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0 }, dt);
        this.computeRates(reactions, &this.temp_state, this.k_liquid[4], this.k_solid[4], this.k_gas[4], &this.k_heat[4]);

        // Stage 6
        this.advanceTemp(state, &.{ 0, 1, 2, 3, 4 }, &.{ 1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0 }, dt);
        this.computeRates(reactions, &this.temp_state, this.k_liquid[5], this.k_solid[5], this.k_gas[5], &this.k_heat[5]);

        const w5 = [_]f64{ 37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0 };
        const w4 = [_]f64{ 2825.0 / 27648.0, 0.0, 18575.0 / 48384.0, 13525.0 / 55296.0, 277.0 / 14336.0, 1.0 / 4.0 };

        this.copyState(&this.trial_state, state);

        var max_error: f64 = 0;

        // Candidate liquid update
        for (0..n_liquids) |p| {
            for (0..n_mol) |m| {
                var sum5: f64 = 0;
                var sum4: f64 = 0;

                inline for (0..6) |s| {
                    sum5 += w5[s] * this.k_liquid[s][p][m];
                    sum4 += w4[s] * this.k_liquid[s][p][m];
                }

                const new_val = @max(0.0, state.liquid_moles[p][m] + dt * sum5);
                const err = @abs(dt * (sum5 - sum4));
                const scale = @max(1e-15, @abs(new_val));

                max_error = @max(max_error, err / scale);
                this.trial_state.liquid_moles[p][m] = new_val;
            }
        }

        // Candidate solid update
        for (0..n_solids) |s_idx| {
            var sum5: f64 = 0;
            var sum4: f64 = 0;

            inline for (0..6) |s| {
                sum5 += w5[s] * this.k_solid[s][s_idx];
                sum4 += w4[s] * this.k_solid[s][s_idx];
            }

            const new_val = @max(0.0, state.solid_moles[s_idx] + dt * sum5);
            const err = @abs(dt * (sum5 - sum4));
            const scale = @max(1e-15, @abs(new_val));

            max_error = @max(max_error, err / scale);
            this.trial_state.solid_moles[s_idx] = new_val;
        }

        // Candidate gas update
        for (0..n_mol) |m| {
            var sum5: f64 = 0;
            var sum4: f64 = 0;

            inline for (0..6) |s| {
                sum5 += w5[s] * this.k_gas[s][m];
                sum4 += w4[s] * this.k_gas[s][m];
            }

            const new_val = @max(0.0, state.gas_moles[m] + dt * sum5);
            const err = @abs(dt * (sum5 - sum4));
            const scale = @max(1e-15, @abs(new_val));

            max_error = @max(max_error, err / scale);
            this.trial_state.gas_moles[m] = new_val;
        }

        // Candidate heat update
        var heat_sum5: f64 = 0.0;
        var heat_sum4: f64 = 0.0;

        inline for (0..6) |s| {
            heat_sum5 += w5[s] * this.k_heat[s];
            heat_sum4 += w4[s] * this.k_heat[s];
        }

        const new_heat = state.heat_energy + dt * heat_sum5;
        const heat_err = @abs(dt * (heat_sum5 - heat_sum4));
        const heat_scale = @max(1.0, @abs(new_heat));

        max_error = @max(max_error, heat_err / heat_scale);
        this.trial_state.heat_energy = new_heat;

        this.updateCachedQuantities(&this.trial_state);
        this.trial_state.updatePH();

        const next_dt = this.adaptStep(dt, max_error);
        const accept = max_error <= this.tolerance or dt <= this.min_dt * (1.0 + 1e-12);

        if (!accept) {
            return .{
                .accepted = false,
                .next_dt = next_dt,
            };
        }

        this.copyState(state, &this.trial_state);

        return .{
            .accepted = true,
            .next_dt = next_dt,
        };
    }

    inline fn copyState(
        this: *const MultiPhaseODESolver,
        dst: *PhaseState,
        src: *const PhaseState,
    ) void {
        _ = this;

        for (dst.liquid_moles, 0..) |phase, i| {
            @memcpy(phase, src.liquid_moles[i]);
        }

        if (dst.solid_moles.len > 0) {
            @memcpy(dst.solid_moles, src.solid_moles);
        }

        @memcpy(dst.gas_moles, src.gas_moles);

        if (dst.liquid_volumes.len > 0) {
            @memcpy(dst.liquid_volumes, src.liquid_volumes);
            @memcpy(dst.liquid_ph, src.liquid_ph);
        }

        if (dst.solid_surface_areas.len > 0) {
            @memcpy(dst.solid_surface_areas, src.solid_surface_areas);
        }

        dst.gas_volume = src.gas_volume;
        dst.heat_energy = src.heat_energy;
        dst.container_heat_capacity = src.container_heat_capacity;
        dst.heat_capacity = src.heat_capacity;
        dst.temperature = src.temperature;
        dst.pressure = src.pressure;
        dst.stirring = src.stirring;
    }

    /// Compute all reaction rates and accumulate into rate arrays
    inline fn computeRates(
        this: *MultiPhaseODESolver,
        reactions: []const MultiPhaseReaction,
        state: *const PhaseState,
        k_liquid: [][]f64,
        k_solid: []f64,
        k_gas: []f64,
        k_heat: *f64,
    ) void {
        // Zero rate arrays
        for (k_liquid) |phase| {
            @memset(phase, 0.0);
        }

        @memset(k_solid, 0.0);
        @memset(k_gas, 0.0);
        k_heat.* = 0.0;

        for (reactions) |rxn| {
            this.computeSingleReaction(&rxn, state, k_liquid, k_solid, k_gas, k_heat);
        }
    }

    /// Check if a reaction should be handled algebraically (not via ODE rates)
    inline fn isAlgebraicEquilibrium(rxn: *const MultiPhaseReaction) bool {
        switch (rxn.locus) {
            .liquid_bulk => {
                switch (rxn.base.kinetics) {
                    .weak_electrolyte, .dissociation, .strong_electrolyte => return true,
                    else => return false,
                }
            },
            else => return false,
        }
    }

    inline fn computeSingleReaction(
        this: *MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        state: *const PhaseState,
        k_liquid: [][]f64,
        k_solid: []f64,
        k_gas: []f64,
        k_heat: *f64,
    ) void {
        if (isAlgebraicEquilibrium(rxn)) {
            return;
        }

        const substrates = rxn.getSubstrates();
        const products = rxn.getProducts();

        switch (rxn.locus) {
            .liquid_bulk => |loc| {
                for (0..state.liquid_moles.len) |phase_idx| {
                    if (this.isEffectivelyDryLiquidPhase(state, phase_idx, loc.solvent)) {
                        continue;
                    }

                    const vol = state.liquid_volumes[phase_idx];

                    if (vol < constants.MIN_VOLUME) {
                        continue;
                    }

                    const v = this.computeVolumetricRate(rxn, state, phase_idx);
                    const molar_rate = v * vol;

                    for (substrates) |s| {
                        k_liquid[phase_idx][@intFromEnum(s.molecule)] -= molar_rate * s.coefficient;
                    }

                    for (products) |p| {
                        k_liquid[phase_idx][@intFromEnum(p.molecule)] += molar_rate * p.coefficient;
                    }

                    this.accumulateReactionHeat(rxn, molar_rate, k_heat);
                }
            },
            .gas_liquid_interface => |loc| {
                for (0..state.liquid_moles.len) |phase_idx| {
                    if (this.isEffectivelyDryLiquidPhase(state, phase_idx, loc.liquid_solvent)) {
                        continue;
                    }

                    const area = state.gas_contact_area;
                    const v = this.computeInterfacialRate(rxn, state, phase_idx, area, substrates);

                    for (substrates) |s| {
                        if (s.molecule.isGas()) {
                            k_gas[@intFromEnum(s.molecule)] -= v * s.coefficient;
                        } else {
                            k_liquid[phase_idx][@intFromEnum(s.molecule)] -= v * s.coefficient;
                        }
                    }

                    for (products) |p| {
                        if (p.molecule.isGas()) {
                            k_gas[@intFromEnum(p.molecule)] += v * p.coefficient;
                        } else {
                            k_liquid[phase_idx][@intFromEnum(p.molecule)] += v * p.coefficient;
                        }
                    }

                    this.accumulateReactionHeat(rxn, v, k_heat);
                }
            },
            .solid_liquid_interface => |loc| {
                for (0..state.liquid_moles.len) |liquid_idx| {
                    if (this.isEffectivelyDryLiquidPhase(state, liquid_idx, loc.liquid_solvent)) {
                        continue;
                    }

                    const liquid_vol_m3 = state.liquid_volumes[liquid_idx] / 1000.0;
                    const max_wetted_area = liquid_vol_m3 / 1e-5;

                    switch (rxn.base.kinetics) {
                        .strong_electrolyte => |se| {
                            const solid_idx = @intFromEnum(se.solid.molecule);
                            const total_area = state.solid_surface_areas[solid_idx];

                            const area = @min(total_area, max_wetted_area);

                            if (area < constants.MIN_AREA) {
                                continue;
                            }

                            const rate = this.computeStrongElectrolyteDissolution(se, state, liquid_idx, area);

                            k_solid[solid_idx] -= rate * se.solid.coefficient;

                            for (se.ions) |ion| {
                                k_liquid[liquid_idx][@intFromEnum(ion.molecule)] += rate * ion.coefficient;
                            }

                            this.accumulateReactionHeat(rxn, rate, k_heat);
                        },
                        .dissolution => |dsl| {
                            const solid_idx = @intFromEnum(dsl.solid.molecule);
                            const total_area = state.solid_surface_areas[solid_idx];

                            const area = @min(total_area, max_wetted_area);

                            const rate = this.computeDissolutionRate(rxn, dsl, state, liquid_idx, solid_idx, area);

                            k_solid[solid_idx] -= rate * dsl.solid.coefficient;

                            for (dsl.ions) |ion| {
                                k_liquid[liquid_idx][@intFromEnum(ion.molecule)] += rate * ion.coefficient;
                            }

                            this.accumulateReactionHeat(rxn, rate, k_heat);
                        },
                        else => {
                            const solid_idx = @intFromEnum(loc.solid);
                            const total_area = state.solid_surface_areas[solid_idx];

                            const area = @min(total_area, max_wetted_area);

                            if (area < constants.MIN_AREA) {
                                continue;
                            }

                            const v = this.computeSolidLiquidRate(rxn, state, liquid_idx, solid_idx, area, substrates);

                            for (substrates) |s| {
                                if (s.molecule == loc.solid) {
                                    k_solid[solid_idx] -= v * s.coefficient;
                                } else {
                                    k_liquid[liquid_idx][@intFromEnum(s.molecule)] -= v * s.coefficient;
                                }
                            }

                            for (products) |p| {
                                k_liquid[liquid_idx][@intFromEnum(p.molecule)] += v * p.coefficient;
                            }

                            this.accumulateReactionHeat(rxn, v, k_heat);
                        },
                    }
                }
            },
            .gas_bulk => {
                if (state.gas_volume < constants.MIN_VOLUME) {
                    return;
                }

                const v = this.computeGasPhaseRate(rxn, state, substrates);

                for (substrates) |s| {
                    k_gas[@intFromEnum(s.molecule)] -= v * s.coefficient;
                }

                for (products) |p| {
                    k_gas[@intFromEnum(p.molecule)] += v * p.coefficient;
                }

                this.accumulateReactionHeat(rxn, v, k_heat);
            },
            .liquid_interface => |loc| {
                for (0..state.liquid_moles.len) |phase_a| {
                    if (this.isEffectivelyDryLiquidPhase(state, phase_a, loc.phase_a_solvent)) {
                        continue;
                    }

                    for (0..state.liquid_moles.len) |phase_b| {
                        if (phase_a == phase_b) {
                            continue;
                        }

                        if (this.isEffectivelyDryLiquidPhase(state, phase_b, loc.phase_b_solvent)) {
                            continue;
                        }

                        const vol_a = state.liquid_volumes[phase_a];
                        const vol_b = state.liquid_volumes[phase_b];
                        const min_vol = @min(vol_a, vol_b);
                        const area = 0.05 * std.math.pow(f64, min_vol * 1000.0, 2.0 / 3.0);

                        const v = this.computeLiquidInterfaceRate(rxn, state, phase_a, phase_b, area, substrates);

                        for (substrates) |s| {
                            const mol_idx = @intFromEnum(s.molecule);
                            const delta_a = @abs(s.molecule.getHildebrand() - this.contents.liquids.items[phase_a].getHildebrand());
                            const delta_b = @abs(s.molecule.getHildebrand() - this.contents.liquids.items[phase_b].getHildebrand());

                            if (delta_a < delta_b) {
                                k_liquid[phase_a][mol_idx] -= v * s.coefficient;
                            } else {
                                k_liquid[phase_b][mol_idx] -= v * s.coefficient;
                            }
                        }

                        for (products) |p| {
                            const mol_idx = @intFromEnum(p.molecule);
                            const delta_a = @abs(p.molecule.getHildebrand() - this.contents.liquids.items[phase_a].getHildebrand());
                            const delta_b = @abs(p.molecule.getHildebrand() - this.contents.liquids.items[phase_b].getHildebrand());

                            if (delta_a < delta_b) {
                                k_liquid[phase_a][mol_idx] += v * p.coefficient;
                            } else {
                                k_liquid[phase_b][mol_idx] += v * p.coefficient;
                            }
                        }

                        this.accumulateReactionHeat(rxn, v, k_heat);
                    }
                }
            },
        }
    }

    /// Compute volumetric reaction rate (mol/L/s) for homogeneous liquid reaction
    /// Compute volumetric reaction rate (mol/L/s) for homogeneous liquid reaction
    inline fn computeVolumetricRate(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        state: *const PhaseState,
        phase_idx: usize,
    ) f64 {
        @setEvalBranchQuota(std.math.maxInt(u32));

        const substrates = rxn.getSubstrates();
        const products = rxn.getProducts();

        var v: f64 = switch (rxn.base.kinetics) {
            .mass_action => |ma| blk: {
                var forward: f64 = ma.k_forward;

                for (substrates) |s| {
                    const c = state.getLiquidConcentration(phase_idx, s.molecule);

                    forward *= std.math.pow(f64, c, s.coefficient);
                }

                var reverse: f64 = 0.0;

                if (rxn.base.reversible and ma.k_reverse > 0) {
                    reverse = ma.k_reverse;

                    for (products) |p| {
                        const c = state.getLiquidConcentration(phase_idx, p.molecule);

                        reverse *= std.math.pow(f64, c, p.coefficient);
                    }
                }

                break :blk forward - reverse;
            },
            .michaelis_menten => |mm| blk: {
                const e = state.getLiquidConcentration(phase_idx, mm.enzyme);
                const s = state.getLiquidConcentration(phase_idx, mm.substrate.molecule);

                if (e < constants.MIN_CONCENTRATION or s < constants.MIN_CONCENTRATION) {
                    break :blk 0.0;
                }

                break :blk mm.kcat * e * s / (mm.km + s);
            },
            .hill => |h| blk: {
                const e = state.getLiquidConcentration(phase_idx, h.enzyme);
                const s = state.getLiquidConcentration(phase_idx, h.substrate.molecule);

                if (e < constants.MIN_CONCENTRATION or s < constants.MIN_CONCENTRATION) {
                    break :blk 0.0;
                }

                const s_n = std.math.pow(f64, s, h.n);
                const k_n = std.math.pow(f64, h.k, h.n);

                break :blk h.kcat * e * s_n / (k_n + s_n);
            },
            .ordered_bi_bi => |ob| blk: {
                // Ordered Bi-Bi: E + A ⇌ EA, EA + B ⇌ EAB → E + P + Q
                // v = kcat * [E] * [A] * [B] / (Km_A * Km_B + Km_B * [A] + Km_A * [B] + [A] * [B])
                // Simplified: assumes rapid equilibrium
                const e = state.getLiquidConcentration(phase_idx, ob.enzyme);
                const a = state.getLiquidConcentration(phase_idx, ob.substrates[0].molecule);
                const b = state.getLiquidConcentration(phase_idx, ob.substrates[1].molecule);

                if (e < constants.MIN_CONCENTRATION) {
                    break :blk 0.0;
                }

                if (a < constants.MIN_CONCENTRATION or b < constants.MIN_CONCENTRATION) {
                    break :blk 0.0;
                }

                const denom = ob.km_a * ob.km_b + ob.km_b * a + ob.km_a * b + a * b;
                break :blk ob.kcat * e * a * b / denom;
            },
            .ping_pong => |pp| blk: {
                // Ping-Pong Bi-Bi: E + A ⇌ F + P, F + B ⇌ E + Q
                // v = kcat * [E] * [A] * [B] / (Km_A * [B] + Km_B * [A] + [A] * [B])
                const e = state.getLiquidConcentration(phase_idx, pp.enzyme);
                const a = state.getLiquidConcentration(phase_idx, pp.substrates[0].molecule);
                const b = state.getLiquidConcentration(phase_idx, pp.substrates[1].molecule);

                if (e < constants.MIN_CONCENTRATION) {
                    break :blk 0.0;
                }

                if (a < constants.MIN_CONCENTRATION or b < constants.MIN_CONCENTRATION) {
                    break :blk 0.0;
                }

                const denom = pp.km_a * b + pp.km_b * a + a * b;
                break :blk pp.kcat * e * a * b / denom;
            },
            .random_bi_bi => |rb| blk: {
                // Random Bi-Bi: substrates can bind in either order
                // v = kcat * [E] * [A] * [B] / (α * Km_A * Km_B + Km_B * [A] + Km_A * [B] + [A] * [B])
                // α < 1: synergistic binding, α > 1: antagonistic binding
                const e = state.getLiquidConcentration(phase_idx, rb.enzyme);
                const a = state.getLiquidConcentration(phase_idx, rb.substrates[0].molecule);
                const b = state.getLiquidConcentration(phase_idx, rb.substrates[1].molecule);

                if (e < constants.MIN_CONCENTRATION) {
                    break :blk 0.0;
                }

                if (a < constants.MIN_CONCENTRATION or b < constants.MIN_CONCENTRATION) {
                    break :blk 0.0;
                }

                const denom = rb.alpha * rb.km_a * rb.km_b + rb.km_b * a + rb.km_a * b + a * b;
                break :blk rb.kcat * e * a * b / denom;
            },
            .first_order => |fo| blk: {
                const s = state.getLiquidConcentration(phase_idx, fo.substrate.molecule);

                break :blk fo.k * s;
            },
            .arrhenius => |arr| blk: {
                const k_forward = arr.A * @exp(-arr.Ea / (constants.R * state.temperature));

                var forward: f64 = k_forward;

                for (substrates) |s| {
                    const c = state.getLiquidConcentration(phase_idx, s.molecule);

                    forward *= std.math.pow(f64, c, s.coefficient);
                }

                var reverse: f64 = 0.0;

                if (arr.A_reverse) |A_rev| {
                    const Ea_rev = arr.Ea_reverse orelse arr.Ea;
                    const k_reverse = A_rev * @exp(-Ea_rev / (constants.R * state.temperature));

                    reverse = k_reverse;

                    for (products) |p| {
                        const c = state.getLiquidConcentration(phase_idx, p.molecule);

                        reverse *= std.math.pow(f64, c, p.coefficient);
                    }
                }

                break :blk forward - reverse;
            },
            .eyring => |eyr| blk: {
                const kB_T_h = constants.k_B * state.temperature / constants.h;
                const k_forward = kB_T_h * eyr.kappa *
                    @exp(eyr.delta_S / constants.R) *
                    @exp(-eyr.delta_H / (constants.R * state.temperature));

                var forward: f64 = k_forward;

                for (substrates) |s| {
                    const c = state.getLiquidConcentration(phase_idx, s.molecule);

                    forward *= std.math.pow(f64, c, s.coefficient);
                }

                var reverse: f64 = 0.0;

                if (eyr.delta_H_reverse) |dH_rev| {
                    const dS_rev = eyr.delta_S_reverse orelse eyr.delta_S;
                    const k_reverse = kB_T_h * eyr.kappa *
                        @exp(dS_rev / constants.R) *
                        @exp(-dH_rev / (constants.R * state.temperature));
                    reverse = k_reverse;

                    for (products) |p| {
                        const c = state.getLiquidConcentration(phase_idx, p.molecule);

                        reverse *= std.math.pow(f64, c, p.coefficient);
                    }
                }

                break :blk forward - reverse;
            },
            .langmuir_hinshelwood => |lh| blk: {
                const c_a = state.getLiquidConcentration(phase_idx, lh.substrates[0].molecule);
                const c_b = state.getLiquidConcentration(phase_idx, lh.substrates[1].molecule);

                const denom = 1.0 + lh.K_ads_a * c_a + lh.K_ads_b * c_b;
                const theta_a = lh.K_ads_a * c_a / denom;
                const theta_b = lh.K_ads_b * c_b / denom;

                var area: f64 = 1.0;

                if (rxn.locus == .solid_liquid_interface) {
                    const idx = @intFromEnum(rxn.locus.solid_liquid_interface.solid);

                    area = state.solid_surface_areas[idx];
                }

                break :blk lh.k_surface * theta_a * theta_b * area;
            },
            .eley_rideal => |er| blk: {
                const c_ads = state.getLiquidConcentration(phase_idx, er.substrates[0].molecule);
                const c_gas = state.getLiquidConcentration(phase_idx, er.substrates[1].molecule);
                const theta = er.K_ads * c_ads / (1.0 + er.K_ads * c_ads);

                break :blk er.k * theta * c_gas;
            },
            .diffusion_limited => |dl| blk: {
                // Smoluchowski with optional viscosity/temperature correction
                var k_diff = 4.0 * std.math.pi * dl.diffusion_coeff * dl.reaction_radius * constants.N_A;

                // Temperature correction: D ∝ T/η
                // η(water) ≈ η(298) * exp(1800 * (1/T - 1/298))
                const eta_ratio = @exp(1800.0 * (1.0 / 298.15 - 1.0 / state.temperature));
                k_diff *= (state.temperature / 298.15) * eta_ratio;

                const c_a = state.getLiquidConcentration(phase_idx, dl.substrates[0].molecule);
                const c_b = state.getLiquidConcentration(phase_idx, dl.substrates[1].molecule);

                // Stirring mildly enhances bulk diffusion-limited reactions (max 5x)
                const stirring_factor = std.math.pow(f64, 5.0, state.stirring);

                break :blk k_diff * c_a * c_b * stirring_factor;
            },
            .autocatalytic => |ac| blk: {
                const c_sub = state.getLiquidConcentration(phase_idx, ac.substrate.molecule);
                const c_cat = state.getLiquidConcentration(phase_idx, ac.catalyst.molecule);

                var rate = ac.k * c_sub * c_cat;

                if (ac.capacity) |cap| {
                    rate *= @max(0.0, 1.0 - c_cat / cap);
                }

                break :blk rate;
            },
            .nth_order => |no| blk: {
                const c = state.getLiquidConcentration(phase_idx, no.substrate.molecule);

                break :blk no.k * std.math.pow(f64, c, no.n);
            },
            .equilibrium => |eq| blk: {
                const delta_H = getReactionEnthalpy(rxn);
                // Van't Hoff temperature correction (reference = 298.15 K)
                const Keq_T = eq.Keq * @exp(-delta_H / constants.R *
                    (1.0 / state.temperature - 1.0 / 298.15));
                const k_reverse = eq.k_forward / Keq_T;

                var forward: f64 = eq.k_forward;

                for (substrates) |s| {
                    const c = state.getLiquidConcentration(phase_idx, s.molecule);
                    forward *= std.math.pow(f64, c, s.coefficient);
                }

                var reverse: f64 = k_reverse;

                for (products) |p| {
                    const c = state.getLiquidConcentration(phase_idx, p.molecule);
                    reverse *= std.math.pow(f64, c, p.coefficient);
                }

                break :blk forward - reverse;
            },
            .dissociation => 0.0,
            .strong_electrolyte => 0.0,
            .weak_electrolyte => 0.0,
            .dissolution => |dsl| blk: {
                // AB(s) ⇌ A⁺(aq) + B⁻(aq)
                // Forward rate depends on solid surface area
                // Reverse rate depends on supersaturation

                // Get solid phase info
                const solid_idx = @intFromEnum(dsl.solid.molecule);
                const surface_area = state.solid_surface_areas[solid_idx];

                if (surface_area < constants.MIN_AREA) {
                    break :blk this.computePrecipitationRate(rxn, dsl, state, phase_idx);
                }

                const delta_H = getReactionEnthalpy(rxn);
                // Temperature-corrected Ksp via van't Hoff
                const Ksp_T = dsl.Ksp * @exp(-delta_H / constants.R *
                    (1.0 / state.temperature - 1.0 / dsl.T_ref));

                // Current ion activity product Q
                const c_cation = state.getLiquidConcentration(phase_idx, dsl.ions[0].molecule);
                const c_anion = state.getLiquidConcentration(phase_idx, dsl.ions[1].molecule);

                const nu_cat = dsl.ions[0].coefficient;
                const nu_an = dsl.ions[1].coefficient;

                var Q = std.math.pow(f64, c_cation, nu_cat) * std.math.pow(f64, c_anion, nu_an);

                // Activity correction
                const ionic_strength = this.calculateIonicStrength(state, phase_idx);
                const gamma_cat = this.calculateSingleIonActivity(ionic_strength, dsl.ions[0].molecule.getFullyProtonatedCharge());
                const gamma_an = this.calculateSingleIonActivity(ionic_strength, dsl.ions[1].molecule.getFullyProtonatedCharge());
                const gamma_product = std.math.pow(f64, gamma_cat, nu_cat) *
                    std.math.pow(f64, gamma_an, nu_an);

                Q *= gamma_product;

                // Saturation ratio: S = Q / Ksp
                // S < 1: undersaturated (dissolution favored)
                // S > 1: supersaturated (precipitation favored)
                // S = 1: equilibrium
                const saturation = Q / Ksp_T;

                // Dissolution rate (positive when dissolving)
                // v_dissolve = k_d * A * (1 - S)  for S < 1
                // v_precipitate = k_p * A * (S - 1)  for S > 1

                var rate: f64 = 0.0;

                if (saturation < 1.0) {
                    // Undersaturated: dissolution
                    // Rate law: v = k * A * (1 - S)^n, n typically 1-2
                    rate = dsl.k_dissolve * surface_area * (1.0 - saturation);
                } else {
                    // Supersaturated: precipitation onto existing solid
                    // Rate law: v = -k * A * (S - 1)^n
                    const k_precipitate = dsl.k_dissolve * 0.1;
                    const supersaturation = saturation - 1.0;
                    rate = -k_precipitate * surface_area * supersaturation;
                }

                break :blk rate;
            },
        };

        // Apply pH effect (bell curve around optimum)
        if (rxn.base.ph_optimum) |opt| {
            const ph = state.getPH(phase_idx);
            const ph_diff = ph - opt;

            v *= @exp(-(ph_diff * ph_diff) / (2.0 * rxn.base.ph_width * rxn.base.ph_width));
        }

        // Q10 temperature correction for kinetics without built-in T-dependence
        switch (rxn.base.kinetics) {
            .mass_action,
            .first_order,
            .michaelis_menten,
            .hill,
            .ordered_bi_bi,
            .ping_pong,
            .random_bi_bi,
            => {
                const temp_factor = std.math.pow(f64, rxn.base.q10, (state.temperature - 310.15) / 10.0);

                v *= temp_factor;
            },
            else => {}, // T-dependence built into kinetics
        }

        // Apply modulators with proper mechanism handling
        v = this.applyModulatorsForPhase(v, &rxn.base, state, phase_idx);

        if (!rxn.base.reversible) {
            v = @max(0.0, v);
        }

        v = std.math.clamp(v, -constants.MAX_RATE, constants.MAX_RATE);

        return v;
    }

    /// Compute dissolution rate for strong electrolyte (mol/s)
    /// Rate = k * A * (1 - C/C_sat) for undersaturated
    /// Rate = -k_cryst * A * (C/C_sat - 1) for supersaturated
    inline fn computeStrongElectrolyteDissolution(
        this: *const MultiPhaseODESolver,
        se: KineticsType.StrongElectrolyte,
        state: *const PhaseState,
        liquid_idx: usize,
        surface_area: f64,
    ) f64 {
        _ = this;

        // Current total salt concentration from ions
        // For 1:1 salt: C_salt = [cation] = [anion]
        // For 1:2 salt: C_salt = [cation] = [anion]/2
        const c_cation = state.getLiquidConcentration(liquid_idx, se.ions[0].molecule);
        const c_anion = state.getLiquidConcentration(liquid_idx, se.ions[1].molecule);

        // Use limiting ion concentration
        const nu_cat = se.ions[0].coefficient;
        const nu_an = se.ions[1].coefficient;

        const c_salt_from_cation = c_cation / nu_cat;
        const c_salt_from_anion = c_anion / nu_an;
        const c_salt = @min(c_salt_from_cation, c_salt_from_anion);

        // Saturation ratio
        const saturation = c_salt / se.saturation;

        // Temperature effect on dissolution kinetics (Arrhenius-like)
        // Dissolution is typically diffusion-controlled, D ∝ T/η
        // η(water) decreases ~2-3% per degree C
        // Simplified: rate doubles roughly every 20-25°C
        const temp_factor = std.math.pow(f64, 2.0, (state.temperature - 298.15) / 25.0);

        if (saturation < 1.0) {
            // Undersaturated: dissolution
            // Rate law: v = k * A * (1 - S)
            const driving_force = 1.0 - saturation;
            // Stirring strongly enhances dissolution by reducing boundary layer (max 10x)
            const stirring_factor = std.math.pow(f64, 10.0, state.stirring);

            return se.k_dissolve * surface_area * driving_force * temp_factor * stirring_factor;
        } else {
            // Supersaturated: crystallization (reverse)
            const k_cryst = se.k_dissolve * 0.1;
            const driving_force = @min(saturation - 1.0, 10.0);
            // Crystallization is less enhanced by stirring than dissolution (max 4x)
            const stirring_factor = std.math.pow(f64, 4.0, state.stirring);

            return -k_cryst * surface_area * driving_force * temp_factor * stirring_factor;
        }
    }

    /// Compute dissolution rate for Ksp-controlled dissolution (mol/s)
    inline fn computeDissolutionRate(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        dsl: KineticsType.Dissolution,
        state: *const PhaseState,
        liquid_idx: usize,
        solid_idx: usize,
        surface_area: f64,
    ) f64 {
        _ = solid_idx;

        const delta_H = getReactionEnthalpy(rxn);

        // Clamp temperature to prevent division by zero and Ksp_T underflow
        const safe_temp = @max(1.0, state.temperature);
        const Ksp_T = @max(1e-300, dsl.Ksp * @exp(-delta_H / constants.R *
            (1.0 / safe_temp - 1.0 / dsl.T_ref)));

        // Current ion activity product Q
        const c_cation = state.getLiquidConcentration(liquid_idx, dsl.ions[0].molecule);
        const c_anion = state.getLiquidConcentration(liquid_idx, dsl.ions[1].molecule);

        const nu_cat = dsl.ions[0].coefficient;
        const nu_an = dsl.ions[1].coefficient;

        var Q = std.math.pow(f64, c_cation, nu_cat) * std.math.pow(f64, c_anion, nu_an);

        // Activity correction
        const ionic_strength = this.calculateIonicStrength(state, liquid_idx);
        const gamma_cat = this.calculateSingleIonActivity(ionic_strength, dsl.ions[0].molecule.getFullyProtonatedCharge());
        const gamma_an = this.calculateSingleIonActivity(ionic_strength, dsl.ions[1].molecule.getFullyProtonatedCharge());
        const gamma_product = std.math.pow(f64, gamma_cat, nu_cat) *
            std.math.pow(f64, gamma_an, nu_an);

        Q *= gamma_product;

        // Saturation ratio
        const saturation = Q / Ksp_T;

        if (saturation < 1.0) {
            // Dissolution - stirring enhances mass transfer (max 10x)
            const stirring_factor = std.math.pow(f64, 10.0, state.stirring);

            return dsl.k_dissolve * surface_area * (1.0 - saturation) * stirring_factor;
        } else {
            // Precipitation - less affected by stirring (max 4x)
            const k_precipitate = dsl.k_dissolve * 0.1;
            const driving_force = @min(saturation - 1.0, 10.0);
            const stirring_factor = std.math.pow(f64, 4.0, state.stirring);

            return -k_precipitate * surface_area * driving_force * stirring_factor;
        }
    }

    /// Calculate ionic strength of solution: I = 0.5 * Σ(c_i * z_i²)
    inline fn calculateIonicStrength(
        this: *const MultiPhaseODESolver,
        state: *const PhaseState,
        phase_idx: usize,
    ) f64 {
        _ = this;

        var ionic_strength: f64 = 0.0;

        for (0..MoleculeId.count()) |i| {
            const mol: MoleculeId = @enumFromInt(i);
            const c = state.getLiquidConcentration(phase_idx, mol);

            if (c < constants.MIN_CONCENTRATION) {
                continue;
            }

            // Get charge at current pH
            const ph = state.getPH(phase_idx);
            const charge = mol.getCharge(@floatCast(ph));

            ionic_strength += c * charge * charge;
        }

        return 0.5 * ionic_strength;
    }

    /// Calculate single-ion activity coefficient using Davies equation
    inline fn calculateSingleIonActivity(
        this: *const MultiPhaseODESolver,
        ionic_strength: f64,
        charge: i8,
    ) f64 {
        _ = this;

        if (ionic_strength < 1e-10) {
            return 1.0;
        }

        // The Davies equation is empirically valid only up to I ≈ 0.5 M.
        // For highly concentrated solutions or nearly dry phases, the linear term (-0.3 * I)
        // causes log_gamma to grow unboundedly positive, resulting in infinite activity.
        // Clamp the effective ionic strength to prevent mathematical explosion.
        const effective_I = @min(ionic_strength, 0.5);

        const A: f64 = 0.509;
        const z: f64 = @floatFromInt(charge);
        const z_sq = z * z;

        const sqrt_I = @sqrt(effective_I);
        const log_gamma = -A * z_sq * (sqrt_I / (1.0 + sqrt_I) - 0.3 * effective_I);

        return std.math.pow(f64, 10.0, log_gamma);
    }

    /// Compute precipitation rate when no solid phase exists yet (homogeneous nucleation)
    inline fn computePrecipitationRate(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        dsl: KineticsType.Dissolution,
        state: *const PhaseState,
        phase_idx: usize,
    ) f64 {
        const c_cation = state.getLiquidConcentration(phase_idx, dsl.ions[0].molecule);
        const c_anion = state.getLiquidConcentration(phase_idx, dsl.ions[1].molecule);

        const nu_cat = dsl.ions[0].coefficient;
        const nu_an = dsl.ions[1].coefficient;

        const Q = std.math.pow(f64, c_cation, nu_cat) * std.math.pow(f64, c_anion, nu_an);

        // Activity correction
        const ionic_strength = this.calculateIonicStrength(state, phase_idx);
        const gamma_cat = this.calculateSingleIonActivity(ionic_strength, dsl.ions[0].molecule.getFullyProtonatedCharge());
        const gamma_an = this.calculateSingleIonActivity(ionic_strength, dsl.ions[1].molecule.getFullyProtonatedCharge());
        const gamma_product = std.math.pow(f64, gamma_cat, nu_cat) *
            std.math.pow(f64, gamma_an, nu_an);

        const Q_activity = Q * gamma_product;

        const delta_H = getReactionEnthalpy(rxn);

        // Clamp temperature to prevent division by zero and Ksp_T underflow
        const safe_temp = @max(1.0, state.temperature);
        const Ksp_T = @max(1e-300, dsl.Ksp * @exp(-delta_H / constants.R *
            (1.0 / safe_temp - 1.0 / dsl.T_ref)));

        const saturation = Q_activity / Ksp_T;

        // Classical nucleation theory: nucleation rate ~ exp(-ΔG*/kT)
        // Simplified: significant nucleation above critical supersaturation
        const S_crit: f64 = 5.0; // Critical supersaturation ratio

        if (saturation < S_crit) {
            return 0.0; // Below nucleation threshold
        }

        // Nucleation rate increases rapidly above S_crit
        // Use empirical power law: rate ~ k * (S - S_crit)^n
        const excess = saturation - S_crit;

        // Homogeneous nucleation much slower
        const k_nucleation = dsl.k_dissolve * 1e-4;

        // Stirring disrupts critical nuclei formation.
        // Higher stirring = slower nucleation. Minimum multiplier is 0.3x (70% reduction).
        const stirring_factor = std.math.pow(f64, 0.3, state.stirring);

        return -k_nucleation * excess * excess * stirring_factor;
    }

    /// Compute interfacial reaction rate (mol/s) for gas-liquid interface
    inline fn computeInterfacialRate(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        state: *const PhaseState,
        liquid_idx: usize,
        area: f64,
        substrates: []const Stoich,
    ) f64 {
        _ = this;

        // For gas-liquid reactions, rate often limited by gas dissolution
        if (rxn.area_rate_constant) |k_area| {
            var driving_force: f64 = 1.0;

            for (substrates) |s| {
                if (s.molecule.isGas()) {
                    const p = state.getPartialPressure(s.molecule);
                    driving_force *= p / constants.P_std;
                } else {
                    const c = state.getLiquidConcentration(liquid_idx, s.molecule);
                    driving_force *= c;
                }
            }

            // Stirring massively enhances gas-liquid mass transfer via bubbles/vortices (max 30x)
            const stirring_factor = std.math.pow(f64, 30.0, state.stirring);

            return k_area * area * driving_force * stirring_factor;
        }

        return 0.0;
    }

    /// Compute solid-liquid interfacial rate (mol/s)
    inline fn computeSolidLiquidRate(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        state: *const PhaseState,
        liquid_idx: usize,
        solid_idx: usize,
        area: f64,
        substrates: []const Stoich,
    ) f64 {
        _ = this;
        _ = solid_idx;

        // Dissolution/reaction rate: v = k * A * (C_sat - C)
        // Or for reaction: v = k * A * [reactant]
        if (rxn.area_rate_constant) |k_area| {
            var driving_force: f64 = 1.0;

            for (substrates) |s| {
                if (!s.molecule.isSolid()) {
                    const c = state.getLiquidConcentration(liquid_idx, s.molecule);
                    driving_force *= std.math.pow(f64, c, s.coefficient);
                }
            }

            // Stirring enhances solid-liquid reactions by reducing the Nernst diffusion layer (max 10x)
            const stirring_factor = std.math.pow(f64, 10.0, state.stirring);

            return k_area * area * driving_force * stirring_factor;
        }

        return 0.0;
    }

    /// Compute gas-phase reaction rate (mol/s)
    inline fn computeGasPhaseRate(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        state: *const PhaseState,
        substrates: []const Stoich,
    ) f64 {
        _ = this;

        // Gas phase kinetics typically use partial pressures
        switch (rxn.base.kinetics) {
            .mass_action => |ma| {
                var forward: f64 = ma.k_forward;

                for (substrates) |s| {
                    const p = state.getPartialPressure(s.molecule);

                    forward *= std.math.pow(f64, p / constants.P_std, s.coefficient);
                }

                return forward * state.gas_volume;
            },
            else => return 0.0,
        }
    }

    /// Compute liquid-liquid interface rate (mol/s)
    inline fn computeLiquidInterfaceRate(
        this: *const MultiPhaseODESolver,
        rxn: *const MultiPhaseReaction,
        state: *const PhaseState,
        phase_a: usize,
        phase_b: usize,
        area: f64,
        substrates: []const Stoich,
    ) f64 {
        _ = this;

        if (rxn.area_rate_constant) |k_area| {
            var driving_force: f64 = 1.0;

            // Use geometric mean of concentrations from both phases
            for (substrates) |s| {
                const c_a = state.getLiquidConcentration(phase_a, s.molecule);
                const c_b = state.getLiquidConcentration(phase_b, s.molecule);
                const c_eff = @sqrt(c_a * c_b + 1e-30);

                driving_force *= std.math.pow(f64, c_eff, s.coefficient);
            }

            // Stirring dramatically increases liquid-liquid interfacial area via emulsification (max 100x)
            const stirring_factor = std.math.pow(f64, 100.0, state.stirring);

            return k_area * area * driving_force * stirring_factor;
        }

        return 0.0;
    }

    /// Apply modulator effects to reaction rate with proper mechanism handling.
    ///
    /// Competitive inhibition: affects Km (approximated via rate reduction)
    /// Uncompetitive inhibition: affects both Km and Vmax equally
    /// Non-competitive inhibition: affects Vmax only (exact)
    /// Mixed inhibition: weighted combination
    /// Allosteric: sigmoidal dose-response
    /// Product inhibition: treated as non-competitive
    /// Irreversible: time-dependent (requires state tracking - approximated here)
    inline fn applyModulatorsForPhase(
        this: *const MultiPhaseODESolver,
        base_rate: f64,
        rxn: *const Reaction,
        state: *const PhaseState,
        phase_idx: usize,
    ) f64 {
        _ = this;

        var v = base_rate;

        for (rxn.modulators) |mod| {
            const i_conc = state.getLiquidConcentration(phase_idx, mod.molecule);

            if (i_conc < constants.MIN_CONCENTRATION) {
                continue;
            }

            const i_ratio = i_conc / mod.ki; // [I] / Ki

            switch (mod.effect) {
                .inhibitor => {
                    const factor: f64 = switch (mod.mechanism) {
                        .competitive => blk: {
                            // True competitive: Km_app = Km * (1 + [I]/Ki)
                            // This changes the shape, not just Vmax.
                            // Approximation: at [S] = Km, rate drops by factor 1/(1 + [I]/Ki)
                            // This underestimates inhibition at low [S], overestimates at high [S]
                            break :blk 1.0 / (1.0 + i_ratio);
                        },
                        .uncompetitive => blk: {
                            // Uncompetitive: binds only to ES complex
                            // Vmax_app = Vmax / (1 + [I]/Ki)
                            // Km_app = Km / (1 + [I]/Ki)
                            // Net effect on v: 1 / (1 + [I]/Ki)
                            break :blk 1.0 / (1.0 + i_ratio);
                        },
                        .non_competitive => blk: {
                            // Pure non-competitive: Vmax_app = Vmax / (1 + [I]/Ki), Km unchanged
                            // Rate reduction is exact
                            break :blk 1.0 / (1.0 + i_ratio);
                        },
                        .mixed => blk: {
                            // Mixed: different affinities for E vs ES
                            // alpha = Kii / Ki (ratio of inhibition constants)
                            // Approximation: use max_effect as effective alpha
                            const alpha = mod.max_effect;
                            // Vmax_app = Vmax / (1 + [I]/Kii) = Vmax / (1 + [I]/(alpha*Ki))
                            // Km_app = Km * (1 + [I]/Ki) / (1 + [I]/Kii)
                            // Simplified: geometric mean of competitive and non-competitive
                            const comp_factor = 1.0 / (1.0 + i_ratio);
                            const noncomp_factor = 1.0 / (1.0 + i_ratio / alpha);

                            break :blk @sqrt(comp_factor * noncomp_factor);
                        },
                        .allosteric => blk: {
                            // Sigmoidal inhibition: Hill-like dose-response
                            // Assumes Hill coefficient of 2 for cooperative binding
                            const hill_n: f64 = 2.0;
                            const i_n = std.math.pow(f64, i_ratio, hill_n);
                            const inhibition = i_n / (1.0 + i_n); // 0 to 1

                            // max_effect is the maximum fractional inhibition (0 to 1)
                            break :blk 1.0 - inhibition * mod.max_effect;
                        },
                        .product => blk: {
                            // Product inhibition: typically competitive or mixed
                            // Treated as non-competitive for simplicity
                            break :blk 1.0 / (1.0 + i_ratio);
                        },
                        .irreversible => blk: {
                            // Irreversible (mechanism-based) inhibition
                            // True modeling requires tracking [E_active] over time
                            // Approximation: rapid inactivation kinetics
                            // k_inact typically 0.01-1 s^-1
                            // Effective [E] = [E]_0 * exp(-k_inact * [I] * t / Ki)
                            // Without time tracking, use steady-state approximation:
                            // At high [I]/Ki, enzyme is mostly inactivated
                            const saturation = i_ratio / (1.0 + i_ratio);

                            // Assume 90% inactivation at saturation * max_effect
                            break :blk 1.0 - 0.9 * saturation * @min(1.0, mod.max_effect);
                        },
                    };
                    v *= @max(0.0, factor);
                },

                .activator => {
                    const factor: f64 = switch (mod.mechanism) {
                        .allosteric => blk: {
                            // Sigmoidal activation with cooperativity
                            const hill_n: f64 = 2.0;
                            const a_n = std.math.pow(f64, i_ratio, hill_n);
                            const activation = a_n / (1.0 + a_n);

                            // max_effect is fold-increase at saturation
                            break :blk 1.0 + activation * mod.max_effect;
                        },
                        else => blk: {
                            // Simple hyperbolic activation
                            const saturation = i_ratio / (1.0 + i_ratio);

                            break :blk 1.0 + saturation * mod.max_effect;
                        },
                    };
                    v *= factor;
                },
            }
        }

        return v;
    }

    /// Helper to advance temp_state for RK stages
    inline fn advanceTemp(
        this: *MultiPhaseODESolver,
        base: *const PhaseState,
        stages: []const usize,
        coeffs: []const f64,
        dt: f64,
    ) void {
        const n_liquids = base.liquid_moles.len;
        const n_solids = base.solid_moles.len;
        const n_mol = MoleculeId.count();

        for (0..n_liquids) |p| {
            for (0..n_mol) |m| {
                var sum: f64 = 0;

                for (stages, coeffs) |s, c| {
                    sum += c * this.k_liquid[s][p][m];
                }

                this.temp_state.liquid_moles[p][m] = @max(0.0, base.liquid_moles[p][m] + dt * sum);
            }
        }

        for (0..n_solids) |s_idx| {
            var sum: f64 = 0;

            for (stages, coeffs) |s, c| {
                sum += c * this.k_solid[s][s_idx];
            }

            this.temp_state.solid_moles[s_idx] = @max(0.0, base.solid_moles[s_idx] + dt * sum);
        }

        for (0..n_mol) |m| {
            var sum: f64 = 0;

            for (stages, coeffs) |s, c| {
                sum += c * this.k_gas[s][m];
            }

            this.temp_state.gas_moles[m] = @max(0.0, base.gas_moles[m] + dt * sum);
        }

        var heat_sum: f64 = 0.0;

        for (stages, coeffs) |s, c| {
            heat_sum += c * this.k_heat[s];
        }

        this.temp_state.heat_energy = base.heat_energy + dt * heat_sum;

        this.updateCachedQuantities(&this.temp_state);
        this.temp_state.updatePH();
    }

    /// Update cached volumes, surface areas, heat capacity, temperature,
    /// and gas volume from current moles. Called at every RK stage.
    inline fn updateCachedQuantities(this: *MultiPhaseODESolver, state: *PhaseState) void {
        const n_mol = MoleculeId.count();

        // 1. Update liquid volumes
        for (state.liquid_volumes, 0..) |*vol, p| {
            var v: f64 = 0;

            for (0..n_mol) |m| {
                const mol: MoleculeId = @enumFromInt(m);
                const moles = state.liquid_moles[p][m];

                if (moles < constants.MIN_MOLES) {
                    continue;
                }

                // Ions don't contribute significantly to volume
                if (mol.isIon()) {
                    continue;
                }

                // V = n * MW / density / 1000 (liters)
                v += moles * mol.getWeight() / mol.getDensity() / 1000.0;
            }

            // Ensure minimum volume if there are any moles at all
            vol.* = @max(v, constants.MIN_VOLUME);
        }

        // 2. Update solid surface areas
        for (0..n_mol) |i| {
            const mol: MoleculeId = @enumFromInt(i);
            const new_moles = state.solid_moles[i];

            var original_moles: f64 = 0;
            var original_area: f64 = 0;

            for (this.contents.solids.items) |solid| {
                if (solid.molecule == mol) {
                    original_moles = solid.moles;
                    original_area = solid.getSurfaceArea();

                    break;
                }
            }

            if (original_moles > 1e-15) {
                const ratio = new_moles / original_moles;

                state.solid_surface_areas[i] = original_area * std.math.pow(f64, ratio, 2.0 / 3.0);
            } else if (new_moles > 1e-15) {
                // V = n * MW / density (in liters -> convert to m^3)
                const volume = new_moles * mol.getWeight() / mol.getDensity() / 1000.0;
                const default_diameter = 1e-6;

                state.solid_surface_areas[i] = 6.0 * volume / default_diameter;
            } else {
                state.solid_surface_areas[i] = 0;
            }
        }

        // 3. Update Thermodynamic State (Cp and T)
        var cp_total = state.container_heat_capacity;

        // Gas heat capacity
        for (0..n_mol) |m| {
            const moles = state.gas_moles[m];

            if (moles > 1e-15) {
                const mol: MoleculeId = @enumFromInt(m);

                cp_total += moles * mol.getHeatCapacity();
            }
        }

        // Liquid heat capacity
        for (state.liquid_moles) |phase| {
            for (0..n_mol) |m| {
                const moles = phase[m];

                if (moles > 1e-15) {
                    const mol: MoleculeId = @enumFromInt(m);

                    cp_total += moles * mol.getHeatCapacity();
                }
            }
        }

        // Solid heat capacity
        for (0..n_mol) |s_idx| {
            const moles = state.solid_moles[s_idx];

            if (moles > 1e-15) {
                const solid_mol: MoleculeId = @enumFromInt(s_idx);
                cp_total += moles * solid_mol.getHeatCapacity();
            }
        }

        for (this.contents.solids.items) |*solid| {
            if (solid.occluded_solution) |*occ| {
                for (0..n_mol) |m| {
                    if (occ.solutes[m] > 1e-15) {
                        const mol: MoleculeId = @enumFromInt(m);

                        cp_total += occ.solutes[m] * mol.getHeatCapacity();
                    }
                }
            }
        }

        // Update cached T and Cp
        state.heat_capacity = cp_total;
        state.temperature = if (cp_total > 0.0)
            state.heat_energy / cp_total
        else
            0.0;

        // 4. Update gas volume via ideal gas law
        var total_gas_moles: f64 = 0;

        for (state.gas_moles) |n| {
            total_gas_moles += n;
        }

        state.gas_volume = total_gas_moles * constants.R * state.temperature / state.pressure * 1000.0;
    }

    inline fn adaptStep(this: *MultiPhaseODESolver, current_dt: f64, err: f64) f64 {
        const safety = 0.9;

        if (err < 1e-15) {
            return @min(current_dt * 2.0, this.max_dt);
        }

        const factor = safety * std.math.pow(f64, this.tolerance / err, 0.2);

        if (factor >= 2.0) {
            return @min(current_dt * 2.0, this.max_dt);
        } else if (factor <= 0.5) {
            return @max(current_dt * 0.5, this.min_dt);
        } else {
            return std.math.clamp(current_dt * factor, this.min_dt, this.max_dt);
        }
    }
};

pub const CACO_WATER_DISSOLUTION = MultiPhaseReaction{
    .base = .{
        .name = "CaCO3 dissolution in water",
        .kinetics = .{
            .dissolution = .{
                .Ksp = 3.3e-9,
                .k_dissolve = 1e-5,
                .solid = Stoich{
                    .molecule = .calcium_carbonate,
                },
                .ions = .{
                    Stoich{
                        .molecule = .calcium_ion,
                    },
                    Stoich{
                        .molecule = .carbonate_ion,
                    },
                },
                .T_ref = 298.15,
            },
        },
        .reversible = true,
        .ph_optimum = 4.0,
        .ph_width = 3.0,
    },
    .locus = .{
        .solid_liquid_interface = .{
            .solid = .calcium_carbonate,
            .liquid_solvent = .oxidane,
        },
    },
    .area_rate_constant = 1e-5,
};

pub const NACL_WATER_DISSOCIATION = MultiPhaseReaction{
    .base = .{
        .name = "NaCl dissociation in water",
        .kinetics = .{
            .strong_electrolyte = .{
                .k_dissolve = 1e-2,
                .saturation = 6.14,
                .solid = Stoich{
                    .molecule = .sodium_chloride,
                },
                .ions = .{
                    Stoich{
                        .molecule = .sodium_ion,
                    },
                    Stoich{
                        .molecule = .chloride_ion,
                    },
                },
            },
        },
    },
    .locus = .{
        .solid_liquid_interface = .{
            .solid = .sodium_chloride,
            .liquid_solvent = .oxidane,
        },
    },
};

pub const ACETIC_ACID_DISSOCIATION = MultiPhaseReaction{
    .base = .{
        .name = "Acetic acid dissocication in water",
        .kinetics = .{
            .weak_electrolyte = .{
                .pKa = 4.76,
                .k_forward = 1e3,
                .parent = Stoich{
                    .molecule = .acetic_acid,
                },
                .products = .{
                    Stoich{
                        .molecule = .acetate_ion,
                    },
                    Stoich{
                        .molecule = .hydrogen_ion,
                    },
                },
            },
        },
    },
    .locus = .{
        .liquid_bulk = .{
            .solvent = .oxidane,
        },
    },
};

const MAX_VOLUME: f64 = 4.0;

test "CaCO3 dissolution" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.addSolid(std.testing.allocator, .calcium_carbonate, MoleculeId.calcium_carbonate.molesPerLiter() * 0.01, 2.5e-5);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(1, contents.solids.items.len);

    contents.stirring = 0.5;

    var solver = MultiPhaseODESolver.init(std.testing.allocator);
    defer solver.deinit();

    try solver.integrate(std.testing.allocator, &.{CACO_WATER_DISSOLUTION}, &contents, 10.0, MAX_VOLUME);
    try contents.settle(std.testing.allocator);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(1, contents.solids.items.len);

    try std.testing.expect(!constants.isNegligible(contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.calcium_ion)]));
    try std.testing.expect(!constants.isNegligible(contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.carbonate_ion)]));
}

test "NaCl dissociation" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.addSolid(std.testing.allocator, .sodium_chloride, MoleculeId.sodium_chloride.molesPerLiter() * 0.01, 2.5e-5);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(1, contents.solids.items.len);

    contents.stirring = 0.5;
    try contents.setTemperature(std.testing.allocator, constants.cToK(60), MAX_VOLUME);

    var solver = MultiPhaseODESolver.init(std.testing.allocator);
    defer solver.deinit();

    try solver.integrate(std.testing.allocator, &.{NACL_WATER_DISSOCIATION}, &contents, 10.0, MAX_VOLUME);
    try contents.settle(std.testing.allocator);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(0, contents.solids.items.len);
}

test "Ions in dry contents" {
    var src = try Contents.init(std.testing.allocator);
    defer src.deinit(std.testing.allocator);

    try src.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter() * 0.025);
    try src.addSolid(std.testing.allocator, .sodium_chloride, MoleculeId.sodium_chloride.molesPerLiter() * 0.025, 2.5e-5);

    helpers.resetStdAir(&src, MAX_VOLUME);
    try src.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    var solver = MultiPhaseODESolver.init(std.testing.allocator);
    defer solver.deinit();

    const REACTIONS = &.{ NACL_WATER_DISSOCIATION, CACO_WATER_DISSOLUTION };

    for (0..3) |_| {
        try src.updatePhaseTransitions(std.testing.allocator, 0.1, MAX_VOLUME);
        try solver.integrate(std.testing.allocator, REACTIONS, &src, 0.1, MAX_VOLUME);
        try src.settle(std.testing.allocator);
    }

    try std.testing.expectEqual(1, src.liquids.items.len);
    try std.testing.expectEqual(1, src.solids.items.len);

    var dst = try Contents.init(std.testing.allocator);
    defer dst.deinit(std.testing.allocator);
    helpers.resetStdAir(&dst, MAX_VOLUME);
    try dst.setTemperature(std.testing.allocator, constants.cToK(25), MAX_VOLUME);

    _ = try src.transferLiquidVolume(std.testing.allocator, &dst, std.math.inf(f64));

    for (0..1) |_| {
        helpers.resetStdAir(&src, MAX_VOLUME);

        try src.updatePhaseTransitions(std.testing.allocator, 0.1, MAX_VOLUME);
        try solver.integrate(std.testing.allocator, REACTIONS, &src, 0.1, MAX_VOLUME);
        try src.settle(std.testing.allocator);

        helpers.resetStdAir(&dst, MAX_VOLUME);

        try dst.updatePhaseTransitions(std.testing.allocator, 0.1, MAX_VOLUME);
        try solver.integrate(std.testing.allocator, REACTIONS, &dst, 0.1, MAX_VOLUME);
        try dst.settle(std.testing.allocator);
    }

    try std.testing.expect(std.math.isFinite(dst.getTemperature()));
}

test "NaCl evaporation" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    const NACL_MOLES = MoleculeId.sodium_chloride.molesPerLiter() * 0.01;

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter() * 0.25);
    try contents.addSolid(std.testing.allocator, .sodium_chloride, NACL_MOLES, 2.5e-5);
    try contents.setTemperature(std.testing.allocator, constants.cToK(60), MAX_VOLUME);

    contents.stirring = 0.5;

    // First dissoaciate
    var solver = MultiPhaseODESolver.init(std.testing.allocator);
    defer solver.deinit();

    try solver.integrate(std.testing.allocator, &.{NACL_WATER_DISSOCIATION}, &contents, 10.0, MAX_VOLUME);
    try contents.settle(std.testing.allocator);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(0, contents.solids.items.len);

    // Then evaporate water

    var secs: f64 = 0.0;

    for (0..90) |_| {
        const dt: f64 = 1.0;

        _ = contents.exchangeHeat(constants.cToK(500), 200, dt);
        try contents.updatePhaseTransitions(std.testing.allocator, dt, MAX_VOLUME);
        try solver.integrate(std.testing.allocator, &.{NACL_WATER_DISSOCIATION}, &contents, dt, MAX_VOLUME);
        try contents.settle(std.testing.allocator);

        secs += dt;
    }

    try std.testing.expectEqual(0, contents.liquids.items.len);
    try std.testing.expectEqual(1, contents.solids.items.len);
    try std.testing.expectApproxEqAbs(NACL_MOLES, contents.solids.items[0].moles, 1e-10);
}

test "Acetic acid dissociation" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerLiter());
    try contents.addLiquid(std.testing.allocator, .acetic_acid, MoleculeId.acetic_acid.molesPerLiter() * 0.1);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(0, contents.solids.items.len);

    var solver = MultiPhaseODESolver.init(std.testing.allocator);
    defer solver.deinit();

    try solver.integrate(std.testing.allocator, &.{ACETIC_ACID_DISSOCIATION}, &contents, 10.0, MAX_VOLUME);
    try contents.settle(std.testing.allocator);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(0, contents.solids.items.len);

    try std.testing.expect(!constants.isNegligible(contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.acetate_ion)]));
    try std.testing.expect(!constants.isNegligible(contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.hydrogen_ion)]));
    // Weak acid, very small dissociation
    try std.testing.expect(!constants.isNegligible(contents.liquids.items[0].solutes[@intFromEnum(MoleculeId.acetic_acid)]));
}

test "Acetic acid freezing" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .acetic_acid, MoleculeId.acetic_acid.molesPerLiter());
    try contents.setTemperature(std.testing.allocator, constants.cToK(15), MAX_VOLUME);

    try std.testing.expectEqual(0, contents.liquids.items.len);
    try std.testing.expectEqual(1, contents.solids.items.len);
}

test "Acetic acid freezing prevention" {
    var contents = try Contents.init(std.testing.allocator);
    defer contents.deinit(std.testing.allocator);

    try contents.addLiquid(std.testing.allocator, .oxidane, MoleculeId.oxidane.molesPerG() * 30);
    try contents.addLiquid(std.testing.allocator, .acetic_acid, MoleculeId.acetic_acid.molesPerG() * 70);

    try contents.setTemperature(std.testing.allocator, constants.cToK(15), MAX_VOLUME);

    try std.testing.expectEqual(1, contents.liquids.items.len);
    try std.testing.expectEqual(0, contents.solids.items.len);
}
