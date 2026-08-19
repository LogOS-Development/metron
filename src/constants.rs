//! Physical constants in typed SI quantities.
//!
//! All constants are expressed as `Quantity<f64, _>` so they carry
//! compile-time unit information.

use crate::aliases::*;
use crate::quantity::Quantity;

/// Speed of light in vacuum (m/s).
pub const C: Velocity = Velocity::new(299_792_458.0);

/// Reduced Planck constant ℏ (J·s).
pub const HBAR: Quantity<f64, crate::dim::Energy> = Quantity::new(1.054571817e-34);

/// Planck constant h = 2πℏ (J·s).
pub const H: Quantity<f64, crate::dim::Energy> = Quantity::new(6.62607015e-34);

/// Boltzmann constant k_B (J/K).
pub const K_B: Quantity<f64, crate::dim::Energy> = Quantity::new(1.380649e-23);

/// Electron mass (kg).
pub const ELECTRON_MASS: Kilograms = Kilograms::new(9.1093837015e-31);

/// Proton mass (kg).
pub const PROTON_MASS: Kilograms = Kilograms::new(1.67262192369e-27);

/// Elementary charge (C).
pub const ELEMENTARY_CHARGE: Charge = Charge::new(1.602176634e-19);

/// Gravitational constant G (m³·kg⁻¹·s⁻²).
pub const G: crate::GConstant = crate::GConstant::new(6.67430e-11);

/// Avogadro number (mol⁻¹).
pub const N_A: Quantity<f64, crate::dim::Wavenumber> = Quantity::new(6.02214076e23);

/// Stefan-Boltzmann constant σ (W·m⁻²·K⁻⁴).
/// σ = 2π⁵k_B⁴ / (15 h³ c²)
pub const SIGMA: Quantity<f64, crate::dim::Power> = Quantity::new(5.670374419e-8);

/// Astronomical unit (m).
pub const AU: Meters = Meters::new(1.495978707e11);

/// Standard gravity g (m/s²).
pub const G_EARTH: crate::Acceleration = crate::Acceleration::new(9.80665);

/// Earth equatorial radius, WGS84 (m).
pub const R_EARTH_EQ: Meters = Meters::new(6_378_137.0);

/// Earth gravitational parameter GM (m³/s²).
pub const GM_EARTH: crate::GravitationalParameter =
    crate::GravitationalParameter::new(3.986004415e14);

/// Sun gravitational parameter GM (m³/s²).
pub const GM_SUN: crate::GravitationalParameter =
    crate::GravitationalParameter::new(1.32712440018e20);

/// Moon gravitational parameter GM (m³/s²).
pub const GM_MOON: crate::GravitationalParameter =
    crate::GravitationalParameter::new(4.902800118e12);

/// Pi (dimensionless).
pub const PI: Dimensionless = Dimensionless::new(core::f64::consts::PI);

/// Two pi (dimensionless).
pub const TAU: Dimensionless = Dimensionless::new(core::f64::consts::TAU);
