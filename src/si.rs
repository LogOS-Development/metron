//! SI base unit constants — 1 of each base unit.
//!
//! Enables natural unit-expression syntax:
//!   use metron::si::*;
//!   5.0 * (m / s)          → Velocity (m/s)
//!   9.8 * (m / pow!(s, 2)) → Acceleration (m/s²)
//!   3.0 * pow!(m, 2)       → Area (m²)
//!   1.0 / s                → Frequency (Hz = s⁻¹)
//!
//! The constants live in a `si` submodule so that `use metron::*` does not
//! pollute the caller's namespace with single-letter names.

use super::aliases::*;

/// 1 metre.
#[allow(non_upper_case_globals)]
pub const m: Meters<f64> = Meters::new(1.0);
/// 1 second.
#[allow(non_upper_case_globals)]
pub const s: Seconds<f64> = Seconds::new(1.0);
/// 1 kilogram.
#[allow(non_upper_case_globals)]
pub const kg: Kilograms<f64> = Kilograms::new(1.0);
/// 1 ampere.
#[allow(non_upper_case_globals)]
pub const a: Amperes<f64> = Amperes::new(1.0);
/// 1 kelvin.
#[allow(non_upper_case_globals)]
pub const k: Kelvins<f64> = Kelvins::new(1.0);
/// 1 mole.
#[allow(non_upper_case_globals)]
pub const mol: Moles<f64> = Moles::new(1.0);
/// 1 candela.
#[allow(non_upper_case_globals)]
pub const cd: Candelas<f64> = Candelas::new(1.0);
