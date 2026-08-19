//! Metron — compile-time SI units for Rust.
//!
//! Zero-cost dimensional analysis with `pow!` macro syntax.
//! Stable Rust, no nightly required.
//!
//! # Quick example
//!
//! ```
//! use metron::si::*;
//! use metron::pow;
//! use metron::{Velocity, Acceleration, Force, Energy, Area};
//!
//! let v: Velocity<f64> = 5.0 * (m / s);
//! let acc: Acceleration<f64> = 9.8 * (m / pow!(s, 2));
//! let area: Area<f64> = 3.0 * pow!(m, 2);
//! ```

#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(feature = "serde")]
extern crate serde;

pub mod aliases;
pub mod constants;
pub mod dim;
pub mod prefix;
pub mod quantity;
pub mod si;
pub mod tensor;
pub mod unit;
pub mod unit_name;
pub mod vector;

// Re-export the core types at the crate root for convenience.
pub use aliases::*;
pub use prefix::SiPrefix;
pub use quantity::{ConvertPrefix, Quantity};
pub use tensor::TensorQuantity;
pub use unit::{Halvable, PowMap, Sqrt, Unit, UnitPow};
pub use unit_name::UnitName;
pub use vector::VectorQuantity;
// `dim` is already a pub mod, no need to re-export

/// Raise a quantity to an integer power.
///
/// `pow!(s, 2)` expands to `s.pow::<2>()`, producing s² at compile time.
/// `pow!(m, 3)` gives m³. `pow!(s, -1)` gives s⁻¹ (frequency).
///
/// Supported range: -10 to +10.
#[macro_export]
macro_rules! pow {
    ($q:expr, $n:literal) => {
        $q.pow::<$n>()
    };
}

/// Assert that a quantity has a specific unit type at compile time.
///
/// `assert_unit!(vel, Velocity)` expands to a type check — if `vel` is not
/// a `Quantity<_, Velocity>`, it's a compile error. This complements
/// `format!("{vel}")` which checks the runtime value + unit symbol.
#[macro_export]
macro_rules! assert_unit {
    ($q:expr, $unit:ty) => {
        let _: &$unit = &$q;
    };
}

#[cfg(test)]
mod tests;
