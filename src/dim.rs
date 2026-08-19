//! Type-level unit tags — zero-sized phantom types encoding SI dimensions.
//!
//! Each alias is a [`Unit<(...)>`](crate::Unit) with type-level integer
//! exponents over `[m, kg, s, A, K, mol, cd]`.  For example, `Meter` is
//! `Unit<(P1, Z0, Z0, Z0, Z0, Z0, Z0)>` — exponent +1 on length, zero
//! elsewhere.
//!
//! These are the **unit tags** only.  The public API wraps them in
//! [`Quantity<T, U>`](crate::Quantity) via aliases like [`Meters`](crate::Meters),
//! [`Velocity`](crate::Velocity), etc.

use super::unit::Unit;
use typenum::consts::*;
use typenum::Z0;

// --- SI base units ---
pub type Meter = Unit<(P1, Z0, Z0, Z0, Z0, Z0, Z0)>;
pub type Kilogram = Unit<(Z0, P1, Z0, Z0, Z0, Z0, Z0)>;
pub type Second = Unit<(Z0, Z0, P1, Z0, Z0, Z0, Z0)>;
pub type Ampere = Unit<(Z0, Z0, Z0, P1, Z0, Z0, Z0)>;
pub type Kelvin = Unit<(Z0, Z0, Z0, Z0, P1, Z0, Z0)>;
pub type Mole = Unit<(Z0, Z0, Z0, Z0, Z0, P1, Z0)>;
pub type Candela = Unit<(Z0, Z0, Z0, Z0, Z0, Z0, P1)>;
pub type Dimensionless = Unit<(Z0, Z0, Z0, Z0, Z0, Z0, Z0)>;

// --- Derived units (m·kg⁻¹·s²·A⁻¹·K⁻¹·mol⁻¹·cd⁻¹) ---
pub type Velocity = Unit<(P1, Z0, N1, Z0, Z0, Z0, Z0)>;
pub type Acceleration = Unit<(P1, Z0, N2, Z0, Z0, Z0, Z0)>;
pub type Force = Unit<(P1, P1, N2, Z0, Z0, Z0, Z0)>;
pub type Energy = Unit<(P2, P1, N2, Z0, Z0, Z0, Z0)>;
pub type Torque = Energy;
pub type Power = Unit<(P2, P1, N3, Z0, Z0, Z0, Z0)>;
pub type Pressure = Unit<(N1, P1, N2, Z0, Z0, Z0, Z0)>;
pub type Area = Unit<(P2, Z0, Z0, Z0, Z0, Z0, Z0)>;
pub type Volume = Unit<(P3, Z0, Z0, Z0, Z0, Z0, Z0)>;
pub type Density = Unit<(N3, P1, Z0, Z0, Z0, Z0, Z0)>;
pub type Frequency = Unit<(Z0, Z0, N1, Z0, Z0, Z0, Z0)>;
pub type AngularVelocity = Frequency;
pub type Charge = Unit<(Z0, Z0, P1, P1, Z0, Z0, Z0)>;
pub type Voltage = Unit<(P2, P1, N3, N1, Z0, Z0, Z0)>;
pub type Resistance = Unit<(P2, P1, N3, N2, Z0, Z0, Z0)>;
pub type Capacitance = Unit<(N2, N1, P3, P2, Z0, Z0, Z0)>;
pub type Inductance = Unit<(P2, P1, N2, N2, Z0, Z0, Z0)>;
pub type MagneticFlux = Unit<(P2, P1, N2, N1, Z0, Z0, Z0)>;
pub type MagneticFluxDensity = Unit<(Z0, P1, N2, N1, Z0, Z0, Z0)>;
pub type GravitationalParameter = Unit<(P3, Z0, N2, Z0, Z0, Z0, Z0)>;
/// Gravitational constant G: m³ kg⁻¹ s⁻².
/// `Kilograms * GConstant = GravitationalParameter`.
pub type GConstant = Unit<(P3, N1, N2, Z0, Z0, Z0, Z0)>;
pub type MomentOfInertia = Unit<(P2, P1, Z0, Z0, Z0, Z0, Z0)>;
pub type SpecificEnergy = Unit<(P2, Z0, N2, Z0, Z0, Z0, Z0)>;
pub type SpringConstant = Unit<(Z0, P1, N2, Z0, Z0, Z0, Z0)>;
pub type DampingCoefficient = Unit<(Z0, P1, N1, Z0, Z0, Z0, Z0)>;
pub type Angle = Dimensionless;
pub type SolidAngle = Dimensionless;
pub type AngularAcceleration = Unit<(Z0, Z0, N2, Z0, Z0, Z0, Z0)>;
pub type MassFlowRate = Unit<(Z0, P1, N1, Z0, Z0, Z0, Z0)>;
pub type SpecificImpulse = Second;
pub type Wavenumber = Unit<(N1, Z0, Z0, Z0, Z0, Z0, Z0)>;
pub type ProbabilityDensity = Wavenumber;
