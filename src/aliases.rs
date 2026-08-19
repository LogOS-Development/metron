//! Quantity type aliases — convenient names for common `Quantity` instantiations.

use num_complex::Complex;

use super::dim;
use super::quantity::Quantity;
use super::vector::VectorQuantity;

// --- Base ---
pub type Meters<T = f64> = Quantity<T, dim::Meter>;
pub type Seconds<T = f64> = Quantity<T, dim::Second>;
pub type Kilograms<T = f64> = Quantity<T, dim::Kilogram>;
pub type Amperes<T = f64> = Quantity<T, dim::Ampere>;
pub type Kelvins<T = f64> = Quantity<T, dim::Kelvin>;
pub type Moles<T = f64> = Quantity<T, dim::Mole>;
pub type Candelas<T = f64> = Quantity<T, dim::Candela>;
pub type Dimensionless<T = f64> = Quantity<T, dim::Dimensionless>;
pub type Radians<T = f64> = Quantity<T, dim::Angle>;
pub type Steradians<T = f64> = Quantity<T, dim::SolidAngle>;

// --- Derived ---
pub type Velocity<T = f64> = Quantity<T, dim::Velocity>;
pub type Acceleration<T = f64> = Quantity<T, dim::Acceleration>;
pub type Force<T = f64> = Quantity<T, dim::Force>;
pub type Energy<T = f64> = Quantity<T, dim::Energy>;
pub type Power<T = f64> = Quantity<T, dim::Power>;
pub type Pressure<T = f64> = Quantity<T, dim::Pressure>;
pub type Area<T = f64> = Quantity<T, dim::Area>;
pub type Volume<T = f64> = Quantity<T, dim::Volume>;
pub type Density<T = f64> = Quantity<T, dim::Density>;
pub type Frequency<T = f64> = Quantity<T, dim::Frequency>;
pub type Charge<T = f64> = Quantity<T, dim::Charge>;
pub type Voltage<T = f64> = Quantity<T, dim::Voltage>;
pub type Resistance<T = f64> = Quantity<T, dim::Resistance>;
pub type Capacitance<T = f64> = Quantity<T, dim::Capacitance>;
pub type Inductance<T = f64> = Quantity<T, dim::Inductance>;
pub type MagneticFlux<T = f64> = Quantity<T, dim::MagneticFlux>;
pub type MagneticFluxDensity<T = f64> = Quantity<T, dim::MagneticFluxDensity>;
pub type GravitationalParameter<T = f64> = Quantity<T, dim::GravitationalParameter>;
pub type GConstant<T = f64> = Quantity<T, dim::GConstant>;
pub type MomentOfInertia<T = f64> = Quantity<T, dim::MomentOfInertia>;
pub type SpecificEnergy<T = f64> = Quantity<T, dim::SpecificEnergy>;
pub type SpringConstant<T = f64> = Quantity<T, dim::SpringConstant>;
pub type DampingCoefficient<T = f64> = Quantity<T, dim::DampingCoefficient>;
pub type Action<T = f64> = Quantity<T, dim::Action>;
pub type Wavenumber<T = f64> = Quantity<T, dim::Wavenumber>;
pub type ProbabilityDensity<T = f64> = Quantity<T, dim::ProbabilityDensity>;
pub type MassFlowRate<T = f64> = Quantity<T, dim::MassFlowRate>;

// --- 3D vector aliases ---
pub type PositionVector = VectorQuantity<f64, 3, dim::Meter>;
pub type VelocityVector = VectorQuantity<f64, 3, dim::Velocity>;
pub type AccelerationVector = VectorQuantity<f64, 3, dim::Acceleration>;
pub type ForceVector = VectorQuantity<f64, 3, dim::Force>;
pub type DimensionlessVector = VectorQuantity<f64, 3, dim::Dimensionless>;

// --- Prefixed (same unit tags, different semantic names) ---
pub type Kilometers<T = f64> = Quantity<T, dim::Meter>;
pub type Millimeters<T = f64> = Quantity<T, dim::Meter>;
pub type Nanometers<T = f64> = Quantity<T, dim::Meter>;
pub type Centimeters<T = f64> = Quantity<T, dim::Meter>;
pub type Milliseconds<T = f64> = Quantity<T, dim::Second>;
pub type Microseconds<T = f64> = Quantity<T, dim::Second>;
pub type Nanoseconds<T = f64> = Quantity<T, dim::Second>;
pub type Grams<T = f64> = Quantity<T, dim::Kilogram>;
pub type Milligrams<T = f64> = Quantity<T, dim::Kilogram>;
pub type Kilohertz<T = f64> = Quantity<T, dim::Frequency>;
pub type Megahertz<T = f64> = Quantity<T, dim::Frequency>;
pub type Gigahertz<T = f64> = Quantity<T, dim::Frequency>;
pub type Nanoteslas<T = f64> = Quantity<T, dim::MagneticFluxDensity>;
pub type Microteslas<T = f64> = Quantity<T, dim::MagneticFluxDensity>;
pub type Milliteslas<T = f64> = Quantity<T, dim::MagneticFluxDensity>;
pub type Millivolts<T = f64> = Quantity<T, dim::Voltage>;
pub type Kilovolts<T = f64> = Quantity<T, dim::Voltage>;
pub type Megavolts<T = f64> = Quantity<T, dim::Voltage>;
pub type Milliamperes<T = f64> = Quantity<T, dim::Ampere>;
pub type Kiloamperes<T = f64> = Quantity<T, dim::Ampere>;
pub type Kilopascals<T = f64> = Quantity<T, dim::Pressure>;
pub type Megapascals<T = f64> = Quantity<T, dim::Pressure>;
pub type Hectopascals<T = f64> = Quantity<T, dim::Pressure>;
pub type Kilojoules<T = f64> = Quantity<T, dim::Energy>;
pub type Megajoules<T = f64> = Quantity<T, dim::Energy>;
pub type Kilowatts<T = f64> = Quantity<T, dim::Power>;
pub type Megawatts<T = f64> = Quantity<T, dim::Power>;
pub type MilliKelvins<T = f64> = Quantity<T, dim::Kelvin>;
pub type Kilokelvins<T = f64> = Quantity<T, dim::Kelvin>;

// --- Complex ---
pub type ComplexMeters = Quantity<Complex<f64>, dim::Meter>;
pub type ComplexSeconds = Quantity<Complex<f64>, dim::Second>;
pub type ComplexVolts = Quantity<Complex<f64>, dim::Voltage>;
pub type ComplexAmperes = Quantity<Complex<f64>, dim::Ampere>;
pub type ComplexFrequency = Quantity<Complex<f64>, dim::Frequency>;
pub type ComplexWavenumber = Quantity<Complex<f64>, dim::Wavenumber>;
