//! Associates a display name (e.g. `"m"`, `"m/s"`, `"N"`) with a unit tag.

use super::dim;

/// Associates a display name (e.g. `"m"`, `"m/s"`, `"N"`) with a unit tag.
pub trait UnitName {
    /// The SI symbol for this unit.
    const NAME: &'static str;
}

impl UnitName for dim::Meter {
    const NAME: &'static str = "m";
}
impl UnitName for dim::Kilogram {
    const NAME: &'static str = "kg";
}
impl UnitName for dim::Second {
    const NAME: &'static str = "s";
}
impl UnitName for dim::Ampere {
    const NAME: &'static str = "A";
}
impl UnitName for dim::Kelvin {
    const NAME: &'static str = "K";
}
impl UnitName for dim::Mole {
    const NAME: &'static str = "mol";
}
impl UnitName for dim::Candela {
    const NAME: &'static str = "cd";
}
impl UnitName for dim::Dimensionless {
    const NAME: &'static str = "";
}
impl UnitName for dim::Velocity {
    const NAME: &'static str = "m/s";
}
impl UnitName for dim::Acceleration {
    const NAME: &'static str = "m/s²";
}
impl UnitName for dim::Force {
    const NAME: &'static str = "N";
}
impl UnitName for dim::Energy {
    const NAME: &'static str = "J";
}
impl UnitName for dim::Power {
    const NAME: &'static str = "W";
}
impl UnitName for dim::Pressure {
    const NAME: &'static str = "Pa";
}
impl UnitName for dim::Area {
    const NAME: &'static str = "m²";
}
impl UnitName for dim::Volume {
    const NAME: &'static str = "m³";
}
impl UnitName for dim::Density {
    const NAME: &'static str = "kg/m³";
}
impl UnitName for dim::Frequency {
    const NAME: &'static str = "Hz";
}
impl UnitName for dim::Charge {
    const NAME: &'static str = "C";
}
impl UnitName for dim::Voltage {
    const NAME: &'static str = "V";
}
impl UnitName for dim::Resistance {
    const NAME: &'static str = "Ω";
}
impl UnitName for dim::Capacitance {
    const NAME: &'static str = "F";
}
impl UnitName for dim::Inductance {
    const NAME: &'static str = "H";
}
impl UnitName for dim::MagneticFlux {
    const NAME: &'static str = "Wb";
}
impl UnitName for dim::MagneticFluxDensity {
    const NAME: &'static str = "T";
}
impl UnitName for dim::GravitationalParameter {
    const NAME: &'static str = "m³/s²";
}
impl UnitName for dim::GConstant {
    const NAME: &'static str = "m³/(kg·s²)";
}
impl UnitName for dim::MomentOfInertia {
    const NAME: &'static str = "kg·m²";
}
impl UnitName for dim::AngularAcceleration {
    const NAME: &'static str = "rad/s²";
}
impl UnitName for dim::MassFlowRate {
    const NAME: &'static str = "kg/s";
}
impl UnitName for dim::Wavenumber {
    const NAME: &'static str = "1/m";
}
