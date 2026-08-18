//! SI metric prefix (yocto through yotta).

use core::fmt;

/// SI metric prefix (yocto through yotta).
///
/// Prefixes are a **construction/display concern**, not a type-level
/// concern.  All quantities store values in SI base units internally;
/// `in_prefix` converts for display only.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum SiPrefix {
    Yocto,
    Zepto,
    Atto,
    Femto,
    Pico,
    Nano,
    Micro,
    Milli,
    Centi,
    Deci,
    #[default]
    None,
    Deca,
    Hecto,
    Kilo,
    Mega,
    Giga,
    Tera,
    Peta,
    Exa,
    Zetta,
    Yotta,
}

impl SiPrefix {
    /// Scale factors indexed by enum discriminant (10⁻²⁴ … 10²⁴).
    pub const SCALES: [f64; 21] = [
        1.0e-24, 1.0e-21, 1.0e-18, 1.0e-15, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-3, 1.0e-2, 1.0e-1, 1.0,
        1.0e1, 1.0e2, 1.0e3, 1.0e6, 1.0e9, 1.0e12, 1.0e15, 1.0e18, 1.0e21, 1.0e24,
    ];
    /// Returns the multiplicative scale factor (e.g. `Kilo → 1000.0`).
    #[inline]
    #[must_use]
    pub const fn scale(self) -> f64 {
        Self::SCALES[self as usize]
    }
}

impl fmt::Display for SiPrefix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Yocto => "y",
            Self::Zepto => "z",
            Self::Atto => "a",
            Self::Femto => "f",
            Self::Pico => "p",
            Self::Nano => "n",
            Self::Micro => "µ",
            Self::Milli => "m",
            Self::Centi => "c",
            Self::Deci => "d",
            Self::None => "",
            Self::Deca => "da",
            Self::Hecto => "h",
            Self::Kilo => "k",
            Self::Mega => "M",
            Self::Giga => "G",
            Self::Tera => "T",
            Self::Peta => "P",
            Self::Exa => "E",
            Self::Zetta => "Z",
            Self::Yotta => "Y",
        })
    }
}
