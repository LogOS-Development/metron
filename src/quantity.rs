//! Scalar quantity — a value tagged with a compile-time unit.

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, AddAssign, Deref, DerefMut, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use num_traits::{Float, NumAssign, Zero};

use super::dim;
use super::prefix::SiPrefix;
use super::unit::{PowMap, Sqrt, UnitPow};
use super::unit_name::UnitName;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A scalar value tagged with a compile-time unit.
///
/// `: `T` is the scalar type (defaults to `f64`; use `Complex<f64>` for
/// phasor domains).  `U` is a [`Unit`](crate::Unit) type tag from the [`dim`](crate::dim) module.
///
/// Arithmetic operators enforce dimensional correctness at compile time:
/// - `Add`/`Sub` require both operands to have the **same unit** `U`.
/// - `Mul`/`Div` between quantities produce a new quantity whose unit is
///   the product/quotient of the operands' units (type-level exponent
///   arithmetic).
/// - `Mul`/`Div` by a raw `T` scalar preserves the unit unchanged.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
#[cfg_attr(
    feature = "serde",
    serde(bound(serialize = "T: Serialize", deserialize = "T: Deserialize<'de>"))
)]
pub struct Quantity<T, U> {
    /// The raw scalar value in SI base units.
    pub value: T,
    _u: PhantomData<U>,
}

impl<T: NumAssign, U> Zero for Quantity<T, U> {
    #[inline]
    fn zero() -> Self {
        Self {
            value: T::zero(),
            _u: PhantomData,
        }
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
}

impl<T: Default, U> Default for Quantity<T, U> {
    fn default() -> Self {
        Self {
            value: T::default(),
            _u: PhantomData,
        }
    }
}

impl<T, U> Quantity<T, U> {
    /// Creates a quantity from a raw scalar value, inferring the unit type `U`.
    #[inline]
    #[must_use]
    pub const fn new(value: T) -> Self {
        Self {
            value,
            _u: PhantomData,
        }
    }
    /// Borrows the raw scalar value.
    #[inline]
    #[must_use]
    pub const fn value(&self) -> &T {
        &self.value
    }
    /// Consumes the quantity, returning the raw scalar value.
    #[inline]
    #[must_use]
    pub fn into_value(self) -> T {
        self.value
    }
    /// Applies a function to the scalar value, preserving the unit tag `U`.
    #[inline]
    #[must_use]
    pub fn map<F, R>(self, f: F) -> Quantity<R, U>
    where
        F: FnOnce(T) -> R,
    {
        Quantity {
            value: f(self.value),
            _u: PhantomData,
        }
    }

    /// Square root of a scalar quantity.
    ///
    /// Produces a `Quantity<T, SqrtU>` where `SqrtU` is the type-level
    /// square root of `U`.  Only callable when every exponent in `U` is
    /// even — `sqrt(m²) = m`, `sqrt(m²/s²) = m/s`.
    /// Odd exponents (`sqrt(m)`) are rejected at compile time.
    #[inline]
    #[must_use]
    pub fn sqrt(self) -> Quantity<T, <U as Sqrt>::Output>
    where
        T: nalgebra::ComplexField,
        U: Sqrt,
    {
        Quantity::new(self.value.sqrt())
    }

    /// Raise a quantity to an integer power.
    ///
    /// Call with `s.pow::<2>()` for s², `m.pow::<3>()` for m³, `s.pow::<-1>()` for s⁻¹.
    ///
    /// For cleaner syntax, use the [`pow!`](crate::pow!) macro:
    /// `pow!(s, 2)` expands to `s.pow::<2>()`.
    ///
    /// Supported range: -10 to +10. Out-of-range powers produce a compile-time error.
    #[inline]
    #[must_use]
    pub fn pow<const N: i32>(self) -> Quantity<T, <U as UnitPow<<() as PowMap<N>>::T>>::Output>
    where
        T: Float,
        (): PowMap<N>,
        U: UnitPow<<() as PowMap<N>>::T>,
    {
        Quantity::new(self.value.powi(N))
    }
}

impl<T, U> Deref for Quantity<T, U> {
    type Target = T;
    #[inline]
    fn deref(&self) -> &T {
        &self.value
    }
}
impl<T, U> DerefMut for Quantity<T, U> {
    #[inline]
    fn deref_mut(&mut self) -> &mut T {
        &mut self.value
    }
}

/// Prefix conversion for `f64`-backed quantities.
pub trait ConvertPrefix {
    /// Returns the value expressed in the given prefix (e.g. 1000 m → 1.0 km).
    fn in_prefix(&self, prefix: SiPrefix) -> f64;
    /// Returns a new quantity with the value converted to the given prefix.
    fn convert_to(&self, prefix: SiPrefix) -> Self;
}
impl<U> ConvertPrefix for Quantity<f64, U> {
    #[inline]
    fn in_prefix(&self, prefix: SiPrefix) -> f64 {
        self.value / prefix.scale()
    }
    #[inline]
    fn convert_to(&self, _p: SiPrefix) -> Self {
        Self::new(self.value)
    }
}

// --- Add / Sub: same-unit only (enforced at compile time) ---

impl<T: NumAssign, U> Add for Quantity<T, U> {
    type Output = Self;
    #[inline]
    fn add(self, r: Self) -> Self {
        Self::new(self.value + r.value)
    }
}
impl<T: NumAssign, U> AddAssign for Quantity<T, U> {
    #[inline]
    fn add_assign(&mut self, r: Self) {
        self.value += r.value;
    }
}
impl<T: NumAssign, U> Sub for Quantity<T, U> {
    type Output = Self;
    #[inline]
    fn sub(self, r: Self) -> Self {
        Self::new(self.value - r.value)
    }
}
impl<T: NumAssign, U> SubAssign for Quantity<T, U> {
    #[inline]
    fn sub_assign(&mut self, r: Self) {
        self.value -= r.value;
    }
}
impl<T: NumAssign + Neg<Output = T>, U> Neg for Quantity<T, U> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.value)
    }
}

// --- Mul / Div between quantities: unit product/quotient at compile time ---

impl<T: NumAssign, A: Mul<B>, B> Mul<Quantity<T, B>> for Quantity<T, A> {
    type Output = Quantity<T, <A as Mul<B>>::Output>;
    #[inline]
    fn mul(self, r: Quantity<T, B>) -> Self::Output {
        Quantity::new(self.value * r.value)
    }
}
impl<T: NumAssign, A: Div<B>, B> Div<Quantity<T, B>> for Quantity<T, A> {
    type Output = Quantity<T, <A as Div<B>>::Output>;
    #[inline]
    fn div(self, r: Quantity<T, B>) -> Self::Output {
        Quantity::new(self.value / r.value)
    }
}
// --- Mul / Div by raw scalar: unit preserved ---
impl<T: NumAssign, U> Mul<T> for Quantity<T, U> {
    type Output = Self;
    #[inline]
    fn mul(self, r: T) -> Self {
        Self::new(self.value * r)
    }
}
impl<T: NumAssign, U> Div<T> for Quantity<T, U> {
    type Output = Self;
    #[inline]
    fn div(self, r: T) -> Self {
        Self::new(self.value / r)
    }
}
impl<T: NumAssign, U> MulAssign<T> for Quantity<T, U> {
    #[inline]
    fn mul_assign(&mut self, r: T) {
        self.value *= r;
    }
}

// --- Reverse scalar ops: f64 * Quantity, f64 / Quantity ---
impl<U> Mul<Quantity<f64, U>> for f64 {
    type Output = Quantity<f64, U>;
    #[inline]
    fn mul(self, rhs: Quantity<f64, U>) -> Self::Output {
        Quantity::new(self * rhs.value)
    }
}
impl<U> Div<Quantity<f64, U>> for f64
where
    dim::Dimensionless: Div<U>,
{
    type Output = Quantity<f64, <dim::Dimensionless as Div<U>>::Output>;
    #[inline]
    fn div(self, rhs: Quantity<f64, U>) -> Self::Output {
        Quantity::new(self / rhs.value)
    }
}

impl<T: fmt::Display, U: UnitName> fmt::Display for Quantity<T, U> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, U::NAME)
    }
}
