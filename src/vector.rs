//! Vector quantity — an N-component vector tagged with a compile-time unit.

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};

use nalgebra::SVector;
use num_traits::{NumAssign, Zero};

use super::dim;
use super::quantity::Quantity;
use super::unit_name::UnitName;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// An `N`-component vector tagged with a compile-time unit.
///
/// Wraps `nalgebra::SVector<T, N>`.  `U` is the unit tag shared by all
/// components.  Arithmetic follows the same dimensional rules as
/// [`Quantity`](crate::Quantity): `Add`/`Sub` require matching units, `Mul`/`Div` by a
/// scalar preserves units, and `dot`/`cross` produce unit products.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
#[cfg_attr(
    feature = "serde",
    serde(bound(
        serialize = "SVector<T, N>: Serialize",
        deserialize = "SVector<T, N>: Deserialize<'de>"
    ))
)]
pub struct VectorQuantity<T, const N: usize, U> {
    /// The raw nalgebra vector, stored in SI base units.
    pub vector: SVector<T, N>,
    _u: PhantomData<U>,
}

impl<T: Zero + Clone + nalgebra::Scalar, const N: usize, U> Default for VectorQuantity<T, N, U> {
    fn default() -> Self {
        Self {
            vector: SVector::zeros(),
            _u: PhantomData,
        }
    }
}

impl<T, const N: usize, U> core::ops::Deref for VectorQuantity<T, N, U> {
    type Target = SVector<T, N>;
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.vector
    }
}

impl<T, const N: usize, U> core::ops::DerefMut for VectorQuantity<T, N, U> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.vector
    }
}

impl<T, const N: usize, U> VectorQuantity<T, N, U> {
    /// Creates a vector quantity from a raw nalgebra vector.
    #[inline]
    #[must_use]
    pub const fn new(vector: SVector<T, N>) -> Self {
        Self {
            vector,
            _u: PhantomData,
        }
    }
    /// Borrows the raw nalgebra vector.
    #[inline]
    #[must_use]
    pub const fn vector(&self) -> &SVector<T, N> {
        &self.vector
    }
    /// Alias for [`vector`](Self::vector).
    #[inline]
    #[must_use]
    pub const fn value(&self) -> &SVector<T, N> {
        &self.vector
    }
    /// Alias for [`vector`](Self::vector).
    #[inline]
    #[must_use]
    pub const fn raw(&self) -> &SVector<T, N> {
        &self.vector
    }
    /// Consumes the vector quantity, returning the raw nalgebra vector.
    #[inline]
    #[must_use]
    pub fn into_vector(self) -> SVector<T, N> {
        self.vector
    }
}

impl<T: NumAssign + Clone + nalgebra::Scalar + nalgebra::ComplexField, const N: usize, U>
    VectorQuantity<T, N, U>
{
    /// Euclidean norm (magnitude) of the vector.
    #[inline]
    #[must_use]
    pub fn norm(&self) -> Quantity<T::RealField, U> {
        Quantity::new(self.vector.norm())
    }

    /// Squared Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm_squared(&self) -> Quantity<T::RealField, <U as Mul<U>>::Output>
    where
        U: Mul<U>,
    {
        Quantity::new(self.vector.norm_squared())
    }

    /// Unit vector in the same direction (dimensionless).
    #[inline]
    #[must_use]
    pub fn normalize(&self) -> VectorQuantity<T, N, dim::Dimensionless> {
        VectorQuantity::new(self.vector.normalize())
    }

    /// Dot product with another vector.
    #[inline]
    #[must_use]
    pub fn dot<B>(&self, other: &VectorQuantity<T, N, B>) -> Quantity<T, <U as Mul<B>>::Output>
    where
        U: Mul<B>,
    {
        Quantity::new(self.vector.dot(&other.vector))
    }
}

impl<T: NumAssign + Clone + nalgebra::Scalar + nalgebra::ComplexField, U> VectorQuantity<T, 3, U> {
    /// Cross product with another 3D vector.
    #[inline]
    #[must_use]
    pub fn cross<B>(
        &self,
        other: &VectorQuantity<T, 3, B>,
    ) -> VectorQuantity<T, 3, <U as Mul<B>>::Output>
    where
        U: Mul<B>,
    {
        VectorQuantity::new(self.vector.cross(&other.vector))
    }

    /// Construct from `[x, y, z]` components.
    #[inline]
    #[must_use]
    pub fn from_xyz(x: T, y: T, z: T) -> Self {
        Self::new(nalgebra::Vector3::new(x, y, z))
    }
}

impl<T: NumAssign + Clone + nalgebra::Scalar, const N: usize, U> Add for VectorQuantity<T, N, U> {
    type Output = Self;
    #[inline]
    fn add(self, r: Self) -> Self {
        Self::new(self.vector + r.vector)
    }
}
impl<T: NumAssign + Clone + nalgebra::Scalar, const N: usize, U> Sub for VectorQuantity<T, N, U> {
    type Output = Self;
    #[inline]
    fn sub(self, r: Self) -> Self {
        Self::new(self.vector - r.vector)
    }
}
impl<T: NumAssign + Clone + nalgebra::Scalar, const N: usize, U> Neg for VectorQuantity<T, N, U>
where
    SVector<T, N>: Neg<Output = SVector<T, N>>,
{
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.vector)
    }
}
impl<T: NumAssign + Clone + nalgebra::Scalar, const N: usize, U> Mul<T>
    for VectorQuantity<T, N, U>
{
    type Output = Self;
    #[inline]
    fn mul(self, r: T) -> Self {
        Self::new(self.vector * r)
    }
}
impl<T: NumAssign + Clone + nalgebra::Scalar, const N: usize, U> Div<T>
    for VectorQuantity<T, N, U>
{
    type Output = Self;
    #[inline]
    fn div(self, r: T) -> Self {
        Self::new(self.vector / r)
    }
}

impl<T: NumAssign + Clone + nalgebra::Scalar, const N: usize, UA: Mul<UB>, UB>
    Mul<VectorQuantity<T, N, UB>> for Quantity<T, UA>
{
    type Output = VectorQuantity<T, N, <UA as Mul<UB>>::Output>;
    #[inline]
    fn mul(self, r: VectorQuantity<T, N, UB>) -> Self::Output {
        VectorQuantity::new(r.vector * self.value)
    }
}

impl<T, const N: usize, U: UnitName> fmt::Display for VectorQuantity<T, N, U>
where
    T: fmt::Display + Clone + nalgebra::Scalar + PartialEq + fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {}", self.vector, U::NAME)
    }
}
