//! Tensor quantity — an M×N matrix tagged with a compile-time unit.

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Mul, Neg, Sub};

use nalgebra::SMatrix;
use num_traits::{NumAssign, Zero};

use super::unit_name::UnitName;
use super::vector::VectorQuantity;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// An `M`×`N` matrix tagged with a compile-time unit.
///
/// Wraps `nalgebra::SMatrix<T, M, N>`.  Used for inertia tensors, stress
/// tensors, transformation matrices, etc.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "serde", serde(transparent))]
#[cfg_attr(
    feature = "serde",
    serde(bound(
        serialize = "SMatrix<T, M, N>: Serialize",
        deserialize = "SMatrix<T, M, N>: Deserialize<'de>"
    ))
)]
pub struct TensorQuantity<T, const M: usize, const N: usize, U> {
    /// The raw nalgebra matrix, stored in SI base units.
    pub matrix: SMatrix<T, M, N>,
    _u: PhantomData<U>,
}

impl<T: Zero + Clone + nalgebra::Scalar, const M: usize, const N: usize, U> Default
    for TensorQuantity<T, M, N, U>
{
    fn default() -> Self {
        Self {
            matrix: SMatrix::zeros(),
            _u: PhantomData,
        }
    }
}

impl<T, const M: usize, const N: usize, U> TensorQuantity<T, M, N, U> {
    /// Creates a tensor quantity from a raw nalgebra matrix.
    #[inline]
    #[must_use]
    pub const fn new(matrix: SMatrix<T, M, N>) -> Self {
        Self {
            matrix,
            _u: PhantomData,
        }
    }
    /// Borrows the raw nalgebra matrix.
    #[inline]
    #[must_use]
    pub const fn matrix(&self) -> &SMatrix<T, M, N> {
        &self.matrix
    }
    /// Consumes the tensor quantity, returning the raw nalgebra matrix.
    #[inline]
    #[must_use]
    pub fn into_matrix(self) -> SMatrix<T, M, N> {
        self.matrix
    }

    /// Returns the SI symbol for this tensor's unit (e.g. "kg·m²", "Pa").
    #[inline]
    #[must_use]
    pub fn unit_name(&self) -> Option<&'static str>
    where
        U: UnitName,
    {
        Some(U::NAME)
    }
}

impl<
        T: NumAssign + Clone + nalgebra::Scalar + nalgebra::ComplexField,
        const M: usize,
        const N: usize,
        U,
    > TensorQuantity<T, M, N, U>
{
    /// Identity matrix of size `M`×`N`, tagged with unit `U`.
    #[inline]
    #[must_use]
    pub fn identity() -> Self {
        Self::new(SMatrix::identity())
    }
    /// Transpose, swapping dimensions `M`↔`N` and preserving the unit.
    #[inline]
    #[must_use]
    pub fn transpose(&self) -> TensorQuantity<T, N, M, U> {
        TensorQuantity::new(self.matrix.transpose())
    }
}

impl<T: NumAssign + Clone + nalgebra::Scalar, const M: usize, const N: usize, U> Add
    for TensorQuantity<T, M, N, U>
{
    type Output = Self;
    #[inline]
    fn add(self, r: Self) -> Self {
        Self::new(self.matrix + r.matrix)
    }
}
impl<T: NumAssign + Clone + nalgebra::Scalar, const M: usize, const N: usize, U> Sub
    for TensorQuantity<T, M, N, U>
{
    type Output = Self;
    #[inline]
    fn sub(self, r: Self) -> Self {
        Self::new(self.matrix - r.matrix)
    }
}
impl<T: NumAssign + Clone + nalgebra::Scalar, const M: usize, const N: usize, U> Neg
    for TensorQuantity<T, M, N, U>
where
    SMatrix<T, M, N>: Neg<Output = SMatrix<T, M, N>>,
{
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.matrix)
    }
}
impl<T: NumAssign + Clone + nalgebra::Scalar, const M: usize, const N: usize, U> Mul<T>
    for TensorQuantity<T, M, N, U>
{
    type Output = Self;
    #[inline]
    fn mul(self, r: T) -> Self {
        Self::new(self.matrix * r)
    }
}

impl<
        T: NumAssign + Clone + nalgebra::Scalar,
        const M: usize,
        const K: usize,
        const N: usize,
        UA: Mul<UB>,
        UB,
    > Mul<TensorQuantity<T, K, N, UB>> for TensorQuantity<T, M, K, UA>
{
    type Output = TensorQuantity<T, M, N, <UA as Mul<UB>>::Output>;
    #[inline]
    fn mul(self, r: TensorQuantity<T, K, N, UB>) -> Self::Output {
        TensorQuantity::new(self.matrix * r.matrix)
    }
}

impl<T: NumAssign + Clone + nalgebra::Scalar, const M: usize, const N: usize, UA: Mul<UB>, UB>
    Mul<VectorQuantity<T, N, UB>> for TensorQuantity<T, M, N, UA>
{
    type Output = VectorQuantity<T, M, <UA as Mul<UB>>::Output>;
    #[inline]
    fn mul(self, r: VectorQuantity<T, N, UB>) -> Self::Output {
        VectorQuantity::new(self.matrix * r.vector)
    }
}

impl<T, const M: usize, const N: usize, U: UnitName> fmt::Display for TensorQuantity<T, M, N, U>
where
    T: fmt::Display + Clone + nalgebra::Scalar + PartialEq + fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}] {}", self.matrix, U::NAME)
    }
}
