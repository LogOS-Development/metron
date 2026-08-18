//! Unit marker types and type-level arithmetic.
//!
//! A [`Unit`] is a zero-sized phantom type carrying a 7-tuple of type-level
//! signed integer exponents over the SI base units `[m, kg, s, A, K, mol, cd]`.
//! Multiplication and division of units are type-level via `Mul`/`Div` —
//! `Meter * Second` produces the correct exponent tuple at compile time
//! with zero runtime cost.

use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Sub};

use typenum::consts::*;
use typenum::{Diff, Sum, Z0};

// ===========================================================================
// UnitPow — raise a unit tuple to an integer power at the type level
// ===========================================================================

/// Raise a unit to an integer power `N` (a typenum signed integer).
///
/// Each exponent in the 7-tuple is multiplied by `N` using typenum's
/// `Mul` impl.  `pow(P2)` doubles all exponents, `pow(N1)` negates them, etc.
///
/// Used by [`Quantity::pow`](crate::Quantity::pow) to produce the correct
/// derived unit at compile time.
pub trait UnitPow<N> {
    /// The unit raised to the `N`-th power.
    type Output;
}

impl<M, Kg, S, A, K, Mol, Cd, N> UnitPow<N> for Unit<(M, Kg, S, A, K, Mol, Cd)>
where
    N: Copy,
    M: Mul<N>,
    Kg: Mul<N>,
    S: Mul<N>,
    A: Mul<N>,
    K: Mul<N>,
    Mol: Mul<N>,
    Cd: Mul<N>,
{
    type Output = Unit<(
        <M as Mul<N>>::Output,
        <Kg as Mul<N>>::Output,
        <S as Mul<N>>::Output,
        <A as Mul<N>>::Output,
        <K as Mul<N>>::Output,
        <Mol as Mul<N>>::Output,
        <Cd as Mul<N>>::Output,
    )>;
}

// ===========================================================================
// PowMap — map const i32 to typenum type for pow(N) syntax
// ===========================================================================

/// Maps a `const N: i32` to the corresponding typenum integer type.
///
/// This enables `quantity.pow::<2>()` syntax while still
/// computing the result unit at compile time via [`UnitPow`].
///
/// Implemented for -10 to +10 via a macro. Out-of-range values
/// produce a compile-time error (trait bound not satisfied).
pub trait PowMap<const N: i32> {
    /// The typenum integer type corresponding to `N`.
    type T;
}

macro_rules! impl_powmap {
    ($($n:literal => $t:ident),* $(,)?) => {
        $(
            impl PowMap<$n> for () { type T = $t; }
        )*
    };
}

impl_powmap! {
    0 => Z0,
    1 => P1, -1 => N1,
    2 => P2, -2 => N2,
    3 => P3, -3 => N3,
    4 => P4, -4 => N4,
    5 => P5, -5 => N5,
    6 => P6, -6 => N6,
    7 => P7, -7 => N7,
    8 => P8, -8 => N8,
    9 => P9, -9 => N9,
    10 => P10, -10 => N10,
}

// ===========================================================================
// Unit marker struct + type-level Mul / Div
// ===========================================================================

/// A unit is a 7-tuple of type-level signed integer exponents over the SI
/// base units `[m, kg, s, A, K, mol, cd]`.
///
/// `Mul` adds exponents, `Div` subtracts them — both at compile time.
/// `Unit` is zero-sized (`PhantomData`); it exists only in the type system.
///
/// Concrete unit tags are defined as type aliases in the [`dim`](crate::dim)
/// module, e.g. `dim::Meter = Unit<(P1, Z0, Z0, Z0, Z0, Z0, Z0)>`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Unit<T>(pub PhantomData<T>);

impl<M1, Kg1, S1, A1, K1, Mol1, Cd1, M2, Kg2, S2, A2, K2, Mol2, Cd2>
    Mul<Unit<(M2, Kg2, S2, A2, K2, Mol2, Cd2)>> for Unit<(M1, Kg1, S1, A1, K1, Mol1, Cd1)>
where
    M1: Add<M2>,
    Kg1: Add<Kg2>,
    S1: Add<S2>,
    A1: Add<A2>,
    K1: Add<K2>,
    Mol1: Add<Mol2>,
    Cd1: Add<Cd2>,
{
    type Output = Unit<(
        Sum<M1, M2>,
        Sum<Kg1, Kg2>,
        Sum<S1, S2>,
        Sum<A1, A2>,
        Sum<K1, K2>,
        Sum<Mol1, Mol2>,
        Sum<Cd1, Cd2>,
    )>;
    fn mul(self, _rhs: Unit<(M2, Kg2, S2, A2, K2, Mol2, Cd2)>) -> Self::Output {
        Unit(PhantomData)
    }
}

impl<M1, Kg1, S1, A1, K1, Mol1, Cd1, M2, Kg2, S2, A2, K2, Mol2, Cd2>
    Div<Unit<(M2, Kg2, S2, A2, K2, Mol2, Cd2)>> for Unit<(M1, Kg1, S1, A1, K1, Mol1, Cd1)>
where
    M1: Sub<M2>,
    Kg1: Sub<Kg2>,
    S1: Sub<S2>,
    A1: Sub<A2>,
    K1: Sub<K2>,
    Mol1: Sub<Mol2>,
    Cd1: Sub<Cd2>,
{
    type Output = Unit<(
        Diff<M1, M2>,
        Diff<Kg1, Kg2>,
        Diff<S1, S2>,
        Diff<A1, A2>,
        Diff<K1, K2>,
        Diff<Mol1, Mol2>,
        Diff<Cd1, Cd2>,
    )>;
    fn div(self, _rhs: Unit<(M2, Kg2, S2, A2, K2, Mol2, Cd2)>) -> Self::Output {
        Unit(PhantomData)
    }
}

// ===========================================================================
// Sqrt trait removed — not needed for physics. Use pow! and division instead.
// ===========================================================================
