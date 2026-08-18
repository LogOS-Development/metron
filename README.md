# Metron

*μέτρον — measure, proportion*

Compile-time SI units for Rust. Zero-cost dimensional analysis.
## Quick start

```toml
[dependencies]
metron = "0.1"
```

```rust
use metron::si::*;
use metron::pow;
use metron::{Velocity, Acceleration, Force, Energy, Area, Volume, Frequency};

// Natural literal syntax — 5.0 * (m / s) is a velocity
let v: Velocity<f64> = 5.0 * (m / s);

// pow! macro for unit powers — compile-time checked
let a: Acceleration<f64> = 9.8 * (m / pow!(s, 2));
let area: Area<f64> = 3.0 * pow!(m, 2);
let vol: Volume<f64> = 2.0 * pow!(m, 3);
let freq: Frequency<f64> = 1.0 / s;          // Hz = s⁻¹

// Derived units compose naturally
let force: Force<f64> = 10.0 * (kg * m / pow!(s, 2));     // N
let energy: Energy<f64> = 4.2 * (kg * pow!(m, 2) / pow!(s, 2));  // J
```

## What it does

Unit mismatches are compile errors, not runtime panics:

```rust
// This does NOT compile:
let wrong: Force<f64> = mass * velocity;  // error: kg*(m/s) ≠ kg*m/s²
```

`Quantity<T, U>` is a newtype around `f64` with `PhantomData<U>` — zero runtime cost. The unit type `U` is a 7-tuple of type-level signed integer exponents over `[m, kg, s, A, K, mol, cd]`, tracked at compile time via `typenum`.

## Features

- **Compile-time checking** — wrong unit combinations are type errors
- **`pow!` macro** — `pow!(s, 2)` for s², `pow!(s, -1)` for s⁻¹, range -10 to +10
- **Natural literal syntax** — `5.0 * (m / s)`, `1.0 / s`, `9.8 * (m / pow!(s, 2))`
- **Scalar, vector, tensor** — `Quantity<T, U>`, `VectorQuantity<T, N, U>`, `TensorQuantity<T, M, N, U>`
- **`no_std` support** — works in embedded contexts
- **Serde support** — optional, behind a feature flag
- **40+ named SI type aliases** — `Meters`, `Velocity`, `Force`, `Energy`, `Pressure`, etc.
- **Complex scalars** — `ComplexMeters`, `ComplexVolts` for phasor-domain analysis
- **SI prefix conversion** — `SiPrefix::Kilo`, `SiPrefix::Milli`, etc.

## How it works

Each SI dimension is a type-level integer exponent in a 7-tuple:

```
Meter        = Unit<(P1, Z0, Z0, Z0, Z0, Z0, Z0)>   // m
Second       = Unit<(Z0, Z0, P1, Z0, Z0, Z0, Z0)>   // s
Velocity     = Unit<(P1, Z0, N1, Z0, Z0, Z0, Z0)>   // m/s = m¹·s⁻¹
Force        = Unit<(P1, P1, N2, Z0, Z0, Z0, Z0)>   // N  = m¹·kg¹·s⁻²
```

Multiplication adds exponents, division subtracts them, `pow!` scales them — all at compile time. `sqrt` halves exponents (only on even exponents; odd ones are a compile error).

## Comparison

| Feature | Metron | uom | dimensional_quantity | danwi |
|---|---|---|---|---|
| Stable Rust | yes | yes | **nightly** | yes |
| Compile-time checking | yes | yes | yes | yes |
| `pow!` macro | **yes** | no | no | no |
| Power syntax | `pow!(s, 2)` | `powi(P2::new())` | `powi::<2>()` | none |
| `5.0 * (m / s)` syntax | **yes** | no | no | no |
| Zero-cost | yes | yes | yes | yes |
| `no_std` | yes | yes | yes | yes |

## License

MIT OR Apache-2.0
