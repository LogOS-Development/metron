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
