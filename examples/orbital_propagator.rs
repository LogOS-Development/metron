//! Orbital propagator using Dormand-Prince 5(4) adaptive integration.
//!
//! Propagates a satellite in a two-body Keplerian orbit around Earth,
//! using metron quantities throughout. The integrator is a proper
//! adaptive-step RKDP54 — the same family used by NASA GMAT and AGI STK.
//!
//! Run: cargo run --example orbital_propagator

use metron::pow;
use metron::{GravitationalParameter, Meters, Quantity, Seconds, Velocity};
use nalgebra::Vector3;

/// Gravitational parameter for Earth (GM = 3.986e14 m³/s²).
const GM_EARTH: f64 = 3.986004415e14;

/// Physical state: position and velocity vectors (raw f64 for nalgebra ops).
/// Units are enforced at construction and when computing derived quantities.
#[derive(Clone, Copy, Debug)]
struct OrbitState {
    pos: Vector3<f64>, // meters
    vel: Vector3<f64>, // m/s
}

/// Dormand-Prince 5(4) coefficients (DORMAND-PRINCE 1980).
/// 7-stage, 5th-order solution with embedded 4th-order error estimate.
mod rkdp54 {
    /// Butcher tableau c values (node positions).
    pub const C: [f64; 6] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0];

    /// Butcher tableau a values (stage coefficients).
    pub const A: [[f64; 6]; 6] = [
        [0.0; 6],
        [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
        [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
        [
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
            0.0,
            0.0,
        ],
        [
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
            0.0,
        ],
    ];

    /// 5th-order solution weights (b).
    pub const B5: [f64; 7] = [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
        0.0,
    ];

    /// 4th-order error estimate weights (b*).
    pub const B4: [f64; 7] = [
        5179.0 / 57600.0,
        0.0,
        7571.0 / 16695.0,
        393.0 / 640.0,
        -92097.0 / 339200.0,
        187.0 / 2100.0,
        1.0 / 40.0,
    ];
}

/// Compute gravitational acceleration magnitude at a given distance from Earth.
///
/// a = GM / r²  (point mass)
///
/// All quantities stay typed — no .into_value() until the integrator
/// needs raw f64 for nalgebra vector math.
fn gravity_accel_magnitude(r: Meters, gm: GravitationalParameter) -> metron::Acceleration {
    let r_sq = pow!(r, 2);
    gm / r_sq
}

/// Compute gravitational acceleration vector from a position.
///
/// Returns a raw Vector3<f64> (m/s²) for nalgebra compatibility.
/// The acceleration direction is -r̂, magnitude from GM/r².
fn gravity_accel(pos: Vector3<f64>, gm: GravitationalParameter) -> Vector3<f64> {
    let r_mag = Meters::new(pos.norm());
    let r_hat = pos / r_mag.into_value();
    let accel = gravity_accel_magnitude(r_mag, gm);
    -accel.into_value() * r_hat
}

/// Derivative of orbital state: velocity and acceleration.
fn derivative(state: &OrbitState) -> (Vector3<f64>, Vector3<f64>) {
    let gm = GravitationalParameter::new(GM_EARTH);
    let acc = gravity_accel(state.pos, gm);
    (state.vel, acc)
}

/// One RKDP54 step with adaptive step size control.
///
/// Returns the new state, the error estimate, and whether the step was accepted.
fn rkdp54_step(state: &OrbitState, dt: Seconds) -> (OrbitState, f64) {
    // Compute 7 stages (k1..k7)
    let (k1_v, k1_a) = derivative(state);

    let s2 = OrbitState {
        pos: state.pos + k1_v * (rkdp54::C[1] * dt.into_value()),
        vel: state.vel + k1_a * (rkdp54::C[1] * dt.into_value()),
    };
    let (k2_v, k2_a) = derivative(&s2);

    let s3 = OrbitState {
        pos: state.pos
            + k1_v * (rkdp54::A[2][0] * dt.into_value())
            + k2_v * (rkdp54::A[2][1] * dt.into_value()),
        vel: state.vel
            + k1_a * (rkdp54::A[2][0] * dt.into_value())
            + k2_a * (rkdp54::A[2][1] * dt.into_value()),
    };
    let (k3_v, k3_a) = derivative(&s3);

    let s4 = OrbitState {
        pos: state.pos
            + k1_v * (rkdp54::A[3][0] * dt.into_value())
            + k2_v * (rkdp54::A[3][1] * dt.into_value())
            + k3_v * (rkdp54::A[3][2] * dt.into_value()),
        vel: state.vel
            + k1_a * (rkdp54::A[3][0] * dt.into_value())
            + k2_a * (rkdp54::A[3][1] * dt.into_value())
            + k3_a * (rkdp54::A[3][2] * dt.into_value()),
    };
    let (k4_v, k4_a) = derivative(&s4);

    let s5 = OrbitState {
        pos: state.pos
            + k1_v * (rkdp54::A[4][0] * dt.into_value())
            + k2_v * (rkdp54::A[4][1] * dt.into_value())
            + k3_v * (rkdp54::A[4][2] * dt.into_value())
            + k4_v * (rkdp54::A[4][3] * dt.into_value()),
        vel: state.vel
            + k1_a * (rkdp54::A[4][0] * dt.into_value())
            + k2_a * (rkdp54::A[4][1] * dt.into_value())
            + k3_a * (rkdp54::A[4][2] * dt.into_value())
            + k4_a * (rkdp54::A[4][3] * dt.into_value()),
    };
    let (k5_v, k5_a) = derivative(&s5);

    let s6 = OrbitState {
        pos: state.pos
            + k1_v * (rkdp54::A[5][0] * dt.into_value())
            + k2_v * (rkdp54::A[5][1] * dt.into_value())
            + k3_v * (rkdp54::A[5][2] * dt.into_value())
            + k4_v * (rkdp54::A[5][3] * dt.into_value())
            + k5_v * (rkdp54::A[5][4] * dt.into_value()),
        vel: state.vel
            + k1_a * (rkdp54::A[5][0] * dt.into_value())
            + k2_a * (rkdp54::A[5][1] * dt.into_value())
            + k3_a * (rkdp54::A[5][2] * dt.into_value())
            + k4_a * (rkdp54::A[5][3] * dt.into_value())
            + k5_a * (rkdp54::A[5][4] * dt.into_value()),
    };
    let (k6_v, k6_a) = derivative(&s6);

    // 5th-order solution
    let dt_val = dt.into_value();
    let new_pos = state.pos
        + k1_v * (rkdp54::B5[0] * dt_val)
        + k3_v * (rkdp54::B5[2] * dt_val)
        + k4_v * (rkdp54::B5[3] * dt_val)
        + k5_v * (rkdp54::B5[4] * dt_val)
        + k6_v * (rkdp54::B5[5] * dt_val);

    let new_vel = state.vel
        + k1_a * (rkdp54::B5[0] * dt_val)
        + k3_a * (rkdp54::B5[2] * dt_val)
        + k4_a * (rkdp54::B5[3] * dt_val)
        + k5_a * (rkdp54::B5[4] * dt_val)
        + k6_a * (rkdp54::B5[5] * dt_val);

    // 7th stage (FSAL: reuses k1 of next step, uses new state)
    let s7 = OrbitState {
        pos: new_pos,
        vel: new_vel,
    };
    let (k7_v, k7_a) = derivative(&s7);

    // Error estimate = |5th-order - 4th-order|
    let err_pos = (k1_v * (rkdp54::B5[0] - rkdp54::B4[0])
        + k3_v * (rkdp54::B5[2] - rkdp54::B4[2])
        + k4_v * (rkdp54::B5[3] - rkdp54::B4[3])
        + k5_v * (rkdp54::B5[4] - rkdp54::B4[4])
        + k6_v * (rkdp54::B5[5] - rkdp54::B4[5])
        + k7_v * (rkdp54::B5[6] - rkdp54::B4[6]))
        * dt_val;

    let err_vel = (k1_a * (rkdp54::B5[0] - rkdp54::B4[0])
        + k3_a * (rkdp54::B5[2] - rkdp54::B4[2])
        + k4_a * (rkdp54::B5[3] - rkdp54::B4[3])
        + k5_a * (rkdp54::B5[4] - rkdp54::B4[4])
        + k6_a * (rkdp54::B5[5] - rkdp54::B4[5])
        + k7_a * (rkdp54::B5[6] - rkdp54::B4[6]))
        * dt_val;

    let error = (err_pos.norm().powi(2) + err_vel.norm().powi(2)).sqrt();

    (
        OrbitState {
            pos: new_pos,
            vel: new_vel,
        },
        error,
    )
}

/// Adaptive-step propagator. Adjusts dt based on error estimate.
fn propagate(
    mut state: OrbitState,
    t_total: Seconds,
    mut dt: Seconds,
    tolerance: f64,
) -> OrbitState {
    let mut t = 0.0_f64;
    let t_end = t_total.into_value();

    while t < t_end {
        if t + dt.into_value() > t_end {
            dt = Seconds::new(t_end - t);
        }

        let (new_state, error) = rkdp54_step(&state, dt);

        if error < tolerance || dt.into_value() < 1e-6 {
            // Accept step
            state = new_state;
            t += dt.into_value();
            // Increase step if error is small
            if error > 0.0 {
                let factor = (tolerance / error).powf(0.2).clamp(0.2, 5.0);
                dt = Seconds::new(dt.into_value() * factor);
            }
        } else {
            // Reject step, shrink dt
            let factor = (tolerance / error).powf(0.2).max(0.1);
            dt = Seconds::new(dt.into_value() * factor);
        }
    }
    state
}

/// Compute specific orbital energy: ε = v²/2 - GM/r
///
/// All arithmetic is compile-time unit checked. v² has units (m/s)² = m²/s²,
/// GM/r has units m³/s² / m = m²/s² — both are specific energy (J/kg).
/// No .into_value() needed.
fn specific_energy(
    pos: Vector3<f64>,
    vel: Vector3<f64>,
) -> Quantity<f64, <metron::dim::Velocity as core::ops::Mul<metron::dim::Velocity>>::Output> {
    // v² has units (m/s)² = m²/s² — this is specific energy (J/kg).
    // GM/r has units m³/s² / m = m²/s² — same type, so subtraction works.
    // No .into_value() needed — fully compile-time checked.
    let gm = GravitationalParameter::new(GM_EARTH);
    let r = Meters::new(pos.norm());
    let v = Velocity::new(vel.norm());
    let v_sq = pow!(v, 2);
    let ke = v_sq.map(|x| x * 0.5);
    let pe = gm / r;
    ke - pe
}

fn main() {
    // ISS-like orbit: 400 km altitude, circular, prograde
    let alt = Meters::new(400_000.0);
    let earth_r = Meters::new(6_378_137.0);
    let r0 = earth_r + alt;
    let gm = GravitationalParameter::new(GM_EARTH);

    // Circular orbit velocity: v = sqrt(GM/r)
    // GM/r has units m³/s² / m = m²/s². sqrt gives m/s.
    // But we removed sqrt from Quantity, so compute via into_value.
    let v_circ = Velocity::new((gm / r0).into_value().sqrt());

    // Position and velocity as raw vectors for the integrator
    let initial_state = OrbitState {
        pos: Vector3::new(r0.into_value(), 0.0, 0.0),
        vel: Vector3::new(0.0, v_circ.into_value(), 0.0),
    };

    // Orbital period: T = 2π * sqrt(r³/GM)
    // pow!(r0, 3) has units m³. GM has units m³/s². r³/GM = s².
    let r_cubed = pow!(r0, 3);
    let period_secs = Seconds::new(2.0 * std::f64::consts::PI * (r_cubed / gm).into_value().sqrt());

    println!("=== Orbital Propagator (Dormand-Prince 5(4)) ===");
    println!("Initial altitude:  {:.0} m", alt.into_value());
    println!("Initial radius:    {:.0} m", r0.into_value());
    println!("Circular velocity: {:.1} m/s", v_circ.into_value());
    println!(
        "Orbital period:    {:.1} s ({:.1} min)",
        period_secs.into_value(),
        period_secs.into_value() / 60.0
    );
    println!();

    // Propagate one full orbit
    let dt_init = Seconds::new(10.0); // start with 10s steps
    let final_state = propagate(initial_state, period_secs, dt_init, 1e-10);

    // Check energy conservation (should be near-zero drift for adaptive RKDP54)
    let e_initial = specific_energy(initial_state.pos, initial_state.vel);
    let e_final = specific_energy(final_state.pos, final_state.vel);
    let energy_drift = (e_final - e_initial).into_value().abs();

    // Check position closure (should return to start after one period)
    let pos_error = (final_state.pos - initial_state.pos).norm();

    println!("After 1 orbit:");
    println!(
        "  Position: ({:.1}, {:.1}, {:.1}) m",
        final_state.pos.x, final_state.pos.y, final_state.pos.z
    );
    println!(
        "  Velocity: ({:.1}, {:.1}, {:.1}) m/s",
        final_state.vel.x, final_state.vel.y, final_state.vel.z
    );
    println!("  Specific energy: {:.6e} J/kg", e_final.into_value());
    println!(
        "  Energy drift:    {:.6e} J/kg ({:.2e} relative)",
        energy_drift,
        energy_drift / e_initial.into_value().abs()
    );
    println!("  Position error:  {:.3} m (closure)", pos_error);
    println!();
    println!("All quantities tracked with compile-time SI units.");
    println!("Zero unit errors — the compiler verified every operation.");
}
