//! Orbital propagator using Dormand-Prince 5(4) adaptive integration.
//!
//! Propagates a satellite in a two-body Keplerian orbit around Earth,
//! using metron quantities throughout. The integrator is a proper
//! adaptive-step RKDP54 — the same family used by NASA GMAT and AGI STK.
//!
//! Run: cargo run --example orbital_propagator

use metron::constants::{GM_EARTH, R_EARTH_EQ, TAU};
use metron::pow;
use metron::{
    Acceleration, AccelerationVector, Area, DimensionlessVector, GravitationalParameter, Meters,
    PositionVector, Seconds, SpecificEnergy, Velocity, VelocityVector, Volume,
};

/// Physical state: position and velocity as typed vectors.
#[derive(Clone, Copy, Debug)]
struct OrbitState {
    pos: PositionVector,
    vel: VelocityVector,
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
/// Fully typed: gm / r² → (m³/s²) / m² = m/s² = Acceleration.
fn gravity_accel_magnitude(r: Meters, gm: GravitationalParameter) -> Acceleration {
    let r_sq: Area = pow!(r, 2);
    gm / r_sq
}

/// Gravitational acceleration vector.
///
/// Takes typed position, returns typed acceleration.
/// a = -GM/|r|² * r̂
fn gravity_accel(pos: PositionVector, gm: GravitationalParameter) -> AccelerationVector {
    let r_mag: Meters = pos.norm();
    let accel: Acceleration = gravity_accel_magnitude(r_mag, gm);
    let r_hat: DimensionlessVector = pos / r_mag;
    r_hat * (-accel)
}

/// Derivative of orbital state: velocity and acceleration.
fn derivative(state: &OrbitState) -> (VelocityVector, AccelerationVector) {
    let gm: GravitationalParameter = GM_EARTH;
    let acc: AccelerationVector = gravity_accel(state.pos, gm);
    (state.vel, acc)
}

/// One RKDP54 step with adaptive step size control.
///
/// Returns the new state and the error estimate (dimensionless relative).
fn rkdp54_step(state: &OrbitState, dt: Seconds) -> (OrbitState, f64) {
    let (k1_v, k1_a): (VelocityVector, AccelerationVector) = derivative(state);

    // k_v * dt → PositionVector, k_a * dt → VelocityVector
    let s2 = OrbitState {
        pos: state.pos + k1_v * (rkdp54::C[1] * dt),
        vel: state.vel + k1_a * (rkdp54::C[1] * dt),
    };
    let (k2_v, k2_a) = derivative(&s2);

    let s3 = OrbitState {
        pos: state.pos + k1_v * (rkdp54::A[2][0] * dt) + k2_v * (rkdp54::A[2][1] * dt),
        vel: state.vel + k1_a * (rkdp54::A[2][0] * dt) + k2_a * (rkdp54::A[2][1] * dt),
    };
    let (k3_v, k3_a) = derivative(&s3);

    let s4 = OrbitState {
        pos: state.pos
            + k1_v * (rkdp54::A[3][0] * dt)
            + k2_v * (rkdp54::A[3][1] * dt)
            + k3_v * (rkdp54::A[3][2] * dt),
        vel: state.vel
            + k1_a * (rkdp54::A[3][0] * dt)
            + k2_a * (rkdp54::A[3][1] * dt)
            + k3_a * (rkdp54::A[3][2] * dt),
    };
    let (k4_v, k4_a) = derivative(&s4);

    let s5 = OrbitState {
        pos: state.pos
            + k1_v * (rkdp54::A[4][0] * dt)
            + k2_v * (rkdp54::A[4][1] * dt)
            + k3_v * (rkdp54::A[4][2] * dt)
            + k4_v * (rkdp54::A[4][3] * dt),
        vel: state.vel
            + k1_a * (rkdp54::A[4][0] * dt)
            + k2_a * (rkdp54::A[4][1] * dt)
            + k3_a * (rkdp54::A[4][2] * dt)
            + k4_a * (rkdp54::A[4][3] * dt),
    };
    let (k5_v, k5_a) = derivative(&s5);

    let s6 = OrbitState {
        pos: state.pos
            + k1_v * (rkdp54::A[5][0] * dt)
            + k2_v * (rkdp54::A[5][1] * dt)
            + k3_v * (rkdp54::A[5][2] * dt)
            + k4_v * (rkdp54::A[5][3] * dt)
            + k5_v * (rkdp54::A[5][4] * dt),
        vel: state.vel
            + k1_a * (rkdp54::A[5][0] * dt)
            + k2_a * (rkdp54::A[5][1] * dt)
            + k3_a * (rkdp54::A[5][2] * dt)
            + k4_a * (rkdp54::A[5][3] * dt)
            + k5_a * (rkdp54::A[5][4] * dt),
    };
    let (k6_v, k6_a) = derivative(&s6);

    // 5th-order solution
    let new_pos: PositionVector = state.pos
        + k1_v * (rkdp54::B5[0] * dt)
        + k3_v * (rkdp54::B5[2] * dt)
        + k4_v * (rkdp54::B5[3] * dt)
        + k5_v * (rkdp54::B5[4] * dt)
        + k6_v * (rkdp54::B5[5] * dt);

    let new_vel: VelocityVector = state.vel
        + k1_a * (rkdp54::B5[0] * dt)
        + k3_a * (rkdp54::B5[2] * dt)
        + k4_a * (rkdp54::B5[3] * dt)
        + k5_a * (rkdp54::B5[4] * dt)
        + k6_a * (rkdp54::B5[5] * dt);

    // 7th stage (FSAL: reuses k1 of next step, uses new state)
    let s7 = OrbitState {
        pos: new_pos,
        vel: new_vel,
    };
    let (k7_v, k7_a) = derivative(&s7);

    // Error estimate = |5th-order - 4th-order|
    // (b5-b4) is dimensionless, * dt gives position/velocity units
    let err_pos: PositionVector = k1_v * ((rkdp54::B5[0] - rkdp54::B4[0]) * dt)
        + k3_v * ((rkdp54::B5[2] - rkdp54::B4[2]) * dt)
        + k4_v * ((rkdp54::B5[3] - rkdp54::B4[3]) * dt)
        + k5_v * ((rkdp54::B5[4] - rkdp54::B4[4]) * dt)
        + k6_v * ((rkdp54::B5[5] - rkdp54::B4[5]) * dt)
        + k7_v * ((rkdp54::B5[6] - rkdp54::B4[6]) * dt);

    let err_vel: VelocityVector = k1_a * ((rkdp54::B5[0] - rkdp54::B4[0]) * dt)
        + k3_a * ((rkdp54::B5[2] - rkdp54::B4[2]) * dt)
        + k4_a * ((rkdp54::B5[3] - rkdp54::B4[3]) * dt)
        + k5_a * ((rkdp54::B5[4] - rkdp54::B4[4]) * dt)
        + k6_a * ((rkdp54::B5[5] - rkdp54::B4[5]) * dt)
        + k7_a * ((rkdp54::B5[6] - rkdp54::B4[6]) * dt);

    // RMS error norm — combines position error (m) and velocity error (m/s).
    // These have different dimensions, so extract raw values for the scalar
    // tolerance comparison (standard RKDP54 adaptive control).
    let error = {
        let ep: f64 = err_pos.norm().value;
        let ev: f64 = err_vel.norm().value;
        (ep * ep + ev * ev).sqrt()
    };

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
    let mut t: Seconds = Seconds::new(0.0);

    while t < t_total {
        if t + dt > t_total {
            dt = t_total - t;
        }

        let (new_state, error) = rkdp54_step(&state, dt);

        if error < tolerance || dt < Seconds::new(1e-6) {
            // Accept step
            state = new_state;
            t += dt;
            // Increase step if error is small
            if error > 0.0 {
                let factor = (tolerance / error).powf(0.2).clamp(0.2, 5.0);
                dt *= factor;
            }
        } else {
            // Reject step, shrink dt
            let factor = (tolerance / error).powf(0.2).max(0.1);
            dt *= factor;
        }
    }
    state
}

/// Compute specific orbital energy: ε = v²/2 - GM/r
///
/// All arithmetic is compile-time unit checked. v² has units (m/s)² = m²/s²,
/// GM/r has units m³/s² / m = m²/s² — both are specific energy (J/kg).
fn specific_energy(pos: PositionVector, vel: VelocityVector) -> SpecificEnergy {
    let gm: GravitationalParameter = GM_EARTH;
    let r: Meters = pos.norm();
    let v: Velocity = vel.norm();
    let v_sq: SpecificEnergy = pow!(v, 2);
    let ke: SpecificEnergy = v_sq.map(|x: f64| x * 0.5);
    let pe: SpecificEnergy = gm / r;
    ke - pe
}

fn main() {
    // ISS-like orbit: 400 km altitude, circular, prograde
    let alt: Meters = Meters::new(400_000.0);
    let earth_r: Meters = R_EARTH_EQ;
    let r0: Meters = earth_r + alt;
    let gm: GravitationalParameter = GM_EARTH;

    // v = sqrt(GM/r)  →  sqrt(m²/s²) = m/s
    let v_circ: Velocity = (gm / r0).sqrt();

    // Position and velocity as typed vectors
    let initial_state = OrbitState {
        pos: PositionVector::from_xyz(r0.value, 0.0, 0.0),
        vel: VelocityVector::from_xyz(0.0, v_circ.value, 0.0),
    };

    // Orbital period: T = 2π * sqrt(r³/GM)
    let r_cubed: Volume = pow!(r0, 3);
    let period_secs: Seconds = TAU * (r_cubed / gm).sqrt();

    println!("=== Orbital Propagator (Dormand-Prince 5(4)) ===");
    println!("Initial altitude:  {:.0}", alt);
    println!("Initial radius:    {:.0}", r0);
    println!("Circular velocity: {:.1}", v_circ);
    println!(
        "Orbital period:    {:.1} ({:.1} min)",
        period_secs,
        period_secs / Seconds::new(60.0)
    );
    println!();

    // Propagate one full orbit
    let dt_init: Seconds = Seconds::new(10.0); // start with 10s steps
    let final_state: OrbitState = propagate(initial_state, period_secs, dt_init, 1e-10);

    // Check energy conservation
    let e_initial: SpecificEnergy = specific_energy(initial_state.pos, initial_state.vel);
    let e_final: SpecificEnergy = specific_energy(final_state.pos, final_state.vel);
    let energy_drift: SpecificEnergy = (e_final - e_initial).abs();

    // Check position closure (should return to start after one period)
    let pos_delta: PositionVector = final_state.pos - initial_state.pos;
    let pos_error: Meters = pos_delta.norm();

    println!("After 1 orbit:");
    println!(
        "  Position: ({:.1}, {:.1}, {:.1}) m",
        final_state.pos.vector.x, final_state.pos.vector.y, final_state.pos.vector.z
    );
    println!(
        "  Velocity: ({:.1}, {:.1}, {:.1}) m/s",
        final_state.vel.vector.x, final_state.vel.vector.y, final_state.vel.vector.z
    );
    println!("  Specific energy: {:.6e}", e_final);
    println!(
        "  Energy drift:    {:.6e} ({:.2e} relative)",
        energy_drift,
        energy_drift / e_initial.abs()
    );
    println!("  Position error:  {:.3} (closure)", pos_error);
    println!();
    println!("All quantities tracked with compile-time SI units.");
    println!("Zero unit errors — the compiler verified every operation.");
}
