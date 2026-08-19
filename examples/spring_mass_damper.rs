//! Damped spring-mass-damper system using velocity-Verlet integration.
//!
//! m·ẍ + c·ẋ + k·x = 0
//!
//! Run: cargo run --example spring_mass_damper

use metron::constants::TAU;
use metron::pow;
use metron::{
    Acceleration, DampingCoefficient, Energy, Frequency, Kilograms, Meters, Seconds,
    SpringConstant, Velocity,
};

fn main() {
    let mass: Kilograms = Kilograms::new(250.0);
    let omega0: Frequency = Frequency::new(f64::from(TAU));
    let q_factor: f64 = 5.0;

    // k = m·ω₀²  (N/m)
    let spring_k: SpringConstant = mass * pow!(omega0, 2);
    // c = m·ω₀/Q  (N·s/m)
    let damping_c: DampingCoefficient = mass * omega0.map(|w| w / q_factor);

    let x0: Meters = Meters::new(0.10);
    let v0: Velocity = Velocity::new(0.0);

    // a = -(k·x + c·v)/m
    let a0: Acceleration = (-spring_k * x0 - damping_c * v0) / mass;

    let zeta: f64 = 1.0 / (2.0 * q_factor);
    let omega_d: Frequency = omega0.map(|w| w * (1.0 - zeta * zeta).sqrt());

    // E = ½kx² + ½mv²
    let e0: Energy =
        (spring_k * pow!(x0, 2)).map(|e| e * 0.5) + (mass * pow!(v0, 2)).map(|e| e * 0.5);

    println!("=== Spring-Mass-Damper System ===");
    println!("Velocity-Verlet integration (symplectic)");
    println!();
    println!("  Mass:           {:.0}", mass);
    println!(
        "  Natural freq:   {:.4} Hz",
        omega0 / Frequency::new(f64::from(TAU))
    );
    println!("  Quality factor: Q = {:.1}", q_factor);
    println!("  Spring constant k:  {:.2}", spring_k);
    println!("  Damping coeff  c:   {:.2}", damping_c);
    println!(
        "  Damping ratio ζ:    {:.4} ({})",
        zeta,
        if zeta < 1.0 {
            "underdamped"
        } else if (zeta - 1.0).abs() < 1e-10 {
            "critically damped"
        } else {
            "overdamped"
        }
    );
    println!(
        "  Damped freq ω_d:   {:.4} ({:.4} Hz)",
        omega_d,
        omega_d / Frequency::new(f64::from(TAU))
    );
    println!("  Initial energy: {:.4}", e0);
    println!();

    let dt: Seconds = Seconds::new(0.001);
    let t_total: Seconds = Seconds::new(5.0);
    let n_steps: usize = (t_total / dt).value as usize;

    let mut x: Meters = x0;
    let mut v: Velocity = v0;
    let mut a: Acceleration = a0;

    println!(
        "  {:>8}  {:>10}  {:>10}  {:>10}  {:>10}",
        "t (s)", "x (m)", "v (m/s)", "E (J)", "E/E₀"
    );
    let print_every: usize = n_steps / 10;

    for step in 0..=n_steps {
        if step % print_every == 0 {
            let t: Seconds = dt * (step as f64);
            let e: Energy = (spring_k * pow!(x, 2)).map(|val| val * 0.5)
                + (mass * pow!(v, 2)).map(|val| val * 0.5);
            let e_ratio: f64 = (e / e0).into();
            println!(
                "  {:>8.2}  {:>10.6}  {:>10.6}  {:>10.4}  {:>10.6}",
                t, x, v, e, e_ratio
            );
        }

        if step < n_steps {
            // Velocity-Verlet — fully typed
            let x_new: Meters = x + v * dt + (a * dt * dt).map(|val| val * 0.5);
            let v_half: Velocity = v + (a * dt).map(|val| val * 0.5);
            let a_new: Acceleration = (-spring_k * x_new - damping_c * v_half) / mass;
            let v_new: Velocity = v + ((a + a_new) * dt).map(|val| val * 0.5);

            x = x_new;
            v = v_new;
            a = a_new;
        }
    }
    println!();

    // Analytic: x(t) = A·e^(-ζω₀t)·cos(ω_d·t)
    let exp_arg: f64 = f64::from(-omega0 * zeta * t_total).exp();
    let cos_arg: f64 = f64::from(omega_d * t_total);
    let x_analytic: Meters = x0 * exp_arg * cos_arg.cos();
    let x_error: Meters = (x - x_analytic).abs();

    println!("Comparison at t = {:.1}:", t_total);
    println!("  Numerical  x = {:.6e}", x);
    println!("  Analytic   x = {:.6e}", x_analytic);
    println!(
        "  Error      = {:.6e} ({:.2e} relative)",
        x_error,
        <f64 as From<_>>::from(x_error / x0)
    );

    let e_final: Energy =
        (spring_k * pow!(x, 2)).map(|val| val * 0.5) + (mass * pow!(v, 2)).map(|val| val * 0.5);
    let e_dissipated: Energy = e0 - e_final;
    let pct: f64 = f64::from(e_dissipated / e0) * 100.0;
    println!(
        "Energy dissipated: {:.4} ({:.1}% of initial)",
        e_dissipated, pct
    );
}
