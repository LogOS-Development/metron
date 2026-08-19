//! Damped spring-mass-damper system using velocity-Verlet integration.
//!
//! Simulates a mass on a spring with viscous damping:
//!
//!   m·ẍ + c·ẋ + k·x = 0
//!
//! Showcases:
//! - Spring constant k (N/m = kg/s²) — a derived unit easy to get wrong
//! - Damping coefficient c (N·s/m = kg/s)
//! - Natural frequency ω₀ (rad/s), quality factor Q
//! - Energy: ½kx² (spring PE), ½mv² (kinetic) — using pow!
//! - Velocity-Verlet integrator (symplectic, preserves energy structure)
//!
//! Run: cargo run --example spring_mass_damper

use metron::pow;
use metron::{Frequency, Kilograms, Meters, Seconds, Velocity};

fn main() {
    // System parameters: car suspension-like
    // m = 250 kg, ω₀ = 2π rad/s (1 Hz), Q = 5 (underdamped)
    let mass = Kilograms::new(250.0);
    let omega0 = Frequency::new(2.0 * std::f64::consts::PI);
    let q_factor: f64 = 5.0;

    // Spring constant: k = m·ω₀²  →  kg·(1/s)² = kg/s² = N/m
    // The type is inferred: mass has kg, ω² has 1/s², product is kg/s².
    let spring_k = mass * pow!(omega0, 2);

    // Damping: c = m·ω₀/Q  →  kg·(1/s) = kg/s = N·s/m
    let omega_over_q = omega0.map(|w| w / q_factor);
    let damping_c = mass * omega_over_q;

    // Initial conditions: displaced 10 cm, at rest
    let x0 = Meters::new(0.10);
    let v0 = Velocity::new(0.0);

    // Spring force: F = -k·x  →  (kg/s²)·m = kg·m/s² = N
    let f_spring0 = -spring_k * x0;
    // Damping force: F = -c·v  →  (kg/s)·(m/s) = kg·m/s² = N
    let f_damp0 = -damping_c * v0;
    // Total force / mass = acceleration
    let _a0 = (f_spring0 + f_damp0) / mass;

    println!("=== Spring-Mass-Damper System ===");
    println!("Velocity-Verlet integration (symplectic)");
    println!();
    println!("Parameters:");
    println!("  Mass:           {:.0} kg", mass.into_value());
    println!(
        "  Natural freq:   {:.4} Hz",
        omega0.into_value() / (2.0 * std::f64::consts::PI)
    );
    println!("  Quality factor: Q = {:.1}", q_factor);
    println!();
    println!("Derived (all compile-time unit checked):");
    println!("  Spring constant k:  {:.2} N/m", spring_k.into_value());
    println!("  Damping coeff  c:   {:.2} N·s/m", damping_c.into_value());
    println!();

    // Damping ratio: ζ = 1/(2Q)
    let zeta = 1.0 / (2.0 * q_factor);
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

    // Damped frequency: ω_d = ω₀·sqrt(1 - ζ²)
    let omega_d = omega0.into_value() * (1.0 - zeta * zeta).sqrt();
    println!(
        "  Damped freq ω_d:   {:.4} rad/s ({:.4} Hz)",
        omega_d,
        omega_d / (2.0 * std::f64::consts::PI)
    );
    println!();

    // Energy at t=0: E = ½kx² + ½mv²
    // k·x² = (kg/s²)·m² = J, m·v² = kg·(m/s)² = J
    let e0 = {
        let spring_pe = spring_k * pow!(x0, 2);
        let kinetic = mass * pow!(v0, 2);
        spring_pe.map(|e| e * 0.5) + kinetic.map(|e| e * 0.5)
    };
    println!("Initial energy: {:.4} J", e0.into_value());
    println!();

    // Simulate
    let dt = Seconds::new(0.001); // 1 ms
    let t_total = 5.0_f64;
    let n_steps = (t_total / dt.into_value()) as usize;

    let mut x = x0;
    let mut v = v0;

    println!("Propagation (every 0.5 s):");
    println!(
        "  {:>8}  {:>10}  {:>10}  {:>10}  {:>10}",
        "t (s)", "x (m)", "v (m/s)", "E (J)", "E/E₀"
    );

    let print_every = n_steps / 10;

    for step in 0..=n_steps {
        if step % print_every == 0 {
            let t = (step as f64) * dt.into_value();
            // Energy: ½kx² + ½mv² — all typed
            let e = {
                let pe = spring_k * pow!(x, 2);
                let ke = mass * pow!(v, 2);
                pe.map(|val| val * 0.5) + ke.map(|val| val * 0.5)
            };
            let e_ratio = e.into_value() / e0.into_value();
            println!(
                "  {:>8.2}  {:>10.6}  {:>10.6}  {:>10.4}  {:>10.6}",
                t,
                x.into_value(),
                v.into_value(),
                e.into_value(),
                e_ratio
            );
        }

        if step < n_steps {
            // Forces from current state — fully typed
            let f_s = -spring_k * x;
            let f_d = -damping_c * v;
            let a_curr = (f_s + f_d) / mass;

            // Velocity-Verlet: x_new = x + v·dt + ½a·dt²
            let dt_val = dt.into_value();
            let x_new = Meters::new(
                x.into_value()
                    + v.into_value() * dt_val
                    + 0.5 * a_curr.into_value() * dt_val * dt_val,
            );
            // Half-step velocity for damping evaluation
            let v_half = Velocity::new(v.into_value() + 0.5 * a_curr.into_value() * dt_val);
            // New acceleration at new position
            let f_s_new = -spring_k * x_new;
            let f_d_new = -damping_c * v_half;
            let a_new = (f_s_new + f_d_new) / mass;
            // v_new = v + ½(a + a_new)·dt
            let v_new = Velocity::new(
                v.into_value() + 0.5 * (a_curr.into_value() + a_new.into_value()) * dt_val,
            );

            x = x_new;
            v = v_new;
        }
    }
    println!();

    // Analytic solution: x(t) = A·e^(-ζω₀t)·cos(ω_d·t)
    let decay = (-zeta * omega0.into_value() * t_total).exp();
    let x_analytic = x0.into_value() * decay * (omega_d * t_total).cos();
    let x_error = (x.into_value() - x_analytic).abs();

    println!("Comparison at t = {:.1} s:", t_total);
    println!("  Numerical  x = {:.6e} m", x.into_value());
    println!("  Analytic   x = {:.6e} m", x_analytic);
    println!(
        "  Error      = {:.6e} m ({:.2e} relative)",
        x_error,
        x_error / x0.into_value()
    );
    println!();

    // Energy dissipated
    let e_final = {
        let pe = spring_k * pow!(x, 2);
        let ke = mass * pow!(v, 2);
        pe.map(|val| val * 0.5) + ke.map(|val| val * 0.5)
    };
    let e_dissipated = e0.into_value() - e_final.into_value();
    println!(
        "Energy dissipated: {:.4} J ({:.1}% of initial)",
        e_dissipated,
        100.0 * e_dissipated / e0.into_value()
    );
    println!();
    println!("All quantities tracked with compile-time SI units.");
    println!("Spring constant (N/m), damping (N·s/m), energy (J) —");
    println!("the compiler verified every unit combination.");
}
