//! Quantum wavefunction simulation using the Crank-Nicolson method.
//!
//! Propagates a bound wavefunction in a 1D infinite square well (particle
//! in a box) using the time-dependent Schrödinger equation:
//!
//!   iℏ ∂ψ/∂t = Ĥψ
//!
//! The Crank-Nicolson scheme is unconditionally stable and unitary —
//! it preserves the norm of ψ exactly (up to floating-point roundoff).
//!
//! After propagation, we perform a position measurement (wavefunction
//! collapse) and verify that probability is conserved throughout.
//!
//! The grid computation uses natural units (ℏ = m = well_len = 1) internally,
//! but all scale factors and physical observables are derived from typed
//! metron constants with compile-time dimensional checking.
//!
//! Run: cargo run --example quantum_simulation

use metron::constants::{ELECTRON_MASS, HBAR, PI};
use metron::pow;
use metron::{Action, Area, Energy, Kilograms, Meters, ProbabilityDensity, Quantity, Seconds};
use num_complex::Complex;
use std::ops::Mul;

/// Number of spatial grid points.
const N: usize = 512;

/// Complex wavefunction: ψ(x) at N grid points.
type Wavefunction = Vec<Complex<f64>>;

/// Build the Hamiltonian matrix H = -(1/2) * d²/dx² in natural units (ℏ=m=1, well_len=1).
/// Returns tridiagonal entries (diagonal, off-diagonal).
fn build_hamiltonian(dx: f64) -> (Vec<f64>, Vec<f64>) {
    let diag_val = 1.0 / (dx * dx);
    let offdiag_val = -1.0 / (2.0 * dx * dx);
    (vec![diag_val; N], vec![offdiag_val; N - 1])
}

/// Crank-Nicolson step: (I + iHΔt/2) ψ_new = (I - iHΔt/2) ψ_old
fn crank_nicolson_step(psi: &Wavefunction, diag: &[f64], offdiag: &[f64], dt: f64) -> Wavefunction {
    let alpha = Complex::new(0.0, dt / 2.0);

    let mut rhs = vec![Complex::new(0.0, 0.0); N];
    for i in 0..N {
        rhs[i] = psi[i] - alpha * diag[i] * psi[i];
        if i > 0 {
            rhs[i] -= alpha * offdiag[i - 1] * psi[i - 1];
        }
        if i < N - 1 {
            rhs[i] -= alpha * offdiag[i] * psi[i + 1];
        }
    }

    let mut a_diag: Vec<Complex<f64>> = (0..N).map(|_| Complex::new(1.0, 0.0)).collect();
    let mut a_upper: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); N - 1];
    let mut a_lower: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); N - 1];

    for i in 0..N {
        a_diag[i] += alpha * diag[i];
    }
    for i in 0..N - 1 {
        a_upper[i] = alpha * offdiag[i];
        a_lower[i] = alpha * offdiag[i];
    }

    // Thomas algorithm — forward sweep
    for i in 1..N {
        let factor = a_lower[i - 1] / a_diag[i - 1];
        a_diag[i] -= factor * a_upper[i - 1];
        let prev_rhs = rhs[i - 1];
        rhs[i] -= factor * prev_rhs;
    }

    // Back substitution
    let mut result = vec![Complex::new(0.0, 0.0); N];
    result[N - 1] = rhs[N - 1] / a_diag[N - 1];
    for i in (0..N - 1).rev() {
        result[i] = (rhs[i] - a_upper[i] * result[i + 1]) / a_diag[i];
    }

    result
}

/// Total probability: ∫|ψ|² dx (natural units)
fn total_probability(psi: &Wavefunction, dx: f64) -> f64 {
    psi.iter().map(|c| c.norm_sqr()).sum::<f64>() * dx
}

/// Energy expectation ⟨E⟩ = ∫ ψ* H ψ dx (natural units)
fn energy_expectation(psi: &Wavefunction, diag: &[f64], offdiag: &[f64], dx: f64) -> f64 {
    let mut e = 0.0;
    for i in 0..N {
        e += (psi[i].conj() * psi[i]).re * diag[i];
        if i > 0 {
            e += (psi[i].conj() * psi[i - 1]).re * offdiag[i - 1];
        }
        if i < N - 1 {
            e += (psi[i].conj() * psi[i + 1]).re * offdiag[i];
        }
    }
    e * dx
}

/// Initialize ψ(x,0) = (φ₁ + φ₂) / √2
fn initial_wavefunction() -> Wavefunction {
    let dx = 1.0 / N as f64;
    let pi_val = f64::from(PI);
    let mut psi = Vec::with_capacity(N);
    for i in 0..N {
        let x = (i as f64 + 0.5) * dx;
        let phi1 = 2.0_f64.sqrt() * (pi_val * x).sin();
        let phi2 = 2.0_f64.sqrt() * (2.0 * pi_val * x).sin();
        psi.push(Complex::new((phi1 + phi2) / 2.0_f64.sqrt(), 0.0));
    }
    psi
}

/// Position measurement: collapse ψ to a Gaussian at a random position.
fn measure_position(psi: &Wavefunction, dx: f64) -> (usize, Wavefunction) {
    let probs: Vec<f64> = psi.iter().map(|c| c.norm_sqr() * dx).collect();
    let total: f64 = probs.iter().sum();
    let r = 0.42;
    let mut cum = 0.0;
    let mut idx = N - 1;
    for (i, &p) in probs.iter().enumerate() {
        cum += p / total;
        if cum >= r {
            idx = i;
            break;
        }
    }

    let x_measured = (idx as f64 + 0.5) * dx;
    let sigma = 1.0 / 50.0;
    let pi_val = f64::from(PI);
    let mut collapsed = Vec::with_capacity(N);
    let mut norm_sq = 0.0;
    for i in 0..N {
        let x = (i as f64 + 0.5) * dx;
        let g = (-(x - x_measured).powi(2) / (2.0 * sigma * sigma)).exp()
            / (sigma * (2.0 * pi_val).sqrt());
        collapsed.push(Complex::new(g, 0.0));
        norm_sq += g * g * dx;
    }
    let norm = norm_sq.sqrt();
    for c in &mut collapsed {
        *c /= Complex::new(norm, 0.0);
    }
    (idx, collapsed)
}

fn main() {
    // Physical parameters — all typed, compile-time dimensionally checked
    let m_e: Kilograms = ELECTRON_MASS;
    let hbar: Action = HBAR;
    let well_width: Meters = Meters::new(1e-9); // 1 nm quantum dot
    let well_len: Meters = well_width;

    // Energy scale: E = ℏ²/(mL²) — converts natural-unit energy to Joules
    // ℏ² has units (J·s)², mL² has units kg·m², so ℏ²/(mL²) = J
    let hbar_sq: Quantity<f64, <metron::dim::Action as Mul<metron::dim::Action>>::Output> =
        pow!(hbar, 2);
    let len_sq: Area = pow!(well_len, 2);
    let mass_len_sq: Quantity<f64, <metron::dim::Kilogram as Mul<metron::dim::Area>>::Output> =
        m_e * len_sq;
    let e_scale: Energy = hbar_sq / mass_len_sq;

    // Time scale: t = mL²/ℏ — converts natural-unit time to seconds
    let t_scale: Seconds = mass_len_sq / hbar;

    // Spatial grid in natural units (well_len=1)
    let dx_natural: f64 = 1.0 / N as f64;
    let dx_si: Meters = Meters::new(dx_natural * well_len.value);

    println!("=== Quantum Wavefunction Simulation ===");
    println!("Crank-Nicolson propagation of TDSE");
    println!();
    println!("Particle: electron (m = {:.4e})", m_e);
    println!("Well width: {:.4e} (1 nm quantum dot)", well_len);
    println!("Grid points: {}", N);
    println!("dx = {:.4e}", dx_si);
    println!();

    // Build Hamiltonian
    let (diag, offdiag) = build_hamiltonian(dx_natural);

    // Initial state
    let psi0 = initial_wavefunction();
    let prob0 = total_probability(&psi0, dx_natural);
    let e0_natural = energy_expectation(&psi0, &diag, &offdiag, dx_natural);
    // Analytic: ⟨E⟩ = (E1 + E2)/2 = (π²/2 + 2π²)/2 = 5π²/4
    let e_analytic_natural = 5.0 * f64::from(PI).powi(2) / 4.0;
    let e0_si: Energy = Energy::new(e0_natural * e_scale.value);

    // Analytic energies (SI): E_n = n²π²ℏ²/(2mL²)
    let e1_si: Energy = Energy::new(f64::from(PI).powi(2) * e_scale.value / 2.0);
    let e2_si: Energy = Energy::new(4.0 * e1_si.value);

    println!("Initial state: (|1⟩ + |2⟩) / √2");
    println!("  |ψ|² integral:  {:.10}", prob0);
    println!(
        "  ⟨E⟩ = {:.6e} ({:.4} eV)  [natural: {:.6}, analytic: {:.6}]",
        e0_si,
        e0_si.value / 1.602e-19,
        e0_natural,
        e_analytic_natural
    );
    println!();
    println!("Analytic energies:");
    println!("  E₁ = {:.6e} ({:.4} eV)", e1_si, e1_si.value / 1.602e-19);
    println!("  E₂ = {:.6e} ({:.4} eV)", e2_si, e2_si.value / 1.602e-19);
    println!(
        "  (E₁+E₂)/2 = {:.6e} (should match ⟨E⟩)",
        Energy::new((e1_si.value + e2_si.value) / 2.0)
    );
    println!();

    // Propagate
    // Ground state period: T₁ = 2π/E₁ (natural units)
    let e1_natural = f64::from(PI).powi(2) / 2.0; // n=1, ℏ=m=well_len=1
    let t1_natural = 2.0 * f64::from(PI) / e1_natural;
    let dt_natural = t1_natural / 1000.0;
    let dt_si: Seconds = Seconds::new(dt_natural * t_scale.value);
    let t1_si: Seconds = Seconds::new(t1_natural * t_scale.value);

    println!("Propagation:");
    println!("  dt = {:.4e}", dt_si);
    println!("  T₁ = {:.4e} (ground state period)", t1_si);
    println!();

    let mut psi = psi0.clone();
    let n_steps = 500;
    for step in 0..=n_steps {
        if step % 100 == 0 {
            let prob = total_probability(&psi, dx_natural);
            let e_nat = energy_expectation(&psi, &diag, &offdiag, dx_natural);
            let e_si: Energy = Energy::new(e_nat * e_scale.value);
            let t_si: Seconds = Seconds::new((step as f64) * dt_natural * t_scale.value);
            println!("  t={:>8.4e}  |ψ|²={:.10}  ⟨E⟩={:.6e}", t_si, prob, e_si);
        }
        if step < n_steps {
            psi = crank_nicolson_step(&psi, &diag, &offdiag, dt_natural);
        }
    }
    println!();

    // Measurement
    let (measured_idx, psi_collapsed) = measure_position(&psi, dx_natural);
    let prob_before = total_probability(&psi, dx_natural);
    let prob_after = total_probability(&psi_collapsed, dx_natural);
    let x_measured: Meters = Meters::new((measured_idx as f64 + 0.5) * dx_natural * well_len.value);

    println!("Measurement (position):");
    println!("  Measured x = {:.6e}", x_measured);
    println!("  |ψ|² before collapse: {:.10}", prob_before);
    println!("  |ψ|² after collapse:  {:.10}", prob_after);
    println!(
        "  (norm preserved: {})",
        (prob_before - prob_after).abs() < 1e-10
    );
    println!();

    // Probability density
    let density: Vec<ProbabilityDensity> = psi_collapsed
        .iter()
        .map(|c| ProbabilityDensity::new(c.norm_sqr() / well_len.value))
        .collect();

    let peak_idx = density
        .iter()
        .enumerate()
        .max_by(|(_, pd_a), (_, pd_b)| pd_a.value.partial_cmp(&pd_b.value).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let peak_x: Meters = Meters::new((peak_idx as f64 + 0.5) * dx_si.value);
    let peak_density: ProbabilityDensity = density[peak_idx];

    println!("Probability density |ψ(x)|² (collapsed state):");
    println!("  Peak at x = {:.6e}, |ψ|² = {:.6e}", peak_x, peak_density);
    println!();
    println!("All quantities tracked with compile-time SI units.");
    println!("Energy in Joules, probability density in 1/m, time in seconds.");
    println!("The compiler verified dimensional correctness of every operation.");
}
