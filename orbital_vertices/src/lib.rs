use pyo3::prelude::*;
use rand::Rng;
use num_complex::Complex;

fn factorial(n: usize) -> f64 {
    if n == 0 || n == 1 {
        1.0
    } else {
        (1..=n).map(|x| x as f64).product()
    }
}

fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        0.0
    } else {
        factorial(n) / (factorial(k) * factorial(n - k))
    }
}

fn assoc_laguerre(n: usize, k: usize, x: f64) -> f64 {
    (0..=n).map(|m| (-1.0_f64).powi(m as i32) * binom(n + k, n - m) * x.powi(m as i32) / factorial(m)).sum()
}

fn assoc_legendre(l: usize, m: usize, x: f64) -> f64 {
    // Simplified implementation for small l, m
    // For full implementation, use proper recurrence
    match (l, m) {
        (0, 0) => 1.0,
        (1, 0) => x,
        (1, 1) => -(1.0 - x * x).sqrt(),
        (2, 0) => (3.0 * x * x - 1.0) / 2.0,
        (2, 1) => -3.0 * x * (1.0 - x * x).sqrt(),
        (2, 2) => 3.0 * (1.0 - x * x),
        (3, 0) => (5.0 * x * x * x - 3.0 * x) / 2.0,
        (3, 1) => (3.0 / 2.0) * (1.0 - x * x).sqrt() * (5.0 * x * x - 1.0),
        (3, 2) => 15.0 * x * (1.0 - x * x),
        (3, 3) => -15.0 * (1.0 - x * x).powf(1.5),
        _ => 0.0, // Not implemented
    }
}

fn spherical_harmonic(l: usize, m: i32, theta: f64, phi: f64) -> Complex<f64> {
    let m_abs = m.abs() as usize;
    let prefactor = (-1.0_f64).powi(m) * ((2 * l + 1) as f64 / (4.0 * std::f64::consts::PI) * factorial(l - m_abs) / factorial(l + m_abs)).sqrt();
    let legendre = assoc_legendre(l, m_abs, theta.cos());
    let exp_part = Complex::new(0.0, m as f64 * phi).exp();
    prefactor * legendre * exp_part
}

#[pyfunction]
fn generate_orbital_vertices(n: usize, l: usize, m: i32, num_points: usize) -> PyResult<(Vec<(f64, f64, f64)>, Vec<f64>)> {
    let mut rng = rand::thread_rng();
    let mut vertices = Vec::with_capacity(num_points);
    let mut phases = Vec::with_capacity(num_points);
    let radius_scale = 0.5;
    let mut attempts = 0;
    let max_attempts = num_points * 512;

    while vertices.len() < num_points && attempts < max_attempts {
        attempts += 1;
        // Exponential distribution
        let u: f64 = rng.gen();
        let r = - (n as f64 * radius_scale) * u.ln();
        if r > n as f64 * 3.0 * radius_scale {
            continue;
        }
        let theta = (2.0 * rng.gen::<f64>() - 1.0).acos();
        let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();
        let rho = 2.0 * r / (n as f64 * radius_scale);
        if rho < 0.01 {
            continue;
        }
        // Radial wave function
        let laguerre_val = assoc_laguerre(n - l - 1, 2 * l + 1, rho);
        let r_wave = (-rho / 2.0).exp() * rho.powi(l as i32) * laguerre_val;
        // Angular wave function
        let y = spherical_harmonic(l, m, theta, phi);
        let prob = r_wave * r_wave * y.norm_sqr() * r * r;
        let threshold = 0.4 * (n as f64).powi(2);
        if prob > threshold && rng.gen::<f64>() < (prob * 2.0).min(1.0) {
            let x = r * theta.sin() * phi.cos();
            let y_pos = r * theta.sin() * phi.sin();
            let z = r * theta.cos();
            vertices.push((x, y_pos, z));
            let psi = r_wave * y;
            let sign = if psi.re >= 0.0 { 1.0 } else { -1.0 };
            let unorm = 0.5 * sign + 0.5;
            phases.push(unorm);
        }
    }
    // Fill with simple sphere if not enough points
    while vertices.len() < num_points {
        let u: f64 = rng.gen();
        let r = - (n as f64 * radius_scale) * u.ln();
        let theta = (2.0 * rng.gen::<f64>() - 1.0).acos();
        let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();
        let x = r * theta.sin() * phi.cos();
        let y_pos = r * theta.sin() * phi.sin();
        let z = r * theta.cos();
        vertices.push((x, y_pos, z));
        phases.push(0.5);
    }
    vertices.truncate(num_points);
    phases.truncate(num_points);
    Ok((vertices, phases))
}

#[pymodule]
fn orbital_vertices(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(generate_orbital_vertices, m)?)?;
    Ok(())
}