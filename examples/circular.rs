//! Analytic circular-orbit check with degree-zero gravity and no drag.
use numerical_propagator::{engine_free, engine_new, engine_propagate, Inputs, Stats};

fn run_case(direction: f64, drag: bool) -> (f64, f64) {
    let gm: f64 = 3.986004418e14;
    let radius: f64 = 7_000_000.;
    let period = std::f64::consts::TAU * (radius.powi(3) / gm).sqrt();
    let identity = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
    let frames: Vec<f64> = identity.into_iter().cycle().take(18).collect();
    let zeros = [0.; 7];
    let c = [1.];
    let flux = [150.];
    let start = 1_700_000_000.;
    let input = Inputs {
        degree: 0,
        gm,
        radius: 6_378_137.,
        coefficients_c: c.as_ptr(),
        coefficients_s: zeros.as_ptr(),
        frame_start: start - 2. * period,
        frame_step: 4. * period,
        nframe: 2,
        q: frames.as_ptr(),
        polar: frames.as_ptr(),
        dut: zeros.as_ptr(),
        body_start: start - 2. * period,
        body_step: 4. * period,
        nbody: 2,
        sun_pos: zeros.as_ptr(),
        sun_vel: zeros.as_ptr(),
        moon_pos: zeros.as_ptr(),
        moon_vel: zeros.as_ptr(),
        weather_start: start - 2. * period,
        weather_step: 4. * period,
        nweather: 1,
        f107: flux.as_ptr(),
        f107a: flux.as_ptr(),
        ap: zeros.as_ptr(),
        doy: c.as_ptr(),
        correction_start: 0.,
        correction_step: 1.,
        ncorrection: 0,
        log_correction: zeros.as_ptr(),
        correction_height_linear: zeros.as_ptr(),
        correction_height_quadratic: zeros.as_ptr(),
        correction_height_origin: 0.,
        correction_height_scale: 1.,
        correction_clip: 0.,
        beta: 0.02,
        drag: drag as i32,
        third_bodies: 0,
    };
    let initial = [radius, 0., 0., 0., (gm / radius).sqrt(), 0.];
    let end = start + direction * period;
    let samples = [start, end];
    let mut out_t = [0.; 3];
    let mut out_y = [0.; 18];
    let mut stats = Stats::default();
    // The context copies the arrays. Output buffers have nsample+1 rows.
    let code = unsafe {
        let ctx = engine_new(&input);
        assert!(!ctx.is_null());
        if drag {
            #[cfg(feature = "msis21")]
            {
                (*ctx)
                    .load_msis21(
                        std::env::var("RMSIS21_LIBRARY").unwrap(),
                        std::env::var("RMSIS21_PARAMETERS").unwrap(),
                    )
                    .unwrap();
                let mut env = [0.; 8];
                assert_eq!(
                    numerical_propagator::engine_environment(
                        ctx,
                        start,
                        initial.as_ptr(),
                        env.as_mut_ptr()
                    ),
                    0
                );
                assert!(env[3] > 0. && env[3] < 1e-9);
            }
        }
        let result = engine_propagate(
            ctx,
            start,
            end,
            initial.as_ptr(),
            6,
            1e-12,
            30.,
            0.,
            samples.as_ptr(),
            samples.len(),
            out_t.as_mut_ptr(),
            out_y.as_mut_ptr(),
            &mut stats,
        );
        engine_free(ctx);
        result
    };
    assert_eq!(code, 0, "propagation error");
    assert!(stats.nout >= 2);
    let last = &out_y[6 * (stats.nout - 1)..6 * stats.nout];
    let distance = |range: std::ops::Range<usize>| {
        range
            .map(|i| (last[i] - initial[i]).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    (distance(0..3), distance(3..6))
}
fn main() {
    let (position, velocity) = run_case(1., false);
    println!(
        "One orbital period: position error {position:.6} m, velocity error {velocity:.9} m/s"
    );
    assert!(position < 0.01 && velocity < 1e-5);
}
#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "msis21")]
    #[test]
    #[ignore = "requires RMSIS21_LIBRARY and RMSIS21_PARAMETERS"]
    fn native_rmsis_drag_orbit() {
        let (p, v) = run_case(1., true);
        assert!(p.is_finite() && v.is_finite() && p > 0.01);
    }
    #[test]
    fn analytic_orbit_forward() {
        let (p, v) = run_case(1., false);
        assert!(p < 0.01 && v < 1e-5, "{p} m, {v} m/s");
    }
    #[test]
    fn analytic_orbit_backward() {
        let (p, v) = run_case(-1., false);
        assert!(p < 0.01 && v < 1e-5, "{p} m, {v} m/s");
    }
}
