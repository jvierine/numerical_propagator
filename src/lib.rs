//! NASA GMAT compatible numerical orbit propagator written in Rust.
//!
//! Adaptive DOP853 integration, fully normalized spherical harmonics, supplied
//! GCRS frame/ephemeris tables, and optional MSIS density callbacks.
//! See README.md for the scope of GMAT compatibility and validation.
mod tableau;
use std::{
    panic::{catch_unwind, AssertUnwindSafe},
    slice,
};
use tableau::{A, B, C, D, E3, E5};

type V = [f64; 3];
type M = [f64; 9];
type Y = [f64; 9];
// Smoothly limit intermediate coefficient excursions to exp(+-2) relative
// to the interpolation of its endpoints. Endpoint coefficients stay fitted.
fn smooth_beta_bend(u: f64, bend: f64, skew: f64) -> f64 {
    let raw = 16. * (bend + 4. * skew * (2. * u - 1.)) * u * u * (1. - u) * (1. - u);
    2. * (raw / 2.).tanh()
}
fn warp_time(u: f64, w: f64) -> f64 {
    if w == 0. {
        u
    } else {
        u / (u + w.exp() * (1. - u))
    }
}
fn linear_drag_coefficient(knots: &[(f64, f64)], t: f64) -> f64 {
    let j = knots.partition_point(|&(epoch, _)| epoch <= t);
    if j == 0 {
        return knots[0].1;
    }
    if j == knots.len() {
        return knots[j - 1].1;
    }
    let (a, ba) = knots[j - 1];
    let (b, bb) = knots[j];
    let u = (t - a) / (b - a);
    (1. - u) * ba + u * bb
}
// A linear spline in log(beta) is positive and bounded without slope limiters.
fn log_linear_drag_coefficient(knots: &[(f64, f64)], t: f64) -> f64 {
    let j = knots.partition_point(|&(epoch, _)| epoch <= t);
    if j == 0 {
        return knots[0].1;
    }
    if j == knots.len() {
        return knots[j - 1].1;
    }
    let (a, ba) = knots[j - 1];
    let (b, bb) = knots[j];
    let u = (t - a) / (b - a);
    ((1. - u) * ba.ln() + u * bb.ln()).exp()
}
// Shape-preserving cubic Hermite interpolation of log(beta). Each interval
// stays between its endpoint values, so node bounds apply to the whole arc.
fn log_spline_drag_coefficient(knots: &[(f64, f64)], t: f64) -> f64 {
    let j = knots.partition_point(|&(epoch, _)| epoch <= t);
    if j == 0 {
        return knots[0].1;
    }
    if j == knots.len() {
        return knots[j - 1].1;
    }
    let secant = |i: usize| (knots[i + 1].1.ln() - knots[i].1.ln()) / (knots[i + 1].0 - knots[i].0);
    let slope = |i: usize| {
        if i == 0 {
            return secant(0);
        }
        if i + 1 == knots.len() {
            return secant(i - 1);
        }
        let a = secant(i - 1);
        let b = secant(i);
        if a * b <= 0. {
            0.
        } else {
            let h0 = knots[i].0 - knots[i - 1].0;
            let h1 = knots[i + 1].0 - knots[i].0;
            let w0 = 2. * h1 + h0;
            let w1 = h1 + 2. * h0;
            (w0 + w1) / (w0 / a + w1 / b)
        }
    };
    let i = j - 1;
    let h = knots[j].0 - knots[i].0;
    let u = (t - knots[i].0) / h;
    let z = (2. * u * u * u - 3. * u * u + 1.) * knots[i].1.ln()
        + (u * u * u - 2. * u * u + u) * h * slope(i)
        + (-2. * u * u * u + 3. * u * u) * knots[j].1.ln()
        + (u * u * u - u * u) * h * slope(j);
    z.exp()
}
#[cfg(test)]
mod log_drag_spline_tests {
    use super::*;
    #[test]
    fn log_linear_has_geometric_midpoint_and_no_overshoot() {
        let knots = [(0., 0.01), (10., 0.09), (20., 0.01)];
        assert!((log_linear_drag_coefficient(&knots, 5.) - 0.03).abs() < 1e-14);
        for i in 0..1001 {
            let b = log_linear_drag_coefficient(&knots, i as f64 / 50.);
            assert!((0.01 - 1e-14..=0.09 + 1e-14).contains(&b));
        }
    }
    #[test]
    fn preserves_nodes_and_factor_ten_between_nodes() {
        for values in [
            [0.01, 0.10, 0.01, 0.10, 0.01],
            [0.01, 0.02, 0.05, 0.08, 0.10],
            [0.10, 0.10, 0.02, 0.01, 0.01],
        ] {
            let knots: Vec<_> = values
                .into_iter()
                .enumerate()
                .map(|(i, v)| (i as f64 * 17., v))
                .collect();
            for &(t, b) in &knots {
                assert!((log_spline_drag_coefficient(&knots, t) - b).abs() < 1e-12);
            }
            for i in 0..1001 {
                let b = log_spline_drag_coefficient(&knots, 68. * i as f64 / 1000.);
                assert!((0.01 - 1e-12..=0.10 + 1e-12).contains(&b));
            }
        }
    }
}
const DAY: f64 = 86400.;
const OMEGA: f64 = 7.292115e-5;
type Msis64 = fn(f64, f64, f64, f64, f64, f64, f64, [f64; 7]) -> f64;

#[repr(C)]
pub struct Inputs {
    pub degree: usize,
    pub gm: f64,
    pub radius: f64,
    pub coefficients_c: *const f64,
    pub coefficients_s: *const f64,
    pub frame_start: f64,
    pub frame_step: f64,
    pub nframe: usize,
    pub q: *const f64,
    pub polar: *const f64,
    pub dut: *const f64,
    pub body_start: f64,
    pub body_step: f64,
    pub nbody: usize,
    pub sun_pos: *const f64,
    pub sun_vel: *const f64,
    pub moon_pos: *const f64,
    pub moon_vel: *const f64,
    pub weather_start: f64,
    pub weather_step: f64,
    pub nweather: usize,
    pub f107: *const f64,
    pub f107a: *const f64,
    pub ap: *const f64,
    pub doy: *const f64,
    pub correction_start: f64,
    pub correction_step: f64,
    pub ncorrection: usize,
    pub log_correction: *const f64,
    pub correction_height_linear: *const f64,
    pub correction_height_quadratic: *const f64,
    pub correction_height_origin: f64,
    pub correction_height_scale: f64,
    pub correction_clip: f64,
    pub beta: f64,
    pub drag: i32,
    pub third_bodies: i32,
}

#[repr(C)]
#[derive(Default)]
pub struct Stats {
    pub nout: usize,
    pub nfev: usize,
    pub accepted: usize,
    pub rejected: usize,
    pub event: f64,
}

fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: V) -> f64 {
    dot(a, a).sqrt()
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn mv(a: M, v: V) -> V {
    std::array::from_fn(|i| a[3 * i] * v[0] + a[3 * i + 1] * v[1] + a[3 * i + 2] * v[2])
}
fn tmv(a: M, v: V) -> V {
    std::array::from_fn(|i| a[i] * v[0] + a[3 + i] * v[1] + a[6 + i] * v[2])
}
fn mm(a: M, b: M) -> M {
    std::array::from_fn(|i| (0..3).map(|k| a[3 * (i / 3) + k] * b[3 * k + i % 3]).sum())
}
fn pos(y: &Y) -> V {
    [y[0], y[1], y[2]]
}
fn vel(y: &Y) -> V {
    [y[3], y[4], y[5]]
}

fn geodetic(p: V) -> V {
    let rxy = p[0].hypot(p[1]);
    let z = p[2];
    let e2 = 6.69437999014e-3;
    let mut lat = z.atan2(rxy * (1. - e2));
    let mut alt = 0.;
    for _ in 0..6 {
        let normal = 6378137. / (1. - e2 * lat.sin().powi(2)).sqrt();
        alt = rxy / lat.cos() - normal;
        lat = z.atan2(rxy * (1. - e2 * normal / (normal + alt)));
    }
    [p[1].atan2(p[0]).to_degrees(), lat.to_degrees(), alt / 1000.]
}

struct Harmonics {
    n: usize,
    gm: f64,
    radius: f64,
    c: Vec<f64>,
    s: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    powers: Vec<f64>,
}
impl Harmonics {
    unsafe fn new(i: &Inputs) -> Self {
        let d = i.degree + 1;
        let source_c = slice::from_raw_parts(i.coefficients_c, d * d);
        let source_s = slice::from_raw_parts(i.coefficients_s, d * d);
        let mut g = Self {
            n: i.degree,
            gm: i.gm,
            radius: i.radius,
            c: vec![0.; d * d],
            s: vec![0.; d * d],
            a: vec![0.; d * d],
            b: vec![0.; d * d],
            powers: vec![0.; d],
        };
        // Column-major recurrence, no allocation or Legendre matrix per RHS.
        for m in 0..d {
            for n in m..d {
                let k = m * d + n;
                g.c[k] = source_c[n * d + m];
                g.s[k] = source_s[n * d + m];
                if n >= m + 2 {
                    let nf = n as f64;
                    let mf = m as f64;
                    g.a[k] = ((2. * nf - 1.) * (2. * nf + 1.) / ((nf - mf) * (nf + mf))).sqrt();
                    g.b[k] = ((2. * nf + 1.) * (nf + mf - 1.) * (nf - mf - 1.)
                        / ((nf - mf) * (nf + mf) * (2. * nf - 3.)))
                        .sqrt();
                }
            }
        }
        g
    }
    fn base(&mut self, p: V) -> (f64, V) {
        let r = norm(p);
        let u = p[0].hypot(p[1]) / r;
        let z = p[2] / r;
        let (sl, cl) = p[1].atan2(p[0]).sin_cos();
        let d = self.n + 1;
        let mut potential = 0.;
        let mut radial = 0.;
        let mut latitude = 0.;
        let mut longitude = 0.;
        self.powers[0] = 1.;
        for n in 1..d {
            self.powers[n] = self.powers[n - 1] * self.radius / r;
        }
        let mut diagonal = 1.;
        let mut ddiagonal = 0.;
        let mut cm = 1.;
        let mut sm = 0.;
        for m in 0..d {
            if m == 1 {
                diagonal = 3_f64.sqrt() * u;
                ddiagonal = -3_f64.sqrt() * z;
            } else if m > 1 {
                let factor = ((2 * m + 1) as f64 / (2 * m) as f64).sqrt();
                ddiagonal = factor * (-z * diagonal + u * ddiagonal);
                diagonal *= factor * u;
            }
            let mut pp = 0.;
            let mut ddp = 0.;
            let mut pn = diagonal;
            let mut dpn = ddiagonal;
            for n in m..d {
                let k = m * d + n;
                if n >= 1 {
                    let az = self.c[k] * cm + self.s[k] * sm;
                    let term = self.powers[n] * pn * az;
                    potential += term;
                    radial += (n + 1) as f64 * term;
                    latitude += self.powers[n] * dpn * az;
                    longitude +=
                        self.powers[n] * pn * m as f64 * (-self.c[k] * sm + self.s[k] * cm);
                }
                if n < self.n {
                    let next = n + 1;
                    let nk = m * d + next;
                    let (newp, newdp) = if next == m + 1 {
                        let factor = (2 * m + 3) as f64;
                        (factor.sqrt() * z * pn, factor.sqrt() * (u * pn + z * dpn))
                    } else {
                        (
                            self.a[nk] * z * pn - self.b[nk] * pp,
                            self.a[nk] * (u * pn + z * dpn) - self.b[nk] * ddp,
                        )
                    };
                    pp = pn;
                    ddp = dpn;
                    pn = newp;
                    dpn = newdp;
                }
            }
            let nextcm = cm * cl - sm * sl;
            sm = sm * cl + cm * sl;
            cm = nextcm;
        }
        potential += self.c[0];
        radial += self.c[0];
        let factor = self.gm / (r * r);
        let ar = -factor * radial;
        let alat = factor * latitude;
        let alon = factor * longitude / u.max(1e-30);
        (
            self.gm / r * potential,
            [
                ar * u * cl - alat * z * cl - alon * sl,
                ar * u * sl - alat * z * sl + alon * cl,
                ar * z + alat * u,
            ],
        )
    }
    fn evaluate(&mut self, p: V) -> (f64, V) {
        let (pot, mut acc) = self.base(p);
        if p[0].hypot(p[1]) < norm(p) * 1e-9 {
            for j in 0..3 {
                let mut plus = p;
                let mut minus = p;
                plus[j] += 1.;
                minus[j] -= 1.;
                acc[j] = (self.base(plus).0 - self.base(minus).0) / 2.;
            }
        }
        (pot, acc)
    }
}

pub struct Context {
    msis_f64: Option<Msis64>,
    weather_slot_override: Option<usize>,
    gravity: Harmonics,
    frame_start: f64,
    frame_step: f64,
    q: Vec<f64>,
    polar: Vec<f64>,
    dut: Vec<f64>,
    body_start: f64,
    body_step: f64,
    sun_pos: Vec<f64>,
    sun_vel: Vec<f64>,
    moon_pos: Vec<f64>,
    moon_vel: Vec<f64>,
    weather_start: f64,
    weather_step: f64,
    f107: Vec<f64>,
    f107a: Vec<f64>,
    ap: Vec<f64>,
    doy: Vec<f64>,
    correction_start: f64,
    correction_step: f64,
    log_correction: Vec<f64>,
    correction_height_linear: Vec<f64>,
    correction_height_quadratic: Vec<f64>,
    correction_height_origin: f64,
    correction_height_scale: f64,
    correction_clip: f64,
    beta: f64,
    drag: bool,
    // Fragment geometry evolution, independent of atmospheric density.
    ballistic_height_profile: Option<(f64, f64, f64)>,
    ballistic_height_power: f64,
    ballistic_time_profile: Option<(f64, f64, f64, f64, f64)>,
    ballistic_time_warp: f64,
    ballistic_linear_knots: Vec<(f64, f64)>,
    ballistic_log_spline: bool,
    ballistic_log_linear: bool,
    lateral_drag_ratio_knots: Vec<(f64, f64)>,
    ballistic_terminal_control: Option<(f64, f64, f64)>,
    smooth_density_correction: bool,
    third_bodies: bool,
    #[cfg(feature = "msis21")]
    msis21: Option<rmsis::msis21::Model>,
    nfev: usize,
}
impl Context {
    /// Load the native rmsis MSIS2.1 adapter and original parameter file.
    ///
    /// # Safety
    /// The library must be trusted and implement the rmsis MSIS2.1 ABI.
    #[cfg(feature = "msis21")]
    pub unsafe fn load_msis21(
        &mut self,
        library: impl AsRef<std::path::Path>,
        parameters: impl AsRef<std::path::Path>,
    ) -> Result<(), String> {
        self.msis21 = Some(rmsis::msis21::Model::load(library, parameters)?);
        self.msis_f64 = None;
        Ok(())
    }
    fn native_density(
        &self,
        day: f64,
        ut: f64,
        lon: f64,
        lat: f64,
        z: f64,
        f: f64,
        fa: f64,
        ap: [f64; 7],
    ) -> Result<f64, i32> {
        #[cfg(feature = "msis21")]
        {
            let input = rmsis::msis21::Input {
                day_of_year: day,
                seconds_of_day: ut,
                longitude_deg: lon,
                latitude_deg: lat,
                altitude_km: z,
                f107: f,
                f107_average: fa,
                ap,
                storm_time: true,
            };
            return self
                .msis21
                .as_ref()
                .ok_or(-1)?
                .evaluate(&input)
                .map_err(|_| -3)?
                .densities[0]
                .ok_or(-3);
        }
        #[cfg(not(feature = "msis21"))]
        {
            let _ = (day, ut, lon, lat, z, f, fa, ap);
            Err(-1)
        }
    }

    /// Install a double precision density callback. Arguments: day of year,
    /// UTC seconds of day, longitude/latitude (degrees), height (km), F10.7,
    /// mean F10.7, seven Ap values. Return density in kg/m^3.
    /// Serialize calls when the supplied MSIS model is not thread-safe.
    pub fn set_density_callback(
        &mut self,
        evaluate: fn(f64, f64, f64, f64, f64, f64, f64, [f64; 7]) -> f64,
    ) {
        self.msis_f64 = Some(evaluate);
    }
    /// Set positive, time-ordered effective Cd*A/m nodes (m^2/kg), log-linear
    /// between nodes and constant outside them. Times use the input UTC scale.
    pub fn set_drag_nodes(&mut self, knots: &[(f64, f64)]) -> Result<(), &'static str> {
        if knots.is_empty()
            || knots
                .iter()
                .any(|(t, b)| !t.is_finite() || !b.is_finite() || *b <= 0.)
            || knots.windows(2).any(|w| w[1].0 <= w[0].0)
        {
            return Err("nodes must have increasing finite times and positive finite coefficients");
        }
        self.ballistic_linear_knots = knots.to_vec();
        self.ballistic_log_linear = true;
        self.ballistic_log_spline = false;
        Ok(())
    }

    unsafe fn new(i: &Inputs) -> Self {
        let copy = |ptr, n| slice::from_raw_parts(ptr, n).to_vec();
        Self {
            msis_f64: None,
            weather_slot_override: None,
            gravity: Harmonics::new(i),
            frame_start: i.frame_start,
            frame_step: i.frame_step,
            q: copy(i.q, 9 * i.nframe),
            polar: copy(i.polar, 9 * i.nframe),
            dut: copy(i.dut, i.nframe),
            body_start: i.body_start,
            body_step: i.body_step,
            sun_pos: copy(i.sun_pos, 3 * i.nbody),
            sun_vel: copy(i.sun_vel, 3 * i.nbody),
            moon_pos: copy(i.moon_pos, 3 * i.nbody),
            moon_vel: copy(i.moon_vel, 3 * i.nbody),
            weather_start: i.weather_start,
            weather_step: i.weather_step,
            f107: copy(i.f107, i.nweather),
            f107a: copy(i.f107a, i.nweather),
            ap: copy(i.ap, 7 * i.nweather),
            doy: copy(i.doy, i.nweather),
            correction_start: i.correction_start,
            correction_step: i.correction_step,
            log_correction: if i.ncorrection > 0 {
                copy(i.log_correction, i.ncorrection)
            } else {
                vec![]
            },
            beta: i.beta,
            ballistic_height_profile: None,
            ballistic_height_power: 1.,
            ballistic_time_profile: None,
            ballistic_time_warp: 0.,
            ballistic_linear_knots: vec![],
            ballistic_log_spline: false,
            ballistic_log_linear: false,
            lateral_drag_ratio_knots: vec![],
            ballistic_terminal_control: None,
            smooth_density_correction: false,
            correction_height_linear: if i.ncorrection > 0 && !i.correction_height_linear.is_null()
            {
                copy(i.correction_height_linear, i.ncorrection)
            } else {
                vec![]
            },
            correction_height_quadratic: if i.ncorrection > 0
                && !i.correction_height_quadratic.is_null()
            {
                copy(i.correction_height_quadratic, i.ncorrection)
            } else {
                vec![]
            },
            correction_height_origin: i.correction_height_origin,
            correction_height_scale: i.correction_height_scale,
            correction_clip: i.correction_clip,
            drag: i.drag != 0,
            third_bodies: i.third_bodies != 0,
            #[cfg(feature = "msis21")]
            msis21: None,
            nfev: 0,
        }
    }
    fn segment(t: f64, start: f64, step: f64, n: usize) -> Result<(usize, f64), i32> {
        let x = (t - start) / step;
        let j = x.floor();
        if j < 0. || j as usize + 1 >= n {
            return Err(-2);
        }
        Ok((j as usize, x - j))
    }
    fn rotation(&self, t: f64) -> Result<M, i32> {
        let (j, f) = Self::segment(t, self.frame_start, self.frame_step, self.dut.len())?;
        let q = std::array::from_fn(|k| self.q[9 * j + k] * (1. - f) + self.q[9 * (j + 1) + k] * f);
        let p = std::array::from_fn(|k| {
            self.polar[9 * j + k] * (1. - f) + self.polar[9 * (j + 1) + k] * f
        });
        let dut = self.dut[j] * (1. - f) + self.dut[j + 1] * f;
        let days = (t - 946728000. + dut) / DAY;
        let angle =
            std::f64::consts::TAU * (0.7790572732640 + 1.00273781191135448 * days).rem_euclid(1.);
        let (s, c) = angle.sin_cos();
        Ok(mm(p, mm([c, s, 0., -s, c, 0., 0., 0., 1.], q)))
    }
    fn body(&self, t: f64, p: &[f64], v: &[f64]) -> Result<V, i32> {
        let (j, f) = Self::segment(t, self.body_start, self.body_step, p.len() / 3)?;
        Ok(std::array::from_fn(|k| {
            (2. * f.powi(3) - 3. * f * f + 1.) * p[3 * j + k]
                + (f.powi(3) - 2. * f * f + f) * self.body_step * v[3 * j + k]
                + (-2. * f.powi(3) + 3. * f * f) * p[3 * (j + 1) + k]
                + (f.powi(3) - f * f) * self.body_step * v[3 * (j + 1) + k]
        }))
    }
    fn environment(&self, t: f64, y: &Y) -> Result<(M, V, V, f64, f64), i32> {
        let rot = self.rotation(t)?;
        let fixed = mv(rot, pos(y));
        let [lon, lat, alt] = geodetic(fixed);
        let j = self
            .weather_slot_override
            .map(|i| i as f64)
            .unwrap_or_else(|| ((t - self.weather_start) / self.weather_step).floor());
        if j < 0. || j as usize >= self.f107.len() {
            return Err(-2);
        }
        let j = j as usize;
        let mut rho = 0.;
        if self.drag {
            if let Some(evaluate) = self.msis_f64 {
                rho = evaluate(
                    self.doy[j],
                    t.rem_euclid(DAY),
                    lon,
                    lat,
                    alt.max(0.),
                    self.f107[j],
                    self.f107a[j],
                    std::array::from_fn(|k| self.ap[7 * j + k]),
                );
            } else {
                rho = self.native_density(
                    self.doy[j],
                    t.rem_euclid(DAY),
                    lon,
                    lat,
                    alt.max(0.),
                    self.f107[j],
                    self.f107a[j],
                    std::array::from_fn(|k| self.ap[7 * j + k]),
                )?;
            }
            if rho >= 9.9e-38 && rho < 1e-37 {
                return Err(-3);
            }
            if !self.log_correction.is_empty() {
                let (i, f) = Self::segment(
                    t,
                    self.correction_start,
                    self.correction_step,
                    self.log_correction.len(),
                )?;
                let f = if self.smooth_density_correction {
                    f * f * (3. - 2. * f)
                } else {
                    f
                };
                let mut log_factor =
                    (1. - f) * self.log_correction[i] + f * self.log_correction[i + 1];
                if !self.correction_height_linear.is_empty() {
                    let z = (alt - self.correction_height_origin) / self.correction_height_scale;
                    let linear = (1. - f) * self.correction_height_linear[i]
                        + f * self.correction_height_linear[i + 1];
                    let quadratic = (1. - f) * self.correction_height_quadratic[i]
                        + f * self.correction_height_quadratic[i + 1];
                    log_factor += linear * z + quadratic * z * z;
                    if self.correction_clip > 0. {
                        log_factor = log_factor.clamp(-self.correction_clip, self.correction_clip);
                    }
                }
                rho *= log_factor.exp();
            }
        }
        if !rho.is_finite() || rho < 0. {
            return Err(-3);
        }
        let vv = mv(rot, vel(y));
        let omega = cross([0., 0., OMEGA], fixed);
        let relative = std::array::from_fn(|k| vv[k] - omega[k]);
        let speed = norm(relative);
        let unit = tmv(
            rot,
            std::array::from_fn(|k| -0.5 * rho * speed * relative[k]),
        );
        Ok((rot, fixed, unit, rho, alt))
    }
    fn rhs(&mut self, t: f64, y: &Y, n: usize) -> Result<Y, i32> {
        self.nfev += 1;
        let (rot, fixed, unit, _, _) = self.environment(t, y)?;
        let acc_fixed = self.gravity.evaluate(fixed).1;
        let mut acc = tmv(rot, acc_fixed);
        let mut perturb = [0.; 3];
        if self.third_bodies {
            for (body, gm) in [
                (
                    self.body(t, &self.sun_pos, &self.sun_vel)?,
                    1.32712440018e20,
                ),
                (
                    self.body(t, &self.moon_pos, &self.moon_vel)?,
                    4.902800066e12,
                ),
            ] {
                let delta = std::array::from_fn(|k| body[k] - y[k]);
                let nr = norm(delta).powi(3);
                let nb = norm(body).powi(3);
                for k in 0..3 {
                    perturb[k] += gm * (delta[k] / nr - body[k] / nb);
                }
            }
        }
        let beta = self.effective_beta(t, fixed);
        let lateral = self.lateral_aerodynamic_acceleration(t, y, unit, beta);
        for k in 0..3 {
            acc[k] += beta * unit[k] + lateral[k] + perturb[k];
        }
        let mut out = [0.; 9];
        out[..3].copy_from_slice(&y[3..6]);
        out[3..6].copy_from_slice(&acc);
        if n == 9 {
            out[6] = dot(vel(y), unit) + dot(vel(y), lateral) / beta.max(1e-30);
            out[7] = dot(vel(y), perturb);
            out[8] = dot(acc_fixed, cross([0., 0., OMEGA], fixed));
        }
        if !out[..n].iter().all(|v| v.is_finite()) {
            return Err(-3);
        }
        Ok(out)
    }
    fn lateral_aerodynamic_acceleration(&self, t: f64, y: &Y, drag_per_beta: V, beta: f64) -> V {
        if self.lateral_drag_ratio_knots.is_empty() {
            return [0.; 3];
        }
        let d = norm(drag_per_beta);
        if d < 1e-30 {
            return [0.; 3];
        }
        let flow = drag_per_beta.map(|v| v / d);
        let normal = cross(pos(y), vel(y));
        let projected: V = std::array::from_fn(|j| normal[j] - dot(normal, flow) * flow[j]);
        let n = norm(projected);
        if n < 1e-30 {
            return [0.; 3];
        }
        let ratio = linear_drag_coefficient(&self.lateral_drag_ratio_knots, t);
        projected.map(|v| v / n * ratio * beta * d)
    }
    fn effective_beta(&self, t: f64, fixed: V) -> f64 {
        if !self.ballistic_linear_knots.is_empty() {
            return if self.ballistic_log_spline {
                if self.ballistic_log_linear {
                    log_linear_drag_coefficient(&self.ballistic_linear_knots, t)
                } else {
                    log_spline_drag_coefficient(&self.ballistic_linear_knots, t)
                }
            } else {
                linear_drag_coefficient(&self.ballistic_linear_knots, t)
            };
        }
        let beta = if let Some((low, high, log_ratio)) = self.ballistic_height_profile {
            self.beta
                * (((geodetic(fixed)[2] - low) / (high - low))
                    .clamp(0., 1.)
                    .powf(self.ballistic_height_power)
                    * log_ratio)
                    .exp()
        } else {
            self.beta
        };
        let beta = if let Some((start, end, log_initial_ratio, bend, skew)) =
            self.ballistic_time_profile
        {
            let u = warp_time(
                ((t - start) / (end - start)).clamp(0., 1.),
                self.ballistic_time_warp,
            );
            beta * ((1. - u * u * (3. - 2. * u)) * log_initial_ratio
                + smooth_beta_bend(u, bend, skew))
            .exp()
        } else {
            beta
        };
        if let Some((start, end, amplitude)) = self.ballistic_terminal_control {
            beta * smooth_beta_bend(((t - start) / (end - start)).clamp(0., 1.), amplitude, 0.)
                .exp()
        } else {
            beta
        }
    }
    fn altitude(&self, t: f64, y: &Y) -> Result<f64, i32> {
        Ok(geodetic(mv(self.rotation(t)?, pos(y)))[2])
    }
}

fn dense(
    ctx: &mut Context,
    start: f64,
    t: f64,
    h: f64,
    old: Y,
    new: Y,
    k: &mut [Y; 16],
    n: usize,
) -> Result<[Y; 7], i32> {
    for stage in 13..16 {
        let mut yi = old;
        for j in 0..n {
            for l in 0..stage {
                yi[j] += h * A[stage][l] * k[l][j];
            }
        }
        k[stage] = ctx.rhs(start + t + C[stage] * h, &yi, n)?;
    }
    let mut f = [[0.; 9]; 7];
    for j in 0..n {
        let dy = new[j] - old[j];
        f[0][j] = dy;
        f[1][j] = h * k[0][j] - dy;
        f[2][j] = 2. * dy - h * (k[12][j] + k[0][j]);
        for l in 0..4 {
            for s in 0..16 {
                f[3 + l][j] += h * D[l][s] * k[s][j];
            }
        }
    }
    Ok(f)
}
fn interpolate(old: Y, f: &[Y; 7], x: f64, n: usize) -> Y {
    let mut y = [0.; 9];
    for (i, row) in f.iter().rev().enumerate() {
        for j in 0..n {
            y[j] = (y[j] + row[j]) * if i % 2 == 0 { x } else { 1. - x };
        }
    }
    for j in 0..n {
        y[j] += old[j];
    }
    y
}
fn integrate(
    ctx: &mut Context,
    start: f64,
    end: f64,
    initial: &[f64],
    rtol: f64,
    max_step: f64,
    threshold: f64,
    samples: &[f64],
    out_t: &mut [f64],
    out_y: &mut [f64],
    stats: &mut Stats,
) -> Result<(), i32> {
    let result = integrate_inner(
        ctx, start, end, initial, rtol, max_step, threshold, samples, out_t, out_y, stats,
    );
    ctx.weather_slot_override = None;
    result
}
fn integrate_inner(
    ctx: &mut Context,
    start: f64,
    end: f64,
    initial: &[f64],
    rtol: f64,
    max_step: f64,
    threshold: f64,
    samples: &[f64],
    out_t: &mut [f64],
    out_y: &mut [f64],
    stats: &mut Stats,
) -> Result<(), i32> {
    let n = initial.len();
    let mut y = [0.; 9];
    y[..n].copy_from_slice(initial);
    if !matches!(n, 6 | 9) || end == start || max_step <= 0. || rtol <= 0. {
        return Err(-1);
    }
    let direction = (end - start).signum();
    if samples
        .iter()
        .any(|t| direction * (*t - start) < 0. || direction * (*t - end) > 0.)
        || samples.windows(2).any(|p| direction * (p[1] - p[0]) <= 0.)
    {
        return Err(-1);
    }
    let base_nfev = ctx.nfev;
    // Optional optimization-trial budget; unset for unrestricted propagation.
    let fit_rhs_limit = std::env::var("PROPAGATOR_RHS_LIMIT")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    let mut k = [[0.; 9]; 16];
    let slot = |time: f64, direction: f64| -> Result<usize, i32> {
        let z = (time - ctx.weather_start + direction * 1e-6) / ctx.weather_step;
        let i = if direction > 0. {
            z.floor()
        } else {
            z.ceil() - 1.
        };
        if i < 0. || i as usize >= ctx.f107.len() {
            Err(-2)
        } else {
            Ok(i as usize)
        }
    };
    let initial_slot = slot(start, direction)?;
    ctx.weather_slot_override = Some(initial_slot);
    k[0] = ctx.rhs(start, &y, n)?;
    let mut t = 0.;
    let duration = end - start;
    let mut h = direction * max_step.min(0.01);
    let mut sample = 0;
    let mut last_alt = ctx.altitude(start, &y)? - threshold;
    if last_alt <= 0. && direction > 0. {
        return Err(-4);
    }
    let event_sign = if last_alt < 0. { -1. } else { 1. };
    last_alt *= event_sign;
    let mut previously_rejected = false;
    while sample < samples.len() && samples[sample] == start {
        out_t[stats.nout] = start;
        out_y[n * stats.nout..n * (stats.nout + 1)].copy_from_slice(&y[..n]);
        stats.nout += 1;
        sample += 1;
    }
    while direction * (duration - t) > 0. {
        if fit_rhs_limit.is_some_and(|limit| ctx.nfev - base_nfev > limit) {
            return Err(-7);
        }
        h = direction * h.abs().min(max_step).min((duration - t).abs());
        // Space-weather drivers are piecewise constant. Stop exactly at their
        // boundaries and use one-sided forcing for RK endpoint stages. This
        // avoids a parameter-dependent error from stepping across a driver
        // jump, and refreshes the FSAL derivative when the forcing changes.
        let now = start + t;
        let z = (now - ctx.weather_start + direction * 1e-6) / ctx.weather_step;
        let slot = if direction > 0. {
            z.floor()
        } else {
            z.ceil() - 1.
        };
        if slot < 0. || slot as usize >= ctx.f107.len() {
            return Err(-2);
        }
        let slot = slot as usize;
        if ctx.weather_slot_override != Some(slot) {
            ctx.weather_slot_override = Some(slot);
            k[0] = ctx.rhs(now, &y, n)?;
        }
        let boundary = ctx.weather_start
            + (slot as f64 + if direction > 0. { 1. } else { 0. }) * ctx.weather_step;
        let distance = direction * (boundary - now);
        if distance > 1e-7 {
            h = direction * h.abs().min(distance);
        }
        if h.abs() < 1e-10 {
            return Err(-5);
        }
        let old = y;
        let mut new = y;
        for s in 1..12 {
            let mut yi = y;
            for j in 0..n {
                let mut sum = 0.;
                for l in 0..s {
                    sum += A[s][l] * k[l][j];
                }
                yi[j] += h * sum;
            }
            k[s] = ctx.rhs(start + t + C[s] * h, &yi, n)?;
        }
        for j in 0..n {
            let mut sum = 0.;
            for s in 0..12 {
                sum += B[s] * k[s][j];
            }
            new[j] += h * sum;
        }
        k[12] = ctx.rhs(start + t + h, &new, n)?;
        let mut err5 = 0.;
        let mut err3 = 0.;
        for j in 0..n {
            let atol = if j < 3 {
                1e-5
            } else if j < 6 {
                1e-8
            } else {
                0.02
            };
            let scale = atol + rtol * old[j].abs().max(new[j].abs());
            let mut e5 = 0.;
            let mut e3 = 0.;
            for s in 0..13 {
                e5 += E5[s] * k[s][j];
                e3 += E3[s] * k[s][j];
            }
            err5 += (e5 / scale).powi(2);
            err3 += (e3 / scale).powi(2);
        }
        let error = if err5 == 0. && err3 == 0. {
            0.
        } else {
            h.abs() * err5 / ((err5 + 0.01 * err3) * n as f64).sqrt()
        };
        if !error.is_finite() {
            return Err(-3);
        }
        if error >= 1. {
            h *= (0.9 * error.powf(-0.125)).max(0.2);
            stats.rejected += 1;
            previously_rejected = true;
            continue;
        }
        stats.accepted += 1;
        let next = t + h;
        let next_alt = event_sign * (ctx.altitude(start + next, &new)? - threshold);
        let event = last_alt > 0. && next_alt <= 0.;
        let needs_dense = event
            || (sample < samples.len() && direction * (samples[sample] - (start + next)) <= 0.);
        if needs_dense {
            let f = dense(ctx, start, t, h, old, new, &mut k, n)?;
            let mut event_dt = next;
            if event {
                let mut low = 0.;
                let mut high = 1.;
                for _ in 0..50 {
                    let mid = (low + high) / 2.;
                    let value = interpolate(old, &f, mid, n);
                    if event_sign * (ctx.altitude(start + t + mid * h, &value)? - threshold) > 0. {
                        low = mid;
                    } else {
                        high = mid;
                    }
                    if (high - low) * h.abs() < 1e-7 {
                        break;
                    }
                }
                event_dt = t + (low + high) * h / 2.;
                stats.event = start + event_dt;
            }
            while sample < samples.len() && direction * (samples[sample] - (start + event_dt)) <= 0.
            {
                let point = interpolate(old, &f, ((samples[sample] - start) - t) / h, n);
                out_t[stats.nout] = samples[sample];
                out_y[n * stats.nout..n * (stats.nout + 1)].copy_from_slice(&point[..n]);
                stats.nout += 1;
                sample += 1;
            }
            if event {
                let value = interpolate(old, &f, (event_dt - t) / h, n);
                out_t[stats.nout] = stats.event;
                out_y[n * stats.nout..n * (stats.nout + 1)].copy_from_slice(&value[..n]);
                stats.nout += 1;
                stats.nfev = ctx.nfev - base_nfev;
                return Ok(());
            }
        }
        y = new;
        t = next;
        k[0] = k[12];
        last_alt = next_alt;
        let factor = if error == 0. {
            10.
        } else {
            (0.9 * error.powf(-0.125)).clamp(0.2, 10.)
        };
        h *= if previously_rejected {
            factor.min(1.)
        } else {
            factor
        };
        previously_rejected = false;
    }
    stats.nfev = ctx.nfev - base_nfev;
    Ok(())
}

// Inputs and output buffers must satisfy the explicit sizes in Inputs/the call.
// Rust owns copied input arrays; callers own all input and output buffers.
/// Construct a context by copying input arrays.
///
/// # Safety
/// `input` and every array pointer must be aligned, readable, and correctly sized
/// (including non-null pointers for zero-length arrays). See README.md.
#[no_mangle]
pub unsafe extern "C" fn engine_new(input: *const Inputs) -> *mut Context {
    match catch_unwind(AssertUnwindSafe(|| {
        Box::into_raw(Box::new(Context::new(&*input)))
    })) {
        Ok(p) => p,
        Err(_) => std::ptr::null_mut(),
    }
}
/// Release a context.
///
/// # Safety
/// The pointer must be null or a live result of `engine_new`, freed exactly once.
/// No concurrent or subsequent access is allowed.
#[no_mangle]
pub unsafe extern "C" fn engine_free(ctx: *mut Context) {
    if !ctx.is_null() {
        drop(Box::from_raw(ctx));
    }
}
/// Set constant Cd*A/m in m^2/kg. Drag nodes, if installed, take precedence.
///
/// # Safety
/// `ctx` must be a live, exclusively accessed context.
#[no_mangle]
pub unsafe extern "C" fn engine_beta(ctx: *mut Context, beta: f64) {
    (*ctx).beta = beta;
}
/// Evaluate Earth-fixed potential and acceleration (SI units).
///
/// # Safety
/// `ctx` must be live and exclusive; `p` has three readable doubles and `out`
/// has four writable doubles. Buffers must not overlap the context.
#[no_mangle]
pub unsafe extern "C" fn engine_gravity(ctx: *mut Context, p: *const f64, out: *mut f64) -> i32 {
    match catch_unwind(AssertUnwindSafe(|| {
        let p = slice::from_raw_parts(p, 3);
        let (pot, acc) = (*ctx).gravity.evaluate([p[0], p[1], p[2]]);
        let out = slice::from_raw_parts_mut(out, 4);
        out[0] = pot;
        out[1..].copy_from_slice(&acc);
    })) {
        Ok(()) => 0,
        Err(_) => -9,
    }
}
/// Evaluate the equations of motion.
///
/// # Safety
/// `ctx` must be live and exclusive, `n` must be 6 or 9, and the input/output
/// pointers must reference valid nonoverlapping arrays of `n` doubles.
#[no_mangle]
pub unsafe extern "C" fn engine_rhs(
    ctx: *mut Context,
    t: f64,
    y: *const f64,
    n: usize,
    out: *mut f64,
) -> i32 {
    match catch_unwind(AssertUnwindSafe(|| {
        if !matches!(n, 6 | 9) {
            return Err(-1);
        }
        let mut state = [0.; 9];
        state[..n].copy_from_slice(slice::from_raw_parts(y, n));
        let result = (*ctx).rhs(t, &state, n)?;
        slice::from_raw_parts_mut(out, n).copy_from_slice(&result[..n]);
        Ok(())
    })) {
        Ok(Ok(())) => 0,
        Ok(Err(e)) => e,
        Err(_) => -9,
    }
}
/// Evaluate environment: longitude/latitude (deg), height (km), density
/// (kg/m^3), height (km), and three inertial drag accelerations per unit beta.
///
/// # Safety
/// `ctx` must be live and exclusive; `y` has six readable doubles and `out`
/// has eight writable doubles, with no overlapping context memory.
#[no_mangle]
pub unsafe extern "C" fn engine_environment(
    ctx: *mut Context,
    t: f64,
    y: *const f64,
    out: *mut f64,
) -> i32 {
    match catch_unwind(AssertUnwindSafe(|| {
        let mut state = [0.; 9];
        state[..6].copy_from_slice(slice::from_raw_parts(y, 6));
        let (_, fixed, unit, rho, alt) = (*ctx).environment(t, &state)?;
        let a = slice::from_raw_parts_mut(out, 8);
        a[..3].copy_from_slice(&geodetic(fixed));
        a[3] = rho;
        a[4] = alt;
        a[5..].copy_from_slice(&unit);
        Ok(())
    })) {
        Ok(Ok(())) => 0,
        Ok(Err(e)) => e,
        Err(_) => -9,
    }
}
/// Propagate forward or backward to requested UTC epochs or a height event.
///
/// # Safety
/// `ctx` must be live and exclusive. `n` is 6 or 9; `initial` has `n` doubles;
/// `samples` has `nsample` doubles; `out_t` has `nsample+1` writable doubles;
/// `out_y` has `(nsample+1)*n` writable doubles; `stats` is writable.
/// Arrays must be aligned, non-null, and correctly sized, and outputs must
/// not overlap inputs or context memory. See README.md for units and ordering.
#[no_mangle]
pub unsafe extern "C" fn engine_propagate(
    ctx: *mut Context,
    start: f64,
    end: f64,
    initial: *const f64,
    n: usize,
    rtol: f64,
    max_step: f64,
    threshold: f64,
    samples: *const f64,
    nsample: usize,
    out_t: *mut f64,
    out_y: *mut f64,
    stats: *mut Stats,
) -> i32 {
    match catch_unwind(AssertUnwindSafe(|| {
        *stats = Stats {
            event: f64::NAN,
            ..Default::default()
        };
        integrate(
            &mut *ctx,
            start,
            end,
            slice::from_raw_parts(initial, n),
            rtol,
            max_step,
            threshold,
            slice::from_raw_parts(samples, nsample),
            slice::from_raw_parts_mut(out_t, nsample + 1),
            slice::from_raw_parts_mut(out_y, (nsample + 1) * n),
            &mut *stats,
        )
    })) {
        Ok(Ok(())) => 0,
        Ok(Err(e)) => e,
        Err(_) => -9,
    }
}

/// Attach native rmsis MSIS2.1 to a context. Returns zero on success.
///
/// # Safety
/// Context must be live and exclusive; paths must be valid UTF-8 NUL-terminated
/// strings. Library must be trusted and implement the rmsis adapter ABI.
#[cfg(feature = "msis21")]
#[no_mangle]
pub unsafe extern "C" fn engine_load_msis21(
    ctx: *mut Context,
    library: *const std::ffi::c_char,
    parameters: *const std::ffi::c_char,
) -> i32 {
    match catch_unwind(AssertUnwindSafe(|| {
        if ctx.is_null() || library.is_null() || parameters.is_null() {
            return Err(-1);
        }
        let lib = std::ffi::CStr::from_ptr(library).to_str().map_err(|_| -1)?;
        let par = std::ffi::CStr::from_ptr(parameters)
            .to_str()
            .map_err(|_| -1)?;
        (*ctx).load_msis21(lib, par).map_err(|_| -1)
    })) {
        Ok(Ok(())) => 0,
        Ok(Err(e)) => e,
        Err(_) => -9,
    }
}
