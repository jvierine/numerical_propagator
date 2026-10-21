# numerical_propagator

NASA GMAT compatible propagator written in Rust, supporting gravity harmonics and atmospheric drag.

A standalone numerical orbit propagator with adaptive DOP853 integration,
fully normalized spherical-harmonic Earth gravity (including EGM2008),
Sun/Moon third-body gravity, rotating-atmosphere drag, and MSIS density
callbacks. It supports forward and backward propagation, dense output at
requested times, and a geodetic-altitude stopping event.

GMAT compatibility means numerical agreement for matched force models, frames,
epochs, and initial states. This is an independent implementation, not a GMAT
script interpreter or a NASA-maintained product. Matching an integrator does
not establish atmospheric-model or orbit-observation accuracy.

## Model and uncertainty

The equation of motion is

$
\dot{\mathbf r}=\mathbf v,\qquad
\dot{\mathbf v}=\mathbf g(\mathbf r,t)
-\frac{1}{2}\rho_{\mathrm{MSIS}}\,\beta_{\mathrm{eff}}(t)
\lvert\mathbf v_{\mathrm{rel}}\rvert\mathbf v_{\mathrm{rel}},
\qquad \beta_{\mathrm{eff}}=\gamma C_D A/m.
$

Gravity includes EGM2008 to degree and order 180, the Sun and Moon in the
validated re-entry configuration. The gravity degree and order are configurable.
MSIS2.1 uses the observed solar and geomagnetic drivers; the atmosphere
co-rotates with Earth. The relative velocity is measured against this
co-rotating atmosphere.

Here, \(C_D\) is the drag coefficient, \(A/m\) is the area-to-mass ratio,
and \(\gamma\) scales the neutral density. The fitted effective coefficient
\(\beta_{\mathrm{eff}}\) combines these terms, so density error absorbed in it
must not be applied again through an independent density correction. The core
propagates supplied states and force parameters; estimating them and their
uncertainty requires an external fit and sampling procedure. The caller supplies
the gravity coefficients, ephemerides, Earth orientation, weather records, and
MSIS callback as described below.

## Build and run

Only the Rust toolchain is required. The core has no external crate dependencies.

```sh
git clone https://github.com/jvierine/numerical_propagator.git
cd numerical_propagator
cargo test --all-targets
cargo build --release
cargo run --release --example circular
```

The circular-orbit example compares a complete revolution with the analytic
two-body solution. Tests run the same check forward and backward, requiring
position error below 1 cm and velocity error below 0.01 mm/s.

The library builds as both a Rust `rlib` and a C-compatible shared library.
Use `cargo doc --open` for Rust API documentation. C declarations are in
[numerical_propagator.h](numerical_propagator.h).

## Inputs and conventions

The caller supplies the environment. Gravity coefficients, Earth orientation,
Sun/Moon ephemerides, space-weather records, and an MSIS implementation are
not bundled or downloaded. No Python runtime is required.

- State: geocentric GCRS Cartesian position in metres and velocity in m/s.
  Use six state components for ordinary trajectory propagation. The optional
  nine-component mode additionally integrates three internal work diagnostics;
  it is not a state transition matrix or covariance.
- Time: UTC Unix seconds throughout, with UT1−UTC supplied separately in `dut`.
  This is not GMAT's internal time coordinate: convert epochs when comparing.
- Gravity: `gm` in m^3/s^2, reference `radius` in metres. `coefficients_c` and
  `coefficients_s` each contain `(degree+1)^2` doubles, row-major `[n,m]`.
  Supply fully normalized coefficients. The central term is computed
  separately; degree zero gives point-mass gravity.
- Frames: `q` and `polar` each contain `nframe` row-major 3×3 matrices; `dut`
  has `nframe` seconds. The Earth-fixed transform is `polar * R_ERA * q`,
  with the GCRS-to-intermediate matrix `q`, Earth rotation angle from UT1, and
  polar motion. Supply sufficiently dense, consistent external frame tables;
  matrices are interpolated linearly. At least two samples are required.
- Ephemerides: Sun/Moon geocentric GCRS positions and velocities have
  `3*nbody` doubles each (metres and m/s); cubic Hermite interpolation is used.
  Supply at least two samples when third-body gravity is enabled.
- Weather: `f107`, `f107a`, and `doy` have `nweather` values; `ap` has
  `7*nweather` values in the MSIS historical-Ap convention. Flux is in sfu;
  `doy` is day of year. Drivers are piecewise constant. Integrator steps stop
  at weather-bin boundaries.
- Uniform table origins and spacings are the corresponding `*_start` and
  `*_step` fields (seconds). Tables must cover the full propagation, including
  endpoint integration stages. Frame/body/correction interpolation requires
  an additional sample after the final evaluation; do not end coverage exactly
  at the propagation endpoint.
- `beta` is effective drag coefficient `gamma * C_D * A/m` in m^2/kg. Drag
  uses a rigidly corotating atmosphere; neutral winds are not modeled.
  `engine_beta` sets a constant value. Rust callers can use
  `Context::set_drag_nodes` for positive log-linear time-dependent values.
- Optional log-density corrections use `ncorrection` values on a uniform
  time grid; zero disables them. Height correction origins/scales use km.
  When fitting an effective beta that already absorbs density error, avoid
  counting that error again with a separate density correction.

The minimal example supplies placeholder tables for unused forces. Input arrays
are copied by `engine_new`; their original storage can then be released. The
context and callbacks must remain live and must not be accessed concurrently.
The raw ABI is unsafe: invalid pointers or lengths can cause undefined behavior;
return codes are not a replacement for meeting the buffer contracts.

## MSIS connection

The core accepts a callback rather than linking a specific MSIS distribution.
For a native Rust adapter, set `drag=1`, construct the context with
`msis_address=0`, and install `Context::set_density_callback` before propagation.
For example, an adapter can call [rmsis](https://github.com/jvierine/rmsis).

The double-precision callback arguments are day of year, UTC seconds of day,
longitude (degrees east), latitude (degrees), geodetic height **in kilometres**,
F10.7, mean F10.7, and seven Ap values. Return mass density in kg/m^3.

The legacy `msis_address` accepts the single-precision, pointer-based pymsis
wrapper ABI used by the original application. It is a live function address,
not a model identifier or path. Keep the loaded library alive. Serialize calls
if the linked Fortran implementation uses shared global state. The callback
must initialize/configure its model, units, switches, and weather convention
consistently with the comparison run.

## Propagation and output

`engine_propagate` takes `start`, `end`, initial state, relative tolerance,
maximum step **in seconds**, and stopping height **in kilometres**. Supply
strictly ordered sample times in the direction of propagation, inside the
start/end interval. Output capacity is `nsample+1` rows to allow an additional
height-event state; `Stats::nout` is the number of valid rows. The endpoint is
only guaranteed in output when explicitly requested as a sample, or when a
height event terminates the trajectory. `Stats::event` is NaN when no event
occurs. Always check the return code, output coverage, and event status.

Return codes: `0` success; `-1` invalid arguments or missing callback; `-2`
outside table coverage; `-3` invalid/nonfinite model result; `-4` initial state
at/below the forward stopping height; `-5` step underflow; `-7` optional RHS
budget exceeded; `-9` caught Rust panic. These do not catch memory misuse.

## GMAT validation

The existing event regression used GMAT's own EGM2008 gravity, Sun/Moon forces,
and RK89 integrator. A C++ force plugin supplied drag using the same native
MSIS context; therefore the atmospheric force is shared in that comparison.
A separate matched-force control shared the full acceleration evaluation
while retaining independent integrators.

| Archived case | Maximum position difference | Maximum velocity difference |
| --- | ---: | ---: |
| F1, GMAT gravity and third bodies | 0.082684 m | 0.00156769 m/s |
| F2, GMAT gravity and third bodies | 0.096365 m | 0.00177886 m/s |
| F1, tighter GMAT step/tolerance | 0.008151 m | 0.00012110 m/s |
| F1, matched-force control | 0.034852 m | 0.00064255 m/s |

The archived summary is in [validation/gmat_reference.log](validation/gmat_reference.log).
These are historical results from the original application, not a GMAT run
performed by this standalone repository's CI. The event fitting application,
observation archive, GMAT installation, and event-reproduction harness are
outside this core-only repository. CI checks the analytic orbit and drag
interpolation; it does not reproduce the archived GMAT experiment.

## Third-party notices

The DOP853 coefficient table is derived from SciPy and retains its BSD notice
in `SCIPY_LICENSE.txt`. `PYMSIS_LICENSE.txt` retains the notice associated with
the legacy callback integration. MSIS itself and gravity/ephemeris datasets
are separate dependencies with their own terms. No project-wide license has
been selected for the remaining source.
