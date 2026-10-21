/* C ABI for numerical_propagator. Buffer contracts and units: README.md. */
#ifndef NUMERICAL_PROPAGATOR_H
#define NUMERICAL_PROPAGATOR_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct Context Context;
typedef struct Inputs {
    size_t degree;
    double gm;
    double radius;
    const double * coefficients_c;
    const double * coefficients_s;
    double frame_start;
    double frame_step;
    size_t nframe;
    const double * q;
    const double * polar;
    const double * dut;
    double body_start;
    double body_step;
    size_t nbody;
    const double * sun_pos;
    const double * sun_vel;
    const double * moon_pos;
    const double * moon_vel;
    double weather_start;
    double weather_step;
    size_t nweather;
    const double * f107;
    const double * f107a;
    const double * ap;
    const double * doy;
    double correction_start;
    double correction_step;
    size_t ncorrection;
    const double * log_correction;
    const double * correction_height_linear;
    const double * correction_height_quadratic;
    double correction_height_origin;
    double correction_height_scale;
    double correction_clip;
    double beta;
    int32_t drag;
    int32_t third_bodies;
} Inputs;
typedef struct Stats {
    size_t nout, nfev, accepted, rejected;
    double event;
} Stats;
Context *engine_new(const Inputs *input);
void engine_free(Context *ctx);
/* Available in default builds (feature msis21). */
int32_t engine_load_msis21(Context *ctx, const char *library, const char *parameters);
void engine_beta(Context *ctx, double beta);
int32_t engine_gravity(Context *ctx, const double *position, double *output);
int32_t engine_rhs(Context *ctx, double time, const double *state, size_t n, double *output);
int32_t engine_environment(Context *ctx, double time, const double *state, double *output);
int32_t engine_propagate(Context *ctx, double start, double end, const double *initial,
    size_t n, double rtol, double max_step, double threshold, const double *samples,
    size_t nsample, double *out_t, double *out_y, Stats *stats);
#ifdef __cplusplus
}
#endif
#endif
