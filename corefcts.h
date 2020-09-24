#ifndef COREFCTS_H
#define COREFCTS_H

#define ROUNDING MPFR_RNDN

void delta(mpfr_t, mpfr_t, mpfr_t);
void bg_method(mpfr_t, mpfr_t*);
void rm_method(mpfr_t, mpfr_t, mpfr_t*);
void rm_method_cosh(mpfr_t, mpfr_t, mpfr_t*);
void rm_method_cosh_lvl1(mpfr_t, mpfr_t, mpfr_t*);
void transform_only_step1(mpfr_t*, double*);
void transform(mpfr_t, mpfr_t, mpfr_t*, mpfr_t*, double*, double*, double*);
void transform_lvl1(mpfr_t, mpfr_t, mpfr_t*, mpfr_t*, double*, double*, double*);
void deltabar(mpfr_t, mpfr_t);
void deltabar_lvl1(mpfr_t, mpfr_t);
void reconstructed_error(mpfr_t, mpfr_t);

void set_params(double, double, int, double);
void set_time_parms(int, int, int);

#endif

