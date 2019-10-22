#ifndef COREFCTS_H
#define COREFCTS_H

#define ROUNDING MPFR_RNDN

void delta(mpfr_t, mpfr_t, mpfr_t);
void bg_method(mpfr_t, mpfr_t*);
void rm_method(mpfr_t, mpfr_t, mpfr_t*);
void rm_method_cosh(mpfr_t, mpfr_t, mpfr_t*);
void transform(mpfr_t, mpfr_t, mpfr_t*, mpfr_t*, double*, double*, double*);
void deltabar(mpfr_t, mpfr_t);

void set_params(double, double, int);
void set_time_parms(int, int, int);

#endif
