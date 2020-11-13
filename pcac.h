#ifndef PCAC_H
#define PCAC_H

#define ROUNDING MPFR_RNDN

void bootstrap_sample(mpfr_t*, mpfr_t*, mpfr_t*);
void full_sample(mpfr_t*, mpfr_t*, mpfr_t*);
void read_corr_AXIAL();
void read_corr_PP();
#endif
