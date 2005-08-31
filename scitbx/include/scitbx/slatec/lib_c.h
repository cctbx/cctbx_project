#ifndef SCITBX_SLATEC_LIB_H
#define SCITBX_SLATEC_LIB_H

#ifdef __cplusplus
extern "C" {
#endif

const char*
slatec_error(void);

void
slatec_clear_error(void);

double
slatec_dgamma(double x);

double
slatec_dlngam(double x);

double
slatec_dlnrel(double x);

double
slatec_dbinom(unsigned n, unsigned m);

#ifdef __cplusplus
}
#endif

#endif /* SCITBX_SLATEC_LIB_H */
