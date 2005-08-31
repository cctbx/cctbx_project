#ifndef SCITBX_SLATEC_LIB_CPP_H
#define SCITBX_SLATEC_LIB_CPP_H

#include <scitbx/slatec/lib_c.h>
#include <scitbx/error.h>

namespace scitbx { namespace slatec {

  inline
  double
  dgamma(double x)
  {
    SCITBX_ASSERT(!slatec_error());
    double result = slatec_dgamma(x);
    if (slatec_error()) {
      std::string msg = slatec_error();
      slatec_clear_error();
      throw std::runtime_error(msg);
    }
    return result;
  }

  inline
  double
  dlngam(double x)
  {
    SCITBX_ASSERT(!slatec_error());
    double result = slatec_dlngam(x);
    if (slatec_error()) {
      std::string msg = slatec_error();
      slatec_clear_error();
      throw std::runtime_error(msg);
    }
    return result;
  }

  inline
  double
  dlnrel(double x)
  {
    SCITBX_ASSERT(!slatec_error());
    double result = slatec_dlnrel(x);
    if (slatec_error()) {
      std::string msg = slatec_error();
      slatec_clear_error();
      throw std::runtime_error(msg);
    }
    return result;
  }

  inline
  double
  dbinom(unsigned n, unsigned m)
  {
    SCITBX_ASSERT(!slatec_error());
    double result = slatec_dbinom(n, m);
    if (slatec_error()) {
      std::string msg = slatec_error();
      slatec_clear_error();
      throw std::runtime_error(msg);
    }
    return result;
  }

}} // namespace scitbx::slatec

#endif /* SCITBX_SLATEC_LIB_CPP_H */
