#ifndef SCITBX_MATH_ABS_H
#define SCITBX_MATH_ABS_H

#include <complex>

namespace scitbx { namespace math {

  template <class T>
  struct abs_traits {
    typedef T result_type;
  };

  template <class T>
  struct abs_traits< std::complex<T> > {
    typedef T result_type;
  };

}}

#endif
