//! backward compatibility 2008-11-12 (whole file)
/*! This header has been deprecated:
    Please #include <scitbx/math/modulo.h> instead
    and don't forget to replace math::mod_positive or
    cctbx::math::mod_positive with scitbx::math::mod_positive
    and similarly for other functions.
 */

#ifndef CCTBX_MATH_MOD_H
#define CCTBX_MATH_MOD_H

#include <scitbx/math/modulo.h>

namespace cctbx { namespace math {

  using scitbx::math::mod_positive;
  using scitbx::math::mod_short;
  using scitbx::math::fmod_short;

}} // namespace cctbx::math

#endif // CCTBX_MATH_MOD_H
