#ifndef CCTBX_SGTBX_MILLER_OPS_H
#define CCTBX_SGTBX_MILLER_OPS_H

#include <cctbx/miller.h>
#include <cctbx/sgtbx/rot_mx.h>

namespace cctbx { namespace sgtbx {

  template <typename NumType>
  inline
  miller::index<NumType>
  operator*(miller::index<NumType> const& h, rot_mx const& r)
  {
    return miller::index<NumType>(h * r.num());
  }

  template <typename NumType>
  inline
  NumType
  operator*(miller::index<NumType> const& h, tr_vec const& t)
  {
    return h * t.num();
  }

  inline
  int
  ht_mod_1(miller::index<> const& h, tr_vec const& t)
  {
    return scitbx::math::mod_positive(h * t, t.den());
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_MILLER_OPS_H
