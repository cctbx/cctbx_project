/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/sgtbx/miller.h (rwgk)
 */

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
    return math::mod_positive(h * t, t.den());
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_MILLER_OPS_H
