// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/sgtbx/miller.h (rwgk)
 */

#ifndef CCTBX_SGTBX_MILLER_OPS_H
#define CCTBX_SGTBX_MILLER_OPS_H

#include <cctbx/sgtbx/matrix.h>

namespace cctbx { namespace sgtbx {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  inline miller::Index
  operator*(const miller::Index& lhs, const RotMx& rhs)
  {
    return miller::Index(
      lhs[0] * rhs[0] + lhs[1] * rhs[3] + lhs[2] * rhs[6],
      lhs[0] * rhs[1] + lhs[1] * rhs[4] + lhs[2] * rhs[7],
      lhs[0] * rhs[2] + lhs[1] * rhs[5] + lhs[2] * rhs[8]);
  }

  inline int operator*(const miller::Index& lhs, const TrVec& rhs)
  {
    int result = 0;
    for(int i=0;i<3;i++) result += lhs[i] * rhs[i];
    return result;
  }

  inline int HT_mod_1(const miller::Index& H, const TrVec& T)
  {
    return modPositive(H * T, T.BF());
  }

#endif // DOXYGEN_SHOULD_SKIP_THIS

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_MILLER_OPS_H
