// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_MATH_H
#define CCTBX_SGTBX_MATH_H

#include <boost/rational.hpp>
#include <cctbx/sgtbx/basic.h>

namespace sgtbx {

  inline int iModPositive(int ix, int iy) {
    if (iy > 0) {
      ix %= iy;
      if (ix < 0) ix += iy;
    }
    return ix;
  }
  inline int iModShort(int ix, int iy) {
        ix = iModPositive(ix, iy);
    if (ix > iy / 2)
        ix -= iy;
    return ix;
  }

  inline Vec3 CrossProduct(const Vec3& r, const Vec3& s) {
    Vec3 result;
    result[0] = r[1] * s[2] - r[2] * s[1];
    result[1] = r[2] * s[0] - r[0] * s[2];
    result[2] = r[0] * s[1] - r[1] * s[0];
    return result;
  }

  using boost::gcd;
  using boost::lcm;
  struct rational : boost::rational<int> {
    rational(int n = 0) : boost::rational<int>(n) { }
    rational(int n, int d) : boost::rational<int>(n, d) { }
    std::string format(bool Decimal = 0) const;
  };

  int rationalize(double fVal, int& iVal, int BF);

  int iRowEchelonFormT(int *M, int mr, int mc, int *T, int tc);
  int iREBacksubst(const int *M, const int *V,
                   const int nr, const int nc,
                   int *Sol, int *FlagIndep);
  int iRESetIxIndep(const int *REMx, int nr, int nc, int *IxIndep, int mIndep);

  int SmithNormalForm(int *M, int mr, int mc, int *P, int *Q);

} // namespace sgtbx

#endif // CCTBX_SGTBX_MATH_H
