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

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <boost/rational.hpp>
#include <cctbx/sgtbx/basic.h>

namespace cctbx { namespace sgtbx {

  inline int modPositive(int ix, int iy) {
    if (iy > 0) {
      ix %= iy;
      if (ix < 0) ix += iy;
    }
    return ix;
  }
  inline int modShort(int ix, int iy) {
        ix = modPositive(ix, iy);
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
    rational(int n = 0) : boost::rational<int>(n) {}
    rational(int n, int d) : boost::rational<int>(n, d) {}
    rational(const boost::rational<int>& r) : boost::rational<int>(r) {}
    std::string format(bool Decimal = 0) const;
  };

  int rationalize(double fVal, int& iVal, int BF);

  inline int norm_denominator(int numerator, int denominator) {
    return denominator / gcd(numerator, denominator);
  }

  int iRowEchelonFormT(int *M, int mr, int mc, int *T, int tc);
  int iREBacksubst(const int *M, const int *V,
                   const int nr, const int nc,
                   int *Sol, int *FlagIndep);
  int iRESetIxIndep(const int *REMx, int nr, int nc, int *IxIndep, int mIndep);
  void SolveHomRE1(const int REMx[3], const int IxIndep[2], Vec3 Sol[4]);

  int SmithNormalForm(int *M, int mr, int mc, int *P, int *Q);

}} // namespace cctbx::sgtbx

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // CCTBX_SGTBX_MATH_H
