// $Id$

#ifndef CCTBX_SGTBX_MATH_H
#define CCTBX_SGTBX_MATH_H

#include <boost/rational.hpp>

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

  using boost::gcd;
  using boost::lcm;
  using boost::rational_cast;
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

} // namespace sgtbx

#endif // CCTBX_SGTBX_MATH_H
