// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <ctype.h> // cannot use cctype b/o non-conforming compilers
#include <stdio.h>
#include <cctbx/fixes/cstdlib>
#include <cctbx/array_family/misc_functions.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/sgtbx/math.h>

namespace cctbx { namespace sgtbx {

  int ChangeBaseFactor(const int *Old, int OldBF, int *New, int NewBF, int n)
  {
    for(std::size_t i=0;i<n;i++) {
          New[i] = Old[i] * NewBF;
      if (New[i] %  OldBF) return -1;
          New[i] /= OldBF;
    }
    return 0;
  }

  int rationalize(double fVal, int& iVal, int BF)
  {
    if (BF == 0) return -1;
        fVal *= BF;
    if (fVal < 0.) iVal = int(fVal - .5);
    else           iVal = int(fVal + .5);
        fVal -= iVal;
        fVal /= BF;
    if (fVal < 0.) fVal = -fVal;
    if (fVal > .0001) return -1;
    return 0;
  }

  std::string rational::format(bool Decimal) const
  {
    if (numerator() == 0) return std::string("0");
    char buf[32];
    if (Decimal) {
      sprintf(buf, "%.6g", double(numerator()) / denominator());
      char* cp = buf;
      if (*cp == '-') cp++;
      if (*cp == '0') {
        char* cpp = cp + 1; while (*cp) *cp++ = *cpp++;
      }
    }
    else if (denominator() == 1) {
      sprintf(buf, "%d", numerator());
    }
    else {
      sprintf(buf, "%d/%d", numerator(), denominator());
    }
    return std::string(buf);
  }

  int SignHemisphere(const af::int3& v)
  {
    if (v[2] >  0) return  1;
    if (v[2] == 0) {
      if (v[1] >  0) return  1;
      if (v[1] == 0) {
        if (v[0] >  0) return  1;
        if (v[0] == 0)
          return 0;
      }
    }
    return -1;
  }

  bool CmpiVect::operator()(const int *a, const int *b) const {
    int i;
    int n0a = 0; for(i=0;i<m_n;i++) if (a[i] == 0) n0a++;
    int n0b = 0; for(i=0;i<m_n;i++) if (b[i] == 0) n0b++;
    if (n0a > n0b) return true;
    if (n0a < n0b) return false;
    for(i=0;i<m_n;i++) {
      if (a[i] != 0 && b[i] == 0) return true;
      if (a[i] == 0 && b[i] != 0) return false;
    }
    for(i=0;i<m_n;i++) {
      if (fn::absolute(a[i]) < fn::absolute(b[i])) return true;
      if (fn::absolute(a[i]) > fn::absolute(b[i])) return false;
    }
    for(i=0;i<m_n;i++) {
      if (a[i] > b[i]) return true;
      if (a[i] < b[i]) return false;
    }
    return false;
  }

}} // namespace cctbx::sgtbx
