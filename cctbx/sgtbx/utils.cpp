/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cstddef>
#include <cctbx/sgtbx/utils.h>
#include <scitbx/array_family/misc_functions.h>

namespace cctbx { namespace sgtbx { namespace utils {

  int change_denominator(
    const int *old_num, int old_den,
          int *new_num, int new_den, int n)
  {
    for(std::size_t i=0;i<n;i++) {
          new_num[i] = old_num[i] * new_den;
      if (new_num[i] %  old_den) return -1;
          new_num[i] /= old_den;
    }
    return 0;
  }

  bool cmp_i_vec::operator()(const int *a, const int *b) const
  {
    using std::size_t;
    using scitbx::fn::absolute;
    size_t n0a = 0; for(size_t i=0;i<n_;i++) if (a[i] == 0) n0a++;
    size_t n0b = 0; for(size_t i=0;i<n_;i++) if (b[i] == 0) n0b++;
    if (n0a > n0b) return true;
    if (n0a < n0b) return false;
    for(size_t i=0;i<n_;i++) {
      if (a[i] != 0 && b[i] == 0) return true;
      if (a[i] == 0 && b[i] != 0) return false;
    }
    for(size_t i=0;i<n_;i++) {
      if (absolute(a[i]) < absolute(b[i])) return true;
      if (absolute(a[i]) > absolute(b[i])) return false;
    }
    for(size_t i=0;i<n_;i++) {
      if (a[i] > b[i]) return true;
      if (a[i] < b[i]) return false;
    }
    return false;
  }

}}} // namespace cctbx::sgtbx::utils
