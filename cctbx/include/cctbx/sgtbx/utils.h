/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_UTILS_H
#define CCTBX_SGTBX_UTILS_H

namespace cctbx { namespace sgtbx { namespace utils {

  int change_denominator(
    const int *old_num, int old_den,
          int *new_num, int new_den, int n);

  class cmp_i_vec
  {
    public:
      cmp_i_vec(std::size_t n) : n_(n) {}

      bool operator()(const int *a, const int *b) const;

    private:
      std::size_t n_;
  };

}}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_UTILS_H
