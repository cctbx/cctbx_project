/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored parts of sgtbx/matrix.cpp (rwgk)
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <boost/rational.hpp>
#include <cctbx/sgtbx/tr_vec.h>
#include <cctbx/sgtbx/utils.h>

namespace cctbx { namespace sgtbx {

  void throw_unsuitable_tr_vec(const char* file, long line)
  {
    throw error_rational_vector(file, line,
      "Unsuitable value for rational translation vector.");
  }

  tr_vec tr_vec::new_denominator(int new_den) const
  {
    tr_vec result(new_den);
    if (utils::change_denominator(num_.begin(), den(),
                                  result.num_.begin(), new_den,
                                  num_.size()) != 0) {
      throw_unsuitable_tr_vec(__FILE__, __LINE__);
    }
    return result;
  }

  tr_vec operator/(tr_vec const& lhs, int rhs)
  {
    sg_vec3 new_num;
    for(std::size_t i=0;i<3;i++) {
      if (lhs.num_[i] % rhs) throw_unsuitable_tr_vec(__FILE__, __LINE__);
      new_num[i] = lhs.num_[i] / rhs;
    }
    return tr_vec(new_num, lhs.den_);
  }

  tr_vec tr_vec::cancel() const
  {
    int g = den();
    for(std::size_t i=0;i<3;i++) g = boost::gcd(g, num_[i]);
    if (g == 0) return *this;
    return tr_vec(num_ / g, den() / g);
  }

  tr_vec tr_vec::plus(tr_vec const& rhs) const
  {
    tr_vec result(boost::lcm(den(), rhs.den()));
    int l = result.den() / den();
    int r = result.den() / rhs.den();
    for(std::size_t i=0;i<3;i++) result[i] = num_[i] * l + rhs[i] * r;
    return result.cancel();
  }

  tr_vec tr_vec::minus(tr_vec const& rhs) const
  {
    tr_vec result(boost::lcm(den(), rhs.den()));
    int l = result.den() / den();
    int r = result.den() / rhs.den();
    for(std::size_t i=0;i<3;i++) result[i] = num_[i] * l - rhs[i] * r;
    return result.cancel();
  }

}} // namespace cctbx::sgtbx
