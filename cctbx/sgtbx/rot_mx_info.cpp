/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored parts of cctbx/properties.cpp (rwgk)
     2001 Sep: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/rot_mx_info.h>
#include <cctbx/sgtbx/row_echelon_solve.h>

namespace cctbx { namespace sgtbx {

  namespace {

    int sense_of_rotation(rot_mx const& r, int type, sg_vec3 const& ev)
    {
      // M.B. Boisen, Jr. & G.V. Gibbs
      // Mathematical Crystallography, Revised Edition 1990
      // pp. 348-349, 354-356

      int f = 1; if (type < 0) f = -1;
      int trace = f * r.num().trace();
      if (trace == 3 || trace == -1) return 0; /* 1-fold or 2-fold */
      if (ev[1] == 0 && ev[2] == 0) {
        if (ev[0] * f * r[7] > 0)
          return 1;
      }
      else {
        if (f * (r[3] * ev[2] - r[6] * ev[1]) > 0)
          return 1;
      }
      return -1;
    }

  } // namespace <anonymous>

  rot_mx_info::rot_mx_info(rot_mx const& r)
  : type_(r.type()), ev_(0,0,0), sense_(0)
  {
    CCTBX_ASSERT(r.den() == 1);
    if (type_ == 0) {
      throw error("Cannot determine type of rotation matrix.");
    }
    rot_mx proper_r = r;
    int proper_order = type_;
    if (proper_order < 0) {
      proper_order *= -1;
      proper_r = -proper_r;
    }
    if (proper_order > 1) {
      rot_mx rmi = proper_r.minus_unit_mx();
      scitbx::mat_ref<int> re_mx(rmi.num().begin(), 3, 3);
      if (row_echelon::form(re_mx) != 2) {
        throw error("Cannot determine Eigenvector of rotation matrix.");
      }
      ev_ = row_echelon::solve::homog_rank_2(re_mx);
      sense_ = sense_of_rotation(r, type_, ev_);
    }
  }

  rot_mx_info rot_mx::info() const
  {
    return rot_mx_info(*this);
  }

}} // namespace cctbx::sgtbx
