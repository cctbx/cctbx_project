/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#ifndef CCTBX_SGTBX_ROW_ECHELON_SOLVE_H
#define CCTBX_SGTBX_ROW_ECHELON_SOLVE_H

#include <cctbx/sgtbx/row_echelon.h>

namespace cctbx { namespace sgtbx { namespace row_echelon { namespace solve {

  af::tiny<sg_vec3, 4>
  homog_rank_1(scitbx::mat_const_ref<int> const& re_mx,
               row_echelon::independent<int> const& indep);

  af::tiny<sg_vec3, 4>
  homog_rank_1(scitbx::mat_const_ref<int> const& re_mx);

  sg_vec3
  homog_rank_2(scitbx::mat_const_ref<int> const& re_mx);

}}}} // namespace cctbx::sgtbx::row_echelon::solve

#endif // CCTBX_SGTBX_ROW_ECHELON_SOLVE_H
