/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created with fragments from several files (rwgk)
 */

#include <cctbx/sgtbx/row_echelon_solve.h>

namespace cctbx { namespace sgtbx { namespace row_echelon { namespace solve {

  af::tiny<sg_vec3, 4>
  homog_rank_1(scitbx::mat_const_ref<int> const& re_mx,
               row_echelon::independent<int> const& indep)
  {
    CCTBX_ASSERT(re_mx.n_rows() == 1);
    CCTBX_ASSERT(indep.indices.size() == 2);
    const int trial_v[4][2] =
      {{ 1,  0 },
       { 0,  1 },
       { 1,  1 },
       { 1, -1 },
      };
    af::tiny<sg_vec3, 4> sol;
    for(std::size_t i_tv=0;i_tv<4;i_tv++) {
      sol[i_tv].fill(0);
      for(std::size_t i=0;i<2;i++) {
        sol[i_tv][indep.indices[i]] = trial_v[i_tv][i];
      }
      int* n_a = 0;
      CCTBX_ASSERT(
        row_echelon::back_substitution(re_mx, n_a, sol[i_tv].begin()) > 0);
    }
    return sol;
  }

  af::tiny<sg_vec3, 4>
  homog_rank_1(scitbx::mat_const_ref<int> const& re_mx)
  {
    row_echelon::independent<int> indep(re_mx);
    return homog_rank_1(re_mx, indep);
  }

  namespace {

    int sign_hemisphere(sg_vec3 const& v)
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

  } // namespace <anonymous>

  sg_vec3
  homog_rank_2(scitbx::mat_const_ref<int> const& re_mx)
  {
    CCTBX_ASSERT(re_mx.n_rows() == 2);
    row_echelon::independent<int> indep(re_mx);
    CCTBX_ASSERT(indep.indices.size() == 1);
    sg_vec3 ev(0,0,0);
    ev[indep.indices[0]] = 1;
    int* n_a = 0;
    CCTBX_ASSERT(
      row_echelon::back_substitution(re_mx, n_a, ev.begin()) >= 1);
    if (sign_hemisphere(ev) < 0) ev *= -1;
    return ev;
  }

}}}} // namespace cctbx::sgtbx::row_echelon::solve
