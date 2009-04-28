#ifndef CCTBX_SGTBX_ROW_ECHELON_SOLVE_H
#define CCTBX_SGTBX_ROW_ECHELON_SOLVE_H

#include <cctbx/sgtbx/basic.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/matrix/row_echelon.h>

namespace cctbx { namespace sgtbx { namespace row_echelon { namespace solve {

  af::tiny<sg_vec3, 4>
  homog_rank_1(af::const_ref<int, af::mat_grid> const& re_mx,
               scitbx::matrix::row_echelon::independent<int> const& indep);

  af::tiny<sg_vec3, 4>
  homog_rank_1(af::const_ref<int, af::mat_grid> const& re_mx);

  sg_vec3
  homog_rank_2(af::const_ref<int, af::mat_grid> const& re_mx);

}}}} // namespace cctbx::sgtbx::row_echelon::solve

#endif // CCTBX_SGTBX_ROW_ECHELON_SOLVE_H
