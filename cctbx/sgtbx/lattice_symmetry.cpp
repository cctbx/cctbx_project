#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/sgtbx/rot_mx_info.h>
#include <scitbx/array_family/sort.h>

namespace cctbx { namespace sgtbx { namespace lattice_symmetry {

  std::size_t
  group_search::n_potential_axes()
  {
    if (potential_axes_.size() == 0) compute_potential_axes();
    return potential_axes_.size();
  }

  space_group
  group_search::operator()(
    uctbx::unit_cell const& niggli_cell,
    double max_delta)
  {
    if (potential_axes_.size() == 0) compute_potential_axes();
    max_delta *= scitbx::constants::pi_180;
    uc_mat3 const& frac = niggli_cell.fractionalization_matrix();
    uc_mat3 const& orth = niggli_cell.orthogonalization_matrix();
    std::vector<double> deltas;
    std::vector<uc_vec3> ts;
    for(potential_axis_t* axis=&*potential_axes_.begin();
                          axis!=&*potential_axes_.end();axis++) {
      uc_vec3 t = orth * axis->u;
      uc_vec3 tau = axis->h * frac;
      double delta = scitbx::fn::absolute(
        std::atan2(t.cross(tau).length(), axis->abs_uh));
      if (delta < max_delta) {
        deltas.push_back(delta);
        ts.push_back(t);
      }
    }
    af::shared<std::size_t> perm = af::sort_permutation(
      af::const_ref<double>(&*deltas.begin(), deltas.size()));
    af::const_ref<std::size_t> perm_cr = perm.const_ref();
    space_group group;
    for(std::size_t i=0;i<perm_cr.size();i++) {
      uc_mat3 w_cart = two_fold_matrix_from_axis_direction(ts[i]);
      sg_mat3 w_frac = as_integer_matrix(frac*w_cart*orth);
      rt_mx s(rot_mx(w_frac,1));
      space_group expanded_group(group);
      try {
        expanded_group.expand_smx(s);
      }
      catch (error_non_crystallographic_rotation_matrix_encountered const&) {
        break;
      }
      group = expanded_group;
    }
    return group;
  }

  void
  group_search::compute_potential_axes()
  {
    if (potential_axes_.size() != 0) return;
    for(int u=0;u<=modulus_;u++)
    for(int v=-modulus_;v<=modulus_;v++)
    for(int w=-modulus_;w<=modulus_;w++)
    for(int h=0;h<=modulus_;h++)
    for(int k=-modulus_;k<=modulus_;k++)
    for(int l=-modulus_;l<=modulus_;l++) {
      int abs_uh = scitbx::fn::absolute(u*h+v*k+w*l);
      if (abs_uh == 1 || abs_uh == 2) {
        potential_axes_.push_back(potential_axis_t(
          sg_vec3(u,v,w),
          sg_vec3(h,k,l),
          abs_uh));
      }
    }
  }

  double
  find_max_delta(
    uctbx::unit_cell const& niggli_cell,
    space_group const& group,
    int modulus)
  {
    CCTBX_ASSERT(group.n_ltr() == 1);
    CCTBX_ASSERT(group.f_inv() == 1);
    uc_mat3 const& frac = niggli_cell.fractionalization_matrix();
    uc_mat3 const& orth = niggli_cell.orthogonalization_matrix();
    double result = 0;
    for(int i_smx=1;i_smx<group.n_smx();i_smx++) {
      rot_mx_info r_info = group.smx()[i_smx].r().info();
      if (r_info.type() != 2) continue;
      int u = r_info.ev()[0];
      int v = r_info.ev()[1];
      int w = r_info.ev()[2];
      uc_vec3 t = orth * r_info.ev();
      double min_delta = -1;
      for(int h=0;h<=modulus;h++)
      for(int k=-modulus;k<=modulus;k++)
      for(int l=-modulus;l<=modulus;l++) {
        int abs_uh = scitbx::fn::absolute(u*h+v*k+w*l);
        if (abs_uh == 1 || abs_uh == 2) {
          uc_vec3 tau = uc_vec3(h,k,l) * frac;
          double delta = scitbx::fn::absolute(
            std::atan2(t.cross(tau).length(), abs_uh));
          if (min_delta == -1 || min_delta > delta) {
            min_delta = delta;
          }
        }
      }
      CCTBX_ASSERT(min_delta != -1);
      math::update_max(result, min_delta);
    }
    return result/scitbx::constants::pi_180;
  }

}}} // namespace cctbx::sgtbx::lattice_symmetry
