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
    double max_delta,
    bool const& only_test_generators,
    bool const& introspection)
  {
    if (potential_axes_.size() == 0) compute_potential_axes();
    max_delta *= scitbx::constants::pi_180;
    uc_mat3 const& frac = niggli_cell.fractionalization_matrix();
    uc_mat3 const& orth = niggli_cell.orthogonalization_matrix();
    std::vector<double> deltas;
    std::vector<uc_vec3> ts;
    candidates = scitbx::af::shared<evaluated_axis_t>();
    for(potential_axis_t* axis=&*potential_axes_.begin();
                          axis!=&*potential_axes_.end();axis++) {
      uc_vec3 t = orth * axis->u;
      uc_vec3 tau = axis->h * frac;
      double delta = scitbx::fn::absolute(
        std::atan2(t.cross(tau).length(), axis->abs_uh));
      if (delta < max_delta) {
        deltas.push_back(delta);
        ts.push_back(t);
        if (introspection) {
          candidates.push_back(evaluated_axis_t(*axis,t,tau,delta));
          }
      }
    }
    af::shared<std::size_t> perm = af::sort_permutation(
      af::const_ref<double>(&*deltas.begin(), deltas.size()));
    af::const_ref<std::size_t> perm_cr = perm.const_ref();
    space_group group;
    for(std::size_t i=0;i<perm_cr.size();i++) {
      uc_mat3 w_cart = two_fold_matrix_from_axis_direction(ts[perm_cr[i]]);
      sg_mat3 w_frac = as_integer_matrix(frac*w_cart*orth);
      rt_mx s(rot_mx(w_frac,1));
      space_group expanded_group(group);
      try {
        expanded_group.expand_smx(s);
      }
      catch (error_non_crystallographic_rotation_matrix_encountered const&) {
        break;
      }
      if (!only_test_generators && scitbx::constants::pi_180*
          find_max_delta(niggli_cell,expanded_group,2) >= max_delta) {
          continue;
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
    double min_cos_delta = 1;
    for(int i_smx=1;i_smx<group.n_smx();i_smx++) {
      rot_mx const& r = group.smx()[i_smx].r();
      rot_mx_info r_info = r.info();
      if (r_info.type() != 2) continue;
      uc_vec3 t = orth * r_info.ev();
      uc_vec3 tau = r.transpose().info().ev() * frac;
      double numerator = scitbx::fn::absolute(t * tau);
      double denominator = std::sqrt(t.length_sq() * tau.length_sq());
      CCTBX_ASSERT(denominator != 0);
      double cos_delta = numerator / denominator;
      scitbx::math::update_min(min_cos_delta, cos_delta);
    }
    return std::acos(min_cos_delta)/scitbx::constants::pi_180;
  }

}}} // namespace cctbx::sgtbx::lattice_symmetry
