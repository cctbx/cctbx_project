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
    bool enforce_max_delta_for_generated_two_folds)
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
      if (enforce_max_delta_for_generated_two_folds
          && scitbx::constants::pi_180*find_max_delta(
               niggli_cell, expanded_group) > max_delta) {
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
    space_group const& group)
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

  space_group
  group_search_fast(
    uctbx::unit_cell const& niggli_cell,
    double max_delta,
    bool enforce_max_delta_for_generated_two_folds)
  {
    uc_mat3 const& frac = niggli_cell.fractionalization_matrix();
    uc_mat3 const& orth = niggli_cell.orthogonalization_matrix();
    double min_cos_delta = std::cos(max_delta * scitbx::constants::pi_180);
    unsigned n_two_folds = sizeof(reduced_cell_two_folds)
                         / sizeof(reduced_cell_two_fold_info);
    std::vector<unsigned> i_two_folds;
    i_two_folds.reserve(n_two_folds);
    std::vector<double> cos_deltas;
    cos_deltas.reserve(n_two_folds);
    for(unsigned i_two_fold=0;i_two_fold<n_two_folds;i_two_fold++) {
      reduced_cell_two_fold_info const&
        two_fold = reduced_cell_two_folds[i_two_fold];
      uc_vec3 t = orth * two_fold.u;
      uc_vec3 tau = two_fold.h * frac;
      double numerator = scitbx::fn::absolute(t * tau);
      double denominator = std::sqrt(t.length_sq() * tau.length_sq());
      CCTBX_ASSERT(denominator != 0);
      double cos_delta = numerator / denominator;
      if (cos_delta >= min_cos_delta) {
        i_two_folds.push_back(i_two_fold);
        cos_deltas.push_back(numerator / denominator);
      }
    }
    bool reverse = true;
    af::shared<std::size_t> perm_memory = af::sort_permutation(
      af::const_ref<double>(&*cos_deltas.begin(), cos_deltas.size()), reverse);
    af::const_ref<std::size_t> perm = perm_memory.const_ref();
    space_group group;
    for(std::size_t i=0;i<perm.size();i++) {
      reduced_cell_two_fold_info const&
        two_fold = reduced_cell_two_folds[i_two_folds[perm[i]]];
      rt_mx s(rot_mx(two_fold.r,1));
      space_group expanded_group(group);
      try {
        expanded_group.expand_smx(s);
      }
      catch (error_non_crystallographic_rotation_matrix_encountered const&) {
        break;
      }
      if (enforce_max_delta_for_generated_two_folds
          && find_max_delta(niggli_cell, expanded_group) > max_delta) {
        continue;
      }
      group = expanded_group;
    }
    return group;
  }

}}} // namespace cctbx::sgtbx::lattice_symmetry
