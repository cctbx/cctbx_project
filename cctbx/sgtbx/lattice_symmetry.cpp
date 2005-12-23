#include <cctbx/sgtbx/lattice_symmetry.h>
#include <cctbx/sgtbx/rot_mx_info.h>
#include <scitbx/array_family/sort.h>

namespace cctbx { namespace sgtbx { namespace lattice_symmetry {

  double
  find_max_delta(
    uctbx::unit_cell const& reduced_cell,
    sgtbx::space_group const& space_group)
  {
    CCTBX_ASSERT(space_group.n_ltr() == 1);
    CCTBX_ASSERT(space_group.f_inv() == 1);
    uc_mat3 const& frac = reduced_cell.fractionalization_matrix();
    uc_mat3 const& orth = reduced_cell.orthogonalization_matrix();
    double min_cos_delta = 1;
    for(int i_smx=1;i_smx<space_group.n_smx();i_smx++) {
      rot_mx const& r = space_group.smx()[i_smx].r();
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
  group(
    uctbx::unit_cell const& reduced_cell,
    double max_delta,
    bool enforce_max_delta_for_generated_two_folds)
  {
    uc_mat3 const& frac = reduced_cell.fractionalization_matrix();
    uc_mat3 const& orth = reduced_cell.orthogonalization_matrix();
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
          && find_max_delta(reduced_cell, expanded_group) > max_delta) {
        continue;
      }
      group = expanded_group;
    }
    return group;
  }

}}} // namespace cctbx::sgtbx::lattice_symmetry
