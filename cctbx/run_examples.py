from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$B/../exe_dev/cctbx.getting_started",
  "$B/../exe_dev/cctbx.sym_equiv_sites",
  "$D/examples/getting_started.py",
  ["$D/examples/space_group_matrices.py", "P31"],
  "$D/examples/space_subgroups.py",
  "$D/examples/trigonal_r_vs_tetragonal_i.py",
  "$D/examples/quartz_structure.py",
  "$D/examples/fft_map_electron_density_around_atom.py",
  "$D/examples/find_sys_abs_equiv_space_groups.py",
  "$D/examples/find_distances_using_cpp_objects.py",
  "$D/examples/analyze_adp.py",
  "$D/examples/g_exp_i_partial_derivatives.py",
  "$D/examples/tst_exp_i_alpha_derivatives.py",
  "$D/examples/tst_g_exp_i_alpha_derivatives.py",
  "$D/examples/tst_structure_factor_derivatives.py",
  "$D/examples/tst_structure_factor_derivatives_2.py",
  ["$D/examples/tst_structure_factor_derivatives_3.py", "P31"],
  ["$D/examples/tst_structure_factor_derivatives_4.py", "--tag=internal"],
  "$D/examples/structure_factor_calculus/site_derivatives.py",
  "$D/examples/structure_factor_calculus/u_star_derivatives.py",
  "$D/examples/structure_factor_calculus/sites_derivatives.py",
  "$D/examples/structure_factor_calculus/sites_least_squares_derivatives.py",
  ["$D/examples/all_axes.py", "P31"],
  ["$D/examples/tst_phase_o_phrenia.py", "P2"],
  "$D/examples/map_skewness.py",
  "$D/examples/site_symmetry_table.py",
  "$D/examples/site_symmetry_constraints.py",
  "$D/examples/adp_symmetry_constraints.py",
  "$D/examples/unit_cell_refinement.py",
  "$D/examples/miller_common_sets.py",
  "$D/examples/miller_transform.py",
  "$D/examples/change_hand_p31.py",
  "$D/examples/steve_collins.py",
  "$D/examples/cr2o3_primitive_cell.py",
  "$D/examples/cr2o3_consistency_checks.py",
  "$D/examples/reduced_cell_two_folds.py",
  "$D/examples/lebedev_2005_perturbation.py",
  "$D/examples/le_page_1982_vs_lebedev_2005.py",
  "$D/examples/random_puckers.py",
  )

  build_dir = libtbx.env.under_build("cctbx")
  dist_dir = libtbx.env.dist_path("cctbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
