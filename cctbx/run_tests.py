import os, os.path
from scitbx import test_utils

def run():
  tst_list = (
  "$D/math/boost_python/tst_math.py",
  "$D/array_family/boost_python/tst_flex.py",
  "$D/uctbx/boost_python/tst_uctbx.py",
  "$D/sgtbx/boost_python/tst_sgtbx.py",
  "$D/adptbx/boost_python/tst_adptbx.py",
  "$D/miller/boost_python/tst_miller.py",
  "$D/eltbx/boost_python/tst_xray_scattering.py",
  "$D/eltbx/boost_python/tst_henke.py",
  "$D/eltbx/boost_python/tst_icsd_radii.py",
  "$D/eltbx/boost_python/tst_neutron.py",
  "$D/eltbx/boost_python/tst_sasaki.py",
  "$D/eltbx/boost_python/tst_tiny_pse.py",
  "$D/eltbx/boost_python/tst_wavelengths.py",
  "$D/xray/boost_python/tst_xray.py",
  "$D/maptbx/boost_python/tst_maptbx.py",
  "$D/mintbx/boost_python/tst_mintbx.py",
  "$D/dmtbx/boost_python/tst_dmtbx.py",
  "$D/translation_search/boost_python/tst_translation_search.py",
  "$D/cctbx/matrix.py",
  ["$D/cctbx/regression/tst_krivy_gruber.py", "--Quick"],
  "$D/cctbx/regression/tst_sgtbx.py",
  ["$D/cctbx/regression/tst_sgtbx_denominators.py", "P31"],
  "$D/cctbx/regression/tst_sgtbx_subgroups.py",
  "$D/cctbx/regression/tst_sgtbx_lattice_symmetry.py",
  "$D/cctbx/regression/tst_sgtbx_harker.py",
  "$D/cctbx/sgtbx/symbol_confidence.py",
  "$D/cctbx/sgtbx/bravais_types.py",
  "$D/cctbx/regression/tst_crystal.py",
  ["$D/cctbx/regression/tst_xray.py", "I41/acd"],
  ["$D/cctbx/regression/tst_miller.py", "P31"],
  ["$D/cctbx/regression/tst_reciprocal_space_asu.py", "P312"],
  ["$D/cctbx/regression/tst_triplet_generator.py", "P41"],
  ["$D/cctbx/regression/tst_emma.py", "P31"],
  ["$D/cctbx/regression/tst_find_centre_of_inversion.py", "P31"],
  ["$D/cctbx/regression/tst_expand_to_p1.py", "P31"],
  ["$D/cctbx/regression/tst_change_basis.py", "P31"],
  ["$D/cctbx/regression/tst_wilson_plot.py", "P31"],
  ["$D/cctbx/regression/tst_xray_derivatives.py", "P31"],
  ["$D/cctbx/regression/tst_xray_fast_gradients.py", "P31"],
  ["$D/cctbx/regression/tst_xray_minimization.py", "P31"],
  ["$D/cctbx/regression/tst_maptbx_structure_factors.py", "P31"],
  ["$D/cctbx/regression/tst_k_b_scaling.py", "P31"],
  ["$D/cctbx/regression/tst_miller_merge_equivalents.py", "P31"],
  ["$D/cctbx/regression/tst_miller_fft_map.py", "P31"],
  ["$D/cctbx/regression/tst_sampled_model_density.py", "P31"],
  ["$D/cctbx/regression/tst_fast_nv1995.py", "F222"],
  ["$D/cctbx/development/tst_cns_epsilon.py", "P31"],
  ["$D/cctbx/development/tst_cns_hl.py", "P31"],
  ["$D/cctbx/development/run_shelx.py", "P31"],
  )

  build_dir = os.path.join(os.environ["LIBTBX_BUILD"], "cctbx")
  dist_dir = os.environ["CCTBX_DIST"]

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
