from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

tst_list = (
    "$D/nanoBragg/tst_nanoBragg_minimal.py",
    "$D/nanoBragg/tst_nanoBragg_mosaic.py",
    "$D/nanoBragg/tst_gaussian_mosaicity.py",
    "$D/nanoBragg/tst_gaussian_mosaicity2.py",
    "$D/nanoBragg/tst_nanoBragg_cbf_write.py",
    "$D/diffBragg/tests/tst_diffBragg_all_refine.py",
    "$D/diffBragg/tests/tst_diffBragg_change_of_basis.py",
    "$D/diffBragg/tests/tst_diffBragg_update_dxtbx_geoms.py",
    "$D/diffBragg/tests/tst_diffBragg_deriv_rois.py",
    "$D/diffBragg/tests/tst_diffBragg_detdist_derivatives.py",
    "$D/diffBragg/tests/tst_diffBragg_nanoBragg_congruency.py",
    "$D/diffBragg/tests/tst_diffBragg_ncells_property.py",
    "$D/diffBragg/tests/tst_diffBragg_ncells_refine.py",
    "$D/diffBragg/tests/tst_diffBragg_originZ_refine.py",
    "$D/diffBragg/tests/tst_diffBragg_regions_of_interest.py",
    "$D/diffBragg/tests/tst_diffBragg_rotXYZ.py",
    "$D/diffBragg/tests/tst_diffBragg_rotXYZ_deriv.py",
    "$D/diffBragg/tests/tst_diffBragg_rotXYZ_refine.py",
    "$D/diffBragg/tests/tst_diffBragg_rotXYZ_ucell_refine.py",
    ["$D/diffBragg/tests/tst_diffBragg_all_refine.py","--umatrix --bmatrix --ncells --curvatures"],
    ["$D/diffBragg/tests/tst_diffBragg_ucell_refine.py","--crystalsystem monoclinic --curvatures"],
    ["$D/diffBragg/tests/tst_diffBragg_ucell_refine.py","--crystalsystem tetragonal"]
    )


def run():
  build_dir = libtbx.env.under_build("simtbx")
  dist_dir = libtbx.env.dist_path("simtbx")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
