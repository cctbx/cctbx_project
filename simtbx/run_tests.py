from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

tst_list = [
    "$D/nanoBragg/tst_nanoBragg_minimal.py",
    "$D/nanoBragg/tst_nanoBragg_mosaic.py",
    "$D/nanoBragg/tst_gaussian_mosaicity.py",
    "$D/nanoBragg/tst_gaussian_mosaicity2.py",
    "$D/nanoBragg/tst_nanoBragg_cbf_write.py",
    "$D/diffBragg/tests/tst_diffBragg_Fhkl_complex.py",
    "$D/diffBragg/tests/tst_diffBragg_change_of_basis.py",
    "$D/diffBragg/tests/tst_diffBragg_update_dxtbx_geoms.py",
    "$D/diffBragg/tests/tst_diffBragg_deriv_rois.py",
    "$D/diffBragg/tests/tst_diffBragg_detdist_derivatives.py",
    "$D/diffBragg/tests/tst_diffBragg_nanoBragg_congruency.py",
    "$D/diffBragg/tests/tst_diffBragg_ncells_property.py",
    ["$D/diffBragg/tests/tst_diffBragg_ncells_property_anisotropic.py", "--idx 0"],
    ["$D/diffBragg/tests/tst_diffBragg_ncells_property_anisotropic.py", "--idx 1"],
    ["$D/diffBragg/tests/tst_diffBragg_ncells_property_anisotropic.py", "--idx 2"],
    ["$D/diffBragg/tests/tst_diffBragg_ncells_offdiag_property.py", "--idx 0"],
    ["$D/diffBragg/tests/tst_diffBragg_ncells_offdiag_property.py", "--idx 1"],
    ["$D/diffBragg/tests/tst_diffBragg_ncells_offdiag_property.py", "--idx 2"],
    ["$D/diffBragg/tests/tst_diffBragg_lambda_coefficients.py", "--idx 0"],
    ["$D/diffBragg/tests/tst_diffBragg_lambda_coefficients.py", "--idx 1"],
    "$D/diffBragg/tests/tst_diffBragg_ncells_refine.py",
    "$D/diffBragg/tests/tst_diffBragg_ncells_offdiag_refine.py",
    ["$D/diffBragg/tests/tst_diffBragg_ncells_refine.py","--testrangeLower"],
    ["$D/diffBragg/tests/tst_diffBragg_ncells_refine.py","--testrangeUpper"],
    "$D/diffBragg/tests/tst_diffBragg_ncells_aniso_refine.py",
    ["$D/diffBragg/tests/tst_diffBragg_ncells_aniso_refine.py", "--constrain"],
    ["$D/diffBragg/tests/tst_diffBragg_ncells_aniso_refine.py","--testrangeLower"],
    ["$D/diffBragg/tests/tst_diffBragg_ncells_aniso_refine.py","--testrangeUpper"],
    "$D/diffBragg/tests/tst_diffBragg_detdist_refine.py",
    "$D/diffBragg/tests/tst_diffBragg_regions_of_interest.py",
    "$D/diffBragg/tests/tst_diffBragg_rotXYZ.py",
    "$D/diffBragg/tests/tst_diffBragg_refine_spectrum.py",
    ["$D/diffBragg/tests/tst_diffBragg_refine_eta.py", "--finitediff"],
    ["$D/diffBragg/tests/tst_diffBragg_refine_eta.py", "--finitediff", "--curvatures"],
    ["$D/diffBragg/tests/tst_diffBragg_refine_eta.py", "--finitediff", "--curvatures", "--aniso 0"],
    ["$D/diffBragg/tests/tst_diffBragg_refine_eta.py", "--finitediff", "--curvatures", "--aniso 1"],
    ["$D/diffBragg/tests/tst_diffBragg_refine_eta.py", "--finitediff", "--curvatures", "--aniso 2"],
    #["$D/diffBragg/tests/tst_diffBragg_refine_eta.py", "--curvatures", "--aniso 0"],# NOTE, this test is really long, recommended to run seperately on an OMP build
    #"$D/diffBragg/tests/tst_diffBragg_refine_eta.py", # NOTE, this test is really long, recommended to run seperately on an OMP build
    #["$D/diffBragg/tests/tst_diffBragg_refine_eta.py","--testUpperBound" ], # NOTE, this test is really long, recommended to run seperately on an OMP build
    #["$D/diffBragg/tests/tst_diffBragg_refine_eta.py", "--curvatures"], # NOTE, this test is really long, recommended to run seperately on an OMP build
    ["$D/diffBragg/tests/tst_diffBragg_rotXYZ_deriv.py", "--curvatures --rotidx 0"],
    ["$D/diffBragg/tests/tst_diffBragg_rotXYZ_deriv.py", "--curvatures --rotidx 1"],
    ["$D/diffBragg/tests/tst_diffBragg_rotXYZ_deriv.py", "--curvatures --rotidx 2"],
    ["$D/diffBragg/tests/tst_diffBragg_rotXYZ_ucell_refine.py", "--curvatures"],
    ["$D/diffBragg/tests/tst_diffBragg_ucell_refine.py", "--crystalsystem monoclinic --curvatures"],
    ["$D/diffBragg/tests/tst_diffBragg_ucell_refine.py", "--crystalsystem tetragonal"],
    ["$D/diffBragg/tests/tst_diffBragg_Fcell_deriv.py", "--curvatures"],
    ["$D/diffBragg/tests/tst_diffBragg_global_refine.py",
        "--nshots 1 --rescale --spotscale --umatrix --ncells " +
        "--bmatrix --bg --fcell --testbg --testfcell --testUmatrix --maxcalls 50 --sz 18"],
    ["$D/diffBragg/tests/tst_diffBragg_global_refine.py",
        "--nshots 1 --rescale --umatrix --bmatrix --curvatures --testUmatrix"],
    ["$D/diffBragg/tests/tst_diffBragg_global_refine.py",
     "--nshots 1 --rescale --umatrix --testUmatrix"],
    ["$D/diffBragg/tests/tst_diffBragg_global_refine.py", "--nshots 1 --rescale --spectra  --maxcalls 20  --nshots 1 " +
        "--refineSpectra  --testSpectra  --perturbSpectra"],
    ["$D/diffBragg/tests/tst_diffBragg_global_refine_older.py", "--spotscale --umatrix --bmatrix --ncells --curvatures --rescale"],
    "$D/diffBragg/tests/tst_diffBragg_panel_RotOFS_shiftXYZ_refine.py",
    ["$D/diffBragg/tests/tst_diffBragg_panelXY_derivs.py", "--panel x"],
    ["$D/diffBragg/tests/tst_diffBragg_panelXY_derivs.py", "--panel y"],
    ["$D/diffBragg/tests/tst_diffBragg_panelXY_derivs.py", "--panel z"],
    #["$D/diffBragg/tests/tst_diffBragg_blue_sausages.py", "--finitediff"],
    ]

OPT = libtbx.env.build_options
if OPT.enable_cuda:

  tst_list_parallel = [
    ["$D/nanoBragg/tst_gauss_argchk.py","GPU"], # tests CPU+GPU, argchk optimization
    "$D/gpu/tst_exafel_api.py",                 # CPU / GPU, polychromatic beam, monolithic detector
  ]
else:
    tst_list.append(
        ["$D/nanoBragg/tst_gauss_argchk.py","CPU"]
    )

def run():
  build_dir = libtbx.env.under_build("simtbx")
  dist_dir = libtbx.env.dist_path("simtbx")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
