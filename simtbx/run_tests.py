from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import libtbx.load_env

tst_list = [
    "$D/nanoBragg/tst_nanoBragg_minimal.py",
    "$D/nanoBragg/tst_nanoBragg_mosaic.py",
    "$D/nanoBragg/tst_gaussian_mosaicity.py",
    "$D/nanoBragg/tst_gaussian_mosaicity2.py",
    "$D/nanoBragg/tst_nanoBragg_cbf_write.py",
    "$D/nanoBragg/tst_multisource_background.py",
    ]

OPT = libtbx.env.build_options
if OPT.enable_cuda:

  tst_list_parallel = [
    ["$D/nanoBragg/tst_gauss_argchk.py","GPU"], # tests CPU+GPU, argchk optimization
    "$D/gpu/tst_exafel_api.py",                 # CPU / GPU, polychromatic beam, monolithic detector
    "$D/gpu/tst_gpu_multisource_background.py", # CPU / GPU background comparison
    "$D/kokkos/tst_kokkos_lib.py",                  # GPU in kokkos
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
