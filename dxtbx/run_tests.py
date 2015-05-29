from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = (
    "$D/tests/tst_dxtbx.py",
    "$D/tests/tst_imageset.py",
    "$D/tests/tst_datablock.py",
    "$D/tests/command_line/tst_to_xds.py",
    "$D/tests/model/tst_beam.py",
    "$D/tests/model/tst_detector.py",
    "$D/tests/model/tst_detector2.py",
    "$D/tests/model/tst_parallax_correction.py",
    "$D/tests/model/tst_pickle.py",
    "$D/tests/model/tst_pixel_to_millimeter.py",
    "$D/tests/model/tst_ray_intersection.py",
    "$D/tests/model/tst_scan_data.py",
    "$D/tests/model/tst_scan_helpers.py",
    "$D/tests/model/tst_crystal_model.py",
    "$D/tests/model/tst_profile.py",
    "$D/tests/model/experiment/tst_experiment_list.py",
    "$D/tests/serialize/tst_serialize.py",
    "$D/tests/serialize/tst_xds.py",
    "$D/tests/serialize/tst_filename.py",
    "$D/tests/serialize/tst_crystal_model_serialize.py",
    )

def run () :
  build_dir = libtbx.env.under_build("dxtbx")
  dist_dir = libtbx.env.dist_path("dxtbx")
  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
