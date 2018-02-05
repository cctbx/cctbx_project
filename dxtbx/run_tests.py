from __future__ import absolute_import, division, print_function

from libtbx.test_utils.pytest import discover

tst_list = [
    "$D/tests/tst_dxtbx.py",
    "$D/tests/tst_imageset.py",
    "$D/tests/tst_datablock.py",
    "$D/tests/tst_filecache.py",
    "$D/tests/tst_beamline_definitions.py",
    "$D/tests/tstFormatCBFFull.py",
    "$D/tests/tst_dials_226.py",
    "$D/tests/tst_image_readers.py",
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
    "$D/tests/model/tst_to_from_dict.py",
    "$D/tests/model/tst_crystal_model.py",
    "$D/tests/model/tst_profile.py",
    "$D/tests/model/experiment/tst_experiment_list.py",
    "$D/tests/serialize/tst_serialize.py",
    "$D/tests/serialize/tst_xds.py",
    "$D/tests/serialize/tst_filename.py",
    "$D/tests/serialize/tst_crystal_model_serialize.py",
    ] + discover()
