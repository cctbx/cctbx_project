from __future__ import absolute_import, division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from io import StringIO
import glob
import os

import libtbx.load_env
from libtbx.test_utils import open_tmp_file
from libtbx.test_utils import approx_equal
from dxtbx.imageset import ImageSetFactory
from dxtbx.serialize import xds

def exercise_to_xds():
  if not libtbx.env.has_module("dials"):
    print("Skipping test: dials not present")
    return
  if not libtbx.env.has_module("dials_regression"):
    print("Skipping exercise_to_xds(): dials_regression not present")
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/centroid_test_data",
    test=os.path.isdir)
  template = os.path.join(data_dir, "centroid_00*.cbf")
  file_names = glob.glob(template)
  sweep = ImageSetFactory.new(file_names)[0]
  to_xds = xds.to_xds(sweep)
  s1 = StringIO()
  to_xds.XDS_INP(out=s1)
  s2 = StringIO()
  real_space_a = (-5.327642, -39.034747, -4.988286)
  real_space_b = (-35.253495, 7.596265, -22.127661)
  real_space_c = (-22.673623, -1.486119, 35.793463)
  to_xds.xparm_xds(real_space_a, real_space_b, real_space_c, space_group=1, out=s2)
  # run coordinate frame converter on xparm.xds as a sanity check
  f = open_tmp_file(suffix="XPARM.XDS", mode="wb")
  s2.seek(0)
  f.writelines(s2.readlines())
  f.close()
  from rstbx.cftbx import coordinate_frame_helpers
  converter = coordinate_frame_helpers.import_xds_xparm(f.name)
  scan = sweep.get_scan()
  detector = sweep.get_detector()
  goniometer = sweep.get_goniometer()
  beam = sweep.get_beam()
  assert approx_equal(real_space_a, converter.get_real_space_a())
  assert approx_equal(real_space_b, converter.get_real_space_b())
  assert approx_equal(real_space_c, converter.get_real_space_c())
  assert approx_equal(goniometer.get_rotation_axis(),
                      converter.get_rotation_axis())
  assert approx_equal(
    beam.get_direction(), converter.get_sample_to_source().elems)
  assert approx_equal(detector[0].get_fast_axis(), converter.get_detector_fast())
  assert approx_equal(detector[0].get_slow_axis(), converter.get_detector_slow())
  assert approx_equal(detector[0].get_origin(), converter.get_detector_origin())


def run():
  exercise_to_xds()
  print('OK')


if __name__ == '__main__':
  run()
