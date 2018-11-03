from __future__ import absolute_import, division, print_function

from six.moves import StringIO
import glob
import os

from libtbx.test_utils import approx_equal
from dxtbx.imageset import ImageSetFactory
from dxtbx.serialize import xds

def test_to_xds(dials_regression, tmpdir):
  tmpdir.chdir()
  template = os.path.join(dials_regression, 'centroid_test_data', "centroid_00*.cbf")
  file_names = glob.glob(template)
  sweep = ImageSetFactory.new(file_names)[0]
  to_xds = xds.to_xds(sweep)
  s1 = StringIO()
  to_xds.XDS_INP(out=s1)
  empty = StringIO()
  s1a = to_xds.XDS_INP(as_str=True, out=empty)
  assert empty.getvalue() == ''
  assert s1.getvalue().strip() == s1a
  s2 = StringIO()
  real_space_a = (-5.327642, -39.034747, -4.988286)
  real_space_b = (-35.253495, 7.596265, -22.127661)
  real_space_c = (-22.673623, -1.486119, 35.793463)
  to_xds.xparm_xds(real_space_a, real_space_b, real_space_c, space_group=1, out=s2)
  # run coordinate frame converter on xparm.xds as a sanity check
  with open('xparm.xds', mode="wb") as fh:
    s2.seek(0)
    fh.writelines(s2.readlines())
  from rstbx.cftbx import coordinate_frame_helpers
  converter = coordinate_frame_helpers.import_xds_xparm('xparm.xds')
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
