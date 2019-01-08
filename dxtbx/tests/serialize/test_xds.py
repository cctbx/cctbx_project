from __future__ import absolute_import, division, print_function

from six.moves import StringIO
import glob
import os

from libtbx.test_utils import approx_equal
from dxtbx.imageset import ImageSetFactory
from dxtbx.serialize import xds

def test_to_xds(dials_regression, tmpdir):
  tmpdir.chdir()
  template = os.path.join(dials_regression, 'centroid_test_data', 'centroid_00*.cbf')
  file_names = glob.glob(template)
  sweep = ImageSetFactory.new(file_names)[0]
  to_xds = xds.to_xds(sweep)
  s1 = to_xds.XDS_INP()
  expected = '''\
DETECTOR=PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD=495976
SENSOR_THICKNESS= 0.320
!SENSOR_MATERIAL / THICKNESS Si 0.320
!SILICON= 3.960382
DIRECTION_OF_DETECTOR_X-AXIS= 1.00000 0.00000 0.00000
DIRECTION_OF_DETECTOR_Y-AXIS= 0.00000 1.00000 0.00000
NX=2463 NY=2527 QX=0.1720 QY=0.1720
DETECTOR_DISTANCE= 190.180000
ORGX= 1235.84 ORGY= 1279.58
ROTATION_AXIS= 1.00000 0.00000 0.00000
STARTING_ANGLE= 0.000
OSCILLATION_RANGE= 0.200
X-RAY_WAVELENGTH= 0.97950
INCIDENT_BEAM_DIRECTION= -0.000 -0.000 1.021
FRACTION_OF_POLARIZATION= 0.999
POLARIZATION_PLANE_NORMAL= 0.000 1.000 0.000
NAME_TEMPLATE_OF_DATA_FRAMES= %s
TRUSTED_REGION= 0.0 1.41
UNTRUSTED_RECTANGLE= 487 495 0 2528
UNTRUSTED_RECTANGLE= 981 989 0 2528
UNTRUSTED_RECTANGLE= 1475 1483 0 2528
UNTRUSTED_RECTANGLE= 1969 1977 0 2528
UNTRUSTED_RECTANGLE= 0 2464 195 213
UNTRUSTED_RECTANGLE= 0 2464 407 425
UNTRUSTED_RECTANGLE= 0 2464 619 637
UNTRUSTED_RECTANGLE= 0 2464 831 849
UNTRUSTED_RECTANGLE= 0 2464 1043 1061
UNTRUSTED_RECTANGLE= 0 2464 1255 1273
UNTRUSTED_RECTANGLE= 0 2464 1467 1485
UNTRUSTED_RECTANGLE= 0 2464 1679 1697
UNTRUSTED_RECTANGLE= 0 2464 1891 1909
UNTRUSTED_RECTANGLE= 0 2464 2103 2121
UNTRUSTED_RECTANGLE= 0 2464 2315 2333
DATA_RANGE= 1 9
JOB=XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT\
''' % os.path.join(dials_regression, 'centroid_test_data', 'centroid_????.cbf')
  assert s1 == expected
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
