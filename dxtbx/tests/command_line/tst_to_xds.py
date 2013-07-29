from __future__ import division

import glob
import os

import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import open_tmp_file
from libtbx.test_utils import show_diff

from dxtbx.serialize import dump
from dxtbx.imageset import ImageSetFactory

def exercise_to_xds():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_to_xds(): dials_regression not present"
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/centroid_test_data",
    test=os.path.isdir)
  template = os.path.join(data_dir, "centroid_00*.cbf")
  file_names = glob.glob(template)
  cmd = " ".join(["dxtbx.to_xds"] + file_names)
  result = easy_run.fully_buffered(cmd)
  assert not show_diff("\n".join(result.stdout_lines), expected_output)

  # now test reading from a json file
  sweep = ImageSetFactory.new(file_names)[0]
  f = open_tmp_file(suffix="sweep.json", mode="wb")
  dump.imageset(sweep, f)
  f.close()
  cmd = " ".join(["dxtbx.to_xds", f.name])
  result = easy_run.fully_buffered(cmd)
  assert not show_diff("\n".join(result.stdout_lines), expected_output)

expected_output = """\
DETECTOR=PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD=495976
SENSOR_THICKNESS= 0.32
DIRECTION_OF_DETECTOR_X-AXIS= 1.000 0.000 0.000
DIRECTION_OF_DETECTOR_Y-AXIS= 0.000 1.000 -0.000
NX=2463 NY=2527 QX=0.1720 QY=0.1720
DETECTOR_DISTANCE= 190.180
ORGX= 1235.8 ORGY= 1279.6
ROTATION_AXIS= 1.000 0.000 0.000
STARTING_ANGLE= 0.000
OSCILLATION_RANGE= 0.200
X-RAY_WAVELENGTH= 0.97950
INCIDENT_BEAM_DIRECTION= -0.000 0.000 1.000
NAME_TEMPLATE_OF_DATA_FRAMES= /Users/rjgildea/software/phenix/sources/dials_regression/centroid_test_data/centroid_????.cbf
TRUSTED_REGION= 0.0 1.41
UNTRUSTED_RECTANGLE= 487 2 493 2528
UNTRUSTED_RECTANGLE= 981 2 987 2528
UNTRUSTED_RECTANGLE= 1475 2 1481 2528
UNTRUSTED_RECTANGLE= 1969 2 1975 2528
UNTRUSTED_RECTANGLE= 0 197 2462 213
UNTRUSTED_RECTANGLE= 0 409 2462 425
UNTRUSTED_RECTANGLE= 0 621 2462 637
UNTRUSTED_RECTANGLE= 0 833 2462 849
UNTRUSTED_RECTANGLE= 0 1045 2462 1061
UNTRUSTED_RECTANGLE= 0 1257 2462 1273
UNTRUSTED_RECTANGLE= 0 1469 2462 1485
UNTRUSTED_RECTANGLE= 0 1681 2462 1697
UNTRUSTED_RECTANGLE= 0 1893 2462 1909
UNTRUSTED_RECTANGLE= 0 2105 2462 2121
UNTRUSTED_RECTANGLE= 0 2317 2462 2333
DATA_RANGE= 1 9
JOB=XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT\
"""

def run():
  exercise_to_xds()
  print "OK"


if __name__ == '__main__':
  run()
