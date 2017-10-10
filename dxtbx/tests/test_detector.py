from __future__ import absolute_import, division
from __future__ import print_function
#!/usr/bin/env python
# test_detector.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Tests for the detector class.

from builtins import range
from scitbx import matrix

from dxtbx.model.detector import DetectorFactory
from dxtbx.model.detector_helpers import compute_frame_rotation

def test_detector():
  '''A test class for the detector class.'''

  d = DetectorFactory.simple('CCD', 100.0, (45.0, 52.0), '+x', '-y',
                              (0.172, 0.172), (516, 590), (0, 1024), [])
  t = DetectorFactory.two_theta(
      'CCD', 60.0, (35.0, 34.0), '+x', '+y', '+x', 30,
      (0.07, 0.07), (1042, 1042), (0, 1024), [])

  import libtbx.load_env
  import os

  dxtbx_dir = libtbx.env.dist_path('dxtbx')

  image = os.path.join(dxtbx_dir, 'tests', 'phi_scan_001.cbf')
  xparm = os.path.join(dxtbx_dir, 'tests', 'example-xparm.xds')

  c = DetectorFactory.imgCIF(image, 'CCD')
  #x = DetectorFactory.XDS(xparm)

  print('OK')

def work_detector():

  for j in range(10000):
    c = DetectorFactory.imgCIF('phi_scan.cbf')

def work_detector_helpers():
  compute_frame_rotation((matrix.col((1, 0, 0)),
                          matrix.col((0, 1, 0)),
                          matrix.col((0, 0, 1))),
                         (matrix.col((1, 0, 0)),
                          matrix.col((0, 1, 0)),
                          matrix.col((0, 0, 1))))

  compute_frame_rotation((matrix.col((-1, 0, 0)),
                          matrix.col((0, 0, 1)),
                          matrix.col((0, 1, 0))),
                         (matrix.col((1, 0, 0)),
                          matrix.col((0, 1, 0)),
                          matrix.col((0, 0, 1))))

if __name__ == '__main__':

  test_detector()
