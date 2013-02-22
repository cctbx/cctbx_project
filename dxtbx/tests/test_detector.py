#!/usr/bin/env python
# test_detector.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Tests for the detector class.

import os
import sys
from scitbx import matrix

from dxtbx.model.detector import detector
from dxtbx.model.detector import detector_factory
from dxtbx.model.detector_helpers import compute_frame_rotation

def test_detector():
    '''A test class for the detector class.'''

    d = detector_factory.simple('CCD', 100.0, (45.0, 52.0), '+x', '-y',
                                (0.172, 0.172), (516, 590), (0, 1024), [])
    t = detector_factory.two_theta(
        'CCD', 60.0, (35.0, 34.0), '+x', '+y', '+x', 30,
        (0.07, 0.07), (1042, 1042), (0, 1024), [])
    c = detector_factory.imgCIF('phi_scan_001.cbf', 'CCD')
    x = detector_factory.XDS('example-xparm.xds')

    print t
    print x

def work_detector():

    for j in range(10000):
        c = detector_factory.imgCIF('phi_scan.cbf')

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
