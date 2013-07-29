from __future__ import division
#!/usr/bin/env python
# to_xds.py
#
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Print out the contents of the dxtbx understanding of a bunch of images to
# an example XDS.INP file. This should illustrate the usage of the dxtbx
# classes.

import sys
import math
from scitbx import matrix

from dxtbx.model.detector_helpers_types import detector_helpers_types

def xds_detector_name(xia2_name):
    '''Translate from a xia2 name from the detector library to an XDS detector
    name.'''

    if 'pilatus' in xia2_name:
        return 'PILATUS'
    if 'rayonix' in xia2_name:
        return 'CCDCHESS'
    if 'adsc' in xia2_name:
        return 'ADSC'
    if 'saturn' in xia2_name:
        return 'SATURN'
    if 'raxis' in xia2_name:
        return 'RAXIS'

    raise RuntimeError, 'detector %s unknown' % xia2_name

class to_xds(object):
    '''A class to export contents of an XSweep2 as XDS.INP.'''

    def __init__(self, sweep):
        self._sweep = sweep

        return

    def get_detector(self):
        return self._sweep.get_detector()

    def get_goniometer(self):
        return self._sweep.get_goniometer()

    def get_beam(self):
        return self._sweep.get_beam()

    def get_scan(self):
        return self._sweep.get_scan()

    def get_template(self):
        return self._sweep.get_template()

    def XDS(self):

        sensor = self.get_detector().get_type()
        fast, slow = map(int, self.get_detector().get_image_size())
        f, s = self.get_detector().get_pixel_size()
        df = int(1000 * f)
        ds = int(1000 * s)

        # FIXME probably need to rotate by pi about the X axis

        R = matrix.col((1.0, 0.0, 0.0)).axis_and_angle_as_r3_rotation_matrix(
            180.0, deg = True)

        detector = xds_detector_name(
            detector_helpers_types.get(sensor, fast, slow, df, ds))
        trusted = self.get_detector().get_trusted_range()

        print 'DETECTOR=%s MINIMUM_VALID_PIXEL_VALUE=%d OVERLOAD=%d' % \
              (detector, trusted[0] + 1, trusted[1])

        if detector == 'PILATUS':
            print 'SENSOR_THICKNESS= 0.32'

        print 'DIRECTION_OF_DETECTOR_X-AXIS= %.3f %.3f %.3f' % \
              (R * matrix.col(self.get_detector().get_fast_axis())).elems

        print 'DIRECTION_OF_DETECTOR_Y-AXIS= %.3f %.3f %.3f' % \
              (R * matrix.col(self.get_detector().get_slow_axis())).elems

        print 'NX=%d NY=%d QX=%.4f QY=%.4f' % (fast, slow, f, s)

        F = R * matrix.col(self.get_detector().get_fast_axis())
        S = R * matrix.col(self.get_detector().get_slow_axis())
        N = F.cross(S)

        origin = R * matrix.col(self.get_detector().get_origin())
        beam = R * matrix.col(self.get_beam().get_direction()) / \
               math.sqrt(matrix.col(self.get_beam().get_direction()).dot())
        centre = -(origin - origin.dot(N) * N)
        x = centre.dot(F)
        y = centre.dot(S)

        print 'DETECTOR_DISTANCE= %.3f' % origin.dot(N)
        print 'ORGX= %.1f ORGY= %.1f' % (x / f, y / s)
        print 'ROTATION_AXIS= %.3f %.3f %.3f' % \
              (R * matrix.col(self.get_goniometer().get_rotation_axis())).elems
        print 'STARTING_ANGLE= %.3f' % \
              self.get_scan().get_oscillation()[0]
        print 'OSCILLATION_RANGE= %.3f' % \
              self.get_scan().get_oscillation()[1]
        print 'X-RAY_WAVELENGTH= %.5f' % \
              self.get_beam().get_wavelength()
        print 'INCIDENT_BEAM_DIRECTION= %.3f %.3f %.3f' % \
              (- beam).elems

        # FIXME LATER
        if hasattr(self.get_beam(), "get_polatization_fraction"):
            print 'FRACTION_OF_POLARIZATION= %.3f' % \
                self.get_beam().get_polarization_fraction()
            print 'POLARIZATION_PLANE_NORMAL= %.3f %.3f %.3f' % \
                self.get_beam().get_polarization()
        print 'NAME_TEMPLATE_OF_DATA_FRAMES= %s' % \
            self.get_template().replace('#', '?')
        print 'TRUSTED_REGION= 0.0 1.41'
        for f0, f1, s0, s1 in self.get_detector().get_mask():
            print 'UNTRUSTED_RECTANGLE= %d %d %d %d' % \
                  (f0 - 1, f1 + 1, s0 - 1, s1 + 1)

        start_end = self.get_scan().get_image_range()

        if start_end[0] == 0:
            start_end = (1, start_end[1])

        print 'DATA_RANGE= %d %d' % start_end
        print 'JOB=XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT'


if __name__ == '__main__':


    file_names = sys.argv[1:]
    if len(file_names) == 1 and file_names[0].endswith('json'):
        from dxtbx.serialize import load
        sweep = load.imageset(file_names[0])
    else:
        from dxtbx.imageset import ImageSetFactory
        sweep = ImageSetFactory.new(file_names)[0]

    xsx = to_xds(sweep)

    xsx.XDS()
