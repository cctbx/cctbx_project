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

import os
import sys
import math
from scitbx import matrix

from dxtbx.model.detector_helpers_types import detector_helpers_types
from dxtbx.format.Registry import Registry

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

class to_xds:
    '''A class to export contents of an XSweep2 as XDS.INP.'''

    def __init__(self, goniometer, detector, beam, scan):
        self._goniometer = goniometer
        self._detector = detector
        self._beam = beam
        self._scan = scan

        return

    def get_detector(self):
        return self._detector

    def get_goniometer(self):
        return self._goniometer

    def get_beam(self):
        return self._beam

    def get_scan(self):
        return self._scan

    def XDS(self):

        sensor = self.get_detector().type
        fast, slow = map(int, self.get_detector().image_size)
        f, s = self.get_detector().pixel_size
        df = int(1000 * f)
        ds = int(1000 * s)

        # FIXME probably need to rotate by pi about the X axis

        R = matrix.col((1.0, 0.0, 0.0)).axis_and_angle_as_r3_rotation_matrix(
            180.0, deg = True)

        detector = xds_detector_name(
            detector_helpers_types.get(sensor, fast, slow, df, ds))
        trusted = self.get_detector().trusted_range

        print 'DETECTOR=%s MINIMUM_VALID_PIXEL_VALUE=%d OVERLOAD=%d' % \
              (detector, trusted[0] + 1, trusted[1])

        if detector == 'PILATUS':
            print 'SENSOR_THICKNESS= 0.32'

        print 'DIRECTION_OF_DETECTOR_X-AXIS= %.3f %.3f %.3f' % \
              (R * matrix.col(self.get_detector().fast_axis)).elems

        print 'DIRECTION_OF_DETECTOR_Y-AXIS= %.3f %.3f %.3f' % \
              (R * matrix.col(self.get_detector().slow_axis)).elems

        print 'NX=%d NY=%d QX=%.4f QY=%.4f' % (fast, slow, f, s)

        F = R * matrix.col(self.get_detector().fast_axis)
        S = R * matrix.col(self.get_detector().slow_axis)
        N = F.cross(S)

        origin = R * matrix.col(self.get_detector().origin)
        beam = R * matrix.col(self.get_beam().direction) / \
               math.sqrt(matrix.col(self.get_beam().direction).dot())
        centre = -(origin - origin.dot(N) * N)
        x = centre.dot(F)
        y = centre.dot(S)

        print 'DETECTOR_DISTANCE= %.3f' % origin.dot(N)
        print 'ORGX= %.1f ORGY= %.1f' % (x / f, y / s)
        print 'ROTATION_AXIS= %.3f %.3f %.3f' % \
              (R * matrix.col(self.get_goniometer().rotation_axis)).elems
        print 'STARTING_ANGLE= %.3f' % \
              self.get_scan().oscillation[0]
        print 'OSCILLATION_RANGE= %.3f' % \
              self.get_scan().oscillation[1]
        print 'X-RAY_WAVELENGTH= %.5f' % \
              self.get_beam().wavelength
        print 'INCIDENT_BEAM_DIRECTION= %.3f %.3f %.3f' % \
              (- beam).elems
        print 'FRACTION_OF_POLARIZATION= %.3f' % \
              self.get_beam().polarization_fraction
        print 'POLARIZATION_PLANE_NORMAL= %.3f %.3f %.3f' % \
              self.get_beam().polarization
        print 'NAME_TEMPLATE_OF_DATA_FRAMES= %s' % os.path.join(
            self.get_scan().get_directory(),
            self.get_scan().get_template().replace('#', '?'))
        print 'TRUSTED_REGION= 0.0 1.41'
        for f0, f1, s0, s1 in self.get_detector().mask:
            print 'UNTRUSTED_RECTANGLE= %d %d %d %d' % \
                  (f0 - 1, f1 + 1, s0 - 1, s1 + 1)

        start_end = self.get_scan().image_range

        if start_end[0] == 0:
            start_end = (1, start_end[1])

        print 'DATA_RANGE= %d %d' % start_end
        print 'JOB=XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT'

def factory(list_of_images):
    '''Instantiate the data model bits and pieces we need...'''

    for image in list_of_images:
        assert(os.path.exists(image))

    list_of_images.sort()

    format = Registry.find(list_of_images[0])

    # verify that these are all the same format i.e. that they are all
    # understood equally well by the format instance.

    format_score = format.understand(list_of_images[0])

    for image in list_of_images:
        assert(format.understand(image) == format_score)

    i = format(list_of_images[0])

    beam = i.get_beam()
    gonio = i.get_goniometer()
    det = i.get_detector()
    scan = i.get_scan()

    # now verify that they share the same detector position, rotation axis
    # and beam properties.

    scans = [scan]

    for image in list_of_images[1:]:
        i = format(image)
        assert(beam == i.get_beam())
        assert(gonio == i.get_goniometer())
        assert(det == i.get_detector())
        scans.append(i.get_scan())

    for s in sorted(scans)[1:]:
        scan += s

    return gonio, det, beam, scan


if __name__ == '__main__':

    # run some tests

    gonio, det, beam, scan = factory(sys.argv[1:])

    xsx = to_xds(gonio, det, beam, scan)

    xsx.XDS()
