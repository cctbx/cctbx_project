from __future__ import division
#!/usr/bin/env python
# detector.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A model for the detector for the "updated experimental model" project
# documented in internal ticket #1555. This is not designed to be used outside
# of the XSweep classes. N.B. this should probably be generalized for non
# flat detectors, or composite detectors constructed from a number of flat
# elements.

import math
import pycbf
from scitbx import matrix
from dxtbx_model_ext import Panel, Detector, SimplePxMmStrategy, \
    ParallaxCorrectedPxMmStrategy

from detector_helpers import detector_helper_sensors
from detector_helpers import find_undefined_value


class detector_factory:
    '''A factory class for detector objects, which will encapsulate standard
    detector designs to make it a little easier to get started with these. In
    cases where a CBF image is provided a full description can be used, in
    other cases assumptions will be made about the experiment configuration.
    In all cases information is provided in the CBF coordinate frame.'''

    def __init__(self):
        pass

    @staticmethod
    def make_detector(stype, fast_axis, slow_axis, origin,
                      pixel_size, image_size, trusted_range = (0.0, 0.0),
                      px_mm=None, name="Panel"):
        """Ensure all types are correct before creating c++ detector class."""

        if stype == 'SENSOR_PAD':
            px_mm = ParallaxCorrectedPxMmStrategy(0.252500934883)
        else:
            px_mm = SimplePxMmStrategy()

        d = Detector(Panel(
            str(stype),
            str(name),
            tuple(map(float, fast_axis)),
            tuple(map(float, slow_axis)),
            tuple(map(float, origin)),
            tuple(map(float, pixel_size)),
            tuple(map(int, image_size)),
            tuple(map(float, trusted_range)),
            px_mm))

        return d

    @staticmethod
    def simple(sensor, distance, beam_centre, fast_direction, slow_direction,
               pixel_size, image_size, trusted_range = (0.0, 0.0), mask = [],
               px_mm=None):
        '''Construct a simple detector at a given distance from the sample
        along the direct beam presumed to be aligned with -z, offset by the
        beam centre - the directions of which are given by the fast and slow
        directions, which are themselves given as +x, +y, -x, -y. The pixel
        size is given in mm in the fast and slow directions and the image size
        is given in pixels in the same order. Everything else is the same as
        for the main reference frame.'''

        assert(fast_direction in ['-x', '+y', '+x', '-y'])
        assert(slow_direction in ['-x', '+y', '+x', '-y'])

        assert(fast_direction[1] != slow_direction[1])

        direction_map = {
            '+x':(1.0, 0.0, 0.0),
            '-x':(-1.0, 0.0, 0.0),
            '+y':(0.0, 1.0, 0.0),
            '-y':(0.0, -1.0, 0.0)
            }

        fast = matrix.col(direction_map[fast_direction])
        slow = matrix.col(direction_map[slow_direction])

        origin = matrix.col((0, 0, -1)) * distance - \
                 fast * beam_centre[0] - slow * beam_centre[1]

        detector = detector_factory.make_detector(
            detector_factory.sensor(sensor),
            fast, slow, origin, pixel_size, image_size, trusted_range)
        detector.mask = mask
        return detector

    @staticmethod
    def two_theta(sensor, distance, beam_centre, fast_direction,
                  slow_direction, two_theta_direction, two_theta_angle,
                  pixel_size, image_size, trusted_range = (0.0, 0.0),
                  mask = [], px_mm = None):
        '''Construct a simple detector at a given distance from the sample
        along the direct beam presumed to be aligned with -z, offset by the
        beam centre - the directions of which are given by the fast and slow
        directions, which are themselves given as +x, +y, -x, -y. The pixel
        size is given in mm in the fast and slow directions and the image size
        is given in pixels in the same order. Everything else is the same as
        for the main reference frame. Also given are the direction of the
        two-theta axis and the angle in degrees by which the detector is
        moved.'''

        assert(fast_direction in ['-x', '+y', '+x', '-y'])
        assert(slow_direction in ['-x', '+y', '+x', '-y'])
        assert(two_theta_direction in ['-x', '+y', '+x', '-y'])

        assert(fast_direction[1] != slow_direction[1])

        direction_map = {
            '+x':(1.0, 0.0, 0.0),
            '-x':(-1.0, 0.0, 0.0),
            '+y':(0.0, 1.0, 0.0),
            '-y':(0.0, -1.0, 0.0)
            }

        fast = matrix.col(direction_map[fast_direction])
        slow = matrix.col(direction_map[slow_direction])

        origin = matrix.col((0, 0, -1)) * distance - \
                 fast * beam_centre[0] - slow * beam_centre[1]

        two_theta = matrix.col(direction_map[two_theta_direction])

        R = two_theta.axis_and_angle_as_r3_rotation_matrix(two_theta_angle,
                                                           deg = True)

        detector = detector_factory.make_detector(
            detector_factory.sensor(sensor),
            (R * fast), (R * slow), (R * origin), pixel_size,
            image_size, trusted_range)

        detector.mask = mask
        return detector

    @staticmethod
    def complex(sensor, origin, fast, slow, pixel, size,
                trusted_range = (0.0, 0.0), px_mm = None):
        '''A complex detector model, where you know exactly where everything
        is. This is useful for implementation of the Rigaku Saturn header
        format, as that is exactly what is in there. Origin, fast and slow are
        vectors in the CBF reference frame, pixel is the dimensions as a tuple
        as is size.'''

        assert(len(origin) == 3)
        assert(len(fast) == 3)
        assert(len(slow) == 3)
        assert(len(pixel) == 2)
        assert(len(size) == 2)

        return detector_factory.make_detector(
                detector_factory.sensor(sensor),
                fast, slow, origin, pixel, size, trusted_range)

    @staticmethod
    def imgCIF(cif_file, sensor):
        '''Initialize a detector model from an imgCIF file.'''

        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

        cbf_detector = cbf_handle.construct_detector(0)

        pixel = (cbf_detector.get_inferred_pixel_size(1),
                 cbf_detector.get_inferred_pixel_size(2))

        # FIXME can probably simplify the code which follows below by
        # making proper use of cctbx vector calls - should not be as
        # complex as it appears to be...

        origin = tuple(cbf_detector.get_pixel_coordinates(0, 0))
        fast = cbf_detector.get_pixel_coordinates(0, 1)
        slow = cbf_detector.get_pixel_coordinates(1, 0)

        dfast = [fast[j] - origin[j] for j in range(3)]
        dslow = [slow[j] - origin[j] for j in range(3)]

        lfast = math.sqrt(sum([dfast[j] * dfast[j] for j in range(3)]))
        lslow = math.sqrt(sum([dslow[j] * dslow[j] for j in range(3)]))

        fast = tuple([dfast[j] / lfast for j in range(3)])
        slow = tuple([dslow[j] / lslow for j in range(3)])

        size = tuple(reversed(cbf_handle.get_image_size(0)))

        try:
            underload = find_undefined_value(cbf_handle)
            overload = cbf_handle.get_overload(0)
            trusted_range = (underload, overload)
        except: # intentional
            trusted_range = (0.0, 0.0)

        cbf_detector.__swig_destroy__(cbf_detector)
        del(cbf_detector)

        # Get the sensor type
        dtype = detector_factory.sensor(sensor)

        # If the sensor type is PAD then create the detector with a
        # parallax corrected pixel to millimeter function
        if dtype == detector_helper_sensors.SENSOR_PAD:
            px_mm = ParallaxCorrectedPxMmStrategy(0.252500934883)
        else:
            px_mm = SimplePxMmStrategy()

        return detector_factory.make_detector(
                  dtype, fast, slow, origin, pixel, size,
                  trusted_range, px_mm)

    @staticmethod
    def imgCIF_H(cbf_handle, sensor):
        '''Initialize a detector model from an imgCIF file handle, where it
        is assumed that the file has already been read.'''

        cbf_detector = cbf_handle.construct_detector(0)

        pixel = (cbf_detector.get_inferred_pixel_size(1),
                 cbf_detector.get_inferred_pixel_size(2))

        # FIXME can probably simplify the code which follows below by
        # making proper use of cctbx vector calls - should not be as
        # complex as it appears to be...

        origin = tuple(cbf_detector.get_pixel_coordinates(0, 0))
        fast = cbf_detector.get_pixel_coordinates(0, 1)
        slow = cbf_detector.get_pixel_coordinates(1, 0)

        dfast = [fast[j] - origin[j] for j in range(3)]
        dslow = [slow[j] - origin[j] for j in range(3)]

        lfast = math.sqrt(sum([dfast[j] * dfast[j] for j in range(3)]))
        lslow = math.sqrt(sum([dslow[j] * dslow[j] for j in range(3)]))

        fast = tuple([dfast[j] / lfast for j in range(3)])
        slow = tuple([dslow[j] / lslow for j in range(3)])

        size = tuple(reversed(cbf_handle.get_image_size(0)))

        try:
            underload = find_undefined_value(cbf_handle)
            overload = cbf_handle.get_overload(0)
            trusted_range = (underload, overload)
        except: # intentional
            trusted_range = (0.0, 0.0)

        cbf_detector.__swig_destroy__(cbf_detector)
        del(cbf_detector)

        return detector_factory.make_detector(
                      detector_factory.sensor(sensor),
                      fast, slow, origin, pixel, size, trusted_range)

    @staticmethod
    def sensor(name):
        '''Return the correct sensor token for a given name, for example:

        ccd, CCD
        image_plate, IMAGE_PLATE
        pad, PAD

        to the appropriate static token which will be used as a handle
        everywhere else in this. Also allow existing token to be passed in.'''

        if detector_helper_sensors.check_sensor(name):
            return name

        if name.upper() == 'PAD':
            return detector_helper_sensors.SENSOR_PAD
        elif name.upper() == 'CCD':
            return detector_helper_sensors.SENSOR_CCD
        elif name.upper() == 'IMAGE_PLATE':
            return detector_helper_sensors.SENSOR_IMAGE_PLATE
        elif name.upper() == 'UNKNOWN':
            return detector_helper_sensors.SENSOR_UNKNOWN

        raise RuntimeError, 'name %s not known' % name
