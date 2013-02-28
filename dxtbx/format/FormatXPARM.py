#!/usr/bin/env python
# FormatXPARM.py
#   Copyright (C) 2011 Diamond Light Source, James Parkhurst
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Format object for XDS XPARM.XDS files

from __future__ import division

from dxtbx.format.Format import Format
from iotbx.xds import xparm

class FormatXPARM(Format):
    '''An image reading class for XDS XPARM.XDS files'''

    @staticmethod
    def understand(image_file):
        '''Check to see if this looks like an CBF format image, i.e. we can
        make sense of it.'''
        return xparm.reader.is_xparm_file(image_file)

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        Format.__init__(self, image_file)

        assert(self.understand(image_file))
        return

    def _start(self):
        '''Open the image file as a cbf file handle, and keep this somewhere
        safe.'''

        # Open and read the xparm file
        xparm_handle = xparm.reader()
        xparm_handle.read_file(self._image_file)

        # Convert the parameters to cbf conventions
        self._convert_to_cbf_convention(xparm_handle)

    def _convert_to_cbf_convention(self, xparm_handle):
        '''Get the parameters from the XPARM file and convert them to CBF
        conventions.

        Params:
            xparm_handle The handle to the xparm file.

        '''
        from dxtbx.model.detector_helpers import compute_frame_rotation
        from scitbx import matrix
        import math

        # fetch out the interesting things that we need
        detector_distance = xparm_handle.detector_distance
        rotation_axis = matrix.col(xparm_handle.rotation_axis)
        beam_vector = matrix.col(xparm_handle.beam_vector)
        detector_normal = matrix.col(xparm_handle.detector_normal)
        image_size = xparm_handle.detector_size
        pixel_size = xparm_handle.pixel_size
        pixel_origin = xparm_handle.detector_origin
        detector_centre = (pixel_size[0] * pixel_origin[0],
                           pixel_size[1] * pixel_origin[1])
        detector_fast = matrix.col(xparm_handle.detector_x_axis)
        detector_slow = matrix.col(xparm_handle.detector_y_axis)

        # compute the detector origin
        origin_xds = detector_distance * detector_normal
        origin = origin_xds - (detector_centre[0] * detector_fast +
                               detector_centre[1] * detector_slow)

        # then convert directions to unit vectors
        rotation_axis = rotation_axis / math.sqrt(rotation_axis.dot())
        beam_vector = beam_vector / math.sqrt(beam_vector.dot())

        # want to now calculate a rotation which will align the axis with
        # the (1, 0, 0) vector, and the component of the beam vector
        # perpendicular to this with (0, 0, 1) - for convenience then start
        # by defining our current reference frame
        x = rotation_axis
        z = beam_vector - (beam_vector.dot(rotation_axis) * rotation_axis)
        z = z / math.sqrt(z.dot())
        y = z.cross(x)

        # and the target reference frame
        _x = matrix.col([1, 0, 0])
        _y = matrix.col([0, 1, 0])
        _z = matrix.col([0, 0, 1])

        # now compute the rotations which need to be applied to move from
        _m = compute_frame_rotation((x, y, z), (_x, _y, _z))

        # rotate all of the parameters from the XDS to CBF coordinate frame
        self._detector_origin = _m * origin
        self._rotation_axis   = _m * rotation_axis
        self._beam_vector     = _m * beam_vector
        self._fast_axis       = _m * detector_fast
        self._slow_axis       = _m * detector_slow

        # Copy other parameters to class members
        self._image_size = image_size
        self._pixel_size = pixel_size
        self._wavelength = xparm_handle.wavelength
        self._starting_angle = xparm_handle.starting_angle
        self._oscillation_range = xparm_handle.oscillation_range
        self._starting_frame = xparm_handle.starting_frame

    def _goniometer(self):
        '''Return a working goniometer instance.'''
        return self._goniometer_factory.known_axis(self._rotation_axis)

    def _detector(self):
        '''Return a working detector instance.'''
        return self._detector_factory.complex(
            self._detector_factory.sensor('unknown'), self._detector_origin,
            self._fast_axis, self._slow_axis, self._pixel_size,
            self._image_size, (0, 0))

    def _beam(self):
        '''Return a working beam instance.'''
        return self._beam_factory.simple_directional(
            self._beam_vector, self._wavelength)

    def _scan(self):
        '''Return a working scan instance.'''
        import os
        from dxtbx.model.scan_helpers import scan_helper_image_formats

        # Set the scan parameters
        image_range = (self._starting_frame, self._starting_frame)
        oscillation = (self._starting_angle, self._oscillation_range)
        template = '#'
        directory = os.path.dirname(self._image_file)
        format = scan_helper_image_formats.FORMAT_CBF

        # Create the scan object
        return self._scan_factory.make_scan(template, directory, format,
            image_range, 0.0, oscillation, [0], deg=True)

    def get_raw_data(self):
        '''Get the raw image data. For GXPARM.XDS file raise am exception.'''
        raise IOError("GXPARM.XDS does not support image data!")

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatXPARM.understand(arg)
