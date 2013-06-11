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

class FormatXPARM(Format):
    '''An image reading class for XDS XPARM.XDS files'''

    @staticmethod
    def understand(image_file):
        '''Check to see if this looks like an CBF format image, i.e. we can
        make sense of it.'''
        from rstbx.cftbx.coordinate_frame_helpers import is_recognized_file
        return is_recognized_file(image_file)

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        Format.__init__(self, image_file)

        assert(self.understand(image_file))
        return

    def _start(self):
        '''Open the image file as a cbf file handle, and keep this somewhere
        safe.'''

        # Convert the parameters to cbf conventions
        self._convert_to_cbf_convention(self._image_file)

    def _convert_to_cbf_convention(self, xparm_filename):
        '''Get the parameters from the XPARM file and convert them to CBF
        conventions.

        Params:
            xparm_handle The handle to the xparm file.

        '''
        from rstbx.cftbx.coordinate_frame_converter import \
            coordinate_frame_converter
        from scitbx import matrix

        # Create a coordinate frame converter and extract other quantities
        cfc = coordinate_frame_converter(xparm_filename)
        self._detector_origin = cfc.get('detector_origin')
        self._rotation_axis = cfc.get('rotation_axis')
        self._fast_axis = cfc.get('detector_fast')
        self._slow_axis = cfc.get('detector_slow')
        self._wavelength  = cfc.get('wavelength')
        self._image_size  = cfc.get('detector_size_fast_slow')
        self._pixel_size  = cfc.get('detector_pixel_size_fast_slow')
        self._starting_angle = cfc.get('starting_angle')
        self._oscillation_range = cfc.get('oscillation_range')
        self._starting_frame = cfc.get('starting_frame')
        self._divergence = 0.0
        self._sigma_divergence = cfc.get('sigma_divergence')
        sample_vector = cfc.get('sample_to_source')
        self._beam_vector = tuple(matrix.col(sample_vector))

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
        return self._beam_factory.make_beam(
            sample_to_source=self._beam_vector,
            wavelength=self._wavelength,
            divergence=self._divergence,
            sigma_divergence=self._sigma_divergence)

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
        return self._scan_factory.make_scan(image_range, 0.0,
            oscillation, [0], deg=True)

    def get_raw_data(self):
        '''Get the raw image data. For GXPARM.XDS file raise am exception.'''
        raise IOError("GXPARM.XDS does not support image data!")

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatXPARM.understand(arg)
