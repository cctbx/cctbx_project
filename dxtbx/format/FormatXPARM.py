#!/usr/bin/env python
# FormatXPARM.py
#   Copyright (C) 2011 Diamond Light Source, James Parkhurst
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Base implementation of CBF formats - which is just really a place holder
# which will tell you whether something is a CBF file (or no.)

from __future__ import division

from dxtbx.format.Format import Format
from iotbx.xds import xparm

class FormatXPARM(Format):
    '''An image reading class for full CBF format images i.e. those from
    a variety of cameras which support this format.'''

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

        self._xparm_handle = xparm.reader()
        self._xparm_handle.read_file(self._image_file)

    def _goniometer(self):
        '''Return a working goniometer instance.'''
        return self._goniometer_factory.known_axis(
            self._xparm_handle.rotation_axis)

    def _detector(self):
        '''Return a working detector instance.'''
        return self._detector_factory.XDS(self._image_file)

    def _beam(self):
        '''Return a working beam instance.'''
        return self._beam_factory.simple_directional(
            self._xparm_handle.beam_vector,
            self._xparm_handle.wavelength)

    def _scan(self):
        '''Return a working scan instance.'''
        import os
        from dxtbx.model.scan_helpers import scan_helper_image_formats
        starting_frame = self._xparm_handle.starting_frame
        starting_angle = self._xparm_handle.starting_angle
        oscillation_range = self._xparm_handle.oscillation_range
        image_range = (starting_frame, starting_frame)
        oscillation = (starting_angle, oscillation_range)
        template = '#'
        directory = os.path.dirname(self._image_file)
        format = scan_helper_image_formats.FORMAT_CBF
        return self._scan_factory.make_scan(template, directory, format,
            image_range, 0.0, oscillation, [0], deg=True)

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatXPARM.understand(arg)
