#!/usr/bin/env python
# FormatSMVADSCSN457.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for ADSC images. Inherits from
# FormatSMVADSC, customised for example on Australian Synchrotron SN 457
# which has reversed phi.

from __future__ import division

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN

class FormatSMVADSCSN457(FormatSMVADSCSN):
    '''A class for reading SMV format ADSC images, and correctly constructing
    a model for the experiment from this, for instrument number 457.'''

    @staticmethod
    def understand(image_file):
        '''Check to see if this is ADSC SN 457.'''

        # check this is detector serial number 457

        size, header = FormatSMVADSCSN.get_smv_header(image_file)

        if int(header['DETECTOR_SN']) != 457:
            return False

        return True

    def __init__(self, image_file):
        '''Initialise the image structure from the given file, including a
        proper model of the experiment.'''

        assert(self.understand(image_file))

        FormatSMVADSCSN.__init__(self, image_file)

        return

    def _goniometer(self):
        '''Return a model for a simple single-axis goniometer. This should
        probably be checked against the image header.'''

        return self._goniometer_factory.known_axis((-1, 0, 0))

    def _scan(self):
        '''Return the scan information for this image. There may be
        no timestamps in there...'''

        format = self._scan_factory.format('SMV')
        exposure_time = float(self._header_dictionary['TIME'])
        epoch = 0
        osc_start = float(self._header_dictionary['OSC_START'])
        osc_range = float(self._header_dictionary['OSC_RANGE'])

        return self._scan_factory.single(
            self._image_file, format, exposure_time,
            osc_start, osc_range, epoch)

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatSMVADSC.understand(arg)
