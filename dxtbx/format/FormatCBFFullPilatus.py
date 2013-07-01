#!/usr/bin/env python
# FormatCBFFullPilatus.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Pilatus implementation of fullCBF format, for use with Dectris detectors.

from __future__ import division

import pycbf

from dxtbx.format.FormatCBFFull import FormatCBFFull
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

class FormatCBFFullPilatus(FormatCBFFull):
    '''An image reading class for full CBF format images from Pilatus
    detectors.'''

    @staticmethod
    def understand(image_file):
        '''Check to see if this looks like an CBF format image, i.e. we can
        make sense of it.'''

        header = FormatCBFFull.get_cbf_header(image_file)

        for record in header.split('\n'):
            if '_array_data.header_convention' in record and \
                   'PILATUS' in record:
                return True

        return False

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        assert(self.understand(image_file))

        FormatCBFFull.__init__(self, image_file)

        return

    def _start(self):
        '''Open the image file as a cbf file handle, and keep this somewhere
        safe.'''

        FormatCBFFull._start(self)

        self._cbf_handle = pycbf.cbf_handle_struct()
        self._cbf_handle.read_widefile(self._image_file, pycbf.MSG_DIGEST)

        return

    def _detector(self):
        '''Return a working detector instance, with added mask regions.'''

        detector = self._detector_factory.imgCIF_H(self._cbf_handle,
                                                   'PAD')

        for f0, s0, f1, s1 in determine_pilatus_mask(detector):
            detector.add_mask(f0, s0, f1, s1)

        return detector

if __name__ == '__main__':

    import sys

    for arg in sys.argv[1:]:
        print FormatCBFFullPilatus.understand(arg)
