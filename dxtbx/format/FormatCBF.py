#!/usr/bin/env python
# FormatCBF.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Base implementation of CBF formats - which is just really a place holder
# which will tell you whether something is a CBF file (or no.)

from __future__ import division

from dxtbx.format.Format import Format

class FormatCBF(Format):
    '''An image reading class for CBF format images i.e. those from Dectris
    amongst others. This is just a first base class which will be used to
    determine whether this is really a CBF file.'''

    @staticmethod
    def understand(image_file):
        '''Check to see if this looks like an CBF format image, i.e. we can
        make sense of it.'''

        if '###CBF' in FormatCBF.open_file(image_file, 'rb').read(6):
            return True

        return False

    @staticmethod
    def get_cbf_header(image_file):
        '''Obtain the text section of the header, which is assumed to be
        everything before --CIF-BINARY-FORMAT-SECTION-- - N.B. for reasons
        of simplicity will read the file in 1k chunks.'''

        fin = FormatCBF.open_file(image_file, 'rb')

        header = fin.read(1024)

        while not '--CIF-BINARY-FORMAT-SECTION--' in header:
            header += fin.read(1024)

        return header.split('--CIF-BINARY-FORMAT-SECTION--')[0]

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        assert(self.understand(image_file))

        Format.__init__(self, image_file)

        return

    def _start(self):
        '''Open the image file, read the image header, copy it into memory
        for future inspection.'''

        Format._start(self)

        self._cif_header = FormatCBF.get_cbf_header(self._image_file)

        self._mime_header = ''

        for record in FormatCBF.open_file(self._image_file, 'rb'):
            if 'X' in record[:1]:
                self._mime_header += record
            if 'X-Binary-Size-Padding' in record:
                break

        return
