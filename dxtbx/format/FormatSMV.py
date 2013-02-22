#!/usr/bin/env python
# FormatSMV.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Implementation of an ImageFormat class to read SMV format image but not -
# in the first instance - actually provide a full image representation. This
# is simply there to set everything up for the ADSC and Rigaku Saturn image
# readers which really will acquire the full image including header information
# and generate the experimental model representations.

from __future__ import division

from dxtbx.format.Format import Format

class FormatSMV(Format):
    '''An image reading class for SMV format images i.e. those from ADSC and
    Rigaku which start with:

    {
    HEADER_BYTES=  512;

    and contain a list of keyword-value pairs thereafter which define the
    header. The keywords etc. for these will depend on the instrument
    manufacturer and will be interpreted by subclasses of this class. Note
    also that every line is finished with a semicolon.'''

    @staticmethod
    def understand(image_file):
        '''Check to see if this looks like an SMV format image, i.e. we can
        make sense of it.'''

        if FormatSMV.open_file(image_file, 'rb').read(15) == \
               '{\nHEADER_BYTES=':
            return True

        return False

    @staticmethod
    def get_smv_header(image_file):
        header_size = int(FormatSMV.open_file(image_file, 'rb').read(45).split(
            '\n')[1].split('=')[1].replace(';', '').strip())
        header_text = FormatSMV.open_file(image_file, 'rb').read(header_size)

        # check we have the whole header in here... it is contained within { }

        assert('}' in header_text)

        header_dictionary = { }

        for record in header_text.split('\n'):
            if not '=' in record:
                continue

            key, value = record.replace(';', '').split('=')

            header_dictionary[key.strip()] = value.strip()

        return header_size, header_dictionary

    def __init__(self, image_file):
        '''Initialise the image structure from the given file.'''

        assert(self.understand(image_file))

        Format.__init__(self, image_file)

        return

    def _start(self):
        '''Open the image file, read the image header, copy the key / value
        pairs into an internal dictionary self._header_dictionary along with
        the length of the header in bytes self._header_size.'''

        self._header_size, self._header_dictionary = FormatSMV.get_smv_header(
            self._image_file)

        return
