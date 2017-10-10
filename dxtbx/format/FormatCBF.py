#!/usr/bin/env python
# FormatCBF.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Base implementation of CBF formats - which is just really a place holder
# which will tell you whether something is a CBF file (or no.)

from __future__ import absolute_import, division
from __future__ import print_function

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
    of simplicity will read the file in 4k chunks.'''

    fin = FormatCBF.open_file(image_file, 'rb')

    header = fin.read(4096)
    # FIXME this is grim as it is searching over longer and longer
    # files
    while not '--CIF-BINARY-FORMAT-SECTION--' in header:
      add = fin.read(4096)
      if add:
        header += add
      else:
        break
    return header.split('--CIF-BINARY-FORMAT-SECTION--')[0]

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    Format.__init__(self, image_file, **kwargs)

    return

  def _start(self):
    '''Open the image file, read the image header, copy it into memory
    for future inspection.'''

    Format._start(self)

    self._cif_header = FormatCBF.get_cbf_header(self._image_file)

    self._mime_header = ''

    in_binary_format_section = False

    for record in FormatCBF.open_file(self._image_file, 'rb'):
      if '--CIF-BINARY-FORMAT-SECTION--' in record:
        in_binary_format_section = True
      elif in_binary_format_section and record[0] == 'X':
        self._mime_header += record
      if in_binary_format_section and len(record.strip()) == 0:
        # http://sourceforge.net/apps/trac/cbflib/wiki/ARRAY_DATA%20Category
        #    In an imgCIF file, the encoded binary data begins after
        #    the empty line terminating the header.
        break

    return

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatCBF.understand(arg))
