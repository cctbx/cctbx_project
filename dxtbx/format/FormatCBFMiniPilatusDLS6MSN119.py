#!/usr/bin/env python
# FormatCBFMiniPilatusDLS6MSN119.py
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images, from the Pilatus
# 6M SN 119 currently on Diamond I24.

from __future__ import division

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus

import os
if 'P6M_60_PANEL' in os.environ:
  single_panel = False
else:
  single_panel = True

def read_cbf_image(cbf_image):
  from cbflib_adaptbx import uncompress
  import binascii
  #from scitbx.array_family import flex

  start_tag = binascii.unhexlify('0c1a04d5')

  data = open(cbf_image, 'rb').read()
  data_offset = data.find(start_tag) + 4
  cbf_header = data[:data_offset - 4]

  fast = 0
  slow = 0
  length = 0

  for record in cbf_header.split('\n'):
    if 'X-Binary-Size-Fastest-Dimension' in record:
      fast = int(record.split()[-1])
    elif 'X-Binary-Size-Second-Dimension' in record:
      slow = int(record.split()[-1])
    elif 'X-Binary-Number-of-Elements' in record:
      length = int(record.split()[-1])
    elif 'X-Binary-Size:' in record:
      size = int(record.split()[-1])

  assert(length == fast * slow)

  pixel_values = uncompress(packed = data[data_offset:data_offset + size],
                            fast = fast, slow = slow)

  return pixel_values

class FormatCBFMiniPilatusDLS6MSN119(FormatCBFMiniPilatus):
  '''A class for reading mini CBF format Pilatus images for 6M SN 119 @ DLS.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniPilatus.get_cbf_header(image_file)

    year = 0

    for record in header.split('\n'):
      if '# 20' in record:
        year = int(record.replace('-', ' ').replace('/', ' ').split()[1])
        break

    assert (year > 0)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'PILATUS' in record and 'S/N 60-0119' in header:
        if year >= 2015:
          return True
        else:
          return False

    return False

  def __init__(self, image_file):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    FormatCBFMiniPilatus.__init__(self, image_file)

    self._raw_data = None

    return

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header, though for miniCBF
    there are limited options for this.'''

    return self._goniometer_factory.known_axis((0, 1, 0))

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMiniPilatusDLS6MSN119.understand(arg)
