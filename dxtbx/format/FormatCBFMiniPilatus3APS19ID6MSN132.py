#!/usr/bin/env python
# FormatCBFMiniPilatus3AOS19ID6MSN132.py
#   Copyright (C) 2016 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus

class FormatCBFMiniPilatus3AOS19ID6MSN132(FormatCBFMiniPilatus):
  '''A class for reading mini CBF format Pilatus3 images for 6M SN 132 @ APS19ID.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'PILATUS3' in record and 'S/N 60-0132' in header:
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatCBFMiniPilatus.__init__(self, image_file, **kwargs)

    return

  def _goniometer(self):
    '''19ID has reversed goniometer'''

    return self._goniometer_factory.single_axis_reverse()

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatCBFMiniPilatus3AOS19ID6MSN132.understand(arg))
