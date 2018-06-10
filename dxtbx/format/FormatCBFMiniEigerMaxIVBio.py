#!/usr/bin/env python
#  Author: Nick Sauter.
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the FormatCBFMiniEiger image reader for the Eiger16M
# detector at the MaxIV BioMAX beamline, which has a vertical goniometer.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger

class FormatCBFMiniEigerMaxIVBio(FormatCBFMiniEiger):
  '''A class for reading mini CBF format Eiger16M images for S/N E-32-105 MaxIV BioMAX,
  which has a vertical goniometer axis.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Eiger mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniEiger.get_cbf_header(image_file)

    for record in header.split('\n'):
      if 'Detector: Dectris Eiger 16M, S/N E-32-0105' in record:
        return True

    return False

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header, though for miniCBF
    there are limited options for this.'''

    return self._goniometer_factory.known_axis((0, 1, 0))


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatCBFMiniEigerMaxIVBio.understand(arg))
