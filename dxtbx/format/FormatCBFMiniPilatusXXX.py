#!/usr/bin/env python
# FormatCBFMiniPilatusDLS6MSN100.py
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images, from the Pilatus
# 6M SN 100 currently on Diamond I04.

from __future__ import division

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
#from dxtbx.model import ParallaxCorrectedPxMmStrategy
#from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

class FormatCBFMiniPilatusXXX(FormatCBFMiniPilatus):
  '''A class for reading mini CBF format Pilatus images for 6M SN XXX.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'PILATUS' in record and 'S/N XX-XXX' in header:
        return True

    return False

  def __init__(self, image_file):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    FormatCBFMiniPilatus.__init__(self, image_file)

    return

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header, though for miniCBF
    there are limited options for this.'''

    if 'Phi' in self._cif_header_dictionary:
      phi_value = float(self._cif_header_dictionary['Phi'].split()[0])

    return self._goniometer_factory.single_axis_reverse()


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMiniPilatusDLS6MSN100F.understand(arg)
