#!/usr/bin/env python
# FormatSMVADSCDBG.py
#   Copyright (C) 2013 LBNL, Aaron Brewster
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for pseudo "ADSC" images, converted
# from Pilatus images using iotbx debug_write.

from __future__ import absolute_import, division
from __future__ import print_function

from dxtbx.format.FormatSMVADSC import FormatSMVADSC

class FormatSMVADSCDBG(FormatSMVADSC):
  ''' Format class for reading images converted from pilatus to smv adsc'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this has the correct dimensions of a pilatus'''

    size, header = FormatSMVADSC.get_smv_header(image_file)

    if int(header['SIZE1']) == 2463 and int(header['SIZE2']) == 2527:
      return True

    return False

  def _start(self):
    # read the headers, then swap size1 and size2
    FormatSMVADSC._start(self)
    self._header_dictionary['SIZE1'] = '2527'
    self._header_dictionary['SIZE2'] = '2463'

  def detectorbase_start(self):
    from iotbx.detectors import SMVImage
    self.detectorbase = SMVImage(self._image_file)
    self.detectorbase.open_file = self.open_file
    self.detectorbase.readHeader()
    self.detectorbase.parameters['SIZE1'] = 2527
    self.detectorbase.parameters['SIZE2'] = 2463

  def _detector(self):
    '''Return a model for a simple detector, presuming no one has
    one of these on a two-theta stage. Assert that the beam centre is
    provided in the Mosflm coordinate frame.'''

    distance = float(self._header_dictionary['DISTANCE'])
    beam_x = float(self._header_dictionary['BEAM_CENTER_X'])
    beam_y = float(self._header_dictionary['BEAM_CENTER_Y'])
    pixel_size = float(self._header_dictionary['PIXEL_SIZE'])
    # size1 and size2 swapped here
    image_size = (float(self._header_dictionary['SIZE2']),
                  float(self._header_dictionary['SIZE1']))
    overload = 65535
    underload = 0

    return self._detector_factory.simple(
        'CCD', distance, (beam_y, beam_x), '+x', '-y',
        (pixel_size, pixel_size), image_size, (underload, overload), [])

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatSMVADSCDBG.understand(arg))
