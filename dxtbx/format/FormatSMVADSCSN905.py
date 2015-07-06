#!/usr/bin/env python
# FormatSMVADSCSN905.py
#   Copyright (C) 2015 Diamond Light Source, Richard Gildea
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for ADSC images. Inherits from
# FormatSMVADSC, customised for example on ALS beamline 8.2.2 from back in the
# day which had it's own way of recording beam centre.

from __future__ import division

from dxtbx.format.FormatSMVADSCSN import FormatSMVADSCSN

class FormatSMVADSCSN445(FormatSMVADSCSN):
  '''A class for reading SMV format ADSC images, and correctly constructing
  a model for the experiment from this, for instrument number 905.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this is ADSC SN 905.'''

    # check this is detector serial number 905

    size, header = FormatSMVADSCSN.get_smv_header(image_file)

    if int(header['DETECTOR_SN']) != 905:
      return False

    return True

  def __init__(self, image_file):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    FormatSMVADSCSN.__init__(self, image_file)

    return

  def _detector(self):
    '''Return a model for a simple detector, presuming no one has
    one of these on a two-theta stage. Assert that the beam centre is
    provided in the Mosflm coordinate frame.'''

    distance = float(self._header_dictionary['DISTANCE'])
    beam_x = float(self._header_dictionary['DENZO_XBEAM'])
    beam_y = float(self._header_dictionary['DENZO_YBEAM'])
    pixel_size = float(self._header_dictionary['PIXEL_SIZE'])
    image_size = (float(self._header_dictionary['SIZE1']),
                  float(self._header_dictionary['SIZE2']))
    image_pedestal = 40

    overload = 65535 - image_pedestal
    underload = 1 - image_pedestal

    return self._detector_factory.simple(
        'CCD', distance, (beam_y, beam_x), '+x', '-y',
        (pixel_size, pixel_size), image_size, (underload, overload), [])


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatSMVADSC.understand(arg)

