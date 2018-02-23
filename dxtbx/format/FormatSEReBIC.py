#!/usr/bin/env python
# FormatSEReBIC.py
#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Experimental format for TIA .ser files used by FEI microscope at eBIC.

from __future__ import absolute_import, division
from dxtbx.format.FormatSER import FormatSER

class FormatSEReBIC(FormatSER):

  def __init__(self, image_file, **kwargs):

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)
    FormatSER.__init__(self, image_file, **kwargs)

  @staticmethod
  def understand(image_file):

    # check the header looks right
    try:
      fmt = FormatSER(image_file)
    except IncorrectFormatError:
      return False

    # Read metadata and check the image dimension is supported
    d = fmt._read_metadata(image_file)
    image_size = d['ArraySizeX'], d['ArraySizeY']
    if image_size not in [(4096, 4096), (2048, 2048)]: return False

    return True

  def _goniometer(self):
    '''Dummy goniometer, 'vertical' as the images are viewed. Not completely
    sure about the handedness yet'''

    return self._goniometer_factory.known_axis((0,-1,0))

  def _detector(self):
    '''Dummy detector'''

    image_size = (self._header_dictionary['ArraySizeX'],
                  self._header_dictionary['ArraySizeY'])
    if image_size == (4096, 4096):
      pixel_size = 0.014, 0.014
    if image_size == (2048, 2048):
      pixel_size = 0.028, 0.028

    distance = 2000
    trusted_range = (-1, 65535)
    beam_centre = [(p * i) / 2 for p, i in zip(pixel_size, image_size)]
    d = self._detector_factory.simple('PAD', distance, beam_centre, '+x', '-y',
                              pixel_size, image_size, trusted_range)
    # Not sure what the gain is
    #for p in d: p.set_gain(8)
    return d

  def _beam(self):
    '''Dummy unpolarized beam, energy 200 keV'''

    wavelength = 0.02508
    return self._beam_factory.make_polarized_beam(
        sample_to_source=(0.0, 0.0, 1.0),
        wavelength=wavelength,
        polarization=(0, 1, 0),
        polarization_fraction=0.5)

  def _scan(self):
    '''Dummy scan for this stack'''

    nframes = self.get_num_images()
    image_range = (1, nframes)
    exposure_times = 0.0
    oscillation = (0, 0.5)

    # FIXME we do actually have acquisition times, might as well use them
    epochs = [0] * nframes

    return self._scan_factory.make_scan(image_range, exposure_times,
      oscillation, epochs, deg=True)

  def get_raw_data(self, index):

    raw_data = super(FormatSEReBIC, self).get_raw_data(index)

    return raw_data

if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    print FormatSEReBIC.understand(arg)

