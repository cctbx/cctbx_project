#!/usr/bin/env python
# FormatSMVFakeADSC.py
#   Copyright (C) 2017 Diamond Light Source, James Parkhurst
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for fake ADSC images converted from
# electron diffraction images

from __future__ import division

import time
from scitbx import matrix

from dxtbx.format.FormatSMV import FormatSMV

class FormatSMVFakeADSC(FormatSMV):
  '''A class for reading SMV format fake ADSC electron diffraction images'''

  @staticmethod
  def understand(image_file):

    size, header = FormatSMV.get_smv_header(image_file)

    # do not understand JHSim images
    if header.get('BEAMLINE') == 'fake': return False

    # do not understand Timepix_SU images
    if header.get('BEAMLINE') == 'TimePix_SU': return False

    wanted_header_items = [
      "HEADER_BYTES",
      "BEAM_CENTER_X",
      "BEAM_CENTER_Y",
      "BIN",
      "BYTE_ORDER",
      "DATE",
      "DETECTOR_SN",
      "DIM",
      "DISTANCE",
      "OSC_RANGE",
      "PHI",
      "SIZE1",
      "SIZE2",
      "PIXEL_SIZE",
      "TIME",
      "WAVELENGTH",
      "TWOTHETA",
      "TYPE"
        ]

    for header_item in wanted_header_items:
      if not header_item in header:
        return False

    return True

  def __init__(self, image_file, **kwargs):

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatSMV.__init__(self, image_file, **kwargs)

    return

  # def _start(self):

  #   FormatSMV._start(self)

  # def detectorbase_start(self):
  #   from iotbx.detectors.saturn import SaturnImage
  #   self.detectorbase = SaturnImage(self._image_file)
  #   self.detectorbase.readHeader()

  def _goniometer(self):
    return self._goniometer_factory.single_axis()

  def _detector(self):

    beam_centre_x = float(self._header_dictionary['BEAM_CENTER_X'])
    beam_centre_y = float(self._header_dictionary['BEAM_CENTER_Y'])
    distance = float(self._header_dictionary['DISTANCE'])
    image_size_x = int(self._header_dictionary['SIZE1'])
    image_size_y = int(self._header_dictionary['SIZE2'])
    pixel_size = float(self._header_dictionary['PIXEL_SIZE'])

    fast_direction = "+x"
    slow_direction = "+y"

    # Just putting in values
    underload = -1
    overload = 1e6

    return self._detector_factory.simple(
        'CCD',
        distance,
        (beam_centre_x, beam_centre_y),
        fast_direction,
        slow_direction,
        (pixel_size, pixel_size),
        (image_size_x, image_size_y),
        (underload, overload))

  def _beam(self):
    '''Return a simple model for the beam.'''

    wavelength = float(self._header_dictionary['WAVELENGTH'])

    return self._beam_factory.simple(wavelength)

  def _scan(self):
    '''Return the scan information for this image.'''
    import calendar

    osc_range =float(self._header_dictionary['OSC_RANGE'])
    osc_start = float(self._header_dictionary['PHI'])


    return self._scan_factory.single(
        self._image_file, None, 0,
        osc_start, osc_range, 0)

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
    flex array of integers.)'''

    from boost.python import streambuf
    from dxtbx import read_uint16, read_uint16_bs, is_big_endian
    from scitbx.array_family import flex
    assert(len(self.get_detector()) == 1)
    panel = self.get_detector()[0]
    size = panel.get_image_size()
    f = FormatSMV.open_file(self._image_file, 'rb')
    f.read(self._header_size)

    if self._header_dictionary['BYTE_ORDER'] == 'big_endian':
      big_endian = True
    else:
      big_endian = False

    if big_endian == is_big_endian():
      raw_data = read_uint16(streambuf(f), int(size[0] * size[1]))
    else:
      raw_data = read_uint16_bs(streambuf(f), int(size[0] * size[1]))

    image_size = panel.get_image_size()
    raw_data.reshape(flex.grid(image_size[1], image_size[0]))

    return raw_data

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatSMVFakeADSC.understand(arg)

  fmt = FormatSMVFakeADSC(arg)
  print fmt.get_beam()
  print fmt.get_detector()
  print fmt.get_goniometer()
  print fmt.get_scan()
  print fmt.get_raw_data()
