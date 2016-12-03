#!/usr/bin/env python
# FormatSMVADSC.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for ADSC images. Inherits from
# FormatSMV.

from __future__ import division

import time

from dxtbx.format.FormatSMV import FormatSMV

class FormatSMVADSC(FormatSMV):
  '''A class for reading SMV format ADSC images, and correctly constructing
  a model for the experiment from this.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an ADSC SMV format image, i.e. we
    can make sense of it. Essentially that will be if it contains all of
    the keys we are looking for and not some we are not (i.e. that belong
    to a Rigaku Saturn.)'''

    size, header = FormatSMV.get_smv_header(image_file)

    # do not understand JHSim images
    if header.get('BEAMLINE') == 'fake': return False

    # this used to include TIME
    wanted_header_items = ['BEAM_CENTER_X', 'BEAM_CENTER_Y',
                           'DISTANCE', 'WAVELENGTH', 'PIXEL_SIZE',
                           'OSC_START', 'OSC_RANGE', 'SIZE1', 'SIZE2',
                           'BYTE_ORDER']

    for header_item in wanted_header_items:
      if not header_item in header:
        return False

    unwanted_header_items = ['DTREK_DATE_TIME']

    for header_item in unwanted_header_items:
      if header_item in header:
        return False

    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    FormatSMV.__init__(self, image_file, **kwargs)

    return

  def _start(self):

    FormatSMV._start(self)

  def detectorbase_start(self):
    if not hasattr(self, "detectorbase") or self.detectorbase is None:
      from iotbx.detectors import SMVImage
      self.detectorbase = SMVImage(self._image_file)
      self.detectorbase.readHeader()

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header.'''

    return self._goniometer_factory.single_axis()

  def _detector(self):
    '''Return a model for a simple detector, presuming no one has
    one of these on a two-theta stage. Assert that the beam centre is
    provided in the Mosflm coordinate frame.'''

    distance = float(self._header_dictionary['DISTANCE'])
    beam_x = float(self._header_dictionary['BEAM_CENTER_X'])
    beam_y = float(self._header_dictionary['BEAM_CENTER_Y'])
    pixel_size = float(self._header_dictionary['PIXEL_SIZE'])
    image_size = (float(self._header_dictionary['SIZE1']),
                  float(self._header_dictionary['SIZE2']))
    overload = 65535
    underload = 1

    return self._detector_factory.simple(
        'CCD', distance, (beam_y, beam_x), '+x', '-y',
        (pixel_size, pixel_size), image_size, (underload, overload), [])

  def _beam(self):
    '''Return a simple model for the beam.'''

    wavelength = float(self._header_dictionary['WAVELENGTH'])

    return self._beam_factory.simple(wavelength)

  def _scan(self):
    '''Return the scan information for this image.'''
    import calendar

    format = self._scan_factory.format('SMV')
    exposure_time = float(self._header_dictionary['TIME'])
    epoch = None

    # PST, PDT timezones not recognised by default...

    epoch = 0
    try:
      date_str = self._header_dictionary['DATE']
      date_str = date_str.replace('PST', '').replace('PDT', '')
    except KeyError, e:
      date_str = ''
    for format_string in ['%a %b %d %H:%M:%S %Y', '%a %b %d %H:%M:%S %Z %Y']:
      try:
        epoch = calendar.timegm(time.strptime(date_str, format_string))
        break
      except ValueError, e:
        pass

    # assert(epoch)
    osc_start = float(self._header_dictionary['OSC_START'])
    osc_range = float(self._header_dictionary['OSC_RANGE'])

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, epoch)

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
    flex array of integers.)'''

    from boost.python import streambuf
    from dxtbx import read_uint16, read_uint16_bs, is_big_endian
    from scitbx.array_family import flex
    assert(len(self.get_detector()) == 1)
    panel = self.get_detector()[0]
    size = panel.get_image_size()
    f = FormatSMVADSC.open_file(self._image_file, 'rb')
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
    print FormatSMVADSC.understand(arg)
