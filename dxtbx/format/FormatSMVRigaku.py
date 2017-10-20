#!/usr/bin/env python
# FormatSMVRigaku.py
#   Copyright (C) 2013 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for Rigaku images.
# Inherits from FormatSMV.

from __future__ import absolute_import, division

from dxtbx.format.FormatSMV import FormatSMV

class FormatSMVRigaku(FormatSMV):
  '''A class for reading SMV format Rigaku images. 'Abstract' class.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like a Rigaku d*TREK SMV format image,
    i.e. we can make sense of it. Essentially that will be if it contains
    all of the keys we are looking for.'''

    size, header = FormatSMV.get_smv_header(image_file)

    wanted_header_items = [
        'DETECTOR_NUMBER', 'DETECTOR_NAMES',
        'BYTE_ORDER', 'DIM', 'SIZE1', 'SIZE2', 'Data_type'
    ]

    for header_item in wanted_header_items:
      if not header_item in header:
        return False

    # code around CMOS1 on ALS 4.2.2

    detector_prefixes = header['DETECTOR_NAMES'].split()

    if len(detector_prefixes) == 1:
      prefix = detector_prefixes[0]
      det_desc = '%sDETECTOR_DESCRIPTION' % prefix
      if 'CMOS-1' in header.get(det_desc, ''):
        return False

    detector_prefixes = header['DETECTOR_NAMES'].split()
    try:
      detector_number = int(header['DETECTOR_NUMBER'].strip())
    except (KeyError,AttributeError,ValueError),e:
      return False

    if detector_number != len(detector_prefixes):
      return False

    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatSMV.__init__(self, image_file, **kwargs)

  def _goniometer(self):
    '''Overload this method to read the image file however you like so
    long as the result is an goniometer.'''

    raise NotImplementedError('overload me')

  def _detector(self):
    '''Overload this method to read the image file however you like so
    long as the result is an detector.'''

    raise NotImplementedError('overload me')

  def _beam(self):
    '''Overload this method to read the image file however you like so
    long as the result is an beam.'''

    raise NotImplementedError('overload me')

  def _scan(self):
    '''Overload this method to read the image file however you like so
    long as the result is an scan.'''

    raise NotImplementedError('overload me')

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
    flex array of integers.)'''

    from boost.python import streambuf
    from dxtbx import read_uint16, read_uint16_bs, is_big_endian
    from scitbx.array_family import flex
    assert(len(self.get_detector()) == 1)
    panel = self.get_detector()[0]
    size = panel.get_image_size()
    f = FormatSMVRigaku.open_file(self._image_file, 'rb')
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
    print FormatSMVRigaku.understand(arg)
