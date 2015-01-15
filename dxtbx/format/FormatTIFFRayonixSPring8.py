#!/usr/bin/env python
# FormatTIFFRayonixSPring8.py
#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division

import struct

from dxtbx.format.FormatTIFFRayonix import FormatTIFFRayonix

class FormatTIFFRayonixSPring8(FormatTIFFRayonix):
  '''A class for reading TIFF format Rayonix images, and correctly
  constructing a model for the experiment from this, for SPring-8 beamlines
  with a reversed phi axis.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Rayonix TIFF format image,
    i.e. we can make sense of it. This simply checks that records which
    describe the size of the image match with the TIFF records which do
    the same.'''

    width, height, depth, order, bytes = FormatTIFFRayonix.get_tiff_header(
        image_file)

    serial_number = -1

    for record in bytes[2464:2464+512].strip().split('\n'):
      if 'detector serial number' in record.lower():
        serial_number = int(record.split()[-1])

    # only understand a square image
    if width != height: return

    # BL26B2 MX225 with 1X1, 2X2, 3X3 or 4X4 binning
    if serial_number == 24 and width in [6144, 3072, 2046, 1536]:
      return True

    # BL32XU MX225-HE with 1X1, 2X2, 3X3 or 4X4 binning
    if serial_number == 31 and width in [6144, 3072, 2046, 1536]:
      return True

    # BL44XU MX225-HE with 1X1, 2X2, 3X3 or 4X4 binning
    if serial_number == 38 and width in [6144, 3072, 2046, 1536]:
      return True

    # BL44XU MX300-HE with 1X1, 2X2, 3X3 or 4X4 binning
    if serial_number == 42 and width in [8192, 4096, 2728, 2048]:
      return True

    # BL41XU MX225-HE with 1X1, 2X2, 3X3 or 4X4 binning
    if serial_number == 40 and width in [6144, 3072, 2046, 1536]:
      return True

    # BL32XU MX225-HS 1X1, 2X2, 3X3, 4X4, 5X5, 6X6, 8X8 or 10X10 binning
    if serial_number == 106 and width in [5760, 2880, 1920, 1440, 1152, 960,
                                          720, 576]:
      return True

  def __init__(self, image_file):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))
    FormatTIFFRayonix.__init__(self, image_file)

    return

  def _detector(self):
    '''Return a model for a simple detector, which at the moment insists
    that the offsets and rotations are all 0.0.'''

    starts, ends, offset, width = self._get_rayonix_scan_angles()
    rotations = self._get_rayonix_detector_rotations()

    # assert that two-theta offset is 0.0

    assert(starts[0] == 0.0)
    assert(ends[0] == 0.0)

    # assert that the rotations are all 0.0

    assert(rotations[0] == 0.0)
    assert(rotations[1] == 0.0)
    assert(rotations[2] == 0.0)

    distance = self._get_rayonix_distance()
    beam_x, beam_y = self._get_rayonix_beam_xy()
    pixel_size = self._get_rayonix_pixel_size()
    image_size = self._tiff_width, self._tiff_height
    overload = struct.unpack(
        self._i, self._tiff_header_bytes[1128:1132])[0]
    underload = 0

    beam = beam_x * pixel_size[0], beam_y * pixel_size[1]

    return self._detector_factory.simple(
        'CCD', distance, beam, '+x', '-y', pixel_size,
        image_size, (underload, overload), [])

  def _goniometer(self):
    '''Return a model for goniometer corresponding to the values stored
    in the image header. In the first instance assume this is a single
    axis and raise exception otherwise.'''

    starts, ends, offset, width = self._get_rayonix_scan_angles()

    return self._goniometer_factory.single_axis_reverse()


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatTIFFRayonixSPring8BL26B2.understand(arg)
