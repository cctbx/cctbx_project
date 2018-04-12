#!/usr/bin/env python
# FormatTIFFRayonixESRF.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Sub class of FormatTiffRayonix to deal with images who have beam centers
# specified in pixels
#

from __future__ import absolute_import, division, print_function

import struct

from dxtbx.format.FormatTIFFRayonix import FormatTIFFRayonix

class FormatTIFFRayonixESRF(FormatTIFFRayonix):
  '''A class for reading TIFF format Rayonix images, and correctly
  constructing a model for the experiment from this.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Rayonix TIFF format image,
    i.e. we can make sense of it.  Returns true if the beam center is specified
    in pixels.'''

    width, height, depth, order, bytes = FormatTIFFRayonix.get_tiff_header(
        image_file)

    import struct
    from scitbx.matrix import col
    from dxtbx.format.FormatTIFFHelpers import LITTLE_ENDIAN, BIG_ENDIAN
    format = {LITTLE_ENDIAN:'<', BIG_ENDIAN:'>'}[order]
    offset = 1024

    detector_size_pixels = col(struct.unpack(format+'ii',bytes[offset+80:offset+88]))
    detector_center_px   = 0.5 * detector_size_pixels

    detector_pixel_sz_mm =  1.E-6 * col( # convert from nano to milli
                                struct.unpack(format+'ii',bytes[offset+772:offset+780]))

    header_beam_center = 0.001 * col( # Rayonix says this should be pixels
                                struct.unpack(format+'ii',bytes[offset+644:offset+652]))

    disagreement = header_beam_center[0]/detector_center_px[0]
    return header_beam_center[0] > 0 and header_beam_center[1] > 0 \
           and disagreement < 0.5  # if header was in mm, disagreement should be
                                   # approximately the pixel size in mm

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)
    FormatTIFFRayonix.__init__(self, image_file, **kwargs)

    return

  def _goniometer(self):
    '''Return a model for goniometer corresponding to the values stored
    in the image header. In the first instance assume this is a single
    axis annd raise exception otherwise.'''

    starts, ends, offset, width = self._get_rayonix_scan_angles()

    # not testing as this the CLS images are not properly structured...
    # and also don't have a serial number in (FAIL)

    return self._goniometer_factory.single_axis()


  ####################################################################
  #                                                                  #
  # Helper methods to get all of the values out of the TIFF header   #
  # - separated out to assist with code clarity                      #
  #                                                                  #
  ####################################################################

  def _get_rayonix_beam_xy(self):
    '''Get the beam x, y positions which are defined in the standard
    to be in pixels. X and Y are not defined by the documentation, so
    taking as a given that these are horizontal and vertical. N.B.
    the documentation states that the horizontal direction is fast.'''

    beam_x, beam_y = struct.unpack(
        self._ii, self._tiff_header_bytes[1668:1676])[:2]
    pixel_x, pixel_y = struct.unpack(
        self._ii, self._tiff_header_bytes[1796:1804])[:2]

    return beam_x * 1000 / pixel_x, beam_y * 1000 / pixel_y

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatTIFFRayonixESRF.understand(arg))
