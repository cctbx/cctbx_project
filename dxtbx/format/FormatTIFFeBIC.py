#!/usr/bin/env python
# FormatTIFFeBIC.py
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
"""Experimental implementation of a format class to recognise images in TIFF
format recorded on a Talos electron microscope at eBIC. This is apparently
an FEI Ceta 16M detector (https://www.fei.com/accessories/ceta-16m/)."""

from __future__ import division
import os
try:
  from PIL import Image
except ImportError:
  import Image
from dxtbx.format.Format import Format
from dxtbx.format.FormatTIFFHelpers import read_basic_tiff_header
from dxtbx.format.FormatTIFFHelpers import tiff_byte_order
from dxtbx.format.FormatTIFFHelpers import LITTLE_ENDIAN
from dxtbx.format.FormatTIFFHelpers import BIG_ENDIAN
from dxtbx.model import ScanFactory

import struct # for struct.error

class FormatTIFFeBIC(Format):
  '''An image reading class for TIFF images from a Talos electron microscope
  at eBIC. We have limited information about the data format at present.

  The header does not contain useful information about the geometry, therefore
  we will construct dummy objects and expect to override on import using
  site.phil.
'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like a TIFF format image of the type seen
    from a microscope at eBIC.'''

    try:
      byte_order = tiff_byte_order(image_file)
    except RuntimeError:
      return False

    if byte_order != LITTLE_ENDIAN: return False

    # check this is *not* something like a Rayonix TIFF
    try:
      read_basic_tiff_header(image_file)
      return False
    except struct.error:
      pass

    # These images do not conform to what is expected by read_basic_tiff_header
    # in FormatTIFFHelpers.py, so we're going to to use PIL to read them.
    im = Image.open(image_file)

    # The TIFF tags used differ between datasets I have seen. The following
    # tags are common to all datasets:
    ImageWidth = 256
    ImageLength = 257
    BitsPerSample = 258
    Compression = 259
    PhotometricInterpretation = 262
    StripOffsets = 273
    RowsPerStrip = 278
    StripByteCounts = 279

    # These tags were in the first dataset, but not others
    SamplesPerPixel = 277
    PlanarConfiguration = 284

    # These tags were in later datasets, but not the first
    NewSubfileType = 254
    XResolution = 282
    YResolution = 283
    ResolutionUnit = 296
    ColorMap = 320
    Helios_metadata = 34682
    unknown1 = 65450
    unknown2 = 65451
    unknown3 = 65452

    # Checks on common tags
    if im.tag[Compression][0] != 1: return False # No compression
    nx = im.tag[ImageWidth][0]
    ny = im.tag[ImageLength][0]
    if (nx, ny) != (4096, 4096): return False # 4096*4096 image
    bytes_per_pixel = im.tag[BitsPerSample][0] // 8
    if bytes_per_pixel not in (1, 2): return False # 8 or 16 bits per pixel

    # Checks on tags in the first dataset
    if im.tag.has_key(SamplesPerPixel):
      if im.tag[SamplesPerPixel][0] != 1: return False # Grayscale, not RGB
    if im.tag.has_key(PlanarConfiguration):
      if im.tag[PlanarConfiguration][0] != 1: return False

    # Checks on tags in later datasets
    if im.tag.has_key(Helios_metadata):
      if 'Microscope TalosArctica 200 kV D3594 CryoTwin' not in im.tag[
        Helios_metadata]: return False

    # It is not necessary for TIFF to store the whole image contiguously.
    # Here at least ensure the last byte of the image is the expected number of
    # bytes from the first byte
    first_byte = im.tag[StripOffsets][0]
    last_byte = im.tag[StripOffsets][-1] + im.tag[StripByteCounts][-1]
    nbytes = last_byte - first_byte
    if nbytes != 4096 * 4096 * bytes_per_pixel: return False

    # Ok if we got this far
    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    Format.__init__(self, image_file, **kwargs)

    return

  def detectorbase_start(self): pass

  def _start(self):
    '''Get the byte offset to the image data'''

    # XXX it is slow to open the image in PIL each time to get these values.
    # Better to have sub classes.
    im = Image.open(self._image_file)
    self._data_offset = im.tag[273][0]
    self._bytes_per_pixel = im.tag[258][0] // 8

    return

  def get_raw_data(self):
    '''Get the pixel intensities'''

    from boost.python import streambuf
    from scitbx.array_family import flex
    if self._bytes_per_pixel == 2:
      from dxtbx import read_uint16 as read_pixel
    else:
      from dxtbx import read_uint8 as read_pixel

    f = FormatTIFFeBIC.open_file(self._image_file, 'rb')
    f.seek(self._data_offset)

    raw_data = read_pixel(streambuf(f), 4096*4096)
    image_size = (4096,4096)
    raw_data.reshape(flex.grid(image_size[1], image_size[0]))

    return raw_data

  def _goniometer(self):
    '''Dummy goniometer, 'vertical' as the images are viewed. Not completely
    sure about the handedness yet'''

    return self._goniometer_factory.known_axis((0,-1,0))

  def _detector(self):
    '''Dummy detector'''

    pixel_size = 0.014, 0.014
    image_size = 4096,4096
    trusted_range = (-1, 65535)
    beam_centre = [(p * i) / 2 for p, i in zip(pixel_size, image_size)]
    d = self._detector_factory.simple('PAD', 2440, beam_centre, '+x', '-y',
                              pixel_size, image_size, trusted_range)
    # Not sure what the gain is
    #for p in d: p.set_gain(8)
    return d

  def _beam(self):
    '''Dummy beam, energy 200 keV'''

    wavelength = 0.02508
    return self._beam_factory.simple(wavelength)

  def _scan(self):
    '''Dummy scan for this image'''

    fname = os.path.split(self._image_file)[-1]
    index = int(fname.split('_')[-1].split('.')[0])
    return ScanFactory.make_scan(
                (index, index), 0.0, (0,0.5),
                {index:0})

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:

    print FormatTIFFeBIC.understand(arg)


