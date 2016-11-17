#!/usr/bin/env python
# FormatCBFFull.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Base implementation of fullCBF format - as used with Dectris detectors
# amongst others - this will read the header and construct the full model,
# but will allow for extension for specific implementations of CBF.

from __future__ import division

import pycbf

from dxtbx.format.FormatCBF import FormatCBF

class FormatCBFFull(FormatCBF):
  '''An image reading class for full CBF format images i.e. those from
  a variety of cameras which support this format.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    header = FormatCBF.get_cbf_header(image_file)

    if not '_diffrn.id' in header and not '_diffrn_source' in header:
      return False

    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    # It appears Pycbf can not handle unicode filenames (see dials/dials#256)
    image_file = str(image_file)
    assert(self.understand(image_file))

    FormatCBF.__init__(self, image_file, **kwargs)

    return

  def __del__(self):
    self._cbf_handle.__swig_destroy__(self._cbf_handle)

  def _start(self):
    '''Open the image file as a cbf file handle, and keep this somewhere
    safe.'''

    FormatCBF._start(self)

  def _get_cbf_handle(self):
    try:
      return self._cbf_handle
    except AttributeError:
      self._cbf_handle = pycbf.cbf_handle_struct()
      self._cbf_handle.read_widefile(self._image_file, pycbf.MSG_DIGEST)
      return self._cbf_handle

  def _goniometer(self):
    '''Return a working goniometer instance.'''

    return self._goniometer_factory.imgCIF_H(self._get_cbf_handle())

  def _detector(self):
    '''Return a working detector instance.'''

    return self._detector_factory.imgCIF_H(self._get_cbf_handle(), 'unknown')

  def _beam(self):
    '''Return a working beam instance.'''

    return self._beam_factory.imgCIF_H(self._get_cbf_handle())

  def _scan(self):
    '''Return a working scan instance.'''

    return self._scan_factory.imgCIF_H(self._image_file, self._get_cbf_handle())

  def detectorbase_start(self):
    from iotbx.detectors.cbf import CBFImage
    self.detectorbase = CBFImage(self._image_file)
    self.detectorbase.readHeader()

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
    flex array.'''
    self.detectorbase_start()
    try:
      image = self.detectorbase
      image.read()
      raw_data = image.get_raw_data()

      return raw_data
    except Exception:
      return None


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFFull.understand(arg)
