#!/usr/bin/env python
# FormatCBFFromTiff.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Base implementation of fullCBF format - as used with Dectris detectors
# amongst others - this will read the header and construct the full model,
# but will allow for extension for specific implementations of CBF.

from __future__ import absolute_import, division

import pycbf

from dxtbx.format.FormatCBF import FormatCBF

class FormatCBFFromTiff(FormatCBF):
  '''An image reading class for full CBF format images i.e. those from
  a variety of cameras which support this format.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    header = FormatCBF.get_cbf_header(image_file)

    if not 'data_fromtiff' in header and not '_array_data.data' in header:
      return False

    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    self._cbf_handle = None

    # It appears Pycbf can not handle unicode filenames (see dials/dials#256)
    image_file = str(image_file)
    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatCBF.__init__(self, image_file, **kwargs)

    return

  def __del__(self):
    if self._cbf_handle:
      self._cbf_handle.__swig_destroy__(self._cbf_handle)

  def _start(self):
    '''Open the image file as a cbf file handle, and keep this somewhere
    safe.'''

    FormatCBF._start(self)

  def _get_cbf_handle(self):
    if self._cbf_handle is not None:
      return self._cbf_handle
    else:
      self._cbf_handle = pycbf.cbf_handle_struct()
      self._cbf_handle.read_widefile(self._image_file, pycbf.MSG_DIGEST)
      return self._cbf_handle

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
    flex array.'''
    from scitbx.array_family import flex
    import numpy
    import pycbf
    self._get_cbf_handle()
    self._cbf_handle.rewind_datablock()
    self._cbf_handle.select_datablock(0)
    self._cbf_handle.select_category(0)
    self._cbf_handle.select_column(0)
    self._cbf_handle.select_row(0)
    dtype = self._cbf_handle.get_typeofvalue()
    assert dtype.find('bnry') > -1
    image_string = self._cbf_handle.get_integerarray_as_string()
    image = numpy.fromstring(image_string, numpy.int16)
    parameters = self._cbf_handle.get_integerarrayparameters_wdims()
    image_size = (parameters[10], parameters[9])
    image.shape = (image_size)
    image = image.astype(numpy.int32)
    return flex.int(image)

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFFromTiff.understand(arg)
