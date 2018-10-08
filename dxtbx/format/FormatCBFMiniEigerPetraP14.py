#!/usr/bin/env python
# FormatCBFMiniEigerPetraP14.py
#   Copyright (C) 2018 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Eiger images. Inherits from
# FormatCBFMiniEiger.

from __future__ import absolute_import, division, print_function

import os

from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger

class FormatCBFMiniEigerPetraP14(FormatCBFMiniEiger):
  '''A class for reading mini CBF format Eiger images, and correctly
  constructing a model for the experiment from this. This tuned for Petra P14'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Eiger mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniEiger.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# detector' in record.lower() and \
             'eiger' in record.lower() and 'E-32-0107' in record:
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatCBFMiniEiger.__init__(self, image_file, **kwargs)

    self._raw_data = None
    return

  def _goniometer(self):
    if 'Phi' in self._cif_header_dictionary:
      phi_value = float(self._cif_header_dictionary['Phi'].split()[0])

    return self._goniometer_factory.known_axis((0, 1, 0))


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatCBFMiniEigerPetraP14.understand(arg))
