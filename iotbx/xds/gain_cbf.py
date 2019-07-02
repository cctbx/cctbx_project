#!/usr/bin/env libtbx.python
#
# iotbx.xds.gain_cbf.py
#
#   James Parkhurst, Diamond Light Source, 2012/OCT/16
#
#   Class to read the GAIN.CBF files used in XDS
#
from __future__ import absolute_import, division, print_function

from iotbx.xds import xds_cbf
from six.moves import range

class reader(xds_cbf.reader):
  """A class to read the GAIN.CBF files used in XDS"""
  def __init__(self):
    xds_cbf.reader.__init__(self)

  def get_data(self, dim):
    """Get the gain array from the file"""
    from math import ceil
    import numpy

    # Get the image data
    sampled_gain = xds_cbf.reader.get_data(self)
    ny = int(ceil(dim[0] / sampled_gain.shape[0]))
    nx = int(ceil(dim[1] / sampled_gain.shape[1]))

    # Resize the gain image
    i1 = numpy.array([list(range(dim[1]))] * dim[0], dtype=numpy.int32)
    j1 = numpy.array([list(range(dim[0]))] * dim[1], dtype=numpy.int32).transpose()
    i2 = numpy.divide(i1, nx)
    j2 = numpy.divide(j1, ny)
    gain = numpy.zeros(dim, dtype=numpy.float64)
    gain[j1,i1] = sampled_gain[j2,i2] / 1000.0
    return gain

if __name__ == '__main__':
    import sys
    import numpy
    handle = reader()
    handle.read_file(sys.argv[1])
    image = handle.get_data(sys.argv[2])
