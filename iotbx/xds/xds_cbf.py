#!/usr/bin/env libtbx.python
#
# iotbx.xds.xds_cbf.py
#
#   James Parkhurst, Diamond Light Source, 2012/OCT/16
#
#   Class to read the CBF files used in XDS
#
from __future__ import absolute_import, division, print_function

class reader:
  """A class to read the CBF files used in XDS"""
  def __init__(self):
    pass

  def read_file(self, filename):
    """Read the CBF file"""
    import pycbf
    self.cbf_handle = pycbf.cbf_handle_struct()
    self.cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
    self.cbf_handle.rewind_datablock()

  def get_data(self):
    """Get the gain array from the file"""
    import numpy

    # Select the first datablock and rewind all the categories
    self.cbf_handle.select_datablock(0)
    self.cbf_handle.select_category(0)
    self.cbf_handle.select_column(2)
    self.cbf_handle.select_row(0)

    # Check the type of the element to ensure it's a binary
    # otherwise raise an exception
    type = self.cbf_handle.get_typeofvalue()
    if type.find('bnry') > -1:

      # Read the image data into an array
      image_string = self.cbf_handle.get_integerarray_as_string()
      image = numpy.fromstring(image_string, numpy.int32)

      # Get the array parameters
      parameters = self.cbf_handle.get_integerarrayparameters_wdims()
      image_size = (parameters[10], parameters[9])

      # Resize the image
      image.shape = (image_size)

    else:
      raise TypeError('Can\'t find image')

    # Return the image
    return image

if __name__ == '__main__':
    import sys
    import numpy
    handle = reader()
    handle.read_file(sys.argv[1])
    image = handle.get_data()
