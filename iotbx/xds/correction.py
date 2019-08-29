#!/usr/bin/env libtbx.python
#
# iotbx.xds.correction.py
#
#   James Parkhurst, Diamond Light Source, 2012/OCT/16
#
#   Class to read the X/Y-CORRECTIONS.CBF files used in XDS
#
from __future__ import absolute_import, division, print_function
from six.moves import range

class reader:
  """A class to read the X/Y-CORRECTIONS.CBF files used in XDS"""
  def __init__(self):
    pass

  def read_file(self, filename):
    """Read the CBF correction file"""
    import pycbf
    self.cbf_handle = pycbf.cbf_handle_struct()
    self.cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
    self.cbf_handle.rewind_datablock()

  def get_correction_array(self):
    """Get the correction array from the file"""
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

  def get_correction(self, dim):
    """Get the correction at each pixel."""
    import numpy

    # Get the raw array
    raw_array = self.get_correction_array()

    # Ensure out dimensions are ok
    if raw_array.shape[0] * 4 < dim[0] or raw_array.shape[1] * 4 < dim[1]:
      raise ValueError("Dimensions are incompatible")

    # Create the array of the given dimension
    correction = numpy.zeros(dim, dtype=numpy.float64)

    # Loop through all pixels and get the correction
    i1 = numpy.array([list(range(dim[1]))] * dim[0], dtype=numpy.int32)
    j1 = numpy.array([list(range(dim[0]))] * dim[1], dtype=numpy.int32).transpose()
    i2 = numpy.divide(i1, 4)
    j2 = numpy.divide(j1, 4)
    correction[j1,i1] = raw_array[j2,i2] / 10.0

    # Return the array of corrections
    return correction
