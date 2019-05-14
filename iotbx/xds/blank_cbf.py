#!/usr/bin/env libtbx.python
#
# iotbx.xds.xds_blank_cbf.py
#
#   James Parkhurst, Diamond Light Source, 2012/OCT/16
#
#   Class to read the BLANK.CBF files used in XDS
#
from __future__ import absolute_import, division, print_function

from iotbx.xds import xds_cbf

class reader(xds_cbf.reader):
  """A class to read the BLANK.CBF files used in XDS"""
  def __init__(self):
    xds_cbf.reader.__init__(self)

  def get_data(self):
    """Get the dark current array from the file"""
    # XXX is this actually necessary here?
    import numpy # import dependency

    # Get the image data
    return xds_cbf.reader.get_data(self)


if __name__ == '__main__':
    import sys
    handle = reader()
    handle.read_file(sys.argv[1])
    image = handle.get_data()
