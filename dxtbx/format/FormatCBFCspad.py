#!/usr/bin/env python
# FormatCBFCspad.py
#
# Contains class methods specific to interacting with CSPAD images
#
# $Id:
#

from __future__ import division

from dxtbx.format.FormatCBFMultiTileHierarchy import FormatCBFMultiTileHierarchy
from dxtbx.format.FormatStill import FormatStill
import pycbf

class FormatCBFCspad(FormatCBFMultiTileHierarchy, FormatStill):
  '''An image reading class CSPAD CBF files'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CSPD CBF format image, i.e. we can
    make sense of it.'''

    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_widefile(image_file, pycbf.MSG_DIGEST)

    cbf_handle.find_category("diffrn_detector")
    if cbf_handle.count_rows() > 1:
      return False # support 1 detector per file for now

    cbf_handle.find_column("type")

    return cbf_handle.get_value() == "CS PAD"
