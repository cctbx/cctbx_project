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
from scitbx.matrix import col, sqr
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

  def sync_detector_to_cbf(self):
    detector = self.get_detector()
    root = detector.hierarchy()
    assert len(root) == 4

    cbf = self._cbf_handle

    for quad in root:
      cbf.find_category("axis")
      cbf.find_column("id")
      while True:
        if cbf.get_value() in quad.get_name():
          axis_name = cbf.get_value()
          break
        cbf.next_row()
      cbf.find_column("equipment_component")
      assert cbf.get_value() in quad.get_name()

      d_mat = quad.get_local_d_matrix()
      fast = col((d_mat[0],d_mat[3],d_mat[6])).normalize()
      slow = col((d_mat[1],d_mat[4],d_mat[7])).normalize()
      orig = col((d_mat[2],d_mat[5],d_mat[8]))

      cbf.find_column("offset[1]"); cbf.set_value(str(orig[0]))
      cbf.find_column("offset[2]"); cbf.set_value(str(orig[1]))
      cbf.find_column("offset[3]"); cbf.set_value(str(orig[2]))

      v3 = fast.cross(slow).normalize()

      r3 = sqr((fast[0],slow[0],v3[0],
                fast[1],slow[1],v3[1],
                fast[2],slow[2],v3[2]))

      angle, axis = r3.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)
      if angle == 0:
        axis = (0,0,1)

      cbf.find_column("vector[1]"); cbf.set_value(str(axis[0]))
      cbf.find_column("vector[2]"); cbf.set_value(str(axis[1]))
      cbf.find_column("vector[3]"); cbf.set_value(str(axis[2]))

      cbf.set_axis_setting(axis_name, angle, 0)
