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

  def sync_detector_to_cbf(self, detector = None):
    ''' If the dectector object has been changed, due to refinment or manual repositioning
    in a gui, call this function to synchronize these changes to the underlying CBF handle'''

    def recursive_sync(cbf, group, root = False):
      ''' Walks the hierarchy and synchronizes the dxtbx matrices to the CBF axes.'''
      try:
        iterator = iter(group) # root object is panel which is not iteratable
      except TypeError, e:
        is_panel = True
      else:
        is_panel = False

      d_mat = group.get_local_d_matrix()
      fast = col((d_mat[0],d_mat[3],d_mat[6])).normalize()
      slow = col((d_mat[1],d_mat[4],d_mat[7])).normalize()
      orig = col((d_mat[2],d_mat[5],d_mat[8]))

      if is_panel:
        cbf.find_category("diffrn_data_frame")
        cbf.find_column("array_id")
        cbf.find_row(group.get_name())
        cbf.find_column("binary_id")
        array_num = int(cbf.get_value()) - 1
        cbf_detector = cbf.construct_detector(array_num)
        axis0 = cbf_detector.get_detector_surface_axes(0)
        axis1 = cbf_detector.get_detector_surface_axes(1)
        assert cbf.get_axis_depends_on(axis0) == axis1
        name = cbf.get_axis_depends_on(axis1)

        # account for the sign change when moving from dxtbx to ImageCIF
        slow = col((slow[0],-slow[1],slow[2]))

        # account for the offset from the center of the sensor to the center
        # of the ASIC.
        orig -= col(cbf.get_axis_offset(axis1))
      else:
        name = group.get_name()

      v3 = fast.cross(slow).normalize()

      if root:
        # the distance is encoded in a different axis.
        distance = orig[2]
        orig = col((orig[0],orig[1],0))

      cbf.find_category("axis")
      cbf.find_column("id")
      while True:
        if cbf.get_value() == name:
          axis_name = cbf.get_value()
          break
        cbf.next_row()

      cbf.find_column("offset[1]"); cbf.set_value(str(orig[0]))
      cbf.find_column("offset[2]"); cbf.set_value(str(orig[1]))
      cbf.find_column("offset[3]"); cbf.set_value(str(orig[2]))

      r3 = sqr((fast[0],slow[0],v3[0],
                fast[1],slow[1],v3[1],
                fast[2],slow[2],v3[2]))

      angle, axis = r3.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)
      if angle == 0:
        axis = col((0,0,1))

      assert axis.length() > 0

      cbf.find_column("vector[1]"); cbf.set_value(str(axis[0]))
      cbf.find_column("vector[2]"); cbf.set_value(str(axis[1]))
      cbf.find_column("vector[3]"); cbf.set_value(str(axis[2]))

      cbf.set_axis_setting(axis_name, angle, 0)

      if root:
        # synchronize the new distance
        while not "_Z" in axis_name:
          axis_name = cbf.get_axis_depends_on(axis_name)
        assert cbf.get_axis_type(axis_name) == 'translation'
        cbf.set_axis_setting(axis_name, distance, 0)

      if not is_panel:
        for subgroup in iterator:
          recursive_sync(cbf, subgroup)


    if detector is None:
      detector = self.get_detector()

    root = detector.hierarchy()
    assert len(root) == 4

    recursive_sync(self._cbf_handle, root, True)

class FormatCBFCspadInMemory(FormatCBFCspad):
  """ Overrides the Format object's init method to accept a cbf handle instead
      of a file name. Used with XFELs when it is desirable to never write
      a file to disk, but to process it only in memory.
  """

  @staticmethod
  def understand(image_file):
    """ This class shouldn't be found by the dxtbx Registry. Instead it should
        be instantiated directly as needed
    """
    return False

  def __init__(self, cbf_handle):
    """ @param cbf_handle In memory cbf_handle, alredy initialized """
    from dxtbx.model.detector import detector_factory
    from dxtbx.model.beam import beam_factory

    self._goniometer_instance = None
    self._scan_instance = None
    self._detector_instance = None
    self._detector_factory = detector_factory
    self._beam_factory = beam_factory

    self._cbf_handle = cbf_handle
    self._raw_data = None
    try:
      detector_instance = self._detector()
      assert(isinstance(detector_instance, Detector))
      self._detector_instance = detector_instance

      beam_instance = self._beam()
      assert(isinstance(beam_instance, Beam))
      self._beam_instance = beam_instance

    except Exception, e:
    #except exceptions.Exception, e:
      # FIXME ideally should not squash the errors here...
      pass
    finally:
      self._end()

    return
