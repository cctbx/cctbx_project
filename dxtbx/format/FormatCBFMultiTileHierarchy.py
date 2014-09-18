#!/usr/bin/env python
# FormatCBFMultiTileHierarchy.py
#
# Reads a multi-tile CBF image, discovering it's detector geometery
# automatically, and builds a hierarchy if present
#
# $Id:
#

from __future__ import division

import pycbf

from FormatCBFMultiTile import FormatCBFMultiTile
from dxtbx.model.detector import HierarchicalDetector
from scitbx.matrix import sqr, col
from libtbx.utils import Sorry

class FormatCBFMultiTileHierarchy(FormatCBFMultiTile):
  '''An image reading class multi-tile CBF files'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_widefile(image_file, pycbf.MSG_DIGEST)

    #check if multiple arrays
    if cbf_handle.count_elements() <= 1:
      return False

    # we need the optional column equipment_component to build a hierarchy
    try:
      cbf_handle.find_category("axis")
      cbf_handle.find_column("equipment_component")
    except Exception, e:
      if "CBF_NOTFOUND" in str(e):
        return False
      else:
        raise e
    return True

  def __init__(self, image_file):
    '''Initialise the image structure from the given file.'''

    assert(self.understand(image_file))

    FormatCBFMultiTile.__init__(self, image_file)

    return

  def _start(self):
    '''Parent class will open the image file as a cbf file handle, and keep
    the handle somewhere safe.'''
    FormatCBFMultiTile._start(self)

  def _get_change_of_basis(self, axis_id):
    """ Get the 4x4 homogenous coordinate matrix for a given axis.  Assumes
    the cbf handle has been intialized
    @param axis_id axis name of basis to get """
    cbf = self._cbf_handle
    axis_type = cbf.get_axis_type(axis_id)

    offset = col(cbf.get_axis_offset(axis_id))
    vector = col(cbf.get_axis_vector(axis_id)).normalize()
    setting, increment = cbf.get_axis_setting(axis_id)

    # change of basis matrix in homologous coordinates
    cob = None

    if axis_type == "rotation":
      r3 = vector.axis_and_angle_as_r3_rotation_matrix(setting + increment, deg = True)
      cob = sqr((r3[0], r3[1], r3[2], offset[0],
                 r3[3], r3[4], r3[5], offset[1],
                 r3[6], r3[7], r3[8], offset[2],
                 0,     0,     0,         1))
    elif axis_type == "translation":
      translation = offset + vector * (setting + increment)
      cob = sqr((1,0,0,translation[0],
                 0,1,0,translation[1],
                 0,0,1,translation[2],
                 0,0,0,1))
    else:
      raise Sorry("Unrecognized vector type: %d"%axis_type)

    return cob

  def _get_cummulative_change_of_basis(self, axis_id):
    """ Get the 4x4 homogenous coordinate matrix for a given axis, combining it with the change of
    basis matrices of parent axes with the same equipment component as the given axis. Assumes
    the cbf handle has been intialized
    @param axis_id axis name of basis to get
    @return (parent, change of basis matrix), where parent is None if the parent in the cbf file
    is ".".  Parent is the axis that the top level axis in this chain of dependent axis depends on
    """

    cbf = self._cbf_handle
    cob = self._get_change_of_basis(axis_id)

    parent_id = cbf.get_axis_depends_on(axis_id)

    if parent_id == ".":
      return None, cob

    eq_comp = cbf.get_axis_equipment_component(axis_id)
    parent_eq_comp = cbf.get_axis_equipment_component(parent_id)
    if eq_comp == parent_eq_comp:
      non_matching_parent, parent_cob = self._get_cummulative_change_of_basis(parent_id)
      return non_matching_parent, parent_cob * cob

    return parent_id, cob

  def _add_panel_group(self, group_id, d):
    """ Adds a panel group to the detector d.  If the group's parent hasn't been
    added yet, recursively add parents to the detector until the detector itself
    is reached.
    @param group_id name of a cbf axis
    @param d detector object
    """

    # group_id will only be "." if the panel being worked on has the same equipment_component name as the
    # last axis in the hierarchy, which isn't really sensible
    assert group_id != "."

    name = group_id

    for subobj in d.hierarchy().iter_preorder():
      if subobj.get_name() == name:
        return subobj

    parent, cob = self._get_cummulative_change_of_basis(group_id)

    if parent is None:
      pg = d.hierarchy() # root object for the detector
      try:
        pg.get_D_matrix() # test to see if we've initialized the detector basis yet
      except RuntimeError, e:
        assert "DXTBX_ASSERT(D_)" in str(e)
      else:
        assert False # shouldn't be reached.  Detector should be initialized only once.
    else:
      parent_pg = self._add_panel_group(parent, d)
      pg = parent_pg.add_group()

    # set up the dxtbx d matrix.  Note use of homogenous coordinates.
    origin = col((cob * col((0,0,0,1)))[0:3])
    fast   = col((cob * col((1,0,0,1)))[0:3]) - origin
    slow   = col((cob * col((0,1,0,1)))[0:3]) - origin

    pg.set_local_frame(
      fast.elems,
      slow.elems,
      origin.elems)

    pg.set_name(name)
    return pg

  def _detector(self):
    '''Return a working detector instance.'''

    cbf = self._cbf_handle

    d = HierarchicalDetector()

    # find the panel elment names. Either array ids or section ids
    cbf.find_category("array_structure_list")
    try:
      cbf.find_column("array_section_id")
    except Exception, e:
      if "CBF_NOTFOUND" not in str(e): raise e
      cbf.find_column("array_id")

    panel_names = []
    for i in xrange(cbf.count_rows()):
      cbf.select_row(i)
      if cbf.get_typeofvalue() == 'null':
        continue

      val = cbf.get_value()
      if val not in panel_names:
        panel_names.append(val)

    # the cbf detector objects are not guaranteed to be in the same order
    # as this array of panel names. re-iterate, associating root axes of
    # detector objects with panel names
    detector_axes = []

    for i in xrange(len(panel_names)):
      cbf_detector = cbf.construct_detector(i)
      axis0 = cbf_detector.get_detector_surface_axes(0)
      detector_axes.append(axis0)
    panel_names_detectororder = []
    cbf.find_category("array_structure_list")
    for detector_axis in detector_axes:
      cbf.find_column("axis_set_id")
      cbf.find_row(detector_axis)
      try:
        cbf.find_column("array_section_id")
      except Exception, e:
        if "CBF_NOTFOUND" not in str(e): raise e
        cbf.find_column("array_id")
      panel_names_detectororder.append(cbf.get_value())


    for panel_name in panel_names:
      cbf_detector = cbf.construct_detector(panel_names_detectororder.index(panel_name))

      # code adapted below from dxtbx.model.detector.detector_factory.imgCIF_H
      pixel = (cbf_detector.get_inferred_pixel_size(1),
               cbf_detector.get_inferred_pixel_size(2))

      axis0 = cbf_detector.get_detector_surface_axes(0)
      axis1 = cbf_detector.get_detector_surface_axes(1)
      assert cbf.get_axis_depends_on(axis0) == axis1

      try:
        size = tuple(cbf.get_image_size_fs(i))
      except Exception, e:
        if "CBF_NOTFOUND" in str(e):
          # no array data in the file, it's probably just a cbf header.  Get the image size elsewhere
          size = [0,0]
          cbf.find_category("array_structure_list")
          for axis in [axis0, axis1]:
            cbf.find_column("axis_set_id")
            cbf.find_row(axis)
            cbf.find_column("precedence")
            idx = int(cbf.get_value()) - 1
            cbf.find_column("dimension")
            size[idx] = int(cbf.get_value())
          assert size[0] != 0 and size[1] != 0
        else:
          raise e

      parent, cob = self._get_cummulative_change_of_basis(axis0)

      pg = self._add_panel_group(parent, d)

      p = pg.add_panel(d.add_panel())

      fast = cbf.get_axis_vector(axis0)
      slow = cbf.get_axis_vector(axis1)
      origin = (cob * col((0,0,0,1)))[0:3]

      p.set_local_frame(fast, slow, origin)

      try:
        cbf.find_category('array_intensities')
        cbf.find_column('undefined_value')
        underload = cbf.get_doublevalue()
        overload = cbf.get_overload(0)
        trusted_range = (underload, overload)
      except: # intentional
        trusted_range = (0.0, 0.0)

      p.set_pixel_size(tuple(map(float, pixel)))
      p.set_image_size(size)
      p.set_trusted_range(tuple(map(float, trusted_range)))
      p.set_name(panel_name)
      #p.set_px_mm_strategy(px_mm) FIXME

      cbf_detector.__swig_destroy__(cbf_detector)
      del(cbf_detector)

    return d

  def _beam(self):
    '''Return a working beam instance.'''

    return self._beam_factory.imgCIF_H(self._cbf_handle)

  def get_raw_data(self, index=None):
    if self._raw_data is None:
      import numpy
      from scitbx.array_family import flex

      self._raw_data = []

      cbf = self._cbf_handle
      cbf.find_category('array_structure')
      cbf.find_column('encoding_type')
      cbf.select_row(0)
      types = []
      for i in xrange(cbf.count_rows()):
        types.append(cbf.get_value())
        cbf.next_row()
      assert len(types) == cbf.count_rows()

      # read the data
      data = {}
      cbf.find_category("array_data")
      for i in xrange(cbf.count_rows()):
        cbf.find_column("array_id")
        name = cbf.get_value()

        cbf.find_column("data")
        assert cbf.get_typeofvalue().find('bnry') > -1

        if types[i] == 'signed 32-bit integer':
          array_string = cbf.get_integerarray_as_string()
          array = flex.int(numpy.fromstring(array_string, numpy.int32))
          parameters = cbf.get_integerarrayparameters_wdims_fs()
          array_size = (parameters[11], parameters[10], parameters[9])
        elif types[i] == 'signed 64-bit real IEEE':
          array_string = cbf.get_realarray_as_string()
          array = flex.double(numpy.fromstring(array_string, numpy.float))
          parameters = cbf.get_realarrayparameters_wdims_fs()
          array_size = (parameters[7], parameters[6], parameters[5])
        else:
          return None # type not supported

        array.reshape(flex.grid(*array_size))
        data[name] = array
        cbf.next_row()

      # extract the data for each panel
      has_sections = True
      try:
        cbf.find_category("array_structure_list_section")
      except Exception, e:
        if "CBF_NOTFOUND" not in str(e): raise e
        has_sections = False

      if has_sections:
        section_shapes = {}
        for i in xrange(cbf.count_rows()):
          cbf.find_column("id")
          section_name = cbf.get_value()
          if not section_name in section_shapes:
            section_shapes[section_name] = {}
          cbf.find_column("array_id")
          if not "array_id" in section_shapes[section_name]:
            section_shapes[section_name]["array_id"] = cbf.get_value()
          else:
            assert section_shapes[section_name]["array_id"] == cbf.get_value()
          cbf.find_column("index"); axis_index = int(cbf.get_value())-1
          cbf.find_column("start"); axis_start = int(cbf.get_value())-1
          cbf.find_column("end");   axis_end   = int(cbf.get_value())

          section_shapes[section_name][axis_index] = slice(axis_start, axis_end)
          cbf.next_row()

        for section_name in sorted(section_shapes):
          section_shape  = section_shapes[section_name]
          section = data[section_shape["array_id"]][ \
            section_shape[2], section_shape[1], section_shape[0]]
          section.reshape(flex.grid(section.focus()[-2], section.focus()[-1]))
          self._raw_data.append(section)
      else:
        for key in sorted(data):
          data[key].reshape(flex.grid(data[key].focus()[-2],data[key].focus()[-1]))
          self._raw_data.append(data[key])

      d = self.get_detector()
      assert len(d) == len(self._raw_data)

    if index is not None:
      return self._raw_data[index]
    return self._raw_data[0]

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMultiTileHierarchy.understand(arg)
