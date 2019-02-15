from __future__ import division

# TODO:
#  - prompt user for missing symmetry
#  - cached scenes

from crys3d.hklview import jsview_3d
from cctbx.miller import display

from libtbx import object_oriented_patterns as oop
from libtbx.str_utils import format_value
from libtbx.utils import Sorry, Abort, to_str
import libtbx.load_env
from libtbx import group_args
from math import sqrt
import copy
import os, sys


argn = 1
argc = len(sys.argv)

# prompt user for value if it's not on the commandline
def Inputarg(varname):
  global argn
  global argc
  if argc > 1 and argn < argc:
    myvar = sys.argv[argn]
    argn = argn + 1
    print varname + " " + myvar
  else:
    myvar = raw_input(varname)
  return myvar



class settings_window () :

  def set_index_span (self, index_span) :
    self._index_span = index_span

  def update_reflection_info (self, hkl, d_min, value) :
    print hkl, value
    if (hkl is None) :
      self.hkl_info.SetValue("")
      self.d_min_info.SetValue("")
      self.value_info.SetValue("")
    else :
      self.hkl_info.SetValue("%d, %d, %d" % hkl)
      d_min_str = format_value("%.3g", d_min)
      self.d_min_info.SetValue(d_min_str)
      value_str = format_value("%.3g", value, replace_none_with="---")
      self.value_info.SetValue(value_str)

  def clear_reflection_info (self) :
    self.update_reflection_info(None, None, None)

  def update_space_group_choices (self, miller_array) :
    from cctbx.sgtbx.subgroups import subgroups
    from cctbx import sgtbx
    sg_info  = miller_array.space_group_info()
    subgrs = subgroups(sg_info).groups_parent_setting()
    choices = []
    for subgroup in subgrs :
      subgroup_info = sgtbx.space_group_info(group=subgroup)
      choices.append(str(subgroup_info))
    if (str(sg_info) in choices) :
      current = choices.index(str(sg_info))
    else :
      choices.insert(0, str(sg_info))
      current = 0
    self.sg_ctrl.SetItems(choices)
    self.sg_ctrl.SetSelection(current)
    self._last_sg_sel = str(sg_info)

# XXXXXXXXXXXXXXXXXXXXXX
# Kau added starts here

  # def update_column_choices (self, array_info,valid_arrays,sel) :

  def update_column_choices (self,array_info,arrays,sel) :
    choices=[]
        # print("Kau printing array_info from within function", array_info)
    for labels in array_info:
        choices.append(str(labels))
    self.column_ctrl.SetItems(choices)
    # for f in valid_arrays:
    #     print("Kau last array selected is ", valid_arrays)
    current=sel
    #self.column_ctrl.SetItems(choices)
    #self.column_ctrl.SetSelection(current)
    #self._last_column_sel = str(sel)

# Kau added ends here
# XXXXXXXXXXXXXXXXXXXXXX

class HKLViewFrame () :
  def __init__ (self, *args, **kwds) :
    self.miller_array = None
    self.settings = display.settings()
    self.viewer = jsview_3d.hklview_3d(self.settings)
    self.viewer.set_miller_array(self.viewer.miller_array)


  def update_clicked (self, index) :#hkl, d_min=None, value=None) :
    if (index is None) :
      self.settings_panel.clear_reflection_info()
    else :
      hkl, d_min, value = self.viewer.scene.get_reflection_info(index)
      self.settings_panel.update_reflection_info(hkl, d_min, value)

  def process_miller_array (self, array) :
    if (array is None) : return
    if (array.is_hendrickson_lattman_array()) :
      raise Sorry("Hendrickson-Lattman coefficients are not supported.")
    info = array.info()
    if isinstance(info, str) :
      labels = "TEST DATA"
    else :
      labels = info.label_string()
    if (array.unit_cell() is None) or (array.space_group() is None) :
      raise Sorry("No space group info is present in data")

    details = []
    merge = None
    if (not array.is_unique_set_under_symmetry()) :
      merge = Inputarg("The data in the selected array are not symmetry-"+
        "unique, which usually means they are unmerged (but could also be due "+
        "to different indexing conventions).  Do you want to merge equivalent "+
        "observations (preserving anomalous data if present), or view the "+
        "array unmodified?  (Note that if you do not merge the array, the "+
        "options to expand to P1 or generate Friedel pairs will be be disabled"+
        ", and the 2D view will only show indices present in the file, rather "+
        "than a full pseudo-precession view.). Y/N?")
      if (merge.lower() == "y") :
        merge = True
        #array = array.merge_equivalents().array().set_info(info)
        details.append("merged")
        #self.update_settings_for_merged(True)
      else :
        merge = False
        details.append("unmerged data")
        self.update_settings_for_unmerged()
        self.settings.expand_to_p1 = False
        self.settings.expand_anomalous = False
    #else :
      #self.update_settings_for_merged()
    if array.is_complex_array() :
      array = array.amplitudes().set_info(info)
      details.append("as amplitudes")
    from iotbx.reflection_file_utils import looks_like_r_free_flags_info
    if (array.is_integer_array()) and (looks_like_r_free_flags_info(info)) :
      from iotbx.reflection_file_utils import get_r_free_flags_scores
      score_array = get_r_free_flags_scores([array], None)
      test_flag_value = score_array.test_flag_values[0]
      array = array.customized_copy(data=(array.data() == test_flag_value))
      array.set_info(info)
    sg = "%s" % array.space_group_info()
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % array.unit_cell().parameters()
    details_str = ""
    if (len(details) > 0) :
      details_str = "(%s)" % ", ".join(details)
    array_info = group_args(
      labels=labels,
      details_str=details_str,
      merge=merge,
      sg=sg,
      uc=uc)
    return array, array_info

  def set_miller_array (self, array) :
    if (array is None) : return
    array, array_info = self.process_miller_array(array)
    print "Data: %s %s (Space group: %s  Unit Cell: %s" \
      % (array_info.labels, array_info.details_str, array_info.sg,
          array_info.uc)
    self.miller_array = array
    # print("kau this is array ",type(array))
    self.viewer.set_miller_array(array, zoom=True, merge=array_info.merge)


  def update_settings (self, *args, **kwds) :
    if (self.miller_array is None) :
      return False
    self.viewer.update_settings(*args, **kwds)

  def set_space_group (self, space_group_info) :
    # print("kau printing space_group_info ", space_group_info)
    if (self.miller_array is None) :
      raise Sorry("No data loaded!")
    from cctbx import crystal
    symm = crystal.symmetry(
      space_group_info=space_group_info,
      unit_cell=self.miller_array.unit_cell())
    array = self.miller_array.expand_to_p1().customized_copy(
      crystal_symmetry=symm)
    print "MERGING 2"
    array = array.merge_equivalents().array().set_info(self.miller_array.info())
    self.viewer.set_miller_array(array, zoom=False)
    self.viewer.DrawNGLJavaScript()

#kau added starts here
  def set_column (self, column_sel) :
    # print("kau printing column_sel ", column_sel)
    self.set_miller_array(self.valid_arrays[column_sel])
    if (self.miller_array is None) :
      raise Sorry("No data loaded!")
    # from cctbx import crystal
    # symm = crystal.symmetry(
    #   space_group_info=space_group_info,
    #   unit_cell=self.miller_array.unit_cell())
    # array = self.miller_array.expand_to_p1().customized_copy(
    #   crystal_symmetry=symm)
    # print "MERGING 2"
    # array = array.merge_equivalents().array().set_info(self.miller_array.info())
    # array=self.miller_array.column_root_label(column_sel)
    # print("Kau this is array from set_column ", array)
    # settings_panel.array_storage()
    # self.viewer.set_miller_array(array[column_sel], zoom=False)
    self.viewer.DrawNGLJavaScript()
##kau added ends here

  def delete_miller_index (self, hkl) :
    if (self.miller_array is None) :
      raise Sorry("No data loaded!")
    info = self.miller_array.info()
    self.miller_array = self.miller_array.delete_index(hkl).set_info(info)
    self.viewer.set_miller_array(self.miller_array, zoom=True)
    self.viewer.DrawNGLJavaScript()

  def load_reflections_file (self, file_name, set_array=True,
      data_only=False) :
    file_name = to_str(file_name)
    if (file_name != "") :
      from iotbx.reflection_file_reader import any_reflection_file
      from iotbx.gui_tools.reflections import get_array_description
      try :
        hkl_file = any_reflection_file(file_name)
      except Exception, e :
        raise Sorry(to_str(e))
      arrays = hkl_file.as_miller_arrays(merge_equivalents=False,
        )#observation_type_callback=misc_dialogs.get_shelx_file_data_type)
      #arrays = f.file_server.miller_arrays
      valid_arrays = []
      array_info = []
      for array in arrays :
        if array.is_hendrickson_lattman_array() :
          continue
        elif (data_only) :
          if (not array.is_real_array()) and (not array.is_complex_array()) :
            continue
        labels = array.info().label_string()
        desc = get_array_description(array)
        array_info.append("%s (%s)" % (labels, desc))
        valid_arrays.append(array)
      self.valid_arrays = valid_arrays
      if (len(valid_arrays) == 0) :
        msg = "No arrays of the supported types in this file."
        raise Sorry(msg)
      elif (len(valid_arrays) == 1) :
        if (set_array) :
          self.set_miller_array(valid_arrays[0])
        return valid_arrays[0]


