
from __future__ import division

# TODO:
#  - prompt user for missing symmetry
#  - cached scenes



r"""
Load mtz file for viewing reflections in a webbrowser using NGL.
Webbrowser input is a javascript file that together with ngl.js displays reflections
positioned in reciprocal space using a webbrowser. Usage:


from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = r"C:\Users\oeffner\Buser\NGL_HKLviewer\myjstr.js")

myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Work\ANI_TNCS\4PA9\4pa9.tncs.mtz")
myHKLview.SetColumn(0)
myHKLview.SetSphereScale(3)
myHKLview.SetColourColumn(3)
myHKLview.SetRadiusColumn(7)

myHKLview.SetColumn(3)
myHKLview.SetSphereScale(3)
myHKLview.SetColumnBinThresholds("asdf", [0, 0.1, 1, 10])
myHKLview.ExpandToP1(True)
myHKLview.SetOpacity(2, 0.1)


myHKLview.ShowSlice(True, "l", 17)
# get the composed javascript for NGL
myHKLview.GetNGLstring()

myHKLview.ShowSlice(False)
myHKLview.ExpandAnomalous(True)

myHKLview.ShowMissing(True)
myHKLview.SetSpaceGroupChoice(3)

myHKLview.load_reflections_file(r"C:\Users\oeffner\Buser\Phenix\dev-2814-working\modules\phenix_examples\lysozyme-MRSAD\lyso2001_scala1.mtz")
myHKLview.SetCameraType("persp")

myHKLview.load_reflections_file(r"C:\Users\oeffner\Buser\Experiments\CRLF3\DLS20151206CRLF3\5840-F11-X1-Hg-SAD-ONsoak\5840-F11-X1_pk_5_5_1_\xia2\dials-run\DataFiles\mx11235v49_x5840F11X1pk551_free.mtz")
myHKLview.SetColumn(4)
myHKLview.ShowSlice(True, "h", 20)








from cctbx.xray import observation_types
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal

xs = crystal.symmetry(unit_cell=(50,50,40, 90,90,120), space_group_symbol="P3 1")
mi = flex.miller_index([ (1,-2,3), (0,0,-4), (1, 2, 3), (0, 1, 2),
                        (1, 0, 2), (-1, 1, -2), (2, -2, -2),
                        (-2, 1, 0) , (1, 0, -2), (0, 0, 2) ]
)

ma = miller.array( miller.set(xs, mi) )

ma1 = miller.array( miller.set(xs, mi), data=flex.double( [11.205, 6.353, 26.167, 14.94, 2.42, 24.921, 16.185, 11.798, 21.183, 4.98] ),
                   sigmas=flex.double( [13.695, 6.353, 24.921, 6.225, 11.193, 26.167, 8.715, 4.538, 27.413, 21.165] )
                   ).set_observation_type( observation_types.intensity() )

mi2 = flex.miller_index([ (1,-2,3), (0,0,-4), (1, 2, 3), (0, 1, 2),
                        (1, 0, 2), (-1, 1, -2), (2, -2, -2),
                        (-2, 1, 0) , (0, 0, 2) ]
                        )

ma2 = miller.array(miller.set(xs, mi2), flex.complex_double( [ -18.691 +7.47j,
                                             14.94 + 1.21j, -24.921 - 4.98j,
                                             -20.572 + 19.937j, 21.183 - 8.715j,
                                             9.378 + 26.167j, -11.205 + 14.824j,
                                             22.429 - 4.98j, -8.471 + 27.413j
                                            ] )
                   )
mi3 = flex.miller_index([ (1,-2,3), (0,0,-4), (0, 1, 2),
                        (1, 0, 2), (-1, 1, -2),
                        (-2, 1, 0) , (1, 0, -2), (0, 0, 2) ]
                        )

ma3 = miller.array(miller.set(xs, mi3), data=flex.double( [22.429, 28.635, 3.328, 3.738, 24.9, 14.521, 3.738, 19.92] ) )

mi4 = flex.miller_index([ (1,-2,3), (0,0,-4), (1, 2, 3), (0, 1, 2),
                        (1, 0, 2), (-1, 1, -2),   (0, 0, 2) ]
                        )

ma4 = miller.array(miller.set(xs, mi4), data=flex.double( [19.937, 12.45, 11.496, 23.675, 4.98, 1.21, 28.659] ) )



"""



from crys3d.hklview import jsview_3d as view_3d
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


class HKLViewFrame () :
  def __init__ (self, *args, **kwds) :
    self.miller_array = None
    self.spacegroup_choices = None
    self.settings = display.settings()
    kwds['settings'] =self.settings
    self.viewer = view_3d.hklview_3d( **kwds )
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
        "than a full pseudo-precession view.). yes/no?\n")
      if (merge.lower()[0] == "y") :
        merge = True
        #array = array.merge_equivalents().array().set_info(info)
        details.append("merged")
        #self.update_settings_for_merged(True)
      else :
        merge = False
        details.append("unmerged data")
        #self.update_settings_for_unmerged()
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
    self.miller_array = array
    self.update_space_group_choices()
    self.viewer.set_miller_array(array, merge=array_info.merge,
       details=array_info.details_str, valid_arrays=self.valid_arrays)


  def update_settings (self, *args, **kwds) :
    if (self.miller_array is None) :
      print "No miller array has been selected"
      return False
    msg = self.viewer.update_settings(*args, **kwds)
    print msg


  def update_space_group_choices (self) :
    from cctbx.sgtbx.subgroups import subgroups
    from cctbx import sgtbx
    sg_info = self.miller_array.space_group_info()
    subgrs = subgroups(sg_info).groups_parent_setting()
    self.spacegroup_choices = []
    for i,subgroup in enumerate(subgrs) :
      subgroup_info = sgtbx.space_group_info(group=subgroup)
      self.spacegroup_choices.append(subgroup_info)
      print i, subgroup_info.symbol_and_number()
    if (sg_info in self.spacegroup_choices) :
      self.current_spacegroup = self.spacegroup_choices.index(sg_info)
    else :
      self.spacegroup_choices.insert(0, sg_info)
      self.current_spacegroup = sg_info


  def delete_miller_index (self, hkl) :
    if (self.miller_array is None) :
      raise Sorry("No data loaded!")
    info = self.miller_array.info()
    self.miller_array = self.miller_array.delete_index(hkl).set_info(info)
    self.viewer.set_miller_array(self.miller_array)
    self.viewer.DrawNGLJavaScript()


  def LoadReflectionsFile (self, file_name, set_array=True,
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
      for i,e in enumerate(array_info):
        print i, e
      if (len(valid_arrays) == 0) :
        msg = "No arrays of the supported types in this file."
        raise Sorry(msg)
      elif (len(valid_arrays) == 1) :
        if (set_array) :
          self.set_miller_array(valid_arrays[0])
        return valid_arrays[0]


  def GetNGLstring(self):
    return self.viewer.NGLscriptstr


  def SetSphereScale(self, nscale):
    self.settings.scale = nscale
    self.update_settings()


  def SetCameraType(self, camtype):
    if camtype.lower() in "orthographic":
      self.viewer.cameratype = "orthographic"
    if camtype.lower() in "perspective":
      self.viewer.cameratype = "perspective"
    self.viewer.DrawNGLJavaScript()


  def ExpandToP1(self, val):
    self.settings.expand_to_p1 = val
    self.update_settings()


  def ExpandAnomalous(self, val):
    self.settings.expand_anomalous = val
    self.update_settings()


  def ShowOnlyMissing(self, val):
    self.settings.show_only_missing = val
    self.update_settings()


  def ShowMissing(self, val):
    self.settings.show_missing = val
    self.update_settings()


  def ShowDataOverSigma(self, val):
    self.settings.show_data_over_sigma = val
    self.update_settings()


  def ShowSystematicAbsences(self, val):
    self.settings.show_systematic_absences = val
    self.update_settings()


  def ShowSlice(self, val, axis="h", index=0):
    self.settings.slice_mode = val
    self.settings.slice_axis = axis.lower()
    self.settings.slice_index = index
    self.update_settings()


  def SetColumnBinThresholds(self, iarray, binvals):
    self.viewer.iarray = iarray
    self.viewer.binvals = binvals
    self.update_settings()


  def SetOpacity(self, bin, alpha):
    self.viewer.SetOpacity(bin, alpha)


  def SetColumn (self, column_sel) :
    # print("kau printing column_sel ", column_sel)
    self.set_miller_array(self.valid_arrays[column_sel])
    if (self.miller_array is None) :
      raise Sorry("No data loaded!")
    print "Miller array %s runs from hkls: %s to %s" \
     %(self.miller_array.info().label_string(), self.miller_array.index_span().min(),
        self.miller_array.index_span().max() )
    self.viewer.binvals = []
    self.viewer.iarray = column_sel
    self.viewer.icolourcol = column_sel
    self.viewer.iradiicol = column_sel
    self.viewer.DrawNGLJavaScript()


  def SetColourColumn(self, colourcol):
    self.viewer.icolourcol = colourcol
    self.update_settings()


  def SetRadiusColumn(self, radiuscol):
    self.viewer.iradiicol = radiuscol
    self.update_settings()


  def SetSpaceGroupChoice(self, n) :
    if (self.miller_array is None) :
      raise Sorry("No data loaded!")
    self.current_spacegroup = self.spacegroup_choices[n]
    from cctbx import crystal
    symm = crystal.symmetry(
      space_group_info= self.current_spacegroup,
      unit_cell=self.miller_array.unit_cell())
    array = self.miller_array.expand_to_p1().customized_copy(
      crystal_symmetry=symm)
    array = array.merge_equivalents().array().set_info(self.miller_array.info())
    othervalidarrays = []
    for validarray in self.valid_arrays:
      #print "Space group casting ", validarray.info().label_string()
      arr = validarray.expand_to_p1().customized_copy(crystal_symmetry=symm)
      arr = arr.merge_equivalents().array().set_info(validarray.info())
      othervalidarrays.append( arr  )
    print "MERGING 2"
    self.viewer.set_miller_array(array, valid_arrays=othervalidarrays)
    self.viewer.DrawNGLJavaScript()

