# LIBTBX_SET_DISPATCHER_NAME cctbx.multiplicity_viewer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import division
from crys3d.hklview.frames import *
from wxtbx import icons
import wxtbx.app
import iotbx.phil
import libtbx.load_env
from cctbx import crystal
import wx
import sys


master_phil = iotbx.phil.parse("""
  include scope cctbx.miller.display.master_phil
  unit_cell = None
    .type = unit_cell
  space_group = None
    .type = space_group
""", process_includes=True)


def run (args) :
  ma = None
  show_2d = False
  if ("--2d" in args) :
    show_2d = True
    args.remove("--2d")
  pcl = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    reflection_file_def="data",
    pdb_file_def="symmetry_file",
    usage_string="%s f_obs.mtz [options]" %libtbx.env.dispatcher_name)
  settings = pcl.work.extract()
  a = wxtbx.app.CCTBXApp(0)
  app_icon = wx.EmptyIcon()
  app_icon.CopyFromBitmap(icons.hklview_3d.GetBitmap())
  if (wx.VERSION >= (2,9)) :
    tb_icon = wx.TaskBarIcon(wx.TBI_DOCK)
  else :
    tb_icon = wx.TaskBarIcon()
  tb_icon.SetIcon(app_icon, "CCTBX multiplicity viewer")
  a.hklview_settings = settings
  viewer_class = MultiplicityViewFrame
  if (show_2d) :
    viewer_class = MultiplicityViewFrame2D
  f = viewer_class(None, -1, "Reflection multiplicity viewer", size=(1024,768))
  f.Show()
  if (ma is not None) :
    f.set_miller_array(ma)
  elif (settings.data is not None) :
    f.load_reflections_file(settings.data)
  else :
    f.OnLoadFile(None)
  a.SetTopWindow(f)
  a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), f)
  a.MainLoop()


class MultiplicityViewFrame(HKLViewFrame):

  def create_settings_panel(self):
    self.settings_panel = multiplicity_settings_window(self, -1, style=wx.RAISED_BORDER)
    self.settings.scale_radii_multiplicity = True
    self.settings.scale_colors_multiplicity = True

  def update_settings_for_merged(self, enable_multiplicity=False):
    pass

  def set_miller_array (self, array) :
    if (array is None) : return
    array, array_info = self.process_miller_array(array)
    self.statusbar.SetStatusText("Data: %s %s (Space group: %s  Unit Cell: %s)"
      % (array_info.labels, array_info.details_str, array_info.sg,
          array_info.uc))
    self.settings_panel.d_min_ctrl.SetValue(array.d_min())
    self.settings_panel.d_min_ctrl.SetRange(array.d_min(), 20.0)
    self.settings_panel.set_index_span(array.index_span())
    #self.settings_panel.update_space_group_choices(array)
    self.settings.spheres = False
    self.settings_panel.spheres_ctrl.SetValue(False)
    self.miller_array = array
    self.viewer.set_miller_array(array, zoom=True, merge=array_info.merge)
    self.viewer.Refresh()
    if (self.view_2d is not None) :
      self.view_2d.set_miller_array(array)

  def load_reflections_file (self, file_name, set_array=True,
      data_only=True) :
    if (isinstance(file_name, unicode)) :
      file_name = str(file_name)
    if (file_name != "") :
      from iotbx.reflection_file_reader import any_reflection_file
      from iotbx.gui_tools.reflections import get_array_description
      try :
        hkl_file = any_reflection_file(file_name)
      except Exception, e :
        raise Sorry(str(e))
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
      if (len(valid_arrays) == 0) :
        msg = "No arrays of the supported types in this file."
        raise Sorry(msg)
      elif (len(valid_arrays) >= 1) :
        if (set_array) :
          self.set_miller_array(valid_arrays[0])
        return valid_arrays[0]
    raise Abort()

  def process_miller_array (self, array) :
    if (array is None) : return
    if (array.is_hendrickson_lattman_array()) :
      raise Sorry("Hendrickson-Lattman coefficients are not supported.")
    info = array.info()
    if isinstance(info, str) :
      labels = "TEST DATA"
    else :
      labels = info.label_string()
    if self.settings.unit_cell is not None:
      symm = crystal.symmetry(unit_cell=self.settings.unit_cell,
                              space_group=array.space_group())
      array = array.customized_copy(crystal_symmetry=symm).set_info(info)
    if self.settings.space_group is not None:
      symm = crystal.symmetry(unit_cell=array.unit_cell(),
                              space_group=self.settings.space_group.group())
      array = array.customized_copy(crystal_symmetry=symm).set_info(info)
    if (array.unit_cell() is None) or (array.space_group() is None) :
      dlg = wxtbx.symmetry_dialog.SymmetryDialog(self, -1, "Enter symmetry")
      dlg.SetUnitCell(array.unit_cell())
      dlg.SetSpaceGroup(array.space_group_info())
      if (dlg.ShowModal() == wx.ID_OK) :
        symm = dlg.GetSymmetry()
        array = array.customized_copy(crystal_symmetry=symm).set_info(info)
      wx.CallAfter(dlg.Destroy)
    details = []
    merge = True
    details.append("merged")
    self.update_settings_for_merged(True)
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

  def OnShow2D (self, evt) :
    if (self.view_2d is None) :
      self.view_2d = MultiplicityViewFrame2D(self, -1, "2D data viewer")
      self.view_2d.set_miller_array(self.viewer.miller_array)
      self.view_2d.Show()
      if (wxtbx.MAC_OS_X_MAVERICKS) :
        wx.GetApp().SafeYield(self.view_2d, True)
    self.view_2d.Raise()


class multiplicity_settings_window(settings_window):

  def add_controls (self) :
    self._index_span = None
    self._last_sg_sel = None
    # d_min control
    self.d_min_ctrl = floatspin.FloatSpin(parent=self, increment=0.05, digits=2)
    self.d_min_ctrl.Bind(wx.EVT_SET_FOCUS, lambda evt: None)
    if (wx.VERSION >= (2,9)) : # XXX FloatSpin bug in 2.9.2/wxOSX_Cocoa
      self.d_min_ctrl.SetBackgroundColour(self.GetBackgroundColour())
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    label = wx.StaticText(self,-1,"High resolution:")
    box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(self.d_min_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(floatspin.EVT_FLOATSPIN, self.OnChangeResolution, self.d_min_ctrl)
    #
    ctrls = self.create_controls(
      setting="black_background",
      label="Black background")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_axes",
      label="Show h,k,l axes")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    if (not self.is_3d_view) :
      ctrls = self.create_controls(
        setting="uniform_size",
        label="Use same radius for all points")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    if (self.is_3d_view) :
      ctrls = self.create_controls(
        setting="expand_to_p1",
        label="Expand data to P1")
      ctrls2 = self.create_controls(
        setting="expand_anomalous",
        label="show Friedel pairs")
      box = wx.BoxSizer(wx.HORIZONTAL)
      self.panel_sizer.Add(box)
      box.Add(ctrls[0], 0, wx.ALL, 5)
      box.Add(ctrls2[0], 0, wx.ALL, 5)
      ctrls = self.create_controls(
        setting="spheres",
        label="Display reflections as spheres")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
      self.spheres_ctrl = ctrls[0]
    else :
      self.spheres_ctrl = oop.null()
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    txt = wx.StaticText(self.panel, -1, "Color scheme:")
    box.Add(txt, 0, wx.TOP|wx.BOTTOM|wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
    self.color_ctrl = wx.Choice(self.panel, -1,
      choices=["rainbow","heatmap","redblue","grayscale","monochrome"])
    self.color_ctrl.SetStringSelection(self.settings.color_scheme)
    box.Add(self.color_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_CHOICE, self.OnChangeColor, self.color_ctrl)
    ctrls = self.create_controls(
      setting="show_missing",
      label="Show missing reflections")
    ctrls2 = self.create_controls(
      setting="show_only_missing",
      label="only")
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    box.Add(ctrls[0], 0, wx.ALL, 5)
    box.Add(ctrls2[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_systematic_absences",
      label="Show systematic absences")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    if (self.is_3d_view) :
      ctrls = self.create_controls(
        setting="sphere_detail",
        label="Sphere detail level",
        min=4,
        max=20)
      box = wx.BoxSizer(wx.HORIZONTAL)
      box.Add(ctrls[0], 0, wx.TOP|wx.BOTTOM|wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
      box.Add(ctrls[1], 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self.panel_sizer.Add(box)
    else:
      box2 = wx.BoxSizer(wx.HORIZONTAL)
      box2.Add(wx.StaticText(self.panel, -1, "View slice:"), 0,
        wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      ctrls = self.create_controls
      self.hkl_choice = wx.Choice(self.panel, -1, choices=["h","k","l"])
      self.hkl_choice.SetStringSelection(self.settings.slice_axis)
      box2.Add(self.hkl_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      box2.Add(wx.StaticText(self.panel, -1, "="), 0,
        wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self.slice_index = wx.SpinCtrl(self.panel, -1)
      self.slice_index.SetValue(self.settings.slice_index)
      box2.Add(self.slice_index, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self.panel_sizer.Add(box2)
      self.Bind(wx.EVT_CHOICE, self.OnSetSlice, self.hkl_choice)
      self.Bind(wx.EVT_SPINCTRL, self.OnSetSlice, self.slice_index)


class multiplicity_settings_window_2d (multiplicity_settings_window) :
  is_3d_view = False


class MultiplicityViewFrame2D(MultiplicityViewFrame, HKLViewFrame2D):

  def create_settings_panel (self) :
    self.settings.scale_radii_multiplicity = True
    self.settings.scale_colors_multiplicity = True
    self.settings.expand_to_p1 = True
    self.settings.expand_anomalous = True
    self.settings.slice_mode = True
    #self.settings.black_background = False
    self.settings_panel = multiplicity_settings_window_2d(self, -1, style=wx.RAISED_BORDER)

  def OnShow3D (self, evt) :
    if (self.view_3d is None) :
      self.view_3d = MultiplicityViewFrame(self, -1, "3D data viewer")
      self.view_3d.Show()
      if (self.viewer.miller_array is not None) :
        self.view_3d.set_miller_array(self.viewer.miller_array)
    self.view_3d.Raise()


if (__name__ == "__main__") :
  run(sys.argv[1:])
