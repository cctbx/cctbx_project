
# TODO:
#  - prompt user for missing symmetry
#  - cached scenes

from crys3d.hklview import view_2d, view_3d
from cctbx.miller.display import settings
from wxtbx import icons
import wxtbx.symmetry_dialog
import wxtbx.utils
import wx.glcanvas
from wx.lib.agw import floatspin
import wx
from libtbx import object_oriented_patterns as oop
from libtbx.str_utils import format_value
from libtbx.utils import Sorry
import libtbx.load_env
from math import sqrt
import copy
import os

class settings_window (wxtbx.utils.SettingsPanel) :
  is_3d_view = True
  def __init__ (self, *args, **kwds) :
    wxtbx.utils.SettingsPanel.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_CHAR, self.OnChar)

  def OnChar (self, event) :
    self.GetParent().viewer.OnChar(event)

  def add_controls (self) :
    self._index_span = None
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
    ctrls = self.create_controls(
      setting="black_background",
      label="Black background")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_axes",
      label="Show h,k,l axes")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_data_over_sigma",
      label="Use I or F over sigma")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    if (not self.is_3d_view) :
      ctrls = self.create_controls(
        setting="uniform_size",
        label="Use same radius for all points")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sqrt_scale_radii",
      label="Scale radii to sqrt(value)")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sqrt_scale_colors",
      label="Scale colors to sqrt(value)")
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
      ctrls = self.create_controls(
        setting="slice_mode",
        label="Show only a slice through reciprocal space")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
      self.slice_ctrl = ctrls[0]
      ctrls = self.create_controls(
        setting="keep_constant_scale",
        label="Keep scale constant across all slices")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
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
    # reflection info box
    box = wx.StaticBox(self.panel, -1, "Reflection info")
    box_szr = wx.StaticBoxSizer(box, wx.VERTICAL)
    self.panel_sizer.Add((1,10))
    self.panel_sizer.Add(box_szr, 0, wx.EXPAND|wx.ALL)
    grid_szr = wx.FlexGridSizer(rows=3, cols=2)
    box_szr.Add(grid_szr, 0, wx.EXPAND|wx.ALL)
    grid_szr.Add(wx.StaticText(self.panel, -1, "Clicked:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.hkl_info = wx.TextCtrl(self.panel, -1, size=(80,-1),
      style=wx.TE_READONLY)
    grid_szr.Add(self.hkl_info, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid_szr.Add(wx.StaticText(self.panel, -1, "Resolution:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.d_min_info = wx.TextCtrl(self.panel, -1, size=(80,-1),
      style=wx.TE_READONLY)
    grid_szr.Add(self.d_min_info, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    grid_szr.Add(wx.StaticText(self.panel, -1, "Value:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.value_info = wx.TextCtrl(self.panel, -1, size=(80,-1),
      style=wx.TE_READONLY)
    grid_szr.Add(self.value_info, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

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

  def OnSetSlice (self, event) :
    self.settings.slice_axis = self.hkl_choice.GetStringSelection()
    axis_index = ["h","k","l"].index(self.settings.slice_axis)
    min_value = self._index_span.min()[axis_index]
    max_value = self._index_span.max()[axis_index]
    self.settings.slice_index = self.slice_index.GetValue()
    self.slice_index.SetRange(min_value, max_value)
    if (self.settings.slice_index > max_value) :
      self.settings.slice_index = max_value
      self.slice_index.SetValue(max_value)
    elif (self.settings.slice_index < min_value) :
      self.settings.slice_index = min_value
      self.slice_index.SetValue(min_value)
    if (not self.is_3d_view) or (self.slice_ctrl.GetValue()) :
      try :
        self.parent.update_settings()
      except ValueError, e : # TODO set limits
        raise Sorry(str(e))

  def OnChangeResolution (self, event) :
    self.settings.d_min = self.d_min_ctrl.GetValue()
    self.parent.update_settings()

  def OnChangeColor (self, event) :
    self.settings.color_scheme = self.color_ctrl.GetStringSelection()
    self.parent.update_settings()

class HKLViewFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.parent = self.GetParent()
    self.view_2d = None
    self.view_3d = None # used by 2D subclass
    self.statusbar = self.CreateStatusBar()
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    self.toolbar.SetToolBitmapSize((32,32))
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)
    app = wx.GetApp()
    if (getattr(app, "hklview_settings", None) is not None) :
      # XXX copying the initial settings avoids awkward interactions when
      # multiple viewer windows are opened
      self.settings = copy.deepcopy(app.hklview_settings)
    else :
      self.settings = settings()
    self.create_settings_panel()
    self.sizer.Add(self.settings_panel, 0, wx.EXPAND)
    self.create_viewer_panel()
    self.sizer.Add(self.viewer, 1, wx.EXPAND|wx.ALL)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Load file",
      bitmap=icons.hkl_file.GetBitmap(),
      shortHelp="Load file",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnLoadFile, btn)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Save image",
      bitmap=icons.save_all.GetBitmap(),
      shortHelp="Save image",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnSave, btn)
    menubar = wx.MenuBar(-1)
    self.file_menu = wx.Menu()
    menubar.Append(self.file_menu, "File")
    item = wx.MenuItem(self.file_menu, -1, "Load data...\tCtrl-O")
    self.Bind(wx.EVT_MENU, self.OnLoadFile, item)
    self.file_menu.AppendItem(item)
    if (libtbx.env.has_module("phenix")) :
      phenix_dir = os.path.dirname(libtbx.env.dist_path("phenix"))
      examples_dir = os.path.join(phenix_dir, "phenix_examples")
      if (os.path.isdir(examples_dir)) :
        submenu = wx.Menu()
        self.file_menu.AppendSubMenu(submenu, "Load example data")
        examples_and_data = [
          ("p9-sad", "p9.sca"),
          ("sec17-sad", "sec17.sca"),
          ("rnase-s", "rnase25.mtz"),
          ("porin-twin", "porin.mtz"),
        ]
        for subdir, data in examples_and_data :
          example_file = os.path.join(examples_dir, subdir, data)
          item = wx.MenuItem(submenu, -1, subdir)
          submenu.AppendItem(item)
          self.Bind(wx.EVT_MENU,
            lambda evt, f=example_file: self.load_reflections_file(f), item)
    self.add_view_specific_functions()
    self.SetMenuBar(menubar)
    self.toolbar.Realize()
    self.SetSizer(self.sizer)
    self.sizer.SetSizeHints(self)
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self.Bind(wx.EVT_ACTIVATE, self.OnActive)
    self.viewer.SetFocus()

  def OnActive (self, event) :
    if (self.IsShown()) :
      self.viewer.Refresh()

  def create_viewer_panel (self) :
    self.viewer = view_3d.hklview_3d(self, size=(800,600),
      style=wx.glcanvas.WX_GL_DOUBLEBUFFER)

  def create_settings_panel (self) :
    self.settings_panel = settings_window(self, -1, style=wx.RAISED_BORDER)

  def add_view_specific_functions (self) :
    item = wx.MenuItem(self.file_menu, -1, "Show 2D view")
    self.file_menu.AppendItem(item)
    self.Bind(wx.EVT_MENU, self.OnShow2D, item)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Show 2D view",
      bitmap=icons.hklview_2d.GetBitmap(),
      shortHelp="Show 2D view",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnShow2D, btn)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Clear labels",
      bitmap=icons.clear_left.GetBitmap(),
      shortHelp="Clear labels",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnClearLabels, btn)

  def update_clicked (self, hkl, d_min=None, value=None) :
    self.settings_panel.update_reflection_info(hkl, d_min, value)

  def update_settings_for_unmerged (self) :
    self.settings.expand_to_p1 = False
    self.settings.expand_anomalous = False
    self.settings_panel.get_control("expand_to_p1").SetValue(False)
    self.settings_panel.get_control("expand_to_p1").Enable(False)
    self.settings_panel.get_control("expand_anomalous").SetValue(False)
    self.settings_panel.get_control("expand_anomalous").Enable(False)

  def update_settings_for_merged (self) :
    self.settings_panel.get_control("expand_to_p1").Enable(True)
    self.settings_panel.get_control("expand_anomalous").Enable(True)

  def set_miller_array (self, array) :
    if (array.is_hendrickson_lattman_array()) :
      raise Sorry("Hendrickson-Lattman coefficients are not supported.")
    info = array.info()
    if isinstance(info, str) :
      labels = "TEST DATA"
    else :
      labels = info.label_string()
    if (array.unit_cell() is None) or (array.space_group() is None) :
      dlg = wxtbx.symmetry_dialog.SymmetryDialog(self, -1, "Enter symmetry")
      dlg.SetUnitCell(array.unit_cell())
      dlg.SetSpaceGroup(array.space_group_info())
      if (dlg.ShowModal() == wx.ID_OK) :
        symm = dlg.GetSymmetry()
        array = array.customized_copy(crystal_symmetry=symm).set_info(info)
      wx.CallAfter(dlg.Destroy)
    details = []
    if (not array.is_unique_set_under_symmetry()) :
      details.append("unmerged data")
      self.update_settings_for_unmerged()
    else :
      self.update_settings_for_merged()
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
    self.statusbar.SetStatusText("Data: %s %s (Space group: %s  Unit Cell: %s)"
      % (labels, details_str, sg, uc))
    self.settings_panel.d_min_ctrl.SetValue(array.d_min())
    self.settings_panel.d_min_ctrl.SetRange(array.d_min(), 20.0)
    self.settings_panel.set_index_span(array.index_span())
    if (type(self).__name__ == "HKLViewFrame") :
      if (array.indices().size() > 100000) :
        if (self.settings.spheres) :
          cnf = wx.MessageBox(message="Warning: this is a lot of reflections; "+
            "unless you have a very powerful graphics card, displaying "+
            "spheres may be slow and/or unstable, especially if data are "+
            "expanded to P1.  Do you want to switch to a faster rendering "+
            "style?", style=wx.YES|wx.NO)
          if (cnf == wx.YES) :
            self.settings.spheres = False
            self.settings_panel.spheres_ctrl.SetValue(False)
    self.viewer.set_miller_array(array)
    self.viewer.Refresh()
    self.viewer.fit_into_viewport()
    if (self.view_2d is not None) :
      self.view_2d.set_miller_array(array)

  def update_settings (self, *args, **kwds) :
    self.viewer.update_settings(*args, **kwds)

  def load_reflections_file (self, file_name) :
    if (file_name != "") :
      from iotbx.reflection_file_reader import any_reflection_file
      from iotbx.gui_tools.reflections import get_array_description
      try :
        hkl_file = any_reflection_file(file_name)
      except Exception, e :
        raise Sorry(str(e))
      arrays = hkl_file.as_miller_arrays(merge_equivalents=False)
      #arrays = f.file_server.miller_arrays
      valid_arrays = []
      array_info = []
      for array in arrays :
        if array.is_complex_array() or array.is_hendrickson_lattman_array() :
          continue
        labels = array.info().label_string()
        desc = get_array_description(array)
        array_info.append("%s (%s)" % (labels, desc))
        valid_arrays.append(array)
      if (len(valid_arrays) == 0) :
        raise Sorry("No arrays of the supported types in this file.")
      elif (len(valid_arrays) == 1) :
        self.set_miller_array(valid_arrays[0])
      else :
        #dlg = SelectArrayDialog(self, -1, "Select data")
        dlg = wx.SingleChoiceDialog(parent=self,
          message="Please select the data you wish to view:",
          caption="Select data",
          choices=array_info)
        if (dlg.ShowModal() == wx.ID_OK) :
          sel = dlg.GetSelection()
          self.set_miller_array(valid_arrays[sel])
        wx.CallAfter(dlg.Destroy)

  def OnLoadFile (self, evt) :
    file_name = wx.FileSelector("Reflections file",
      wildcard="Reflection files (*.mtz, *.sca, *.hkl)|*.mtz;*.sca;*.hkl",
      default_path="",
      flags=wx.OPEN)
    self.load_reflections_file(file_name)

  def OnSave (self, evt) :
    output_file = wx.FileSelector("Save image as:",
      default_filename="hklview.png",
      wildcard="PNG image (*.png)|*.png",
      flags=wx.SAVE)
    if (output_file != "") :
      self.viewer.save_screen_shot(file_name=output_file,
        extensions=["png"])

  def OnShow2D (self, evt) :
    if (self.view_2d is None) :
      self.view_2d = HKLViewFrame2D(self, -1, "2D data viewer")
      self.view_2d.set_miller_array(self.viewer.miller_array)
      self.view_2d.Show()
    self.view_2d.Raise()

  def OnClearLabels (self, evt) :
    self.viewer.clear_labels()

  def OnClose (self, event) :
    self.Destroy()
    event.Skip()

  def OnDestroy (self, event) :
    if (self.parent is not None) :
      self.parent.view_3d = None
    event.Skip()

class settings_window_2d (settings_window) :
  is_3d_view = False

class HKLViewFrame2D (HKLViewFrame) :
  def create_viewer_panel (self) :
    self.viewer = view_2d.hklview_2d(self, -1, size=(640,640))
    self.viewer.SetMinSize((640,640))

  def create_settings_panel (self) :
    self.settings.expand_to_p1 = True
    self.settings.expand_anomalous = True
    self.settings.slice_mode = True
    #self.settings.black_background = False
    self.settings_panel = settings_window_2d(self, -1, style=wx.RAISED_BORDER)

  def add_view_specific_functions (self) :
    item = wx.MenuItem(self.file_menu, -1, "Show 3D view")
    self.file_menu.AppendItem(item)
    self.Bind(wx.EVT_MENU, self.OnShow3D, item)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Show 3D view",
      bitmap=icons.hklview_3d.GetBitmap(),
      shortHelp="Show 3D view",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnShow3D, btn)

  def update_settings_for_merged (self) :
    self.settings.expand_to_p1 = True
    self.settings.expand_anomalous = True

  def OnClose (self, evt) :
    self.Destroy()
    evt.Skip()

  def OnDestroy (self, event) :
    if (self.parent is not None) :
      self.parent.view_2d = None
    event.Skip()

  def OnShow2D (self, evt) :
    pass

  def OnShow3D (self, evt) :
    if (self.view_3d is None) :
      self.view_3d = HKLViewFrame(self, -1, "3D data viewer")
      self.view_3d.Show()
      if (self.viewer.miller_array is not None) :
        self.view_3d.set_miller_array(self.viewer.miller_array)
    self.view_3d.Raise()
