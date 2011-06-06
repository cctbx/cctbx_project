
# TODO:
#  - cached scenes

from crys3d.hklview import settings
from crys3d.hklview import view_2d, view_3d
from wxtbx import icons
import wxtbx.utils
from wx.lib.agw import floatspin
import wx
from libtbx.utils import Sorry
from math import sqrt

class settings_window (wxtbx.utils.SettingsPanel) :
  is_3d_view = True
  def add_controls (self) :
    self.d_min_ctrl = floatspin.FloatSpin(parent=self, increment=0.05, digits=2)
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
    ctrls = self.create_controls(
      setting="uniform_size",
      label="Use same radius for all points")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sqrt_scale_radii",
      label="Scale radii to sqrt(I)")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sqrt_scale_colors",
      label="Scale colors to sqrt(I)")
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
        setting="display_as_spheres",
        label="Display reflections as spheres")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="color_scheme",
      label="Color scheme",
      captions=["Rainbow","Grayscale","Monochrome"])
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    box.Add(ctrls[0], 0, wx.TOP|wx.BOTTOM|wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(ctrls[1], 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    ctrls = self.create_controls(
      setting="show_missing_reflections",
      label="Show missing reflections")
    ctrls2 = self.create_controls(
      setting="show_only_missing",
      label="only")
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    box.Add(ctrls[0], 0, wx.ALL, 5)
    box.Add(ctrls2[0], 0, wx.ALL, 5)
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
    self.hkl_choice.SetSelection(self.settings.slice_axis)
    box2.Add(self.hkl_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box2.Add(wx.StaticText(self.panel, -1, "="), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.slice_index = wx.SpinCtrl(self.panel, -1)
    self.slice_index.SetValue(self.settings.slice_index)
    box2.Add(self.slice_index)
    self.panel_sizer.Add(box2)
    self.Bind(wx.EVT_CHOICE, self.OnSetSlice, self.hkl_choice)
    self.Bind(wx.EVT_SPINCTRL, self.OnSetSlice, self.slice_index)

  def OnSetSlice (self, event) :
    self.settings.slice_axis = self.hkl_choice.GetSelection()
    self.settings.slice_index = self.slice_index.GetValue()
    if (not self.is_3d_view) or (self.slice_ctrl.GetValue()) :
      try :
        self.parent.update_settings()
      except ValueError, e : # TODO set limits
        raise Sorry(str(e))

  def OnChangeResolution (self, event) :
    self.settings.d_min = self.d_min_ctrl.GetValue()
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
    self.settings = settings()
    self.create_settings_panel()
    self.sizer.Add(self.settings_panel, 0, wx.EXPAND)
    self.create_viewer_panel()
    self.sizer.Add(self.viewer, 1, wx.EXPAND)
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
    self.add_view_specific_functions()
    self.SetMenuBar(menubar)
    self.toolbar.Realize()
    self.SetSizer(self.sizer)
    self.sizer.SetSizeHints(self)
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)

  def create_viewer_panel (self) :
    self.viewer = view_3d.hklview_3d(self, size=(800,600))

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
    self.statusbar.SetFieldsCount(2)

  def update_clicked (self, hkl, resolution=None) :
    if (hkl is not None) :
      self.statusbar.SetStatusText("CLICKED: %d,%d,%d (d = %g)" % (hkl[0],
        hkl[1], hkl[2], resolution), 1)
    else :
      self.statusbar.SetStatusText("", 1)

  def set_miller_array (self, array) :
    if array.is_complex_array() or array.is_hendrickson_lattman_array() :
      raise Sorry("Complex (map) data and HL coefficients not supported.")
    info = array.info()
    if isinstance(info, str) :
      labels = "TEST DATA"
    else :
      labels = info.label_string()
    details = ""
    if (not array.is_unique_set_under_symmetry()) :
      details = "(unmerged data)"
    sg = "%s" % array.space_group_info()
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % array.unit_cell().parameters()
    self.statusbar.SetStatusText("Data: %s %s (Space group: %s  Unit Cell: %s)"
      % (labels, details, sg, uc))
    self.settings_panel.d_min_ctrl.SetValue(array.d_min())
    self.settings_panel.d_min_ctrl.SetRange(array.d_min(), 20.0)
    self.viewer.set_miller_array(array)
    self.viewer.Refresh()
    self.viewer.fit_into_viewport()
    if (self.view_2d is not None) :
      self.view_2d.set_miller_array(array)

  def update_settings (self, *args, **kwds) :
    self.viewer.update_settings(*args, **kwds)

  def load_reflections_file (self, file_name) :
    if (file_name != "") :
      from iotbx import file_reader
      from iotbx.gui_tools.reflections import get_array_description
      f = file_reader.any_file(file_name, force_type="hkl")
      f.assert_file_type("hkl")
      arrays = f.file_object.as_miller_arrays(merge_equivalents=False)
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

  def OnDestroy (self, event) :
    if (self.parent is not None) :
      self.parent.view_3d = None

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

  def OnClose (self, evt) :
    self.Destroy()

  def OnDestroy (self, evt) :
    if (self.parent is not None) :
      self.parent.view_2d = None

  def OnShow2D (self, evt) :
    pass

  def OnShow3D (self, evt) :
    if (self.view_3d is None) :
      self.view_3d = HKLViewFrame(self, -1, "3D data viewer")
      self.view_3d.set_miller_array(self.viewer.miller_array)
      self.view_3d.Show()
    self.view_3d.Raise()
