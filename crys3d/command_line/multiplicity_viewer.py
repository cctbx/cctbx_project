# LIBTBX_SET_DISPATCHER_NAME cctbx.multiplicity_viewer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import absolute_import, division, print_function
from crys3d.hklview.frames import *
from wxtbx import icons
import wxtbx.app
import iotbx.phil
import libtbx.load_env
from cctbx import crystal
import wx
if wx.VERSION >= (4, 0):
  import wx.adv

import sys
from six.moves import range
from six.moves import zip


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
  if wx.VERSION >= (4,0):
    app_icon = wx.Icon()
  else:
    app_icon = wx.EmptyIcon()
  app_icon.CopyFromBitmap(icons.hklview_3d.GetBitmap())
  if wx.VERSION >= (4,0):
    tb_icon = wx.adv.TaskBarIcon(wx.adv.TBI_DOCK)
  else:
    tb_icon = wx.TaskBarIcon(wx.TBI_DOCK)
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
  def __init__ (self, *args, **kwds) :
    HKLViewFrame.__init__(self, *args, **kwds)
    self.create_colour_bar_panel()
    self.sizer.Add(self.colour_bar, 0, wx.EXPAND)
    self.SetSizer(self.sizer)
    self.sizer.SetSizeHints(self)
    self.viewer.SetFocus()

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
    self.colour_bar.Refresh()
    if (self.view_2d is not None) :
      self.view_2d.set_miller_array(array)

  def load_reflections_file (self, file_name, set_array=True,
      data_only=True) :
    if (file_name != "") :
      from iotbx.reflection_file_reader import any_reflection_file
      from iotbx.gui_tools.reflections import get_array_description
      try :
        hkl_file = any_reflection_file(file_name)
      except Exception as e :
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

  def create_colour_bar_panel(self):
    self.colour_bar = ColourBar(parent=self, size=(100,640))
    self.colour_bar.SetMinSize((100,640))
    self.colour_bar.viewer = self.viewer

  def OnActive (self, event) :
    HKLViewFrame.OnActive(self, event)


  def update_settings (self, *args, **kwds) :
    if (self.miller_array is None) :
      return False
    self.viewer.update_settings(*args, **kwds)
    self.colour_bar.update_settings(*args, **kwds)


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
    ctrls = self.create_controls(
      setting="show_anomalous_pairs",
      label="Show anomalous pairs")
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

  def update_reflection_info (self, hkl, d_min, value) :
    return


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


class ColourBar(wx.Panel):
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    font = wx.Font(14, wx.MODERN, wx.NORMAL, wx.NORMAL)
    self.SetFont(font)
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    #self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftClick, self)
    #self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp, self)
    #self.Bind(wx.EVT_MOTION, self.OnMouseMotion, self)
    self.viewer = None
    self.parent = self.GetParent()
    self.settings = self.parent.settings
    self.setup_colors()

  @property
  def scene(self):
    if self.viewer is not None:
      return self.viewer.scene

  def paint(self, gc):
    font = self.GetFont()
    font.SetFamily(wx.FONTFAMILY_MODERN)
    if (self.settings.black_background) :
      gc.SetFont(gc.CreateFont(font, (255,255,255)))
    else :
      gc.SetFont(gc.CreateFont(font, (0,0,0)))
    self.render(gc)

  def render (self, canvas) :
    from scitbx.array_family import flex
    from libtbx.utils import frange
    import math
    size = self.GetSize()
    border = 10
    i_rows = flex.double_range(border, size[1]-border)
    scene = self.scene
    if self.scene.settings.scale_colors_multiplicity:
      data = self.scene.multiplicities.data()
    else:
      data = self.scene.data
      if self.settings.sqrt_scale_colors:
        data = flex.sqrt(data)
    min_data = flex.min(data)
    max_data = flex.max(data)
    data_for_colors = flex.double(frange(
      max_data, min_data, -(max_data-min_data)/len(i_rows)))
    tick_step = int(math.ceil((max_data-min_data)/10))
    i_row_ticks = []
    tick_text = []
    start_tick = math.floor(max_data)
    i_tick = 0
    for i in range(len(data_for_colors)-1):
      tick_d = start_tick - tick_step * i_tick
      if abs(data_for_colors[i]-tick_d) < abs(data_for_colors[i+1]-tick_d):
        i_row_ticks.append(i_rows[i])
        tick_text.append(str(int(tick_d)))
        i_tick += 1
    tick_d = start_tick - tick_step * i_tick
    if tick_d == min_data:
      i_row_ticks.append(i_rows[-1])
      tick_text.append(str(int(tick_d)))

    from scitbx import graphics_utils
    if (self.settings.color_scheme in ["rainbow", "heatmap", "redblue"]) :
      colors = graphics_utils.color_by_property(
        properties=data_for_colors,
        selection=flex.bool(data_for_colors.size(), True),
        color_all=False,
        gradient_type=self.settings.color_scheme)
    elif (self.settings.color_scheme == "grayscale") :
      colors = graphics_utils.grayscale_by_property(
        properties=data_for_colors,
        selection=flex.bool(data_for_colors.size(), True),
        shade_all=False,
        invert=self.settings.black_background)
    else :
      if (self.settings.black_background) :
        base_color = (1.0,1.0,1.0)
      else :
        base_color = (0.0,0.0,0.0)
      colors = flex.vec3_double(data_for_colors.size(), base_color)

    l_padding = border
    r_padding = 4 * border

    for i_row, color in zip(i_rows, colors):
      self.draw_line(canvas, l_padding, i_row, size[0]-r_padding, i_row, color=color)

    for i_row, text in zip(i_row_ticks, tick_text):
      self.draw_text(canvas, text, size[0]-0.8*r_padding, i_row-5)
      self.draw_line(canvas, size[0]-r_padding-10, i_row, size[0]-r_padding, i_row)

  def setup_colors (self) :
    if (self.settings.black_background) :
      self._background = (0.,0.,0.)
      self._foreground = (0.95,0.95,0.95)
      if (self.settings.color_scheme == "heatmap") :
        self._missing = (0.,1.,0.)
      elif (not self.settings.color_scheme in ["rainbow", "redblue"]) :
        self._missing = (1.,0.,0.)
      else :
        self._missing = (1.,1.,1.)
    else :
      self._background = (1.,1.,1.)
      self._foreground = (0.,0.,0.)
      if (self.settings.color_scheme == "heatmap") :
        self._missing = (0.,1.,0.)
      elif (not self.settings.color_scheme in ["rainbow", "redblue"]) :
        self._missing = (1.,0.,0.)
      else :
        self._missing = (0.,0.,0.)

  def get_color (self, c) :
    return (int(c[0]*255), int(c[1]*255), int(c[2]*255))

  def draw_line (self, canvas, x1, y1, x2, y2, color=None) :
    gc = canvas
    x_axis = gc.CreatePath()
    x_axis.MoveToPoint(x1, y1)
    x_axis.AddLineToPoint(x2, y2)
    x_axis.CloseSubpath()
    if (color is None) :
      color = self._foreground
    gc.SetPen(wx.Pen(self.get_color(color)))
    gc.PushState()
    gc.StrokePath(x_axis)
    gc.PopState()

  def draw_text (self, canvas, text, x, y) :
    gc = canvas
    gc.SetPen(wx.Pen(self.get_color(self._foreground)))
    gc.DrawText(text, x, y)

  def OnPaint (self, event) :
    if (self.scene is None) :
      return
    if (self.settings.black_background) :
      self.SetBackgroundColour((0,0,0))
    else :
      self.SetBackgroundColour((255,255,255))
    dc = wx.AutoBufferedPaintDCFactory(self)
    if (self.settings.black_background) :
      dc.SetBackground(wx.BLACK_BRUSH)
    else :
      dc.SetBackground(wx.WHITE_BRUSH)
    dc.Clear()
    gc = wx.GraphicsContext.Create(dc)
    self.paint(gc)

  def update_settings (self) :
    self.Refresh()

if (__name__ == "__main__") :
  run(sys.argv[1:])
