
from __future__ import absolute_import, division, print_function
from mmtbx.validation import ramalyze
from mmtbx.validation import rotalyze
from mmtbx.validation import graphics
import mmtbx.validation.utils
from wxtbx import utils
import wxtbx.plots
import wx
from math import sqrt, floor
from six.moves import range

class rotarama_plot(wxtbx.plots.plot_container):
  hit_test_radius = 3.0
  hit_test_minimum_difference = 0.5
  def __init__(self, parent, figure_size=(8,8), xyz_shift=None):
    wxtbx.plots.plot_container.__init__(self,
      parent=parent,
      figure_size=figure_size,
      handle_left_click=True)
    self.xyz_shift = xyz_shift

  def show_plot(self, *args, **kwds):
    if self.disabled : return
    self.draw_plot(*args, **kwds)
    self.parent.Refresh()

  def process_mouse_click(self, mpl_event):
    (xdata, ydata) = (mpl_event.xdata, mpl_event.ydata)
    if xdata is None or ydata is None :
      return False
    min_dist = self.hit_test_radius
    closest_point = None
    for i, (x, y) in enumerate(self._points):
      dist = sqrt((xdata - x)**2 + (ydata - y)**2)
      if dist < min_dist :
        closest_point = i
        min_dist = dist
    if closest_point is not None :
      xyz = self._xyz[closest_point]
      if (self.xyz_shift is not None):
        xyz = (xyz[0] + self.xyz_shift[0],
               xyz[1] + self.xyz_shift[1],
               xyz[2] + self.xyz_shift[2])
      self.parent.zoom_callback(xyz=xyz)

class ramalyze_plot(rotarama_plot,
                     ramalyze.ramachandran_plot_mixin):
  def __init__(self, *args, **kwds):
    rotarama_plot.__init__(self, *args, **kwds)
    ramalyze.ramachandran_plot_mixin.__init__(self)

class rotalyze_plot(rotarama_plot, rotalyze.rotamer_plot_mixin):
  def __init__(self, *args, **kwds):
    rotarama_plot.__init__(self, *args, **kwds)
    rotalyze.rotamer_plot_mixin.__init__(self)

class rotarama_frame(wxtbx.plots.plot_frame):
  frame_name = "rotarama_frame"
  show_controls_default = True
  def __init__(self, parent, title, validation):
    self._validation = validation
    wxtbx.plots.plot_frame.__init__(self,
      parent=parent,
      title=title,
      style=wx.DEFAULT_FRAME_STYLE)
    self._map_cache = {}
    self._point_cache = {}
    self._xyz_cache = {}
    self.OnUpdatePlot(None)

  def OnDestroy(self, event):
    if self and hasattr(self.GetParent(), self.frame_name):
      setattr(self.GetParent(), self.frame_name, None)

  def OnUpdatePlot(self, event):
    pass

  def zoom_callback(self, **kwds):
    if getattr(self.GetParent(), "main_window", None) is not None :
      self.GetParent().main_window.show_gfx_selection(**kwds)

class ramalyze_frame(rotarama_frame):
  frame_name = "rama_frame"
  def draw_top_panel(self):
    self.top_panel = wx.Panel(self)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.top_panel.SetSizer(szr)
    grid = wx.FlexGridSizer(cols=4)
    szr.Add(grid, 0, wx.ALL)
    grid.Add(utils.bold_text(self.top_panel, "Position type:"), 0,
      utils.std_sizer_flags, 5)
    pos_choice = wx.Choice(parent=self.top_panel,
      choices=ramalyze.res_type_labels)
    self.Bind(wx.EVT_CHOICE, self.OnUpdatePlot, pos_choice)
    pos_choice.SetSelection(ramalyze.RAMA_GENERAL)
    grid.Add(pos_choice, 0, utils.std_sizer_flags, 5)
    grid.Add(utils.bold_text(self.top_panel, "Residue name:"), 0,
      utils.std_sizer_flags, 5)
    aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'LEU',
          'LYS', 'MET', 'PHE', 'SER', 'THR', 'TRP', 'TYR',]
    res_choice = wx.Choice(parent=self.top_panel,
      choices=['*']+aa)
    res_choice.SetSelection(0)
    self.Bind(wx.EVT_CHOICE, self.OnUpdatePlot, res_choice)
    res_choice.Disable()
    grid.Add(res_choice, 0, utils.std_sizer_flags, 5)
    grid.Add(utils.bold_text(self.top_panel, "Show data points:"), 0,
      utils.std_sizer_flags, 5)
    default = ramalyze.RAMALYZE_ANY
    pt_choice = wx.Choice(parent=self.top_panel,
      choices=ramalyze.rama_type_labels)
    pt_choice.SetSelection(default)
    self.Bind(wx.EVT_CHOICE, self.OnUpdatePlot, pt_choice)
    grid.Add(pt_choice, 0, utils.std_sizer_flags, 5)
    grid.Add(utils.bold_text(self.top_panel, "Color scheme:"), 0,
      utils.std_sizer_flags, 5)
    cm_choice = wx.Choice(parent=self.top_panel,
      choices=wxtbx.plots.colormap_names)
    cm_choice.SetStringSelection("Blue")
    self.Bind(wx.EVT_CHOICE, self.OnUpdatePlot, cm_choice)
    grid.Add(cm_choice, 0, utils.std_sizer_flags, 5)
    self.choices = [pos_choice, res_choice, pt_choice, cm_choice]
    self.show_labels = wx.CheckBox(self.top_panel, label="Show labels")
    self.show_labels.SetValue(True)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdatePlot, self.show_labels)
    grid.Add(self.show_labels, 0, utils.std_sizer_flags, 5)

  def get_current_state(self):
    pos_type = self.choices[0].GetSelection()
    if (pos_type != ramalyze.RAMA_GENERAL):
      self.choices[1].SetSelection(0)
      self.choices[1].Disable()
    else :
      self.choices[1].Enable()
    res_type = self.choices[1].GetStringSelection()
    pt_type = self.choices[2].GetSelection()
    cm_name = self.choices[3].GetStringSelection()
    show_labels = self.show_labels.GetValue()
    return (pos_type, res_type, pt_type, cm_name, show_labels)

  def create_plot_panel(self):
    panel = ramalyze_plot(self)
    return panel

  def set_plot_type(self,
      pos_type,
      res_type="*",
      pt_type=ramalyze.RAMALYZE_ANY,
      cm_name='jet',
      show_labels=True):
    if not pos_type in self._map_cache :
      z_data = mmtbx.validation.utils.get_rotarama_data(
        pos_type=ramalyze.res_types[pos_type],
        convert_to_numpy_array=True)
      self._map_cache[pos_type] = z_data
    if not res_type in self._point_cache :
      pass
    title = ramalyze.format_ramachandran_plot_title(
      position_type=pos_type,
      residue_type=res_type)
    points = None
    coords = None
    if not (pos_type, res_type, pt_type) in self._point_cache :
      (points, coords) = self._validation.get_plot_data(
        position_type=pos_type,
        residue_name=res_type,
        point_type=pt_type)
      self._point_cache[(pos_type, res_type, pt_type)] = points
      self._xyz_cache[(pos_type, res_type, pt_type)] = coords
    else :
      points = self._point_cache[(pos_type, res_type, pt_type)]
      coords = self._xyz_cache[(pos_type, res_type, pt_type)]
    contours = ramalyze.get_contours(pos_type)
    self.plot_panel.show_plot(
      stats=self._map_cache[pos_type],
      title=title,
      points=points,
      show_labels=show_labels,
      colormap=wxtbx.plots.colormap_id_dict[cm_name],
      contours=contours,
      xyz=coords)

  def OnUpdatePlot(self, event):
    (pos_type,res_type,pt_type,cm_name,show_labels) = self.get_current_state()
    self.set_plot_type(pos_type, res_type, pt_type, cm_name, show_labels)

class rotalyze_frame(rotarama_frame):
  frame_name = "rota_frame"
  residues = ["Arg", "Asn", "Asp", "Gln", "Glu", "His", "Ile", "Leu", "Lys",
              "Met", "Phe", "Trp", "Tyr"]

  def draw_top_panel(self):
    self.top_panel = wx.Panel(self)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.top_panel.SetSizer(szr)
    grid = wx.FlexGridSizer(cols=4)
    szr.Add(grid, 0, wx.ALL)
    grid.Add(utils.bold_text(self.top_panel, "Residue name:"), 0,
      utils.std_sizer_flags, 5)
    res_choice = wx.Choice(parent=self.top_panel, choices=self.residues)
    res_choice.SetSelection(0)
    self.Bind(wx.EVT_CHOICE, self.OnUpdatePlot, res_choice)
    grid.Add(res_choice, 0, utils.std_sizer_flags, 5)
    grid.Add(utils.bold_text(self.top_panel, "Show data points:"), 0,
      utils.std_sizer_flags, 5)
    pt_names = ["All", "None", "Outlier"]
    pt_choice = wx.Choice(self.top_panel, choices=pt_names)
    pt_choice.SetSelection(0)
    self.Bind(wx.EVT_CHOICE, self.OnUpdatePlot, pt_choice)
    grid.Add(pt_choice, 0, utils.std_sizer_flags, 5)
    grid.Add(utils.bold_text(self.top_panel, "Color scheme:"), 0,
      utils.std_sizer_flags, 5)
    cm_choice = wx.Choice(parent=self.top_panel,
      choices=wxtbx.plots.colormap_names)
    cm_choice.SetSelection(5)
    self.Bind(wx.EVT_CHOICE, self.OnUpdatePlot, cm_choice)
    grid.Add(cm_choice, 0, utils.std_sizer_flags, 5)
    self.choices = [res_choice, pt_choice, cm_choice]
    self.show_labels = wx.CheckBox(self.top_panel, label="Show labels")
    self.show_labels.SetValue(False)
    self.Bind(wx.EVT_CHECKBOX, self.OnUpdatePlot, self.show_labels)
    grid.Add(self.show_labels, 0, utils.std_sizer_flags, 5)

  def get_current_state(self):
    res_type = self.choices[0].GetStringSelection()
    pt_type = self.choices[1].GetStringSelection()
    cm_name = self.choices[2].GetStringSelection()
    show_labels = self.show_labels.GetValue()
    return (res_type, pt_type, cm_name, show_labels)

  def create_plot_panel(self):
    panel = rotalyze_plot(self)
    return panel

  def set_plot_type(self, res_type, pt_type="All", cm_name='jet',
      show_labels=False):
    if not res_type in self._map_cache :
      z_data = mmtbx.validation.utils.get_rotarama_data(
        residue_type=res_type,
        db="rota",
        convert_to_numpy_array=True)
      self._map_cache[res_type] = z_data
    title = "Chi1-Chi2 plot for %s" % res_type
    points = coords = None
    if (pt_type != "None"):
      if not (res_type, pt_type) in self._point_cache :
        (points, coords) = self._validation.get_plot_data(
          residue_name=res_type.upper(),
          point_type=pt_type)
        # shift chi2 values by 180 to fit in contours
        if (res_type.lower() in ["asp", "phe", "tyr"]):
          for i in range(len(points)):
            if (points[i][1] > 180.0):
              point = list(points[i])
              point[1] -= 180.0
              points[i] = tuple(point)
        self._point_cache[(res_type, pt_type)] = points
        self._xyz_cache[(res_type, pt_type)] = coords
      else :
        points = self._point_cache[(res_type, pt_type)]
        coords = self._xyz_cache[(res_type, pt_type)]
    extent = y_marks = None
    if (res_type.lower() in ["asp", "phe", "tyr"]):
      extent = [0, 360, 0, 180]
      y_marks = [90]
    self.plot_panel.show_plot(
      stats=self._map_cache[res_type],
      title=title,
      points=points,
      show_labels=show_labels,
      colormap=wxtbx.plots.colormap_id_dict[cm_name],
      contours=None,
      xyz=coords,
      extent=extent,
      y_marks=y_marks)

  def OnUpdatePlot(self, event):
    (res_type,pt_type,cm_name,show_labels) = self.get_current_state()
    self.set_plot_type(res_type, pt_type, cm_name, show_labels)

class multi_criterion_plot(wxtbx.plots.plot_container,
    graphics.multi_criterion_plot_mixin):
  def __init__(self, parent, binner, y_limits):
    self._reset = False
    graphics.multi_criterion_plot_mixin.__init__(self,
      binner=binner,
      y_limits=y_limits)
    wxtbx.plots.plot_container.__init__(self,
      parent=parent,
      figure_size=(16,9),
      handle_left_click=True)
    self.figure.subplots_adjust(hspace=0.001)
    self.canvas.mpl_connect('motion_notify_event', self.OnHover)

  def show_plot(self, *args, **kwds):
    if self.disabled : return
    self.draw_plot(*args, **kwds)
    self.parent.Refresh()

  def OnPickRange(self, event):
    bin = self.choose_range.GetSelection()
    self.plot_range(bin)

  def OnHover(self, mpl_event):
    (xdata, ydata) = (mpl_event.xdata, mpl_event.ydata)
    if xdata is None or ydata is None :
      if self._reset :
        self.parent.residue_status.SetValue("")
      self._reset = False
      return False
    idx = int(floor(xdata))
    try :
      residue = self._current_bin.get_selected(index=idx)
    except IndexError :
      pass
    else :
      if (residue is not None):
        self.parent.residue_status.SetValue(residue.id_str())
        self._reset = True

  def process_mouse_click(self, mpl_event):
    (xdata, ydata) = (mpl_event.xdata, mpl_event.ydata)
    if xdata is None or ydata is None :
      return False
    idx = int(floor(xdata))
    selection_string = None
    xyz = None
    try :
      residue = self._current_bin.get_selected(idx)
      if (residue is not None):
        selection_string = residue.id_str()
        xyz = residue.xyz
    except IndexError :
      pass
    else :
      if ( (selection_string is not None) and (xyz is not None) ):
        self.parent.zoom_callback(selection_string=selection_string,
                                  xyz=xyz)

class multi_criterion_frame(wxtbx.plots.plot_frame):
  show_controls_default = True
  def __init__(self, parent, title, validation):
    self._binner = validation.binned_data()
    self._y_limits = validation.get_y_limits()
    wxtbx.plots.plot_frame.__init__(self,
      parent=parent,
      title=title,
      style=wx.DEFAULT_FRAME_STYLE)
    self._map_cache = {}
    self._point_cache = {}
    self._xyz_cache = {}
    self.OnUpdatePlot(None)

  def create_plot_panel(self):
    panel = multi_criterion_plot(self,
      binner=self._binner,
      y_limits=self._y_limits)
    return panel

  def draw_top_panel(self):
    self.top_panel = wx.Panel(self)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.top_panel.SetSizer(szr)
    grid = wx.FlexGridSizer(cols=4)
    szr.Add(grid, 0, wx.ALL)
    grid.Add(utils.bold_text(self.top_panel, "Show residues:"), 0,
      utils.std_sizer_flags, 5)
    self.range_choice = wx.Choice(parent=self.top_panel,
      choices=self._binner.get_ranges())
    self.range_choice.SetSelection(0)
    self.Bind(wx.EVT_CHOICE, self.OnUpdatePlot, self.range_choice)
    grid.Add(self.range_choice, 0, utils.std_sizer_flags, 5)
    grid.Add(utils.bold_text(self.top_panel, "Click to zoom residue:"), 0,
      utils.std_sizer_flags, 5)
    self.residue_status = wx.TextCtrl(self.top_panel, size=(160,-1),
      style=wx.TE_READONLY)
    grid.Add(self.residue_status, 0, utils.std_sizer_flags, 5)

  def OnUpdatePlot(self, event):
    bin = self.range_choice.GetSelection()
    self.plot_panel.plot_range(bin)

  def zoom_callback(self, **kwds):
    if getattr(self.GetParent(), "main_window", None) is not None :
      self.GetParent().main_window.show_gfx_selection(**kwds)
