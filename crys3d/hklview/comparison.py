from __future__ import absolute_import, division, print_function

from crys3d.hklview.frames import *
from crys3d.hklview import view_2d
import wx
import os

class twin_settings_window (settings_window) :
  def add_value_widgets (self, sizer) :
    sizer.SetRows(4)
    sizer.Add(wx.StaticText(self.panel, -1, "Value 1:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.value_info_1 = wx.TextCtrl(self.panel, -1, size=(80,-1),
      style=wx.TE_READONLY)
    sizer.Add(self.value_info_1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    sizer.Add(wx.StaticText(self.panel, -1, "Value 2:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.value_info_2 = wx.TextCtrl(self.panel, -1, size=(80,-1),
      style=wx.TE_READONLY)
    sizer.Add(self.value_info_2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def update_reflection_info (self, hkl, d_min, value_1, value_2) :
    print(hkl, value_1, value_2)
    if (hkl is None) :
      self.hkl_info.SetValue("")
      self.d_min_info.SetValue("")
      self.value_info.SetValue("")
    else :
      self.hkl_info.SetValue("%d, %d, %d" % hkl)
      d_min_str = format_value("%.3g", d_min)
      self.d_min_info.SetValue(d_min_str)
      value_str_1 = format_value("%.3g", value_1, replace_none_with="---")
      self.value_info_1.SetValue(value_str_1)
      value_str_2 = format_value("%.3g", value_2, replace_none_with="---")
      self.value_info_2.SetValue(value_str_2)

class twin_viewer_panel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self.settings = self.GetParent().settings
    szr = wx.BoxSizer(wx.HORIZONTAL)
    self.SetSizer(szr)
    self.view_1 = view_2d.hklview_2d(self, -1, size=(480,480))
    self.view_2 = view_2d.hklview_2d(self, -1, size=(480,480))
    self.view_1.SetMinSize((480,480))
    self.view_2.SetMinSize((480,480))
    szr.Add(self.view_1, 1, wx.EXPAND)
    szr.Add(self.view_2, 1, wx.EXPAND)
    self.SetMinSize((960,480))
    szr.Fit(self)

  def add_view_specific_functions (self) :
    pass

  def _propagate_action (self, name, *args, **kwds) :
    for viewer in [self.view_1, self.view_2] :
      getattr(viewer, name)(*args, **kwds)

  def clear_labels (self) :
    self._propagate_action("clear_labels")

  def set_miller_arrays (self, array1, array2) :
    self.view_1.set_miller_array(array1)
    self.view_2.set_miller_array(array2)
    self.Refresh()

  def update_settings (self, *args, **kwds) :
    self.view_1.update_settings(*args, **kwds)
    self.view_2.update_settings(*args, **kwds)
    self.Refresh()

  def update_clicked (self, index) :
    hkl_1, d_min_1, value_1 = self.view_1.scene.get_reflection_info(index)
    hkl_2, d_min_2, value_2 = self.view_2.scene.get_reflection_info(index)
    assert (hkl_1 == hkl_2)
    self.GetParent().update_reflection_info(
      hkl=hkl_1,
      d_min=d_min_1,
      value_1=value_1,
      value_2=value_2)

  def save_screen_shot (self, file_name, extensions=None) :
    base, ext = os.path.splitext(file_name)
    file_1 = base + "_1" + ext
    file_2 = base + "_2" + ext
    self.view_1.save_screen_shot(file_1)
    self.view_2.save_screen_shot(file_2)

class ComparisonFrame (HKLViewFrame) :
  def __init__ (self, *args, **kwds) :
    HKLViewFrame.__init__(self, *args, **kwds)
    self._array1 = None
    self._array2 = None

  def create_viewer_panel (self) :
    self.viewer = twin_viewer_panel(self)
    self.viewer.SetMinSize((960, 480))

  def update_reflection_info (self, *args, **kwds) :
    self.settings_panel.update_reflection_info(*args, **kwds)

  def create_settings_panel (self) :
    self.settings.expand_to_p1 = True
    self.settings.expand_anomalous = True
    self.settings.slice_mode = True
    #self.settings.black_background = False
    self.settings_panel = twin_settings_window(self, style=wx.RAISED_BORDER)

  def update_settings_for_merged (self) :
    self.settings.expand_to_p1 = True
    self.settings.expand_anomalous = True

  def update_settings (self, *args, **kwds) :
    if (None in [self._array1, self._array2]) :
      return False
    self.viewer.update_settings(*args, **kwds)

  def SetupMenus (self) :
    self.menubar = wx.MenuBar(-1)
    self.file_menu = wx.Menu()
    self.menubar.Append(self.file_menu, "File")
    item = wx.MenuItem(self.file_menu, -1, "Load data...\tCtrl-O")
    self.Bind(wx.EVT_MENU, self.OnLoadFile, item)
    self.file_menu.AppendItem(item)

  def OnLoadFile (self, evt) :
    file_name_1 = wx.FileSelector("Reflections file 1",
      wildcard="Reflection files (*.mtz, *.sca, *.hkl)|*.mtz;*.sca;*.hkl",
      default_path="",
      flags=wx.FD_OPEN)
    file_name_2 = wx.FileSelector("Reflections file 2",
      wildcard="Reflection files (*.mtz, *.sca, *.hkl)|*.mtz;*.sca;*.hkl",
      default_path="",
      flags=wx.FD_OPEN)
    self.load_files(file_name_1, file_name_2)

  def load_files (self, file_name_1, file_name_2) :
    array1 = self.load_reflections_file(file_name_1, set_array=False,
      data_only=True)
    array2 = self.load_reflections_file(file_name_2, set_array=False,
      data_only=True)
    symm1 = array1.crystal_symmetry()
    symm2 = array2.crystal_symmetry()
    if (symm1 is None) :
      raise Sorry(("No crystal symmetry found in %s!  Please convert to a "+
        "more appropriate format.") % file_name_1)
    if (symm2 is None) :
      raise Sorry(("No crystal symmetry found in %s!  Please convert to a "+
        "more appropriate format.") % file_name_2)
    if (symm1.unit_cell() is None) :
      symm1 = symm1.customized_copy(unit_cell=symm2.unit_cell())
      array1 = array1.customized_copy(crystal_symmetry=symm1)
    if (symm2.unit_cell() is None) :
      symm2 = symm2.customized_copy(unit_cell=symm1.unit_cell())
      array2 = array2.customized_copy(crystal_symmetry=symm2)
    if (not array1.is_similar_symmetry(array2)) :
      from cctbx import crystal
      space_group_1 = array1.space_group_info()
      space_group_2 = array2.space_group_info()
      if (str(space_group_1) != str(space_group_2)) :
        # TODO need to figure out if these are really incompatible!
        confirm = wx.MessageBox(("The space groups for the two datasets are "+
          "different (%s versus %s).  The space group from the first dataset "+
          "will be used for both.") % (space_group_1, space_group_2),
          style=wx.OK)
      unit_cell_1 = array1.unit_cell()
      unit_cell_2 = array2.unit_cell()
      if (not unit_cell_1.is_similar_to(unit_cell_2)) :
        uc_str_1 = "%g %g %g %g %g %g" % unit_cell_1.parameters()
        uc_str_2 = "%g %g %g %g %g %g" % unit_cell_2.parameters()
        confirm = wx.MessageBox(("The unit cells for the two datasets are "+
          "different (%s versus %s).  The unit cell from the first dataset "+
          "will be used to calculated the resolution of clicked reflections.")%
            (uc_str_1, uc_str_2),
          style=wx.OK)
      symm = crystal.symmetry(
        space_group_info=space_group_1,
        unit_cell=unit_cell_1)
      array2 = array2.customized_copy(crystal_symmetry=symm)
    if (array1.anomalous_flag() != array2.anomalous_flag()) :
      wx.MessageBox("Warning: only one array contains anomalous data; to "+
        "allow comparison, Bijvoet mates will be generated for the "+
        "non-anomalous array.")
      if (not array1.anomalous_flag()) :
        array1 = array1.generate_bijvoet_mates()
      else :
        array2 = array2.generate_bijvoet_mates()
    array1 = array1.common_set(other=array2)
    array2 = array2.common_set(other=array1)
    is_intensities = [ array1.is_xray_intensity_array(),
                       array2.is_xray_intensity_array() ]
    if (len(set(is_intensities)) == 2) :
      convert = wx.MessageBox("You appear to be comparing intensities with "+
        "another type of data.  Do you want to convert the intensities to "+
        "amplitudes?  If you leave them as intensities the program will "+
        "still run, but the scale of the data may be much different.",
        style=wx.YES_NO)
      if (convert == wx.YES) :
        if (array1.is_xray_intensity_array()) :
          array1 = array1.f_sq_as_f()
        else :
          array2 = array2.f_sq_as_f()
    self._array1 = array1
    self._array2 = array2
    self.settings_panel.d_min_ctrl.SetValue(array1.d_min())
    self.settings_panel.d_min_ctrl.SetRange(array1.d_min(), 20.0)
    self.settings_panel.set_index_span(array1.index_span())
    self.settings_panel.update_space_group_choices(array1)
    self.viewer.set_miller_arrays(array1, array2)

  def add_view_specific_functions (self) :
    pass

  def delete_miller_index (self, hkl) :
    self._array1 = self._array1.delete_index(hkl)
    self._array2 = self._array2.delete_index(hkl)
    self.viewer.set_miller_arrays(self._array1, self._array2)
