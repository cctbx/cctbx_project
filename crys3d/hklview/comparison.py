
from crys3d.hklview.frames import *
from crys3d.hklview import view_2d
import wx
import os

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

  def update_clicked (self, *args, **kwds) :
    self.GetParent().update_clicked(*args, **kwds)

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

  def create_settings_panel (self) :
    self.settings.expand_to_p1 = True
    self.settings.expand_anomalous = True
    self.settings.slice_mode = True
    #self.settings.black_background = False
    self.settings_panel = settings_window_2d(self, -1, style=wx.RAISED_BORDER)

  def update_settings_for_merged (self) :
    self.settings.expand_to_p1 = True
    self.settings.expand_anomalous = True

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
      flags=wx.OPEN)
    file_name_2 = wx.FileSelector("Reflections file 2",
      wildcard="Reflection files (*.mtz, *.sca, *.hkl)|*.mtz;*.sca;*.hkl",
      default_path="",
      flags=wx.OPEN)
    self.load_files(file_name_1, file_name_2)

  def load_files (self, file_name_1, file_name_2) :
    array1 = self.load_reflections_file(file_name_1, set_array=False)
    array2 = self.load_reflections_file(file_name_2, set_array=False)
    array1 = array1.common_set(other=array2)
    array2 = array2.common_set(other=array1)
    self._array1 = array1
    self._array2 = array2
    self.settings_panel.d_min_ctrl.SetValue(array1.d_min())
    self.settings_panel.d_min_ctrl.SetRange(array1.d_min(), 20.0)
    self.settings_panel.set_index_span(array1.index_span())
    self.settings_panel.update_space_group_choices(array1)
    self.viewer.set_miller_arrays(array1, array2)

  def add_view_specific_functions (self) :
    pass
