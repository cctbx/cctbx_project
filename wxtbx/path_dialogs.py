
from __future__ import absolute_import, division, print_function
from wxtbx.utils import add_ok_cancel_buttons, std_sizer_flags
from wxtbx.phil_controls import choice_multi, path
import wx
from libtbx.utils import Abort, Sorry, to_unicode
import libtbx.path
import os.path
from six.moves import zip

class manager(object):
  """
  Because wxPython is inconsistent (if not outright buggy on some platforms)
  in its selection of the starting directory for the file/folder dialogs,
  this manager remembers where it was last opened.  This is not the same as
  using the wx.FD_CHANGE_DIR flag, which actually changes the working
  directory.  Instead, a manager can be associated with a specific window (or
  group of windows), so that any changes in path are localized.
  """
  def __init__(self, start_dir=None):
    if (start_dir is None):
      start_dir = os.getcwd()
    assert os.path.isdir(start_dir), start_dir
    self.start_dir = start_dir
    self.last_dir = start_dir

  def reset_dir(self):
    self.last_dir = self.start_dir

  def set_current_dir(self, dir_name):
    assert os.path.isdir(dir_name)
    self.last_dir = dir_name

  def select_file(self,
                   parent,
                   message,
                   style=wx.FD_OPEN,
                   wildcard="All files (*.*)|*.*",
                   current_file=None,
                   multiple=False,
                   save=None):
    if (save):
      style = wx.FD_SAVE
    default_dir = self.last_dir
    default_file = ""
    if (current_file is not None) and (current_file != ""):
      if (os.path.isabs(current_file)):
        default_dir, default_file = os.path.split(current_file)
      else :
        default_file = current_file
    default_dir = to_unicode(default_dir)
    default_file = to_unicode(default_file)
    dlg = wx.FileDialog(
      parent=parent,
      message=message,
      defaultDir=default_dir,
      defaultFile=default_file,
      style=style,
      wildcard=wildcard)
    if (dlg.ShowModal() == wx.ID_OK):
      new_path = None
      if (multiple):
        new_path = dlg.GetPaths()
        if (new_path is not None) and (len(new_path) > 0):
          new_dir = os.path.dirname(new_path[0])
          if (os.path.isdir(new_dir)):
            self.last_dir = new_dir
      else :
        new_path = dlg.GetPath()
        if (new_path != ""):
          self.last_dir = os.path.dirname(new_path)
        else :
          new_path = None
      wx.CallAfter(dlg.Destroy)
      return new_path
    else :
      wx.CallAfter(dlg.Destroy)
      raise Abort()

  def select_directory(self,
                        parent,
                        message,
                        current_path,
                        style):
    default_path = self.last_dir
    if (current_path is not None) and (current_path != ""):
      default_path = current_path
    default_path = to_unicode(default_path)
    if wx.version().startswith('3'):
      new_path = wx.DirSelector(
        message=message,
        defaultPath=default_path,
        style=style|wx.DD_NEW_DIR_BUTTON,
        parent=parent)
    else:
      new_path = wx.DirSelector(
        message=message,
        default_path=default_path,
        style=style|wx.DD_NEW_DIR_BUTTON,
        parent=parent)
    if (new_path != ""):
      return new_path
    return None

class DirectoryCleanupDialog(wx.Dialog):
  def __init__(self, *args, **kwds):
    wx.Dialog.__init__(self, *args, **kwds)
    outer_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(outer_sizer)
    szr = wx.BoxSizer(wx.VERTICAL)
    outer_sizer.Add(szr, 0, wx.ALL, 5)
    txt = wx.StaticText(self, label="Please select the file and/or directory "+
      "type(s) you want to delete.  You will be prompted to confirm this "+
      "action once a list of targeted paths has been collected.")
    txt.Wrap(500)
    szr.Add(txt, 0, wx.ALL, 5)
    grid = wx.FlexGridSizer(cols=2, rows=2)
    szr.Add(grid, 0, wx.ALL, 5)
    grid.Add(wx.StaticText(self, label="Directory:"), 0, std_sizer_flags, 5)
    self.dir_ctrl = path.PathCtrl(parent=self,
      name="Directory",
      style=path.WXTBX_PHIL_PATH_DIRECTORY|path.WXTBX_PHIL_PATH_DEFAULT_CWD)
    self.dir_ctrl.SetOptional(False)
    grid.Add(self.dir_ctrl, 0, std_sizer_flags, 5)
    grid.Add(wx.StaticText(self, label="Items to remove:"), 0,
      std_sizer_flags, 5)
    self.remove_ctrl = choice_multi.MultiChoiceCtrl(
      parent=self,
      name="Items to remove")
    self.remove_ctrl.SetCols(2)
    items = [ "kin", "geo", "map", "probe", "temp", ]
    labels = [ "Kinemage files (.kin)", "Geometry files (.geo)",
      "Map files (.ccp4, .xplor)", "Probe files (probe.txt)",
      "Temporary folders (TEMP0)", ]
    for choice, label in zip(items, labels):
      self.remove_ctrl.AddChoice("*" + choice, label)
    self.remove_ctrl.Realize()
    grid.Add(self.remove_ctrl, 0, std_sizer_flags, 5)
    add_ok_cancel_buttons(self, szr)
    self.Fit()
    self.Centre(wx.BOTH)

  def SetDirectory(self, dir_name):
    assert os.path.isdir(dir_name)
    self.dir_ctrl.SetValue(dir_name)

  def GetDirectory(self):
    return self.dir_ctrl.GetPhilValue()

  def GetChoices(self):
    return self.remove_ctrl.GetPhilValue()

class ConfirmCleanupDialog(wx.Dialog):
  def __init__(self, parent, cleanup_obj):
    wx.Dialog.__init__(self, parent=parent, title="Confirm delete")
    outer_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(outer_sizer)
    szr = wx.BoxSizer(wx.VERTICAL)
    outer_sizer.Add(szr, 0, wx.ALL, 5)
    txt = wx.StaticText(self, label=("%d file(s) and %d directories have "+
      "been found matching the specified criteria; this will save %s of "+
      "disk space.  Their paths are listed below.  Click OK to "+
      "confirm, or Cancel to abort.  You will not be able to recover these "+
      "paths once they have been deleted!") % (cleanup_obj.n_files,
        cleanup_obj.n_dirs, cleanup_obj.get_freed_space()))
    txt.Wrap(480)
    szr.Add(txt, 0, std_sizer_flags, 5)
    list_font = wx.Font(10, wx.FONTFAMILY_MODERN, wx.NORMAL, wx.NORMAL)
    if (cleanup_obj.n_dirs > 0):
      szr.Add(wx.StaticText(self, label="Directories to remove:"), 0,
        std_sizer_flags, 5)
      dir_box = wx.TextCtrl(self,
        style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_DONTWRAP,
        size=(480, 200))
      dir_box.SetFont(list_font)
      szr.Add(dir_box, 0, std_sizer_flags, 5)
      for dir_name in cleanup_obj.dir_paths :
        dir_box.AppendText(dir_name + "\n")
    if (cleanup_obj.n_files > 0):
      szr.Add(wx.StaticText(self, label="Files to remove:"), 0,
        std_sizer_flags, 5)
      file_box = wx.TextCtrl(self,
        style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_DONTWRAP,
        size=(480, 200))
      file_box.SetFont(list_font)
      szr.Add(file_box, 0, std_sizer_flags, 5)
      for file_name in cleanup_obj.file_paths :
        file_box.AppendText(file_name + "\n")
    add_ok_cancel_buttons(self, szr)
    self.Fit()
    self.Centre(wx.BOTH)

def clean_out_directory(parent=None, dir_name=None):
  if (dir_name is None):
    dir_name = os.getcwd()
  dlg1 = DirectoryCleanupDialog(
    parent=parent,
    title="Clean up directory",
    style=wx.CAPTION|wx.CLOSE_BOX|wx.RAISED_BORDER|\
      wx.WS_EX_VALIDATE_RECURSIVELY)
  dlg1.SetDirectory(dir_name)
  dlg2 = None
  try :
    if (dlg1.ShowModal() == wx.ID_OK):
      dir_name = dlg1.GetDirectory()
      choices = dlg1.GetChoices()
      cleanup_obj = libtbx.path.clean_out_directory(
        path_name=dir_name,
        delete_kin_files="kin" in choices,
        delete_geo_files="geo" in choices,
        delete_map_files="map" in choices,
        delete_probe_files="probe" in choices,
        delete_temp_dirs="temp" in choices)
      if (cleanup_obj.n_total == 0):
        raise Sorry("No files or directories matching the specified criteria.")
      else :
        dlg2 = ConfirmCleanupDialog(
          parent=parent,
          cleanup_obj=cleanup_obj)
        if (dlg2.ShowModal() == wx.ID_OK):
          cleanup_obj.run()
          wx.MessageBox(message="Action complete; %s of disk space freed." %
            cleanup_obj.get_freed_space())
          return cleanup_obj
        else :
          raise Abort()
    else :
      raise Abort()
  finally :
    dlg1.Destroy()
    if (dlg2 is not None) : dlg2.Destroy()

if (__name__ == "__main__"):
  app = wx.App(0)
  mgr = manager()
  fn = mgr.select_file(parent=None,
    message="Select a file")
  print(fn)
