
import wx
import os.path
from libtbx.utils import Abort

class manager (object) :
  """
  Because wxPython is inconsistent (if not outright buggy on some platforms)
  in its selection of the starting directory for the file/folder dialogs,
  this manager remembers where it was last opened.  This is not the same as
  using the wx.FD_CHANGE_DIR flag, which actually changes the working
  directory.  Instead, a manager can be associated with a specific window (or
  group of windows), so that any changes in path are localized.
  """
  def __init__ (self, start_dir=None) :
    if (start_dir is None) :
      start_dir = os.getcwd()
    assert os.path.isdir(start_dir)
    self.start_dir = start_dir
    self.last_dir = start_dir

  def reset_dir (self) :
    self.last_dir = self.start_dir

  def set_current_dir (self, dir_name) :
    assert os.path.isdir(dir_name)
    self.last_dir = dir_name

  def select_file (self,
                   parent,
                   message,
                   style=wx.OPEN,
                   wildcard="All files (*.*)|*.*",
                   current_file=None,
                   multiple=False) :
    default_dir = self.last_dir
    default_file = ""
    if (current_file is not None) and (current_file != "") :
      if (os.path.isabs(current_file)) :
        default_dir, default_file = os.path.split(current_file)
      else :
        default_file = current_file
    dlg = wx.FileDialog(
      parent=parent,
      message=message,
      defaultDir=default_dir,
      defaultFile=default_file,
      style=style,
      wildcard=wildcard)
    if (dlg.ShowModal() == wx.ID_OK) :
      new_path = None
      if (multiple) :
        new_path = dlg.GetPaths()
        if (new_path is not None) and (len(new_path) > 0) :
          new_dir = os.path.dirname(new_path[0])
          assert os.path.isdir(new_dir)
          self.last_dir = new_dir
      else :
        new_path = dlg.GetPath()
        if (new_path != "") :
          self.last_dir = os.path.dirname(new_path)
        else :
          new_path = None
      wx.CallAfter(dlg.Destroy)
      return new_path
    else :
      wx.CallAfter(dlg.Destroy)
      raise Abort()

  def select_directory (self,
                        parent,
                        message,
                        current_path,
                        style) :
    default_path = self.last_dir
    if (current_path is not None) and (current_path != "") :
      default_path = current_path
    new_path = wx.DirSelector(
      message=message,
      defaultPath=default_path,
      style=style|wx.DD_NEW_DIR_BUTTON,
      parent=parent)
    if (new_path != "") :
      return new_path
    return None
