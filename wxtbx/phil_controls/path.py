from __future__ import absolute_import, division, print_function
from wxtbx.phil_controls import text_base
from wxtbx import phil_controls
import wxtbx.icons
from libtbx.utils import Sorry, to_unicode, to_str
from libtbx import Auto
import wx
import os
import sys

WXTBX_PHIL_PATH_SAVE = 1
WXTBX_PHIL_PATH_DIRECTORY = 2
WXTBX_PHIL_PATH_VIEW_BUTTON = 4
WXTBX_PHIL_PATH_NARROW = 8
WXTBX_PHIL_PATH_UPDATE_ON_KILL_FOCUS = 16
WXTBX_PHIL_PATH_DEFAULT_CWD = 32

class PathCtrl(wx.Panel, phil_controls.PhilCtrl):
  def __init__(self, *args, **kwds):
    phil_controls.PhilCtrl.__init__(self)
    self.SetOptional(True) # this will be overridden elsewhere if necessary
    wx.SystemOptions.SetOption("osx.openfiledialog.always-show-types", "1")
    kwds = dict(kwds)
    self._path_style = kwds.get("style", WXTBX_PHIL_PATH_VIEW_BUTTON)
    assert ((self._path_style & WXTBX_PHIL_PATH_DIRECTORY) or
            (not self._path_style & WXTBX_PHIL_PATH_DEFAULT_CWD))
    value = kwds.pop("value", None)
    kwds['style'] = wx.NO_BORDER
    self._formats = ()
    if (kwds.get("size", None) == wx.DefaultSize):
      kwds.pop("size")
    if (self._path_style & WXTBX_PHIL_PATH_NARROW):
      szr_type = wx.VERTICAL
      path_size = kwds.get("size", (300,-1))
      szr2_pad = wx.TOP
    else :
      szr_type = wx.HORIZONTAL
      path_size = kwds.get("size", (400, -1))
      szr2_pad = wx.LEFT
    kwds['size'] = wx.DefaultSize
    wx.Panel.__init__(self, *args, **kwds)
    self._formats = ()
    szr = wx.BoxSizer(szr_type)
    self.SetSizer(szr)
    self._path_text = wx.TextCtrl(self, -1, size=path_size,
      style=wx.TE_PROCESS_ENTER,
      name=self.GetName())
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter, self._path_text)
    if (self._path_style & WXTBX_PHIL_PATH_UPDATE_ON_KILL_FOCUS):
      self._path_text.Bind(wx.EVT_KILL_FOCUS, self.OnEnter)
    self._path_text.SetValidator(PathValidator())
    self._path_text.GetValidator().Bind(wx.EVT_TEXT_ENTER, self.OnEnter)
    szr.Add(self._path_text, 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT)
    browse_btn = wx.Button(self, -1, "Browse...")
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    szr.Add(szr2, 0, wx.ALIGN_CENTER_VERTICAL|szr2_pad, 5)
    szr2.Add(browse_btn, 0, wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnBrowse, browse_btn)
    self.browse_btn = browse_btn
    if (self._path_style & WXTBX_PHIL_PATH_VIEW_BUTTON):
      bitmap = wxtbx.icons.viewmag.GetBitmap()
      view_btn = wx.BitmapButton(self, -1, bitmap)
      self.Bind(wx.EVT_BUTTON, self.OnDisplayFile, view_btn)
      view_btn.Bind(wx.EVT_CONTEXT_MENU, self.OnDisplayFileMenu)
      #self.Bind(wx.EVT_RIGHT_DOWN, self.OnDisplayFileMenu, view_btn)
      szr2.Add(view_btn, 0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT, 5)
    szr.Layout()
    szr.Fit(self)
    drop_target = PathDropTarget(self)
    self.SetDropTarget(drop_target)
    self.SetValue(value)
    self._pathmgr = None
    self._cached_paths = None

  def SetFormats(self, formats):
    if (isinstance(formats, str)):
      formats = formats.split(",")
    else :
      assert (isinstance(formats, list) or isisntance(formats, tuple))
    self._formats = formats

  def Clear(self):
    self._path_text.Clear()

  def GetPathStyle(self):
    return self._path_style

  def GetValue(self):
    val = self._path_text.GetValue().strip()
    # use unicode check to avoid bytes in Python 3
    check_type = bytes
    if sys.version_info.major == 2:
      check_type = unicode
    if (isinstance(val, check_type)) and wxtbx.is_unicode_build():
      return to_str(val)
    else :
      assert isinstance(val, str)
      return val

  def SetValue(self, value):
    if (value is None) or (value is Auto):
      self._path_text.SetValue("")
    else :
      value = to_unicode(value)
      self._path_text.SetValue(value)

  def SetBackgroundColour(self, *args, **kwds):
    self._path_text.SetBackgroundColour(*args, **kwds)

  def GetPhilValue(self):
    ok = self._path_text.GetValidator().Validate(self)
    if ok:
      val_str = self.GetValue()
    else:
      val_str = ""
    assert isinstance(val_str, str)
    if (val_str == ""):
      return self.ReturnNoneIfOptional()
    return val_str

  def GetStringValue(self):
    return str(self.GetPhilValue())

  def FormatValue(self, value):
    if (value != ""):
      if (not os.path.isabs(value)):
        print("ABSPATH")
        return os.path.abspath(value)
      else :
        return value
    return value

  def Validate(self):
    return self._path_text.GetValidator().Validate(self)

  def OnBrowse(self, event):
    flags = 0
    if (self._path_style & WXTBX_PHIL_PATH_SAVE):
      flags |= wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT
    else :
      flags |= wx.FD_OPEN
    path_manager = self.GetPathManager()
    new_path = None
    if (self._path_style & WXTBX_PHIL_PATH_DIRECTORY):
      new_path = path_manager.select_directory(
        message="Choose a directory: %s" % self.GetName(),
        current_path=to_unicode(self.GetValue()),
        style=flags|wx.DD_NEW_DIR_BUTTON,
        parent=self)
    else :
      from iotbx import file_reader
      wildcard = file_reader.get_wildcard_strings(self._formats)
      new_path = path_manager.select_file(
        parent=self,
        message="Choose a file: %s" % self.GetName(),
        current_file=to_unicode(self.GetValue()),
        style=flags,
        wildcard=wildcard)
    if (new_path is not None):
      if ('}' in new_path) or ('{' in new_path):
        raise Sorry("Curly brackets ({}) are not allowed in pathnames.")
      self.SetValue(new_path)
      self.DoSendEvent()

  def GetPathManager(self):
    if (self._pathmgr is None):
      main_window = self.GetTopLevelParent()
      if hasattr(main_window, "get_path_manager"):
        self._pathmgr = main_window.get_path_manager()
    if (self._pathmgr is None):
      from wxtbx import path_dialogs
      self._pathmgr = path_dialogs.manager()
    return self._pathmgr

  def OnDisplayFile(self, event):
    file_name = self.GetValue()
    if (file_name == ""):
      return
    elif (not os.path.exists(file_name)):
      raise Sorry("Path does not exist.")
    tlp = self.GetTopLevelParent()
    if (hasattr(tlp, "display_file")):
      tlp.display_file(file_name)
    else :
      print("NotImplemented")

  def OnDisplayFileMenu(self, event):
    file_name = self.GetValue()
    if (file_name == ""):
      return
    elif (not os.path.exists(file_name)):
      raise Sorry("Path does not exist.")
    tlp = self.GetTopLevelParent()
    if (hasattr(tlp, "display_file_menu")):
      tlp.display_file_menu(file_name)
    else :
      print("NotImplemented")

  def OnEnter(self, event):
    if (self._path_style & WXTBX_PHIL_PATH_DEFAULT_CWD):
      value = self.GetValue()
      if (value == ""):
        from libtbx.utils import getcwd_safe
        self.SetValue(getcwd_safe())
    self.Validate()
    self.DoSendEvent()

  def SetCachedPaths(self, paths):
    self._cached_paths = paths
    # XXX on Mac, right-clicking on the TextCtrl is already handled by the OS,
    # so only the browse button will work for this; need to check other OSes
    self.browse_btn.Bind(wx.EVT_RIGHT_DOWN, self.OnRightClick, self.browse_btn)

  def OnRightClick(self, event):
    assert (self._cached_paths is not None)
    menu = wx.Menu()
    current_path_name = self.GetValue()
    n_items = 0
    def set_path(pn):
      self.SetValue(pn)
    if (current_path_name not in [None, ""]):
      item = menu.Append(-1, current_path_name)
      self.Bind(wx.EVT_MENU, lambda evt, pn=current_path_name: set_path(pn),
        item)
      n_items +=1
    for path_name in self._cached_paths :
      if (path_name != current_path_name):
        item = menu.Append(-1, path_name)
        self.Bind(wx.EVT_MENU, lambda evt, pn=path_name: set_path(pn),
          item)
        n_items +=1
    if (n_items > 0):
      self.PopupMenu(menu)
    menu.Destroy()

class PathValidator(text_base.TextCtrlValidator):
  def CheckFormat(self, value):
    style = self.GetWindow().GetParent().GetPathStyle()
    if ('}' in value) or ('{' in value):
      raise ValueError("Curly brackets ({}) are not allowed in pathnames.")
    elif (value == ""):
      return ""
    elif (os.path.isfile(value)):
      if (style & WXTBX_PHIL_PATH_DIRECTORY):
        raise ValueError("file specified, but this parameter requires "+
          "a directory.")
      if (not os.path.isabs(value)):
        print("ABSPATH")
        return os.path.abspath(value)
      else :
        return value
    elif (os.path.isdir(value)):
      if (not style & WXTBX_PHIL_PATH_DIRECTORY):
        # XXX hack to allow app and package bundles on OS X
        if ((sys.platform == "darwin") and (value.endswith(".app") or
            (value.endswith(".pkg")))):
          pass
        else :
          raise ValueError("directory specified, but this parameter requires "+
            "a file.")
      if (not os.path.isabs(value)):
        print("ABSPATH")
        return os.path.abspath(value)
      else :
        return value
    else :
      if (style & WXTBX_PHIL_PATH_SAVE):
        return os.path.abspath(value)
      else :
        raise ValueError("path does not exist")

class PathDropTarget(wx.FileDropTarget):
  def __init__(self, window):
    wx.FileDropTarget.__init__(self)
    self.window = window

  def OnDropFiles(self, x, y, filenames):
    self.window.SetValue(filenames[-1])
    self.window.Validate()
    self.window.DoSendEvent()
    return True

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Path test")
  panel = wx.Panel(frame, -1, size=(800,480))
  txt1 = wx.StaticText(panel, -1, "Reflections file", pos=(20,100))
  txt2 = wx.StaticText(panel, -1, "Output PDB file", pos=(20,160))
  txt3 = wx.StaticText(panel, -1, "Output directory", pos=(20,220))
  txt4 = wx.StaticText(panel, -1, "Default parameters", pos=(20,280))
  path1 = PathCtrl(panel, -1, pos=(200,100), name="Reflections file")
  path1.SetOptional(False)
  path1.SetFormats("hkl")
  path2 = PathCtrl(panel, -1, pos=(200,160), name="Output PDB file",
    style=WXTBX_PHIL_PATH_SAVE|WXTBX_PHIL_PATH_VIEW_BUTTON)
  path2.SetFormats("pdb")
  path3 = PathCtrl(panel, -1, pos=(200,220), name="Output directory",
    size=(300,-1),
    style=WXTBX_PHIL_PATH_DIRECTORY|WXTBX_PHIL_PATH_VIEW_BUTTON|
      WXTBX_PHIL_PATH_UPDATE_ON_KILL_FOCUS|WXTBX_PHIL_PATH_DEFAULT_CWD)
  path3.SetOptional(False)
  path3.SetCachedPaths(["/var/tmp", "/tmp"])
  path4 = PathCtrl(panel, -1, pos=(20,300), name="Default parameters",
    style=WXTBX_PHIL_PATH_NARROW|WXTBX_PHIL_PATH_VIEW_BUTTON)
  path4.SetFormats("phil")
  pdb_out = os.path.join(os.getcwd(), "model.pdb")
  path2.SetValue(pdb_out)
  button = wx.Button(panel, -1, "Process inputs", pos=(600,400))
  def OnProcess(event):
    print(path1.GetPhilValue())
    print(path2.GetPhilValue())
    print(path3.GetPhilValue())
    print(path4.GetPhilValue())
  def OnPhilEvent(event):
    print("PhilEvent")
  frame.Bind(wx.EVT_BUTTON, OnProcess, button)
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnPhilEvent)
  frame.Fit()
  frame.Show()
  assert (path2.GetPhilValue() == pdb_out)
  assert (path2.GetStringValue() == pdb_out)
  try :
    path1.GetPhilValue()
  except Exception :
    pass
  else :
    raise RuntimeError("Exception expected!")
  app.MainLoop()
