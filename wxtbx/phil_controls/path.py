
from wxtbx.phil_controls import text_base
from wxtbx import phil_controls
import wxtbx.icons
import wx
import os

WXTBX_PHIL_PATH_SAVE = 1
WXTBX_PHIL_PATH_DIRECTORY = 2
WXTBX_PHIL_PATH_VIEW_BUTTON = 4
WXTBX_PHIL_PATH_NARROW = 8

class PathCtrl (wx.PyPanel, phil_controls.PhilCtrl) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    self._path_style = kwds.get("style", WXTBX_PHIL_PATH_VIEW_BUTTON)
    kwds['style'] = wx.NO_BORDER
    wx.PyPanel.__init__(self, *args, **kwds)
    self.SetValidator(PathValidator())
    self._formats = ()
    if (self._path_style & WXTBX_PHIL_PATH_NARROW) :
      szr = wx.BoxSizer(wx.VERTICAL)
      path_size = (400,-1)
      szr2_pad = wx.TOP
    else :
      szr = wx.BoxSizer(wx.HORIZONTAL)
      path_size = (300, -1)
      szr2_pad = wx.LEFT
    self.SetSizer(szr)
    self._path_text = wx.TextCtrl(self, -1, size=path_size,
      style=wx.TE_PROCESS_ENTER)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter, self._path_text)
    self.GetValidator().Bind(wx.EVT_TEXT_ENTER, self.OnEnter)
    szr.Add(self._path_text, 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT)
    browse_btn = wx.Button(self, -1, "Browse...")
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    szr.Add(szr2, 0, wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL|szr2_pad, 5)
    szr2.Add(browse_btn, 0, wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnBrowse, browse_btn)
    if (self._path_style & WXTBX_PHIL_PATH_VIEW_BUTTON) :
      bitmap = wxtbx.icons.viewmag.GetBitmap()
      view_btn = wx.BitmapButton(self, -1, bitmap)
      self.Bind(wx.EVT_BUTTON, self.OnDisplayFile, view_btn)
      self.Bind(wx.EVT_CONTEXT_MENU, self.OnDisplayFileMenu, view_btn)
      szr2.Add(view_btn, 0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT, 5)
    szr.Fit(self)

  def SetWidth (self, width) :
    self._path_text.SetSize((width, -1))
    self.Layout()
    self.GetSizer().Fit(self)

  def SetFormats (self, formats) :
    if (isinstance(formats, str)) :
      formats = formats.split(",")
    else :
      assert (isinstance(formats, list) or isisntance(formats, tuple))
    self._formats = formats

  def GetPathStyle (self) :
    return self._path_style

  def GetValue (self) :
    return self._path_text.GetValue()

  def SetValue (self, value) :
    self._path_text.SetValue(value)

  def SetBackgroundColour (self, *args, **kwds) :
    self._path_text.SetBackgroundColour(*args, **kwds)

  def GetPhilValue (self) :
    self._path_text.Validate()
    val_str = self._path_text.GetValue()
    if (val_str == "") :
      return self.ReturnNoneIfOptional()
    return val_str

  def GetStringValue (self) :
    return str(self.GetPhilValue())

  def FormatValue (self, value) :
    if (value != "") :
      return os.path.abspath(value)
    return value

  def Validate (self) :
    self.GetValidator().Validate(self)

  def OnBrowse (self, event) :
    flags = 0
    if (self._path_style & WXTBX_PHIL_PATH_SAVE) :
      flags |= wx.SAVE|wx.OVERWRITE_PROMPT
    new_path = ""
    if (self._path_style & WXTBX_PHIL_PATH_DIRECTORY) :
      new_path = wx.DirSelector(
        message="Choose a directory: %s" % self.GetName(),
        defaultPath=self._path_text.GetValue(),
        style=flags,
        parent=self)
    else :
      from iotbx import file_reader
      wildcard = file_reader.get_wildcard_strings(self._formats)
      current_path = self._path_text.GetValue()
      defaultDir = defaultFile = ""
      if (current_path != "") :
        defaultDir, defaultFile = os.path.split(current_path)
      dlg = wx.FileDialog(
        parent=self,
        message="Choose a file: %s" % self.GetName(),
        defaultDir=defaultDir,
        defaultFile=defaultFile,
        style=flags,
        wildcard="All files (*.*)|*.*")#wildcard)
      if (dlg.ShowModal() == wx.ID_OK) :
        new_path = dlg.GetPath()
      wx.CallAfter(dlg.Destroy)
    if (new_path != "") :
      self.SetValue(new_path)
      self.DoSendEvent()

  def OnDisplayFile (self, event) :
    print "NotImplemented"

  def OnDisplayFileMenu (self, event) :
    print "NotImplemented"

  def OnEnter (self, event) :
    self.Validate()
    self.DoSendEvent()

class PathValidator (text_base.TextCtrlValidator) :
  def CheckFormat (self, value) :
    style = self.GetWindow().GetPathStyle()
    if (value == "") :
      return ""
    elif (os.path.isfile(value)) :
      if (style & WXTBX_PHIL_PATH_DIRECTORY) :
        raise ValueError("file specified, but this parameter requires "+
          "a directory.")
      return os.path.abspath(value)
    elif (os.path.isdir(value)) :
      if (not style & WXTBX_PHIL_PATH_DIRECTORY) :
        raise ValueError("directory specified, but this parameter requires "+
          "a file.")
      return os.path.abspath(value)
    else :
      if (style & WXTBX_PHIL_PATH_SAVE) :
        return os.path.abspath(value)
      else :
        raise ValueError("path does not exist")

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Path test")
  panel = wx.Panel(frame, -1, size=(720,480))
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
    style=WXTBX_PHIL_PATH_DIRECTORY|WXTBX_PHIL_PATH_VIEW_BUTTON)
  path4 = PathCtrl(panel, -1, pos=(20,300), name="Default parameters",
    style=WXTBX_PHIL_PATH_NARROW|WXTBX_PHIL_PATH_VIEW_BUTTON)
  path3.SetFormats("phil")
  pdb_out = os.path.join(os.getcwd(), "model.pdb")
  path2.SetValue(pdb_out)
  button = wx.Button(panel, -1, "Process inputs", pos=(600,400))
  def OnProcess (event) :
    print path1.GetPhilValue()
    print path2.GetPhilValue()
    print path3.GetPhilValue()
    print path4.GetPhilValue()
  def OnPhilEvent (event) :
    print "PhilEvent"
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
