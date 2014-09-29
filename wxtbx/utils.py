
from __future__ import division
from libtbx import object_oriented_patterns as oop
import wx.lib.colourselect
import wx
import sys

UNICODE_BUILD = (wx.PlatformInfo[2] == 'unicode')

class SettingsToolBase (object) :
  def __init__ (self, *args, **kwds) :
    self.parent = self.GetParent()
    self.settings = self.parent.settings
    self._controls = {}
    assert hasattr(self.parent, "update_settings")

  def add_controls (self) :
    raise NotImplementedError()

  def get_control (self, name) :
    return self._controls.get(name, oop.null())

  def create_controls (self,
                       setting,
                       label,
                       captions=None,
                       min=-sys.maxint,
                       max=sys.maxint) :
    panel = self.panel
    value = getattr(self.settings, setting)
    ctrls = []
    if isinstance(value, bool) :
      box = wx.CheckBox(panel, -1, label)
      self._controls[setting] = box
      box.SetValue(value)
      ctrls.append(box)
      self.Bind(wx.EVT_CHECKBOX,
        lambda evt: self.update_values(evt, setting),
        box)
    elif isinstance(value, int) and (captions is None) :
      txt = wx.StaticText(panel, -1, label)
      ctrls.append(txt)
      spinctrl = wx.SpinCtrl(panel, -1)
      spinctrl.SetValue(value)
      spinctrl.SetRange(min, max)
      self._controls[setting] = spinctrl
      ctrls.append(spinctrl)
      self.Bind(wx.EVT_SPINCTRL,
        lambda evt: self.update_values(evt, setting),
        spinctrl)
    elif isinstance(value, int) and isinstance(captions, list) :
      txt = wx.StaticText(panel, -1, label)
      ctrls.append(txt)
      choice = wx.Choice(panel, -1, choices=captions)
      choice.SetSelection(value)
      self._controls[setting] = choice
      ctrls.append(choice)
      self.Bind(wx.EVT_CHOICE,
        lambda evt: self.update_values(evt, setting),
        choice)
    else :
      assert 0
    return ctrls

  def update_values (self, evt, setting) :
    assert hasattr(self.settings, setting)
    ctrl = evt.GetEventObject()
    if isinstance(ctrl, wx.Choice) :
      value = ctrl.GetSelection()
    elif (type(ctrl).__name__ in ["CheckBox", "SpinCtrl", "Slider"]) :
      value = ctrl.GetValue()
    elif isinstance(ctrl, wx.lib.colourselect.ColourSelect) :
      value = [ x / 255.0 for x in ctrl.GetValue() ]
    assert (type(value) == type(getattr(self.settings, setting)))
    setattr(self.settings, setting, value)
    self.parent.update_settings()

class SettingsFrame (wx.MiniFrame, SettingsToolBase) :
  def __init__ (self, *args, **kwds) :
    wx.MiniFrame.__init__(self, *args, **kwds)
    SettingsToolBase.__init__(self, *args, **kwds)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.panel = wx.Panel(self, -1)
    self.panel_sizer = wx.BoxSizer(wx.VERTICAL)
    self.panel.SetSizer(self.panel_sizer)
    self.sizer.Add(self.panel, 1, wx.EXPAND)
    self.add_controls()
    self.sizer.Fit(self.panel)
    self.Fit()
    self.Bind(wx.EVT_CLOSE, self.OnClose, self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)

  def OnClose (self, event) :
    self.Destroy()

  def OnDestroy (self, event) :
    self.parent.settings_window = None

class SettingsPanel (wx.Panel, SettingsToolBase) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    SettingsToolBase.__init__(self, *args, **kwds)
    self.panel = self
    self.panel_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.panel_sizer)
    self.add_controls()
    self.panel_sizer.Layout()

def bold_text (parent, label) :
  txt = wx.StaticText(parent, label=label)
  font = txt.GetFont()
  font.SetWeight(wx.FONTWEIGHT_BOLD)
  txt.SetFont(font)
  return txt

std_sizer_flags = wx.ALL|wx.ALIGN_CENTER_VERTICAL

def add_ok_cancel_buttons (self, sizer) :
  assert isinstance(self, wx.Dialog)
  ok_btn = wx.Button(self, wx.ID_OK)
  cancel_btn = wx.Button(self, wx.ID_CANCEL)
  btn_szr = wx.StdDialogButtonSizer()
  btn_szr.Add(cancel_btn, 0, wx.ALL, 5)
  btn_szr.Add(ok_btn, 0, wx.ALL, 5)
  ok_btn.SetDefault()
  btn_szr.Realize()
  sizer.Add(btn_szr, 0, wx.ALL|wx.ALIGN_RIGHT, 5)
  return ok_btn, cancel_btn

class LogViewer (wx.TextCtrl) :
  font_size = 12
  if (wx.Platform == '__WXGTK__') :
    font_size = 11
  if (wx.Platform == '__WXMSW__') :
    font_size = 9
  def __init__ (self, *args, **kwds) :
    kwds['style'] = wx.TE_MULTILINE|wx.TE_WORDWRAP|wx.TE_READONLY
    wx.TextCtrl.__init__(self, *args, **kwds)
    #self.SetFont(wx.Font(self.font_size, wx.MODERN, wx.NORMAL, wx.NORMAL))

  def Clear (self) :
    wx.TextCtrl.Clear(self)
    self.SetFont(wx.Font(self.font_size, wx.MODERN, wx.NORMAL, wx.NORMAL))

  def WriteText (self, text) :
    if isinstance(text, str) and UNICODE_BUILD :
      try :
        text = text.decode("utf8")
      except UnicodeDecodeError, e :
        print "warning: %s (string = '%s')" % (str(e), text)
        return
    self.SetFont(wx.Font(self.font_size, wx.MODERN, wx.NORMAL, wx.NORMAL))
    wx.TextCtrl.WriteText(self, text)

  def AppendText (self, text) :
    if isinstance(text, str) and UNICODE_BUILD :
      try :
        text = text.decode("utf8")
      except UnicodeDecodeError, e :
        print "warning: %s (string = '%s')" % (str(e), text)
        return
    wx.TextCtrl.AppendText(self, text)

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, title="Test frame")
  panel = wx.Panel(frame, -1, size=(600,400))
  log = LogViewer(panel, size=(580,380), pos=(10,10))
  frame.Show()
  log.WriteText("This is a log line\n")
  angstrom = u"\u00C5".encode("utf-8", "strict").strip()
  log.WriteText("This is line containing special characters: %s" % angstrom)
  app.MainLoop()
