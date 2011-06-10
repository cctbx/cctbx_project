
from libtbx import object_oriented_patterns as oop
import wx.lib.colourselect
import wx
import sys

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
