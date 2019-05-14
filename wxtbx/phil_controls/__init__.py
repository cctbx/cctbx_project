from __future__ import absolute_import, division, print_function

from libtbx.utils import Sorry
from libtbx import Auto
import wx

class PhilCtrl(object):
  def __init__(self):
    self.phil_name = None
    self.optional = False
    self._none_is_auto = False
    self._wxtbx_style = 0

  def SetUseAuto(self, enable=True):
    self._none_is_auto = enable

  def UseAuto(self):
    return getattr(self, "_none_is_auto", False)

  def ReturnNoneIfOptional(self):
    if (self.IsOptional()):
      if (self.UseAuto()):
        return Auto
      return None
    else :
      raise Sorry("Value required for '%s'." % self.GetName())

  def SetOptional(self, optional=True):
    if (optional is None) : optional = True
    self.optional = optional

  def IsOptional(self):
    return getattr(self, "optional", True)

  def SetPhilName(self, name):
    self.phil_name = name

  def GetPhilName(self):
    return getattr(self, "phil_name", None)

  def GetStringValue(self):
    raise NotImplementedError()

  def GetPhil(self, full_path=True, indent=0):
    assert (self.phil_name is not None)
    value = self.GetStringValue()
    if (full_path):
      phil_name = self.phil_name
    else :
      phil_name = self.phil_name.split(".")[-1]
    format = "%s%s = %s"
    return format % (" "*indent, phil_name, value)

  def __str__(self):
    return type(self).__name__ + (" (%s)" % self.phil_name)

  def DoSendEvent(self, original_window=None):
    event = PhilCtrlEvent(wxEVT_PHIL_CONTROL, self.GetId())
    event.SetEventObject(self)
    if (original_window is not None):
      event.SetOriginalWindow(original_window)
    self.GetEventHandler().ProcessEvent(event)

  def GetWxtbxStyle(self):
    return self._wxtbx_style

  def SetWxtbxStyle(self, style):
    assert isinstance(style, int)
    self._wxtbx_style = style
    return self

wxEVT_PHIL_CONTROL = wx.NewEventType()
EVT_PHIL_CONTROL = wx.PyEventBinder(wxEVT_PHIL_CONTROL, 1)
class PhilCtrlEvent(wx.PyCommandEvent):
  def __init__(self, eventType, eventId):
    wx.PyCommandEvent.__init__(self, eventType, eventId)
    self._original_window = None

  def SetOriginalWindow(self, window):
    assert isinstance(window, wx.Window)
    self._original_window = window

  def GetOriginalWindow(self):
    return self._original_window
