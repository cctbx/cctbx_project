
from __future__ import absolute_import, division, print_function
from wxtbx import phil_controls
import wxtbx
from libtbx.utils import Abort, to_unicode, to_str
from libtbx import Auto
import wx
import sys

class ValidatedTextCtrl(wx.TextCtrl, phil_controls.PhilCtrl):
  def __init__(self, *args, **kwds):
    saved_value = None
    if (kwds.get('value', "") != ""):
      saved_value = kwds['value']
      kwds['value'] = ""
    super(ValidatedTextCtrl, self).__init__(*args, **kwds)
    font = wx.Font(wxtbx.default_font_size, wx.MODERN, wx.NORMAL, wx.NORMAL)
    self.SetFont(font)
    style = self.GetWindowStyle()
    if (not style & wx.TE_PROCESS_ENTER):
      style |= wx.TE_PROCESS_ENTER
      self.SetWindowStyle(style)
    self.SetValidator(self.CreateValidator())
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter, self)
    self.Bind(wx.EVT_KILL_FOCUS, self.OnFocusLost, self)
    if saved_value is not None:
      if (type(saved_value) == str):
        save_value = to_unicode(saved_value)
      self.SetValue(saved_value)

  def GetValue(self):
    try:
      val = wx.TextCtrl.GetValue(self)
    except Exception as e:
      val = "" # Probably C++ object deleted
    if wxtbx.is_unicode_build():
      return to_str(val)
    else :
      assert isinstance(val, str)
      return val

  def OnEnter(self, evt=None):
    #self.Validate()
    self.DoSendEvent()

  def OnFocusLost(self, event):
    self.DoSendEvent()
    event.Skip()

  def CreateValidator(self):
    raise NotImplementedError()

  def Validate(self):
    # XXX why doesn't self.Validate() work?
    if self.GetValidator().Validate(self.GetParent()):
      return True
    else :
      raise Abort()

  def FormatValue(self, value):
    raise NotImplementedError()

  def GetPhilValue(self):
    raise NotImplementedError()

  def GetStringValue(self):
    value = self.GetPhilValue()
    if (value is not None) and (value is not Auto):
      return self.FormatValue(value)
    elif (self.UseAuto()) or (value is Auto):
      return Auto
    return None

  def Enable(self, enable=True):
    wx.TextCtrl.Enable(self, enable)
    if enable :
      self.SetBackgroundColour((255,255,255))
    else :
      self.SetBackgroundColour((200,200,200))

Validator = wx.Validator
if wx.VERSION < (4, 0):
  Validator = wx.PyValidator

class TextCtrlValidator(Validator):
  def __init__(self):
    Validator.__init__(self)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter)

  def Clone(self):
    return self.__class__()

  def TransferToWindow(self):
    return True

  def TransferFromWindow(self):
    return True

  def CheckFormat(self, value):
    raise NotImplementedError()

  def Validate(self, win):
    ctrl = self.GetWindow()
    try :
      value = to_unicode(ctrl.GetValue())
      # if isinstance(value, str):
      #   value = value.decode("utf-8")
      if (value == ""):
        ctrl.SetBackgroundColour(
          wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
        return True
      reformatted = self.CheckFormat(value)
      if isinstance(reformatted, str):
        reformatted = to_unicode(reformatted)
      ctrl.SetValue(reformatted)
      ctrl.SetBackgroundColour(
        wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
      #ctrl.SetFocus()
      ctrl.Refresh()
      return True
    except NotImplementedError :
      raise
    except Exception as e :
      ctrl_name = str(ctrl.GetName())
      msg = "Inappropriate value given for \"%s\": %s" %(ctrl_name,str(e))
      if (type(e).__name__ == "UnicodeEncodeError"):
        msg = ("You have entered characters which cannot be converted to "+
          "Latin characters in the control '%s'; due to limitations of the "+
          "underlying code, only the standard UTF-8 character set is "+
          "allowed.") % ctrl_name
      wx.MessageBox(caption="Format error", message=msg)
      ctrl.SetBackgroundColour("red")
      # Don't set focus on Windows since messagebox is modal and thus
      # would automatically recapture focus leading to an endless UI loop
      if (sys.platform != 'win32'):
        ctrl.SetFocus()
      ctrl.Refresh()
      return False

  def OnEnter(self, event):
    #self.Validate(None)
    ctrl = self.GetWindow()
    ctrl.DoSendEvent()
