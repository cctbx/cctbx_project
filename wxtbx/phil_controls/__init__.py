
from libtbx.utils import Abort
import wx

class ValidatedTextCtrl (wx.TextCtrl) :
  def __init__ (self, *args, **kwds) :
    super(ValidatedTextCtrl, self).__init__(*args, **kwds)
    style = self.GetWindowStyle()
    style = self.GetWindowStyle()
    if (not style & wx.TE_PROCESS_ENTER) :
      style |= wx.TE_PROCESS_ENTER
      self.SetWindowStyle(style)
    self.SetValidator(self.CreateValidator())
    self.Bind(wx.EVT_TEXT_ENTER, lambda evt: self.Validate(), self)

  def CreateValidator (self) :
    raise NotImplementedError()

  def Validate (self) :
    # XXX why doesn't self.Validate() work?
    if self.GetValidator().Validate(self.GetParent()) :
      return True
    else :
      raise Abort()

class TextCtrlValidator (wx.PyValidator) :
  def __init__ (self) :
    wx.PyValidator.__init__(self)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter)

  def Clone (self) :
    return self.__class__()

  def TransferToWindow (self) :
    return True

  def TransferFromWindow (self) :
    return True

  def CheckFormat (self, value) :
    raise NotImplementedError()

  def Validate (self, win) :
    ctrl = self.GetWindow()
    value_str = str(ctrl.GetValue())
    try :
      reformatted = self.CheckFormat(value_str)
      ctrl.SetValue(reformatted)
      ctrl.SetBackgroundColour(
        wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
      #ctrl.SetFocus()
      ctrl.Refresh()
      return True
    except NotImplementedError :
      raise
    except Exception, e :
      ctrl_name = ctrl.GetName()
      wx.MessageBox(caption="Format error",
        message="Inappropriate value given for \"%s\": %s" %(ctrl_name,str(e)))
      ctrl.SetBackgroundColour("red")
      ctrl.SetFocus()
      ctrl.Refresh()
      return False

  def OnEnter (self, event) :
    self.Validate(None)
