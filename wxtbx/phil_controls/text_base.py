
from wxtbx import phil_controls
from libtbx.utils import Abort
import wx

class ValidatedTextCtrl (wx.TextCtrl, phil_controls.PhilCtrl) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    saved_value = None
    if (kwds.get('value', "") != "") :
      saved_value = kwds['value']
      kwds['value'] = ""
    super(ValidatedTextCtrl, self).__init__(*args, **kwds)
    style = self.GetWindowStyle()
    style = self.GetWindowStyle()
    if (not style & wx.TE_PROCESS_ENTER) :
      style |= wx.TE_PROCESS_ENTER
      self.SetWindowStyle(style)
    self.SetValidator(self.CreateValidator())
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter, self)
    if (saved_value is not None) :
      self.SetValue(saved_value)

  def GetValue (self) :
    val = wx.TextCtrl.GetValue(self)
    if (isinstance(val, unicode)) :
      return val.encode('utf8')
    else :
      assert isinstance(val, str)
      return val

  def OnEnter (self) :
    self.Validate()
    self.DoSendEvent()

  def CreateValidator (self) :
    raise NotImplementedError()

  def Validate (self) :
    # XXX why doesn't self.Validate() work?
    if self.GetValidator().Validate(self.GetParent()) :
      return True
    else :
      raise Abort()

  def FormatValue (self, value) :
    raise NotImplementedError()

  def GetPhilValue (self) :
    raise NotImplementedError()

  def GetStringValue (self) :
    value = self.GetPhilValue()
    if (value is not None) :
      return self.FormatValue(value)
    return None

  def Enable (self, enable=True) :
    wx.TextCtrl.Enable(self, enable)
    if enable :
      self.SetBackgroundColour((255,255,255))
    else :
      self.SetBackgroundColour((200,200,200))

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
    if (value_str == "") :
      ctrl.SetBackgroundColour(
        wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
      return True
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
      ctrl_name = str(ctrl.GetName())
      wx.MessageBox(caption="Format error",
        message="Inappropriate value given for \"%s\": %s" %(ctrl_name,str(e)))
      ctrl.SetBackgroundColour("red")
      ctrl.SetFocus()
      ctrl.Refresh()
      return False

  def OnEnter (self, event) :
    self.Validate(None)
    ctrl = self.GetWindow()
    ctrl.DoSendEvent()
