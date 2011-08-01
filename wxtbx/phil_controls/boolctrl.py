
from wxtbx import phil_controls
import wx
from libtbx import Auto

WXTBX_PHIL_BOOL_TRIBOOL = 1
WXTBX_PHIL_BOOL_AUTO = 2
#WXTBX_PHIL_BOOL

class BoolCtrl (wx.CheckBox, phil_controls.PhilCtrl) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    self._bool_style = kwds.get("style", 0)
    kwds['style'] = 0
    if ((self._bool_style & WXTBX_PHIL_BOOL_TRIBOOL) or
        (self._bool_style & WXTBX_PHIL_BOOL_AUTO)) :
      kwds['style'] |= wx.CHK_ALLOW_3RD_STATE_FOR_USER|wx.CHK_3STATE
    else :
      kwds['style'] |= wx.CHK_3STATE # wx.CHK_ALLOW_3RD_STATE_FOR_USER?
    wx.CheckBox.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_CHECKBOX, lambda evt: self.DoSendEvent())

  def SetValue (self, value) :
    if (value is None) or (value is Auto) :
      assert (self.Is3State())
      self.Set3StateValue(wx.CHK_UNDETERMINED)
    else :
      if (self.Is3State()) :
        if (value == True) :
          self.Set3StateValue(wx.CHK_CHECKED)
        else :
          self.Set3StateValue(wx.CHK_UNCHECKED)
      else :
        wx.CheckBox.SetValue(self, value)

  def GetValue (self) :
    if (self.Is3State()) :
      value = self.Get3StateValue()
      if (value == wx.CHK_UNDETERMINED) :
        if (self._bool_style & WXTBX_PHIL_BOOL_AUTO) :
          return Auto
        else :
          return None
      else :
        return (value == wx.CHK_CHECKED)
    else :
      return wx.CheckBox.GetValue(self)

  def GetPhilValue (self) :
    return self.GetValue()

  def GetStringValue (self) :
    return str(self.GetValue())

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "PHIL bool test")
  panel = wx.Panel(frame, -1, size=(600,400))
  box1 = BoolCtrl(panel, label="Use NCS restraints", pos=(100,100))
  box2 = BoolCtrl(panel, label="Find NCS groups automatically", pos=(100,150))
  box3 = BoolCtrl(panel, label="Fast search mode", pos=(100,200),
    style=WXTBX_PHIL_BOOL_AUTO)
  box1.SetValue(False)
  box2.SetValue(None)
  box3.SetValue(Auto)
  assert (box1.GetValue() == box1.GetPhilValue() == False)
  assert (box2.GetValue() is None)
  assert (box3.GetValue() is Auto)
  assert (box2.GetStringValue() == "None")
  assert (box3.GetStringValue() == "Auto")
  box3.SetValue(False)
  assert (box3.GetStringValue() == "False")
  box1.SetValue(True)
  assert (box1.GetStringValue() == "True")
  def OnChange (event) :
    print event.GetEventObject().GetPhilValue()
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnChange)
  frame.Show()
  app.MainLoop()
