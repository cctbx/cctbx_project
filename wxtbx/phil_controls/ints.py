
from wxtbx.phil_controls.numbers import NumbersCtrlBase, NumbersValidator
import wx

class IntsCtrl (NumbersCtrlBase) :
  def CheckType (self, value) :
    return isinstance(value, int)

  def CreateValidator (self) :
    return IntsValidator()

  def SetInts (self, value) :
    if (value is None) :
      NumbersCtrlBase.SetValue(self, "")
    elif (isinstance(value, int) or isinstance(value, float)) :
      NumbersCtrlBase.SetValue(self, str(int(value)))
    elif (isinstance(value, list) or isinstance(value, tuple)) :
      NumbersCtrlBase.SetValue(self, self.FormatValue(value))
    else :
      raise TypeError("Type '%s' not allowed!" % type(value).__name__)

  def SetValue (self, value) :
    self.SetInts(value)

  def GetPhilValue (self) :
    self.Validate()
    val_str = str(NumbersCtrlBase.GetValue(self))
    if (val_str == "") :
      return self.ReturnNoneIfOptional()
    return [ int(field) for field in val_str.split() ]

  def FormatValue (self, value) :
    format = " ".join([ "%d" for x in range(len(value)) ])
    return format % tuple(value)

class IntsValidator (NumbersValidator) :
  def ConvertValue (self, value) :
    return int(value)

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Integer list test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Run SA on cycle numbers:", pos=(100,180))
  ints_ctrl = IntsCtrl(panel, -1, pos=(300,180), size=(200,-1),
    value=[1,2,5,10],
    name="Run SA on cycle numbers")
  ints_ctrl.SetMin(1)
  txt2 = wx.StaticText(panel, -1, "Window color:", pos=(100,240))
  ints_ctrl2 = IntsCtrl(panel, -1, pos=(300,240), size=(200,-1),
    name="Window color")
  ints_ctrl2.SetSize(3)
  ints_ctrl2.SetMin(0)
  ints_ctrl2.SetMax(255)
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay (evt) :
    ints = ints_ctrl.GetPhilValue()
    ints2 = ints_ctrl2.GetPhilValue()
    print type(ints).__name__, str(ints)
    print type(ints2).__name__, str(ints2)
  def OnChange (evt) :
    print "OnChange"
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  import wxtbx.phil_controls
  frame.Bind(wxtbx.phil_controls.EVT_PHIL_CONTROL, OnChange)
  frame.Fit()
  frame.Show()
  app.MainLoop()
