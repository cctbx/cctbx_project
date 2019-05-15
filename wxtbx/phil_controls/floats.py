from __future__ import absolute_import, division, print_function

from wxtbx.phil_controls.numbers import NumbersCtrlBase, NumbersValidator
import wx
from six.moves import range

class FloatsCtrl(NumbersCtrlBase):
  def CheckType(self, value):
    return isinstance(value, int) or isinstance(value, float)

  def CreateValidator(self):
    return FloatsValidator()

  def SetFloats(self, value):
    if (value is None):
      NumbersCtrlBase.SetValue(self, "")
    elif (isinstance(value, int) or isinstance(value, float)):
      NumbersCtrlBase.SetValue(self, str(float(value)))
    elif (isinstance(value, list) or isinstance(value, tuple)):
      NumbersCtrlBase.SetValue(self, self.FormatValue(value))

  def SetValue(self, value):
    self.SetFloats(value)

  def GetPhilValue(self):
    self.Validate()
    val_str = str(NumbersCtrlBase.GetValue(self))
    if (val_str == ""):
      return self.ReturnNoneIfOptional()
    return [ float(field) for field in val_str.split() ]

  def FormatValue(self, value):
    format = " ".join([ "%g" for x in range(len(value)) ])
    return format % tuple(value)

class FloatsValidator(NumbersValidator):
  def ConvertValue(self, value):
    return float(value)

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Floating-point list test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Search model RMSDs:", pos=(100,180))
  floats_ctrl = FloatsCtrl(panel, -1, pos=(300,180), size=(200,-1),
    value=[1.0,2.0,3.5,4.7],
    name="Search model RMSDs")
  floats_ctrl.SetMin(0.5)
  floats_ctrl.SetOptional(False)
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay(evt):
    floats = floats_ctrl.GetPhilValue()
    print(type(floats).__name__, str(floats))
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
