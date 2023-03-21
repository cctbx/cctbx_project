from __future__ import absolute_import, division, print_function

from wxtbx.phil_controls.text_base import ValidatedTextCtrl, TextCtrlValidator
from libtbx.utils import Sorry
from libtbx import Auto
import wx
import sys

class FloatCtrl(ValidatedTextCtrl):
  def __init__(self, *args, **kwds):
    super(FloatCtrl, self).__init__(*args, **kwds)
    self.min = float(-sys.maxsize)
    self.max = float(sys.maxsize)

  def SetMin(self, min):
    assert isinstance(min, float) or isinstance(min, int)
    self.min = float(min)

  def SetMax(self, max):
    assert isinstance(max, float) or isinstance(max, int)
    self.max = max

  def GetMin(self):
    return self.min

  def GetMax(self):
    return self.max

  def CreateValidator(self):
    return FloatValidator()

  def SetFloat(self, value):
    if (value is None) or (value is Auto):
      ValidatedTextCtrl.SetValue(self, "")
    elif (isinstance(value, float) or isinstance(value, int)):
      ValidatedTextCtrl.SetValue(self, str(value))
    elif (isinstance(value, str)):
      try : # TODO remove this after testing
        value = float(value)
      except ValueError :
        raise Sorry("Inappropriate value '%s' for %s." % (value,
          self.GetName()))
      else :
        ValidatedTextCtrl.SetValue(self, str(value))
    else :
      raise TypeError("Type '%s' not allowed!" % type(value).__name__)

  def SetValue(self, value):
    self.SetFloat(value)

  def GetPhilValue(self):
    try:
      self.Validate()
    except Exception as e:
      return self.ReturnNoneIfOptional() # probably C++ object deleted

    val_str = ValidatedTextCtrl.GetValue(self)
    if (val_str == ""):
      return self.ReturnNoneIfOptional()
    return float(val_str)

  def FormatValue(self, value):
    return str(value)

class FloatValidator(TextCtrlValidator):
  def CheckFormat(self, value):
    value = float(value)
    window = self.GetWindow()
    if (value > window.GetMax()):
      raise ValueError("Value exceeds maximum allowed (%d)." % window.GetMax())
    elif (value < window.GetMin()):
      raise ValueError("Value is less than minimum allowed (%d)." %
        window.GetMin())
    return value #return window.FormatValue(value)

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Float control test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "High resolution:", pos=(100,180))
  float_ctrl = FloatCtrl(panel, -1, pos=(300,180), size=(80,-1),
    value=3,
    name="High resolution")
  float_ctrl.SetMax(20)
  float_ctrl.SetMin(1.5)
  float_ctrl.SetOptional(False)
  txt2 = wx.StaticText(panel, -1, "mFo-DFc level", pos=(100,240))
  float_ctrl2 = FloatCtrl(panel, -1, pos=(300,240), size=(80,-1),
    name="mFo-DFc level")
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay(evt):
    float1 = float_ctrl.GetPhilValue()
    float2 = float_ctrl2.GetPhilValue()
    print(type(float1).__name__, str(float1))
    print(type(float2).__name__, str(float2))
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
