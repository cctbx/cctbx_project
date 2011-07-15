
from wxtbx import phil_controls
from wxtbx.phil_controls import ValidatedTextCtrl, TextCtrlValidator
import wx
import re

class FloatsCtrl (ValidatedTextCtrl, phil_controls.PhilCtrl) :
  def __init__ (self, *args, **kwds) :
    super(FloatsCtrl, self).__init__(*args, **kwds)
    self.elems = None

  def SetElems (self, elems) :
    assert (elems is None) or (isinstance(elems, int))
    self.elems = elems

  def CreateValidator (self) :
    return FloatsValidator()

  def SetFloats (self, value) :
    if (value is None) :
      ValidatedTextCtrl.SetValue(self, "")
    elif (isinstance(value, int) or isinstance(value, float)) :
      ValidatedTextCtrl.SetValue(self, str(float(value)))
    elif (isinstance(value, list) or isinstance(value, tuple)) :
      ValidatedTextCtrl.SetValue(self, self.FormatValue(value))

  def SetValue (self, value) :
    self.SetFloats(value)

  def GetPhilValue (self) :
    self.Validate()
    val_str = str(ValidatedTextCtrl.GetValue(self))
    if (val_str == "") :
      return None
    return [ float(field) for field in val_str.split() ]

  def FormatValue (self, value) :
    format = " ".join([ "%g" for x in range(len(value)) ])
    return format % tuple(value)

class FloatsValidator (TextCtrlValidator) :
  def CheckFormat (self, value) :
    if ("," in value) :
      value = re.sub(",", " ", value)
    floats_list = [ float(field) for field in value.split() ]
    if (self.GetWindow().elems is not None) :
      if (self.GetWindow().elems != len(floats_list)) :
        raise RuntimeError(("Wrong number of items - %d values required, "+
          "but %d are entered.") % (self.GetWindow().elems, len(floats_list)))
    return self.GetWindow().FormatValue(floats_list)

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Floating-point list test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Search model RMSDs:", pos=(100,180))
  floats_ctrl = FloatsCtrl(panel, -1, pos=(300,180), size=(200,-1),
    value=[1.0,2.0,3.5,4.7],
    name="Search model RMSDs")
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay (evt) :
    floats = floats_ctrl.GetPhilValue()
    print type(floats).__name__, str(floats)
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
