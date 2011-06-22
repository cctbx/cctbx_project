
from wxtbx.phil_controls import ValidatedTextCtrl, TextCtrlValidator
import wx
import re

def format_int_list (value) :
  format = " ".join([ "%d" for x in range(len(value)) ])
  return format % tuple(value)

class IntsControl (ValidatedTextCtrl) :
  def __init__ (self, *args, **kwds) :
    super(IntsControl, self).__init__(*args, **kwds)

  def CreateValidator (self) :
    return IntsValidator()

  def SetInts (self, value) :
    if (value is None) :
      ValidatedTextCtrl.SetValue(self, "")
    elif (isinstance(value, int) or isinstance(value, float)) :
      ValidatedTextCtrl.SetValue(self, str(int(value)))
    elif (isinstance(value, list) or isinstance(value, tuple)) :
      ValidatedTextCtrl.SetValue(self, format_int_list(value))

  def SetValue (self, value) :
    self.SetInts(value)

  def GetPhilValue (self) :
    self.Validate()
    val_str = str(ValidatedTextCtrl.GetValue(self))
    if (val_str == "") :
      return None
    return [ int(field) for field in val_str.split() ]

class IntsValidator (TextCtrlValidator) :
  def CheckFormat (self, value) :
    if ("," in value) :
      value = re.sub(",", " ", value)
    ints_list = [ int(field) for field in value.split() ]
    print ints_list
    return format_int_list(ints_list)

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Integer list test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Run SA on cycle numbers:", pos=(100,180))
  ints_ctrl = IntsControl(panel, -1, pos=(300,180), size=(200,-1),
    name="Run SA on cycle numbers")
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay (evt) :
    ints = ints_ctrl.GetPhilValue()
    print type(ints).__name__, str(ints)
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
