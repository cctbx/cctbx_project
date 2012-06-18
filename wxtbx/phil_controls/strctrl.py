
from wxtbx.phil_controls.text_base import ValidatedTextCtrl, TextCtrlValidator
import wx

class StrCtrl (ValidatedTextCtrl) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    if (kwds.get("size", wx.DefaultSize) == wx.DefaultSize) :
      kwds['size'] = (200,-1)
    super(StrCtrl, self).__init__(*args, **kwds)

  def CreateValidator (self) :
    return StrValidator()

  def SetValue (self, value) :
    if (value is None) :
      ValidatedTextCtrl.SetValue(self, "")
    else :
      assert isinstance(value, str)
      ValidatedTextCtrl.SetValue(self, value)

  def GetPhilValue (self) :
    self.Validate()
    val_str = self.GetValue()
    if (val_str in ["", "none", "None"]) :
      return self.ReturnNoneIfOptional()
    return val_str

  def GetStringValue (self) :
    value = self.GetPhilValue()
    if (value is None) :
      return "None"
    else :
      return '"""%s"""' % value

  def FormatValue (self, value) :
    return str(value)

class StrValidator (TextCtrlValidator) :
  def CheckFormat (self, value) :
    if (";" in value) :
      raise ValueError("Semicolons are not allowed in text input.")
    return value # XXX does anything else need to be done here?

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "String parameter test")
  panel = wx.Panel(frame, -1, size=(720,400))
  txt1 = wx.StaticText(panel, -1, "Job title:", pos=(10,100))
  ctrl1 = StrCtrl(panel, -1, value=None, pos=(160, 100), size=(400,-1),
    name="Job title")
  txt2 = wx.StaticText(panel, -1, "Output file prefix:", pos=(10,200))
  ctrl2 = StrCtrl(panel, -1, value="refine", pos=(160,200),
    name="Output file prefix")
  ctrl2.SetOptional(False)
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay (evt) :
    value1 = ctrl1.GetPhilValue()
    value2 = ctrl2.GetPhilValue()
    print value1
    print value2
  assert (ctrl1.GetPhilValue() is None)
  assert (ctrl1.GetStringValue() == "None")
  assert (ctrl2.GetPhilValue() == "refine")
  assert (ctrl2.GetStringValue() == '"""refine"""')
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  import wxtbx.phil_controls
  def OnChange (evt) :
    pass
  frame.Bind(wxtbx.phil_controls.EVT_PHIL_CONTROL, OnChange)
  frame.Fit()
  frame.Show()
  app.MainLoop()
