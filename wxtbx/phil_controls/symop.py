
from wxtbx.phil_controls.text_base import ValidatedTextCtrl, TextCtrlValidator
from libtbx.utils import Sorry
from libtbx import Auto
import wx

class SymopCtrl (ValidatedTextCtrl) :
  def CreateValidator (self) :
    return SymopValidator()

  def SetSymop (self, value) :
    if (value is None) or (value is Auto) :
      ValidatedTextCtrl.SetValue(self, "")
    elif (isinstance(value, str)) :
      try :
        from cctbx import sgtbx
        rt_mx = sgtbx.rt_mx(symbol=value)
      except ValueError :
        raise Sorry("Inappropriate value '%s' for %s." % (value,
          sself.GetName()))
      else :
        ValidatedTextCtrl.SetValue(self, str(value))
    else :
      raise TypeError("Type '%s' not allowed!" % type(value).__name__)

  def SetValue (self, value) :
    self.SetSymop(value)

  def GetPhilValue (self) :
    self.Validate()
    val_str = ValidatedTextCtrl.GetValue(self)
    if (val_str == "") :
      return self.ReturnNoneIfOptional()
    return val_str

  def FormatValue (self, value) :
    return str(value)

class SymopValidator (TextCtrlValidator) :
  def CheckFormat (self, value) :
    from cctbx import sgtbx
    rt_mx = sgtbx.rt_mx(symbol=value)
    return value

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Symop control test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Symmetry operator:", pos=(100,180))
  sym_ctrl = SymopCtrl(panel, -1, pos=(300,180), size=(80,-1),
      name="Symmetry operator")
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay (evt) :
    symop = sym_ctrl.GetPhilValue()
    print symop
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
