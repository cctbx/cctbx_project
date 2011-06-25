
from wxtbx.phil_controls import ValidatedTextCtrl, TextCtrlValidator
import wx

class UnitCellControl (ValidatedTextCtrl) :
  def __init__ (self, *args, **kwds) :
    super(UnitCellControl, self).__init__(*args, **kwds)
    self.SetToolTip(wx.ToolTip(
      "If all unit cell edges are perpendicular to each "+
      "other, you only need to enter the edge lengths; the angles will be "+
      "filled in as 90 degrees each.  You can trigger this conversion by "+
      "typing in the cell lengths and pressing 'Enter'."))

  def CreateValidator (self) :
    return UnitCellValidator()

  def SetUnitCell (self, uc) :
    if (type(uc).__name__ == "unit_cell") :
      uc = self.FormatValue(uc)
    elif (uc is None) :
      uc = ""
    assert isinstance(uc, str)
    wx.TextCtrl.SetValue(self, uc)

  def SetValue (self, value) :
    self.SetUnitCell(value)

  def FormatValue (self, value) :
    if (value is None) :
      return "None"
    return "%g %g %g %g %g %g" % value.parameters()

  def GetPhilValue (self) :
    self.Validate()
    val_str = str(wx.TextCtrl.GetValue(self))
    if (val_str == "") :
      return None
    from cctbx import uctbx
    return uctbx.unit_cell(val_str)

class UnitCellValidator (TextCtrlValidator) :
  def CheckFormat (self, value) :
    if (len(value) == 0) :
      return ""
    from cctbx import uctbx
    uc = uctbx.unit_cell(value)
    return self.GetWindow().FormatValue(uc)

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Unit cell test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Unit cell:", pos=(100,180))
  sg_ctrl = UnitCellControl(panel, -1, pos=(200,180), size=(300,-1),
    name="Unit cell")
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay (evt) :
    uc = sg_ctrl.GetPhilValue()
    print type(uc).__name__, str(uc)
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
