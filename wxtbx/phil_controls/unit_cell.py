
import wx

def format_unit_cell (uc) :
  if (uc is None) :
    return ""
  return "%g %g %g %g %g %g" % uc.parameters()

class UnitCellControl (wx.TextCtrl) :
  def __init__ (self, *args, **kwds) :
    super(UnitCellControl, self).__init__(*args, **kwds)
    style = self.GetWindowStyle()
    if (not style & wx.TE_PROCESS_ENTER) :
      style |= wx.TE_PROCESS_ENTER
      self.SetWindowStyle(style)
    self.SetValidator(UnitCellValidator())
    self.Bind(wx.EVT_TEXT_ENTER, lambda evt: self.Validate(), self)

  def SetUnitCell (self, uc) :
    if (type(uc).__name__ == "unit_cell") :
      uc = format_unit_cell(uc)
    elif (uc is None) :
      uc = ""
    assert isinstance(uc, str)
    wx.TextCtrl.SetValue(self, uc)

  def SetValue (self, value) :
    self.SetUnitCell(value)

  def GetPhilValue (self) :
    self.Validate()
    val_str = wx.TextCtrl.GetValue(self)
    if (val_str == "") :
      return None
    from cctbx import uctbx
    return uctbx.unit_cell(val_str)

class UnitCellValidator (wx.PyValidator) :
  def __init__ (self) :
    wx.PyValidator.__init__(self)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnEnter)

  def Clone (self) :
    return self.__class__()

  def TransferToWindow (self) :
    return True

  def TransferFromWindow (self) :
    return True

  def Validate (self, win) :
    ctrl = self.GetWindow()
    uc_str = ctrl.GetValue()
    if (len(uc_str) == 0) :
      return True
    try :
      from cctbx import uctbx
      uc = uctbx.unit_cell(uc_str)
      ctrl.SetValue(format_unit_cell(uc))
      ctrl.SetBackgroundColour(
        wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
      #ctrl.SetFocus()
      ctrl.Refresh()
      return True
    except Exception, e :
      wx.MessageBox(caption="Format error", message=str(e))
      ctrl.SetBackgroundColour("red")
      ctrl.SetFocus()
      ctrl.Refresh()
      return False

  def OnEnter (self, event) :
    self.Validate(None)

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Unit cell test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Unit cell:", pos=(100,180))
  sg_ctrl = UnitCellControl(panel, -1, pos=(200,180), size=(300,-1))
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  def OnOkay (evt) :
    uc = sg_ctrl.GetPhilValue()
    print type(uc).__name__, str(uc)
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
