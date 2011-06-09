
from wxtbx.phil_controls import space_group, unit_cell
from libtbx.utils import Sorry
import wx

class SymmetryDialog (wx.Dialog) :
  def __init__ (self, *args, **kwds) :
    super(SymmetryDialog, self).__init__(*args, **kwds)
    style = self.GetWindowStyle()
    style |= wx.WS_EX_VALIDATE_RECURSIVELY|wx.RAISED_BORDER|wx.CAPTION
    self.SetWindowStyle(style)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    szr.Add(szr2, 1, wx.ALL|wx.EXPAND, 10)
    txt = wx.StaticText(self, -1,
      "Missing or incomplete symmetry information.  Please enter a space "+
      "group and unit cell.")
    txt.Wrap(480)
    szr2.Add(txt, 0, wx.ALL, 5)
    szr3 = wx.FlexGridSizer(rows=2, cols=2)
    txt2 = wx.StaticText(self, -1, "Unit cell:")
    self.unit_cell_ctrl = unit_cell.UnitCellControl(
      parent=self,
      id=-1,
      size=(300,-1))
    txt3 = wx.StaticText(self, -1, "Space group:")
    self.space_group_ctrl = space_group.SpaceGroupControl(
      parent=self,
      id=-1)
    szr3.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr3.Add(self.unit_cell_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr3.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr3.Add(self.space_group_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr2.Add(szr3, 0, wx.ALL, 0)
    cancel_btn = wx.Button(self, wx.ID_CANCEL)
    ok_btn = wx.Button(self, wx.ID_OK)
    ok_btn.SetDefault()
    szr4 = wx.StdDialogButtonSizer()
    szr4.Add(cancel_btn)
    szr4.Add(ok_btn, 0, wx.LEFT, 5)
    szr2.Add(szr4, 0, wx.ALL|wx.ALIGN_RIGHT, 5)
    szr.Layout()
    self.Fit()

  def SetUnitCell (self, uc) :
    self.unit_cell_ctrl.SetValue(uc)

  def SetSpaceGroup (self, sg) :
    self.space_group_ctrl.SetValue(sg)

  def GetSymmetry (self, allow_incomplete=False) :
    uc = self.unit_cell_ctrl.GetPhilValue()
    sg = self.space_group_ctrl.GetPhilValue()
    if (not allow_incomplete) :
      if (uc is None) :
        raise Sorry("Missing unit cell parameters.")
      elif (sg is None) :
        raise Sorry("Missing space group.")
    from cctbx import crystal
    symm = crystal.symmetry(
      unit_cell=uc,
      space_group_info=sg)
    return symm

  def OnOkay (self, event) :
    if (not self.Validate()) :
      pass
    else :
      symm = self.GetSymmetry()
      self.EndModal(wx.ID_OK)

  def OnCancel (self, event) :
    self.EndModal(wx.ID_CANCEL)

if (__name__ == "__main__") :
  app = wx.App(0)
  dlg = SymmetryDialog(None, -1, "Enter symmetry")
  dlg.SetSpaceGroup("P21")
  if (dlg.ShowModal() == wx.ID_OK) :
    symm = dlg.GetSymmetry()
    assert (symm.space_group_info() is not None)
    assert (symm.unit_cell() is not None)
  wx.CallAfter(dlg.Destroy)
