from __future__ import absolute_import, division, print_function

from wxtbx.phil_controls import space_group, unit_cell
from wxtbx.utils import add_ok_cancel_buttons
import wxtbx.icons
from libtbx.utils import Sorry, Abort
import wx
from six.moves import zip

class SymmetryDialog(wx.Dialog):
  def __init__(self, *args, **kwds):
    caption = kwds.get("caption",
      "Missing or incomplete symmetry information.  Please enter a space "+
      "group and unit cell.")
    if ("caption" in kwds):
      del kwds['caption']
    super(SymmetryDialog, self).__init__(*args, **kwds)
    style = self.GetWindowStyle()
    style |= wx.WS_EX_VALIDATE_RECURSIVELY|wx.RAISED_BORDER|wx.CAPTION
    self.SetWindowStyle(style)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    szr.Add(szr2, 0, wx.ALL|wx.EXPAND, 10)
    caption_sizer = wx.BoxSizer(wx.HORIZONTAL)
    bmp = wx.StaticBitmap(self, -1, wxtbx.icons.symmetry.GetBitmap())
    caption_sizer.Add(bmp, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
    txt = wx.StaticText(self, -1, caption)
    txt.Wrap(480)
    caption_sizer.Add(txt, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
    szr2.Add(caption_sizer, 0, wx.ALL, 0)
    szr3 = wx.FlexGridSizer(rows=3, cols=2)
    szr2.Add(szr3, 0, wx.ALL, 0)
    txt2 = wx.StaticText(self, -1, "Unit cell:")
    self.unit_cell_ctrl = unit_cell.UnitCellCtrl(
      parent=self,
      id=-1,
      size=(300,-1),
      name="Unit cell")
    txt3 = wx.StaticText(self, -1, "Space group:")
    self.space_group_ctrl = space_group.SpaceGroupCtrl(
      parent=self,
      id=-1,
      name="Space group")
    szr3.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr3.Add(self.unit_cell_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr3.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr3.Add(self.space_group_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr3.Add((1,1), 0, wx.ALL, 5)
    load_btn = wx.Button(self, -1, "Load symmetry from file...")
    szr3.Add(load_btn, 0, wx.ALL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnLoadSymmetry, load_btn)
    add_ok_cancel_buttons(self, szr2)
    szr.Layout()
    self.Fit()
    self.Centre(wx.BOTH)

  def SetUnitCell(self, uc):
    self.unit_cell_ctrl.SetValue(uc)

  def SetSpaceGroup(self, sg):
    self.space_group_ctrl.SetValue(sg)

  def SetSymmetry(self, symm):
    if (symm is not None):
      self.SetSpaceGroup(symm.space_group_info())
      self.SetUnitCell(symm.unit_cell())

  def GetSymmetry(self, allow_incomplete=False):
    uc = self.unit_cell_ctrl.GetPhilValue()
    sg = self.space_group_ctrl.GetPhilValue()
    if (not allow_incomplete):
      if (uc is None):
        raise Sorry("Missing unit cell parameters.")
      elif (sg is None):
        raise Sorry("Missing space group.")
    from cctbx import crystal
    symm = crystal.symmetry(
      unit_cell=uc,
      space_group_info=sg)
    return symm

  def OnLoadSymmetry(self, event):
    file_name = wx.FileSelector(
      message="Select a reflection or PDB file containing symmetry",
      flags=wx.FD_OPEN)
    if (file_name != ""):
      from iotbx import crystal_symmetry_from_any
      symm = crystal_symmetry_from_any.extract_from(file_name)
      if (symm is not None):
        space_group = symm.space_group_info()
        if (space_group is not None):
          self.space_group_ctrl.SetSpaceGroup(space_group)
        unit_cell = symm.unit_cell()
        if (unit_cell is not None):
          self.unit_cell_ctrl.SetUnitCell(unit_cell)
      else :
        raise Sorry("This file does not contain valid symmetry information.")

  def OnOkay(self, event):
    print(1)
    if (not self.Validate()):
      pass
    else :
      symm = self.GetSymmetry()
      self.EndModal(wx.ID_OK)

  def OnCancel(self, event):
    self.EndModal(wx.ID_CANCEL)

class SymmetryChoiceDialog(wx.Dialog):
  def __init__(self, *args, **kwds):
    default_caption = "The input data specify multiple unique symmetry "+ \
      "and/or unit cell parameters.  Please select the consensus settings "+\
      "from the menus below."
    caption = kwds.pop("caption", default_caption)
    if (caption is None):
      caption = default_caption
    super(SymmetryChoiceDialog, self).__init__(*args, **kwds)
    style = self.GetWindowStyle()
    style |= wx.WS_EX_VALIDATE_RECURSIVELY|wx.RAISED_BORDER|wx.CAPTION
    self.SetWindowStyle(style)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    szr.Add(szr2, 1, wx.ALL|wx.EXPAND, 10)
    caption_sizer = wx.BoxSizer(wx.HORIZONTAL)
    bmp = wx.StaticBitmap(self, -1, wxtbx.icons.symmetry.GetBitmap())
    caption_sizer.Add(bmp, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
    txt = wx.StaticText(self, -1, caption)
    txt.Wrap(480)
    caption_sizer.Add(txt, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
    szr2.Add(caption_sizer, 0, wx.ALL, 0)
    szr3 = wx.BoxSizer(wx.HORIZONTAL)
    txt2 = wx.StaticText(self, -1, "Space group:")
    szr3.Add(txt2, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
    self.space_group_chooser = wx.Choice(self, -1, size=(200,-1))
    szr3.Add(self.space_group_chooser, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
    szr2.Add(szr3, 0, wx.ALL, 5)
    szr4 = wx.BoxSizer(wx.HORIZONTAL)
    txt3 = wx.StaticText(self, -1, "Unit cell:")
    szr4.Add(txt3, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
    self.unit_cell_chooser = wx.Choice(self, -1, size=(400,-1))
    szr4.Add(self.unit_cell_chooser, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
    szr2.Add(szr4, 0, wx.ALL, 5)
    add_ok_cancel_buttons(self, szr2)
    szr.Layout()
    self.Fit()
    self.Centre(wx.BOTH)
    self.unit_cells = []
    self.space_groups = []

  def SetUnitCells(self, unit_cells, source_info=None):
    assert (source_info is None) or (len(source_info) == len(unit_cells))
    self.unit_cells = unit_cells
    if (source_info is None):
      choices = [ "%g %g %g %g %g %g" % uc.parameters() for uc in unit_cells ]
    else :
      choices = [ "%s (%s)" % ("%g %g %g %g %g %g" % uc.parameters(), si)
         for (uc, si) in zip(unit_cells, source_info) ]
    self.unit_cell_chooser.SetItems(choices)

  def SetSpaceGroups(self, space_groups, source_info=None):
    assert (source_info is None) or (len(source_info) == len(space_groups))
    self.space_groups = space_groups
    if (source_info is None):
      choices = [ str(sg) for sg in space_groups ]
    else :
      choices = [ "%s (%s)" % (str(sg), si)
                  for (sg, si) in zip(space_groups, source_info) ]
    self.space_group_chooser.SetItems(choices)

  def GetUnitCell(self):
    if (len(self.unit_cells) == 0) : return None
    sel = self.unit_cell_chooser.GetSelection()
    return self.unit_cells[sel]

  def GetSpaceGroup(self):
    if (len(self.space_groups) == 0) : return None
    sel = self.space_group_chooser.GetSelection()
    return self.space_groups[sel]

def get_unique_symmetry(
    space_groups=(),
    unit_cells=(),
    sg_source_info=None,
    uc_source_info=None,
    parent=None,
    caption=None):
  assert (len(space_groups) > 0) or (len(unit_cells) > 0)
  non_unique = False
  for i_cell, uc in enumerate(unit_cells):
    j_cell = i_cell + 1
    while (j_cell < len(unit_cells)):
      if (not uc.is_similar_to(unit_cells[j_cell])):
        non_unique = True
        break
      j_cell += 1
  for i_sg, sg in enumerate(space_groups):
    j_sg = i_sg + 1
    sg_number = sg.group().type().number()
    while (j_sg < len(space_groups)):
      other_number = space_groups[j_sg].group().type().number()
      if (sg_number != other_number):
        non_unique = True
        break
      j_sg += 1
  final_cell = final_group = None
  if (non_unique):
    dlg = SymmetryChoiceDialog(parent=parent,
      title="Select unit cell and/or space group",
      caption=caption)
    dlg.SetUnitCells(unit_cells=unit_cells, source_info=uc_source_info)
    dlg.SetSpaceGroups(space_groups=space_groups, source_info=sg_source_info)
    if (dlg.ShowModal() == wx.ID_OK):
      final_cell = dlg.GetUnitCell()
      final_group = dlg.GetSpaceGroup()
    wx.CallAfter(dlg.Destroy)
    if (final_cell is None) and (final_group is None):
      raise Abort()
  else :
    final_cell = unit_cells[0]
    final_group = space_groups[0]
  from cctbx import crystal
  return crystal.symmetry(
    unit_cell=final_cell,
    space_group_info=final_group)

if (__name__ == "__main__"):
  app = wx.App(0)
  dlg = SymmetryDialog(None, -1, "Enter symmetry")
  dlg.SetSpaceGroup("P21")
  if (dlg.ShowModal() == wx.ID_OK):
    symm = dlg.GetSymmetry()
    assert (symm.space_group_info() is not None)
    assert (symm.unit_cell() is not None)
  wx.CallAfter(dlg.Destroy)
  from cctbx import uctbx, sgtbx
  unit_cells = [
    uctbx.unit_cell((10,20,30,90,90,90)),
    uctbx.unit_cell((11,21,29,90.05,90.2,89.8)),
  ]
  source_info = [ "File 1", "File 2" ]
  space_groups = [sgtbx.space_group_info("P21"), sgtbx.space_group_info("P1")]
  symm = get_unique_symmetry(
    space_groups=space_groups,
    unit_cells=unit_cells,
    sg_source_info=source_info,
    uc_source_info=source_info)
  symm.show_summary()
