
"""
Interface for editing custom geometry restraints (as defined in
mmtbx.monomer_library.pdb_interpretation).  Currently only bonds and angles
are supported.  The parameters are specific to phenix.refine right now but
they could theoretically be used in other programs.
"""

import wxtbx.bitmaps
from wxtbx.phil_controls import simple_dialogs
import wx
import sys

#-----------------------------------------------------------------------
# Basic panel

class CustomRestraintsPanel (wx.Panel) :
  n_atom_selections = -1
  list_label = None
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self._default_label = "---"
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self.sizer = szr
    label = wx.StaticText(self, -1, self.list_label)
    font = label.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    label.SetFont(font)
    szr.Add(label, 0, wx.ALL, 5)
    self.CreateList()
    self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.OnSelect, self.lc)
    self.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnDeSelect, self.lc)
    self.sizer.Add(self.lc, 1, wx.EXPAND|wx.ALL, 5)
    self.buttons = wx.BoxSizer(wx.HORIZONTAL)
    add_btn = self.AddControlButton(
      label="Add",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "edit_add", 16))
    del_btn = self.AddControlButton(
      label="Delete",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "cancel", 16))
    del_btn = self.AddControlButton(
      label="Clear all",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "editdelete", 16))
    update_btn = self.AddControlButton(
      label="Update selections",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "recur", 16))
    options_btn = self.AddControlButton(
      label="Other options",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "run", 16))
    self.Bind(wx.EVT_BUTTON, self.OnAdd, add_btn)
    self.Bind(wx.EVT_BUTTON, self.OnDelete, del_btn)
    self.Bind(wx.EVT_BUTTON, self.OnDeleteAll, del_btn)
    self.Bind(wx.EVT_BUTTON, self.OnUpdate, update_btn)
    self.Bind(wx.EVT_BUTTON, self.OnOptionsMenu, options_btn)
    szr.Add(self.buttons, 0, wx.LEFT|wx.BOTTOM|wx.RIGHT, 5)
    self.sizer = szr
    self.CreateSelectionFields()
    szr.Layout()

  def CreateList (self) :
    raise NotImplementedError()

  def CreateSelectionFields (self) :
    raise NotImplementedError()

  def FillSelections (self) :
    selections = self.lc.GetSelections()
    assert (len(selections) == self.n_atom_selections)
    for i, sele in enumerate(selections) :
      if (sele is not None) :
        field = getattr(self, "selection_%d" % (i+1))
        field.SetValue(sele)
      else :
        field.SetValue("")

  def ClearSelections (self) :
    for i in range(self.n_atom_selections) :
      field = getattr(self, "selection_%d" % (i+1))
      field.SetValue("")

  def OnSelect (self, event) :
    self.FillSelections()

  def OnDeSelect (self, event) :
    self.ClearSelections()

  def OnDeleteAll (self, evt) :
    if (self.lc.GetItemCount() == 0) :
      return False
    confirm = wx.MessageBox(caption="Confirm delete",
      message="Are you sure you want to delete all items in the list?")
    if (confirm == wx.OK) :
      self.lc.DeleteAllItems()
      self.lc.ClearSelections()

  def __getattr__ (self, name) :
    return getattr(self.lc, name)

  def OnOptionsMenu (self, event) :
    raise NotImplementedError()

  def OnAdd (self, event) :
    self.lc.AddRestraint()

  def OnDelete (self, event) :
    self.lc.DeleteRestraint()

  def OnUpdate (self, event) :
    selections = []
    for i in range(self.n_atom_selections) :
      field = getattr(self, "selection_%d" % (i+1))
      value = field.GetValue()
      if (value.isspace()) :
        value = None
      selections.append(value)
    self.lc.SetSelections(selections)

class CustomBondPanel (CustomRestraintsPanel) :
  n_atom_selections = 2
  list_label = "Custom bond restraints:"
  def CreateList (self) :
    self.lc = BondRestraintsList(self)

  def CreateSelectionFields (self) :
    edit_szr = wx.FlexGridSizer(cols=2)
    self.sizer.Add(edit_szr, 0, wx.ALL)
    edit_szr.Add(wx.StaticText(self, -1, "Atom selections:"))
    for i in range(2) :
      edit_field = wx.TextCtrl(self, -1, size=(400,-1))
      setattr(self, "selection_%d" % (i+1), edit_field)
      edit_szr.Add(edit_field, 0, wx.ALL, 5)
      if (i == 1) : edit_szr.Add((1,1))

class CustomAnglePanel (CustomRestraintsPanel) :
  n_atom_selections = 3
  list_label = "Custom angle restraints:"
  def CreateList (self) :
    self.lc = AngleRestraintsList(self, size=(800,200))

  def CreateSelectionFields (self) :
    edit_szr = wx.FlexGridSizer(cols=2)
    self.sizer.Add(edit_szr, 0, wx.ALL)
    edit_szr.Add(wx.StaticText(self, -1, "Atom selections:"))
    for i in range(3) :
      edit_field = wx.TextCtrl(self, -1, size=(400,-1))
      setattr(self, "selection_%d" % (i+1), edit_field)
      edit_szr.Add(edit_field, 0, wx.ALL, 5)
      if (i == 0) : edit_szr.Add((1,1))

#-----------------------------------------------------------------------
# ListCtrl-derived classes

class RestraintsListBase (wx.ListCtrl) :
  restraint_name = None
  n_columns = -1
  def __init__ (self, *args, **kwds) :
    kwds['style'] = wx.LC_REPORT|wx.LC_SINGLE_SEL
    wx.ListCtrl.__init__(self, *args, **kwds)
    self.CreateColumns()
    self._params = []
    self._prefix = None

  def GetPhilPath (self) :
    if (self._prefix is not None) :
      return self._prefix + ".geometry_restraints.edits." + self.restraint_name
    else :
      return "geometry_restraints.edits." + self.restraint_name

  def ClearRestraints (self) :
    self._params = []

  def AddRestraint (self) :
    item = self.InsertStringItem(sys.maxint, "---")
    for i in range(self.n_columns - 1) :
      self.SetStringItem(item, i+1, "---")
    new_params = self._index.get_template_copy(self.GetPhilPath())
    self._params.append(new_params)
    self.SelectItem(item)

  def DeleteRestraint (self) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      assert (item < len(self._params))
      self.DeleteItem(item)
      del self._params[item]

  def CreateColumns (self) :
    raise NotImplementedError()

  def SetPhilIndex (self, index) :
    self._index = index

  def SetParams (self, params) :
    self._params = params
    self.PopulateList()

  def UpdateSelections (self, selections) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      for i in range(len(selections)) :
        sel_attr = "atom_selection_%d" % (i+1)
        assert hasattr(self._params[item], sel_attr)
        setattr(self._params[item], sel_attr, selections[i])

  def GetParams (self) :
    return self._params

  def GetPhil (self) :
    return self._index.master_phil.format(python_object=self._params)

class BondRestraintsList (RestraintsListBase) :
  restraint_name = "bond"
  def CreateColumns (self) :
    self.InsertColumn(0, "Selections", width=400)
    self.InsertColumn(1, "Sym. exp.", width=120)
    self.InsertColumn(2, "Distance", width=100)
    self.InsertColumn(3, "Sigma", width=80)
    self.InsertColumn(4, "Slack", width=80)

  def GetSelections (self) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      bond = self._params[item]
      return [bond.atom_selection_1, bond.atom_selection_2]
    return [""] * 2

  def SetSelections (self, selections) :
    assert (len(selections) == 2)
    item = self.GetFirstSelected()
    if (item >= 0) :
      bond = self._params[item]
      bond.atom_selection_1 = selections[0]
      bond.atom_selection_2 = selections[1]

  def OnSetSymexp (self, event) :
    pass # TODO

  def OnSetDistance (self, event) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      dlg = simple_dialogs.FloatDialog(
        parent=self,
        title="Ideal distance",
        label="Bond length",
        value=self._params[item].distance_ideal,
        caption="Please enter an ideal distance (in Angstroms) for the "+
          "selected bond; it will usually be a decimal number between 1 "+
          "and 4.")
      dlg.SetMin(0.1)
      dlg.SetOptional(False)
      if (dlg.ShowModal() == wx.ID_OK) :
        distance_ideal = dlg.GetPhilValue()
        self._params[item].distance_ideal = distance_ideal
        self.SetStringItem(item, 2, "%g" % distance_ideal)

  def OnSetSigma (self, event) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      dlg = simple_dialogs.FloatDialog(
        parent=self,
        title="Restraint sigma",
        label="Bond sigma",
        value=self._params[item].sigma,
        caption="Please enter a sigma value (in Angstroms) for the selected "+
          "bond; it should be a decimal number between 0.01 and 1.")
      dlg.SetMin(0.001)
      dlg.SetOptional(False)
      if (dlg.ShowModal() == wx.ID_OK) :
        sigma = dlg.GetPhilValue()
        self._params[item].sigma = sigma
        self.SetStringItem(item, 3, "%g" % sigma)

  def OnSetSlack (self, event) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      dlg = simple_dialogs.FloatDialog(
        parent=self,
        title="Restraint slack",
        label="Bond slack",
        value=self._params[item].sigma,
        caption="Please enter a slack value (in Angstroms) for the selected "+
          "bond; this is the distance the bond is allowed to deviate from "+
          "ideal before the restraint will take effect.  By default this is "+
          "undefined (or zero).")
      dlg.SetMin(0)
      dlg.SetOptional(True)
      if (dlg.ShowModal() == wx.ID_OK) :
        slack = dlg.GetPhilValue()
        self._params[item].slack = slack
        if (slack is not None) :
          self.SetStringItem(item, 4, "%g" % slack)
        else :
          self.SetStringItem(item, 4, "---")

class AngleRestraintsList (RestraintsListBase) :
  restraint_name = "angle"
  def CreateColumns (RestraintsListBase) :
    self.InsertColumn(0, "Selections", width=600)
    self.InsertColumn(1, "Angle", width=100)
    self.InsertColumn(2, "Sigma", width=80)

  def GetSelections (self) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      angle = self._params[item]
      return [angle.atom_selection_1, angle.atom_selection_2,
              angle.atom_selection_3]
    return [""] * 3

  def SetSelections (self, selections) :
    assert (len(selections) == 3)
    item = self.GetFirstSelected()
    if (item >= 0) :
      angle = self._params[item]
      angle.atom_selection_1 = selections[0]
      angle.atom_selection_2 = selections[1]
      angle.atom_selection_3 = selections[2]

  def OnSetAngle (self, event) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      dlg = simple_dialogs.FloatDialog(
        parent=self,
        title="Ideal angle",
        label="Bond angle",
        value=self._params[item].sigma,
        caption="Please enter an angle value (in degrees) for the selected "+
          "angle; it usually be a decimal number between 1 and 180.")
      dlg.SetMin(1)
      dlg.SetMax(180)
      dlg.SetOptional(False)
      if (dlg.ShowModal() == wx.ID_OK) :
        angle = dlg.GetPhilValue()
        self._params[item].angle = angle
        self.SetStringItem(item, 1, "%g" % angle)

  def OnSetSigma (self, event) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      dlg = simple_dialogs.FloatDialog(
        parent=self,
        title="Restraint sigma",
        label="Bond sigma",
        value=self._params[item].sigma,
        caption="Please enter a sigma value (in degrees) for the selected "+
          "angle; it should usually be a decimal number between 1 and 5.")
      dlg.SetMin(0.1)
      dlg.SetOptional(False)
      if (dlg.ShowModal() == wx.ID_OK) :
        sigma = dlg.GetPhilValue()
        self._params[item].sigma = sigma
        self.SetStringItem(item, 2, "%g" % sigma)

#-----------------------------------------------------------------------
# Frame class

class RestraintsFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.SetSizer(self.sizer)

#-----------------------------------------------------------------------
# REGRESSION TESTING
if (__name__ == "__main__") :
  pass
