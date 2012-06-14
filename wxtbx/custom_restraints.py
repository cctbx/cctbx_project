
"""
Interface for editing custom geometry restraints (as defined in
mmtbx.monomer_library.pdb_interpretation).  Currently only bonds and angles
are supported.  The parameters are specific to phenix.refine right now but
they could theoretically be used in other programs.
"""

import wxtbx.bitmaps
from wxtbx import metallicbutton
from wxtbx.phil_controls import simple_dialogs
from wx.lib.agw import flatnotebook
import wx
from libtbx import str_utils
import sys

ATOM_SELECTION_BUTTONS = 1

#-----------------------------------------------------------------------
# Basic panel

class CustomRestraintsPanel (wx.Panel) :
  n_atom_selections = -1
  list_label = None
  def __init__ (self, *args, **kwds) :
    style = kwds.get('style', 0)
    self.flag_atom_selection_buttons = (style & ATOM_SELECTION_BUTTONS)
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
    clear_btn = self.AddControlButton(
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
    self.Bind(wx.EVT_BUTTON, self.OnDeleteAll, clear_btn)
    self.Bind(wx.EVT_BUTTON, self.OnUpdate, update_btn)
    self.Bind(wx.EVT_BUTTON, self.OnOptionsMenu, options_btn)
    szr.Add(self.buttons, 0, wx.LEFT|wx.BOTTOM|wx.RIGHT, 5)
    self.sizer = szr
    self.CreateSelectionFields()
    szr.Layout()

  def AddControlButton (self, label, bitmap) :
    btn = metallicbutton.MetallicButton(
      parent=self,
      label=label,
      bmp=bitmap,
      highlight_color=(200,220,240))
    self.buttons.Add(btn, 0, wx.RIGHT, 5)
    return btn

  def CreateList (self) :
    raise NotImplementedError()

  def CreateSelectionFields (self) :
    raise NotImplementedError()

  def FillSelections (self) :
    selections = self.lc.GetSelections()
    assert (len(selections) == self.n_atom_selections)
    for i, sele in enumerate(selections) :
      field = self.GetSelectionControl(i)
      if (sele is not None) :
        field.SetValue(sele)
      else :
        field.SetValue("")

  def GetSelectionControl (self, i) :
    return getattr(self, "selection_%d" % (i+1))

  def ClearSelections (self) :
    for i in range(self.n_atom_selections) :
      field = self.GetSelectionControl(i)
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
      self.lc.ClearRestraints()
      self.ClearSelections()

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
      field = self.GetSelectionControl(i)
      value = field.GetValue()
      if (value.isspace()) or (value == "") :
        value = None
      selections.append(value)
    self.lc.SetSelections(selections)

  def CreateAtomSelectionButton (self, selection_ctrl) :
    btn = AtomSelectionButton(
      parent=self,
      selection_ctrl=selection_ctrl)
    self.Bind(wx.EVT_BUTTON, self.GetTopLevelParent().OnView, btn)
    return btn

class AtomSelectionButton (metallicbutton.MetallicButton) :
  def __init__ (self, parent, selection_ctrl) :
    metallicbutton.MetallicButton.__init__(self,
      parent=parent,
      label="Select...",
      bmp=wxtbx.bitmaps.fetch_icon_bitmap("actions", "select", 16),
      highlight_color=(200,220,240))
    self.selection_ctrl = selection_ctrl

class CustomBondPanel (CustomRestraintsPanel) :
  n_atom_selections = 2
  list_label = "Custom bond restraints:"
  def CreateList (self) :
    self.lc = BondRestraintsList(self)

  def CreateSelectionFields (self) :
    edit_szr = wx.FlexGridSizer(cols=3)
    self.sizer.Add(edit_szr, 0, wx.ALL, 5)
    edit_szr.Add(wx.StaticText(self, -1, "Atom selections:"))
    for i in range(2) :
      edit_field = wx.TextCtrl(self, -1, size=(540,-1),
        style=wx.TE_PROCESS_ENTER, name="bond")
      setattr(self, "selection_%d" % (i+1), edit_field)
      edit_szr.Add(edit_field, 0, wx.ALL, 5)
      select_btn = (1,1)
      if (self.flag_atom_selection_buttons) :
        select_btn = self.CreateAtomSelectionButton(edit_field)
      edit_szr.Add(select_btn, 0, wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
      if (i == 0) : edit_szr.Add((1,1))
      self.Bind(wx.EVT_TEXT_ENTER, self.OnUpdate, edit_field)

  def OnOptionsMenu (self, event) :
    menu = wx.Menu()
    item1 = menu.Append(-1, "Set ideal distance")
    item2 = menu.Append(-1, "Set distance sigma")
    item3 = menu.Append(-1, "Set bond slack")
    item4 = menu.Append(-1, "Set symmetry operator")
    self.Bind(wx.EVT_MENU, self.lc.OnSetDistance, item1)
    self.Bind(wx.EVT_MENU, self.lc.OnSetSigma, item2)
    self.Bind(wx.EVT_MENU, self.lc.OnSetSlack, item3)
    self.Bind(wx.EVT_MENU, self.lc.OnSetSymop, item4)
    event.GetEventObject().PopupMenu(menu)
    menu.Destroy()

class CustomAnglePanel (CustomRestraintsPanel) :
  n_atom_selections = 3
  list_label = "Custom angle restraints:"
  def CreateList (self) :
    self.lc = AngleRestraintsList(self, size=(800,200))

  def CreateSelectionFields (self) :
    edit_szr = wx.FlexGridSizer(cols=3)
    self.sizer.Add(edit_szr, 0, wx.ALL, 5)
    edit_szr.Add(wx.StaticText(self, -1, "Atom selections:"))
    for i in range(3) :
      edit_field = wx.TextCtrl(self, -1, size=(540,-1),
        style=wx.TE_PROCESS_ENTER, name="angle")
      setattr(self, "selection_%d" % (i+1), edit_field)
      edit_szr.Add(edit_field, 0, wx.ALL, 5)
      select_btn = (1,1)
      if (self.flag_atom_selection_buttons) :
        select_btn = self.CreateAtomSelectionButton(edit_field)
      edit_szr.Add(select_btn, 0, wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
      if (i < 2) : edit_szr.Add((1,1))
      self.Bind(wx.EVT_TEXT_ENTER, self.OnUpdate, edit_field)

  def OnOptionsMenu (self, event) :
    menu = wx.Menu()
    item1 = menu.Append(-1, "Set ideal angle")
    item2 = menu.Append(-1, "Set angle sigma")
    self.Bind(wx.EVT_MENU, self.lc.OnSetAngle, item1)
    self.Bind(wx.EVT_MENU, self.lc.OnSetSigma, item2)
    event.GetEventObject().PopupMenu(menu)
    menu.Destroy()

class CustomPlanarityPanel (CustomRestraintsPanel) :
  n_atom_selections = 1
  list_label = "Custom planarity restraints:"
  def CreateList (self) :
    self.lc = PlanarityRestraintsList(self, size=(800,200))

  def CreateSelectionFields (self) :
    edit_szr = wx.FlexGridSizer(cols=3)
    self.sizer.Add(edit_szr, 0, wx.ALL, 5)
    edit_szr.Add(wx.StaticText(self, -1, "Atom selection:"))
    self.selection_1 = wx.TextCtrl(self, -1, size=(540,64),
      style=wx.TE_PROCESS_ENTER, name="plane")
    edit_szr.Add(self.selection_1, 0, wx.ALL, 5)
    select_btn = (1,1)
    if (self.flag_atom_selection_buttons) :
      select_btn = self.CreateAtomSelectionButton(self.selection_1)
    edit_szr.Add(select_btn, 0, wx.TOP|wx.BOTTOM|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnUpdate, self.selection_1)

  def OnOptionsMenu (self, event) :
    menu = wx.Menu()
    item1 = menu.Append(-1, "Set planarity sigma")
    self.Bind(wx.EVT_MENU, self.lc.OnSetSigma, item1)
    event.GetEventObject().PopupMenu(menu)
    menu.Destroy()

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

  def SetPhilPrefix (self, prefix) :
    self._prefix = prefix

  def ClearRestraints (self) :
    self._params = []

  def AddRestraint (self) :
    item = self.InsertStringItem(sys.maxint, "---")
    for i in range(self.n_columns - 1) :
      self.SetStringItem(item, i+1, "---")
    new_params = self._index.get_template_copy(self.GetPhilPath()).extract()
    self._params.append(new_params)
    self.Select(item)

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
    if (self.GetItemCount() > 0) :
      self.DeleteAllItems()
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

def fv (fs, val) :
  return str_utils.format_value(fs, val, replace_none_with="---")

class BondRestraintsList (RestraintsListBase) :
  restraint_name = "bond"
  def CreateColumns (self) :
    self.InsertColumn(0, "Selections", width=400)
    self.InsertColumn(1, "Distance", width=100)
    self.InsertColumn(2, "Sigma", width=80)
    self.InsertColumn(3, "Slack", width=80)
    self.InsertColumn(4, "Sym. exp.", width=120)

  def PopulateList (self) :
    for bond in self._params :
      selections = [bond.atom_selection_1, bond.atom_selection_2]
      sele_str = "; ".join([ fv("%s", s) for s in selections ])
      item = self.InsertStringItem(sys.maxint, sele_str)
      self.SetStringItem(item, 1, fv("%g", bond.distance_ideal))
      self.SetStringItem(item, 2, fv("%g", bond.sigma))
      self.SetStringItem(item, 3, fv("%g", bond.slack))
      self.SetStringItem(item, 4, fv("%s", bond.symmetry_operation))

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
      sele_text = "; ".join([ str(s) for s in selections ])
      self.SetStringItem(item, 0, sele_text)

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
        self.SetStringItem(item, 1, "%g" % distance_ideal)
      dlg.Destroy()

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
        self.SetStringItem(item, 2, fv("%g", sigma))
      dlg.Destroy()

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
        self.SetStringItem(item, 3, fv("%g", slack))
      dlg.Destroy()

  def OnSetSymop (self, event) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      dlg = simple_dialogs.SymopDialog(
        parent=self,
        title="Symmetry operator",
        label="Symmetry operator",
        value=self._params[item].symmetry_operation,
        caption="You may specify a bond across symmetry mates; in this case "+
          "the symmetry operation will be applied to the coordinates of the "+
          "second atom when calculating the distance.  The operation should "+
          "be in a format like 'x,y,-z'.")
      if (dlg.ShowModal() == wx.ID_OK) :
        symop = dlg.GetPhilValue()
        self._params[item].symmetry_operation = symop
        self.SetStringItem(item, 4, fv("%s", symop))

class AngleRestraintsList (RestraintsListBase) :
  restraint_name = "angle"
  def CreateColumns (self) :
    self.InsertColumn(0, "Selections", width=600)
    self.InsertColumn(1, "Angle", width=100)
    self.InsertColumn(2, "Sigma", width=80)

  def PopulateList (self) :
    for angle in self._params :
      selections = [angle.atom_selection_1, angle.atom_selection_2,
                    angle.atom_selection_3]
      if (selections == [None, None]) :
        sele_str = "---"
      else :
        sele_str = "; ".join([ str(s) for s in selections ])
      item = self.InsertStringItem(sys.maxint, sele_str)
      self.SetStringItem(item, 1, fv("%g", angle.angle_ideal))
      self.SetStringItem(item, 2, fv("%g", angle.sigma))

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
      sele_text = "; ".join(selections)
      self.SetStringItem(item, 0, sele_text)

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
        self._params[item].angle_ideal = angle
        self.SetStringItem(item, 1, "%g" % angle)
      dlg.Destroy()

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
      dlg.Destroy()

class PlanarityRestraintsList (RestraintsListBase) :
  restraint_name = "planarity"
  def CreateColumns (self) :
    self.InsertColumn(0, "Atom selection", width=680)
    self.InsertColumn(1, "Sigma", width=100)

  def PopulateList (self) :
    for plane in self._params :
      item = self.InsertStringItem(sys.maxint, fv("%s", plane.atom_selection))
      self.SetStringItem(item, 1, fv("%g", plane.sigma))

  def GetSelections (self) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      plane = self._params[item]
      return plane.atom_selection
    return [""]

  def SetSelections (self, selections) :
    assert (len(selections) == 1)
    item = self.GetFirstSelected()
    if (item >= 0) :
      plane = self._params[item]
      plane.atom_selection = selections[0]
      self.SetStringItem(item, 0, fv("%s", plane.atom_selection))

  def OnSetSigma (self, event) :
    item = self.GetFirstSelected()
    if (item >= 0) :
      dlg = simple_dialogs.FloatDialog(
        parent=self,
        title="Restraint sigma",
        label="Plane sigma",
        value=self._params[item].sigma,
        caption="Please enter a sigma value for the selected planarity "+
          "restraint; the sigma for all planes defined by the CCP4 monomer "+
          "library is 0.02.")
      dlg.SetMin(0.01)
      dlg.SetOptional(False)
      if (dlg.ShowModal() == wx.ID_OK) :
        sigma = dlg.GetPhilValue()
        self._params[item].sigma = sigma
        self.SetStringItem(item, 1, "%g" % sigma)
      dlg.Destroy()

#-----------------------------------------------------------------------
# Frame class

# XXX fix for true division bug
flatnotebook.PageContainer.OnMouseWheel = lambda win, evt: False
class RestraintsFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    kwds['style'] = wx.DEFAULT_FRAME_STYLE
    wx.Frame.__init__(self, *args, **kwds)
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.SetSizer(self.sizer)
    tb = wx.ToolBar(self, style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    tb.SetToolBitmapSize((32,32))
    self.SetToolBar(tb)
    btn1 = tb.AddLabelTool(-1,
      label="Update and exit",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "button_ok"),
      kind=wx.ITEM_NORMAL)
    btn2 = tb.AddLabelTool(-1,
      label="Update settings",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "redo"),
      kind=wx.ITEM_NORMAL)
    btn3 = tb.AddLabelTool(-1,
      label="Cancel",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "cancel"),
      kind=wx.ITEM_NORMAL)
    btn4 = tb.AddLabelTool(-1,
      label="Delete all",
      bitmap=wxtbx.bitmaps.fetch_icon_bitmap("actions", "editdelete"),
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnExit, btn1)
    self.Bind(wx.EVT_MENU, self.OnUpdate, btn2)
    self.Bind(wx.EVT_MENU, self.OnCancel, btn3)
    self.Bind(wx.EVT_MENU, self.OnClearAll, btn4)
    tb.Realize()
    self.statusbar = self.CreateStatusBar()
    self.statusbar.SetStatusText("No restraints loaded.")
    self.nb = flatnotebook.FlatNotebook(
      parent=self,
      agwStyle=flatnotebook.FNB_TABS_BORDER_SIMPLE|flatnotebook.FNB_NODRAG|
        flatnotebook.FNB_NO_X_BUTTON|flatnotebook.FNB_NO_NAV_BUTTONS)
    self.bonds_panel = CustomBondPanel(self.nb, style=ATOM_SELECTION_BUTTONS)
    self.angles_panel = CustomAnglePanel(self.nb, style=ATOM_SELECTION_BUTTONS)
    self.planes_panel = CustomPlanarityPanel(self.nb,
      style=ATOM_SELECTION_BUTTONS)
    self.nb.AddPage(self.bonds_panel, "Bonds")
    self.nb.AddPage(self.angles_panel, "Angles")
    self.nb.AddPage(self.planes_panel, "Planes")
    self.sizer.Add(self.nb, 1, wx.EXPAND|wx.ALL)
    self.sizer.Fit(self.nb)
    self.Fit()
    self.Bind(wx.EVT_CLOSE, lambda evt: self.Destroy(), self)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self._index = self._params = self._prefix = None

  def SetPhilIndex (self, index) :
    self._index = index
    self.bonds_panel.SetPhilIndex(index)
    self.angles_panel.SetPhilIndex(index)
    self.planes_panel.SetPhilIndex(index)

  def SetParams (self, params, prefix=None) :
    self._params = params
    self._prefix = prefix
    self.bonds_panel.SetPhilPrefix(prefix)
    self.angles_panel.SetPhilPrefix(prefix)
    self.planes_panel.SetPhilPrefix(prefix)
    self.SetPanelParams()
    geo_params = self.GetGeoParams()
    self.statusbar.SetStatusText(
      "%d custom bond, %d custom angle, and %d custom plane restraints loaded."
      % (len(geo_params.bond),len(geo_params.angle),len(geo_params.planarity)))

  def SetPanelParams (self) :
    geo_params = self.GetGeoParams()
    self.bonds_panel.SetParams(geo_params.bond)
    self.angles_panel.SetParams(geo_params.angle)
    self.planes_panel.SetParams(geo_params.planarity)

  def GetGeoParams (self) :
    geo_params_parent = self._params
    if (self._prefix is not None) :
      geo_params_parent = getattr(geo_params_parent, self._prefix)
    return geo_params_parent.geometry_restraints.edits

  def GetPhilObject (self) :
    geo_params = self.GetGeoParams()
    geo_params.bond = self.bonds_panel.GetParams()
    geo_params.angle = self.angles_panel.GetParams()
    geo_params.planarity = self.planes_panel.GetParams()
    self.Validate()
    new_phil = self._index.master_phil.format(python_object=self._params)
    return new_phil

  def Validate (self) :
    from mmtbx.monomer_library import pdb_interpretation
    pdb_interpretation.validate_geometry_edits_params(self.GetGeoParams())

  def OnCancel (self, event) :
    self.Close()

  def OnUpdate (self, event) :
    self.GetPhilObject()

  def OnExit (self, event) :
    self.OnUpdate(event)
    self.Close()

  def OnClearAll (self, event) :
    confirm = wx.MessageBox("Are you sure you want to delete all custom "+
      "restraints?  This action cannot be undone from the restraints editor, "+
      " but clicking the Cancel "+
      "button on the toolbar will revert to the original parameters.")
    if (confirm == wx.YES) :
      geo_params = self.GetGeoParams()
      geo_params.bond = []
      geo_params.angle = []
      geo_params.planarity = []
      self.SetPanelParams()

  def OnDestroy (self, event) :
    pass

  def OnView (self, event) :
    pass

#-----------------------------------------------------------------------
# REGRESSION TESTING
if (__name__ == "__main__") :
  from mmtbx.monomer_library import pdb_interpretation
  import iotbx.phil
  import libtbx.phil.interface
  app = wx.App(0)
  master_phil = iotbx.phil.parse(pdb_interpretation.grand_master_phil_str,
    process_includes=True)
  index = libtbx.phil.interface.index(master_phil)
  frame = RestraintsFrame(None, -1, "Custom restraints editor")
  frame.SetPhilIndex(index)
  params = index.get_python_object()
  frame.SetParams(index.get_python_object())
  frame.Show()
  app.MainLoop()
