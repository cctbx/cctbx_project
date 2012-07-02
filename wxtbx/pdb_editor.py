
import wxtbx.app
import wxtbx.bitmaps
from wxtbx.phil_controls import simple_dialogs
from wxtbx import path_dialogs
from libtbx.utils import Abort
from wx.lib.agw import customtreectrl
import wx
import sys

def format_residue_group (rg) :
  return "residue %s (resseq='%s', icode='%s') [%s]" % \
    (rg.resid().strip(), rg.resseq, rg.icode,
      ",".join([ ag.resname for ag in rg.atom_groups() ]))

def format_atom_group (ag) :
  return "atom group %s (altloc='%s')" % (ag.resname, ag.altloc)

def format_atom (atom) :
  return "atom '%s' (x=%.3f, y=%.3f, z=%.3f, b=%.2f, occ=%.2f, elem='%s', charge='%s')" % (atom.name, atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.b, atom.occ,
           atom.element, atom.charge)
  #return "atom '%s'" % atom.name

class PDBTree (customtreectrl.CustomTreeCtrl) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    kwds['agwStyle'] = wx.TR_HAS_VARIABLE_ROW_HEIGHT|wx.TR_HAS_BUTTONS| \
      wx.TR_TWIST_BUTTONS|wx.TR_HIDE_ROOT|wx.TR_SINGLE
    customtreectrl.CustomTreeCtrl.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_TREE_KEY_DOWN, self.OnChar)
    self.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnRightClick)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick)
    self.DeleteAllItems()
    self.flag_multiple_selections = False

  def DeleteAllItems (self) :
    customtreectrl.CustomTreeCtrl.DeleteAllItems(self)
    self.AddRoot("pdb_hierarchy")
    self._node_lookup = {}
    self._state = []
    self._state_tmp = []
    self._changes_made = False

  def SetHierarchy (self, pdb_hierarchy) :
    self.DeleteAllItems()
    self._state.append(pdb_hierarchy)
    self.PopulateTree(pdb_hierarchy)

  def PopulateTree (self, pdb_hierarchy) :
    root_node = self.GetRootItem()
    for model in pdb_hierarchy.models() :
      model_node = self.AppendItem(root_node, "model '%s'" % model.id)
      self._node_lookup[model_node] = model
      for chain in model.chains() :
        chain_node = self.AppendItem(model_node, "chain '%s'" % chain.id)
        self._node_lookup[chain_node] = chain
        for rg in chain.residue_groups() :
          rg_node = self.AppendItem(chain_node, format_residue_group(rg))
          self._node_lookup[rg_node] = rg
          for ag in rg.atom_groups() :
            ag_node = self.AppendItem(rg_node, format_atom_group(ag))
            self._node_lookup[ag_node] = ag
            for atom in ag.atoms() :
              atom_node = self.AppendItem(ag_node, format_atom(atom))
              self._node_lookup[atom_node] = atom

  def PropagateAtomChanges (self, node) :
    child, cookie = self.GetFirstChild(node)
    if (child is None) :
      pdb_object = self._node_lookup.get(node, None)
      if (type(pdb_object).__name__ == 'atom') :
        self.SetItemText(node, format_atom(pdb_object))
      else :
        print "Can't modify object %s'" % type(pdb_object).__name__
    else :
      while (child is not None) :
        self.PropagateAtomChanges(child)
        child, cookie = self.GetNextChild(node, cookie)

  def HaveUnsavedChanges (self) :
    return self._changes_made

  def SaveChanges (self) :
    self._changes_made = False

  def OnChar (self, event) :
    pass

  def OnRightClick (self, event) :
    self.ActionsForSelection()

  def OnDoubleClick (self, event) :
    pass
  #  self.ActionsForSelection()

  def ActionsForSelection (self) :
    if (self.flag_multiple_selections) :
      pass
    else :
      sel = self.GetSelection()
      pdb_object = self._node_lookup.get(sel, None)
      if (pdb_object is not None) :
        pdb_type = type(pdb_object).__name__
        if (pdb_type == "atom") :
          self.ShowAtomMenu(pdb_object)
        elif (pdb_type == "atom_group") :
          self.ShowAtomGroupMenu(pdb_object)
        elif (pdb_type == "residue_group") :
          self.ShowResidueGroupMenu()
        elif (pdb_type == "chain") :
          self.ShowChainMenu()
        elif (pdb_type == "model") :
          self.ShowModelMenu()
        else :
          raise RuntimeError("Unrecognized object type '%s'" % pdb_type)

  def ShowAtomMenu (self, atom) :
    labels_and_actions = [
      ("Delete atom", self.OnDeleteObject),
      ("Set name...", self.OnSetName),
      ("Set occupancy...", self.OnSetOccupancy),
      ("Set B-factor...", self.OnSetBfactor),
      ("Set element...", self.OnSetElement),
      ("Set charge...", self.OnSetCharge),
    ]
    self.ShowMenu(labels_and_actions)

  def ShowAtomGroupMenu (self, atom_group) :
    labels_and_actions = [
      ("Delete atom group", self.OnDeleteObject),
      ("Set altloc...", self.OnSetAltloc),
      ("Set occupancy...", self.OnSetOccupancy),
      ("Set B-factor...", self.OnSetBfactor),
      ("Set residue name...", self.OnSetResname),
    ]
    if (atom_group.resname == "MET") :
      labels_and_actions.append(("Convert to SeMet", self.OnConvertMet))
    elif (atom_group.resname == "MSE") :
      labels_and_actions.append(("Convert to Met", self.OnConvertSeMet))
    self.ShowMenu(labels_and_actions)

  def ShowResidueGroupMenu (self) :
    pass

  def ShowChainMenu (self) :
    pass

  def ShowModelMenu (self) :
    pass

  def ShowMenu (self, items) :
    menu = wx.Menu()
    for label, action in items :
      item = menu.Append(-1, label)
      self.Bind(wx.EVT_MENU, action, item)
    self.PopupMenu(menu)
    menu.Destroy()

  def GetSelectedObject (self, object_type=None) :
    item = self.GetSelection()
    pdb_object = self._node_lookup.get(item, None)
    if (object_type is not None) :
      assert (type(pdb_object).__name__ == object_type)
    return item, pdb_object

  #---------------------------------------------------------------------
  # EDITING ACTIONS
  def OnDeleteObject (self, event) :
    pass

  # atom
  def OnSetName (self, event) :
    pass

  # atom
  def OnSetElement (self, event) :
    item, atom = self.GetSelectedObject('atom')
    new_elem = self.GetNewElement(atom.element)
    assert (new_elem is None) or (len(new_elem) <= 2)
    self._changes_made = True
    if (new_elem is None) :
      atom.element = '  '
    elif (new_elem.strip() == '') :
      new_elem = new_elem.strip()
      atom.element = '%2d' % new_elem
    # TODO validate element symbol
    self.SetItemText(item, format_atom(atom))

  # atom
  def OnSetCharge (self, event) :
    item, atom = self.GetSelectedObject('atom')
    new_charge = self.GetNewCharge(atom.charge)
    assert (new_charge is None) or (-9 <= new_charge <= 9)
    self._changes_made = True
    if (new_charge in [0, None]) :
      atom.charge = '  '
    elif (new_charge < 0) :
      atom.charge = '%d-' % abs(new_charge)
    else :
      atom.charge = '%d+' % new_charge
    self.SetItemText(item, format_atom(atom))

  # atom_group
  def OnSetAltloc (self, event) :
    item, atom_group = self.GetSelectedObject('atom_group')
    new_altloc = self.GetNewAltloc(atom_group.altloc)
    assert (new_altloc is None) or (len(new_altloc) in [0,1])
    if (new_altloc in [None, '']) :
      rg = atom_group.parent()
      if (len(rg.atom_groups()) > 1) :
        self.Confirm("You have specified a blank altloc ID for this atom "+
          "group, but it is part of a residue containing multiple "+
          "conformations.  Are you sure this is what you want to do?")
      atom_group.altloc = ''
    else :
      atom_group.altloc = new_altloc
    self._changes_made = True
    self.SetItemText(item, format_atom_group(atom_group))

  # atom_group
  def OnSetResname (self, event) :
    item, atom_group = self.GetSelectedObject('atom_group')
    new_resname = self.GetNewResname(atom_group.resname)
    assert (new_resname is not None) and (len(new_resname) in [1,2,3])
    self._changes_made = True
    atom_group.resname = new_resname
    self.SetItemText(item, format_atom_group(atom_group))
    rg_item = self.GetItemParent(item)
    self.SetItemText(rg_item, format_residue_group(atom_group.parent()))

  # atom_group (resname == MET)
  def OnConvertMet (self, event) :
    item, atom_group = self.GetSelectedObject('atom_group')
    assert (atom_group.resname == "MET")
    self._changes_made = True
    atom_group.resname = "MSE"
    for atom in atom_group.atoms() :
      if (atom.name == ' SD ') :
        atom.name = ' SE '
        atom.element = 'Se'
        break
    self.SetItemText(item, format_atom_group(atom_group))
    rg_item = self.GetItemParent(item)
    self.SetItemText(rg_item, format_residue_group(atom_group.parent()))
    self.PropagateAtomChanges(item)

  # atom_group (resname == MSE)
  def OnConvertSeMet (self, event) :
    item, atom_group = self.GetSelectedObject('atom_group')
    assert (atom_group.resname == "MSE")
    self._changes_made = True
    atom_group.resname = "MET"
    for atom in atom_group.atoms() :
      if (atom.name == ' SE ') :
        atom.name = ' SD '
        atom.element = ' S'
        break
    self.SetItemText(item, format_atom_group(atom_group))
    rg_item = self.GetItemParent(item)
    self.SetItemText(rg_item, format_residue_group(atom_group.parent()))
    self.PropagateAtomChanges(item)

  # residue_group
  def OnSetResseq (self, event) :
    pass

  # residue_group
  def OnSetIcode (self, event) :
    item, residue_group = self.GetSelectedObject('residue_group')
    new_icode = self.GetNewIcode(residue_group.icode)
    assert (new_icode is None) or (len(new_icode) == 1)
    self._changes_made = True

  # chain
  def OnSetChainID (self, event) :
    pass

  # model
  def OnSetModelID (self, event) :
    pass

  # all
  def OnSetOccupancy (self, event) :
    if (self.flag_multiple_selections) :
      assert 0
    else :
      item = self.GetSelection()
      pdb_object = self._node_lookup.get(item, None)
      pdb_type = type(pdb_object).__name__
      occ = None
      if (pdb_type == 'atom') :
        occ = pdb_object.occ
      else :
        all_occ = pdb_object.atoms().extract_occ()
        if (all_occ.all_eq(all_occ[0])) :
          occ = all_occ[0]
      new_occ = self.GetNewOccupancy(occ)
      print new_occ
      assert (0 <= occ <= 1.0)
      self._changes_made = True
      if (pdb_type == 'atom') :
        pdb_object.occ = new_occ
        self.SetItemText(item, format_atom(pdb_object))
      else :
        for atom in pdb_object.atoms() :
          atom.occ = new_occ
        self.PropagateAtomChanges(item)

  # all
  def OnSetBfactor (self, event) :
    if (self.flag_multiple_selections) :
      assert 0
    else :
      item = self.GetSelection()
      pdb_object = self._node_lookup.get(item, None)
      pdb_type = type(pdb_object).__name__
      b_iso = None
      if (pdb_type == 'atom') :
        b_iso = pdb_object.b
      else :
        all_b = pdb_object.atoms().extract_b()
        if (all_b.all_eq(all_b[0])) :
          b_iso = all_b[0]
      new_b = self.GetNewBiso(b_iso)
      assert (0 < new_b < 1000)
      self._changes_made = True
      if (pdb_type == 'atom') :
        pdb_object.b = new_b
        self.SetItemText(item, format_atom(pdb_object))
      else :
        for atom in pdb_object.atoms() :
          atom.b = new_b
        self.PropagateAtomChanges(item)

  #---------------------------------------------------------------------
  # USER INPUT
  def Confirm (self, msg) :
    confirm = wx.MessageBox(
      message=msg,
      style=wx.YES_NO)
    if (confirm == wx.NO) :
      raise Abort()
    return True

  def GetNewOccupancy (self, occ=None) :
    dlg = simple_dialogs.FloatDialog(
      parent=self,
      title="Set new occupancy",
      label="New occupancy",
      caption="Please specify the occupancy for the selected atom(s); this "+
        "value represents the fraction of unit cells in which the atom(s) "+
        "is/are present, and must be a value between 0 (no contribution to "+
        "F_calc) and 1.0; note that atoms with multiple conformers should "+
        "always have a sum of occupancies of 1.0.  The precision will be "+
        "truncated to two digits after the decimal point.",
      value=occ)
    dlg.SetMin(0.0)
    dlg.SetMax(1.0)
    dlg.SetOptional(False)
    new_occ = simple_dialogs.get_phil_value_from_dialog(dlg)
    return new_occ

  def GetNewBiso (self, b=None) :
    dlg = simple_dialogs.FloatDialog(
      parent=self,
      title="Set new isotropic B-factor",
      label="New B_iso",
      caption="Please specify the isotropic B-factor for the selected "+
        "atom(s).  This should be a value between 0.01 and 999.9; if you "+
        "are not sure what value to use but will be performing refinement "+
        "on the modified PDB file, 20 is a good guess.  The precision will "+
        "be truncated to two digits after the decimal point.",
      value=b)
    dlg.SetMin(0.01)
    dlg.SetMax(999.99)
    dlg.SetOptional(False)
    new_b = simple_dialogs.get_phil_value_from_dialog(dlg)
    return new_b

  def GetNewCharge (self, charge=None) :
    if (charge.isspace()) :
      charge = None
    elif (charge is not None) :
      charge = int(charge)
    dlg = simple_dialogs.IntegerDialog(
      parent=self,
      title="Set new atomic charge",
      label="New charge",
      caption="Please specify the atomic charge; this should be a number "+
        "between -9 and 9, but usually between -2 and 2.  If you want to "+
        "set the charge to zero you can simply leave the field blank.",
      value=charge)
    dlg.SetMin(-9)
    dlg.SetMax(9)
    dlg.SetOptional(True)
    new_charge = simple_dialogs.get_phil_value_from_dialog(dlg)
    return new_charge

  def GetNewElement (self, elem=None) :
    if (elem.isspace()) :
      elem = None
    dlg = simple_dialogs.StringDialog(
      parent=self,
      title="Set atomic element",
      label="Element symbol",
      caption="The atomic element field is optional if it can be guessed "+
        "based on the atom name, but it is best to specify explicitly.  It "+
        "must be one or two characters in length.",
      value=elem)
    dlg.SetMaxLength(2)
    dlg.SetOptional(True)
    new_elem = simple_dialogs.get_phil_value_from_dialog(dlg)
    return new_elem

  def GetNewIcode (self, icode=None) :
    if (icode.isspace()) :
      icode = None
    dlg = simple_dialogs.StringDialog(
      parent=self,
      title="Set residue insertion code",
      label="Insertion code",
      caption="The insertion code is a single character used to disambiguate "+
        "between multiple residues with identical numbering, such as "+
        "proteins where specific numbering is desired which may not match "+
        "the linear sequence (antibodies, proteases, proteins with "+
        "engineered insertions, etc.).",
      value=icode)
    dlg.SetMaxLength(1)
    dlg.SetOptional(True)
    new_icode = simple_dialogs.get_phil_value_from_dialog(dlg)
    return new_icode

  def GetNewChainID (self, chain_id=None) :
    if (chain_id.isspace()) :
      chain_id = None
    dlg = simple_dialogs.StringDialog(
      parent=self,
      title="Set chain ID",
      label="New ID",
      caption="Please specify a chain ID.  If you have only a single chain "+
        "in your structure the chain ID may be left blank, but we recommend "+
        "labeling it 'A' instead.  Otherwise the chain ID may be up to two "+
        "characters in Phenix, Coot, and CCP4, but the official limit for "+
        "the format is one character.",
      chain_id=chain_id)
    #dlg.SetMinLength(1)
    dlg.SetMaxLength(2)
    new_id = simple_dialogs.get_phil_value_from_dialog(dlg)
    return new_id

########################################################################
class PDBTreeFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.statusbar = self.CreateStatusBar()
    self.statusbar.SetStatusText("Right-click any item for a list of editing actions")
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    bmp = wxtbx.bitmaps.fetch_custom_icon_bitmap("phenix.pdbtools")
    btn = self.toolbar.AddLabelTool(-1, "Load file", bmp,
      shortHelp="Load file", kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnOpen, btn)
    self.toolbar.Realize()
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self.panel = wx.Panel(self, -1)
    szr.Add(self.panel, 1, wx.EXPAND)
    pszr = wx.BoxSizer(wx.VERTICAL)
    self.panel.SetSizer(pszr)
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    szr2.Add(wx.StaticText(self.panel, -1, "PDB file:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._path_field = wx.TextCtrl(self.panel, -1, size=(540,-1))
    szr2.Add(self._path_field, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    pszr.Add(szr2, 0, wx.EXPAND)
    szr3 = wx.BoxSizer(wx.HORIZONTAL)
    self._header_box = wx.CheckBox(self.panel, -1,
      "Preserve header records when saving file")
    self._header_box.SetValue(True)
    szr3.Add(self._header_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    pszr.Add(szr3, 0, wx.EXPAND)
    self._tree = PDBTree(self.panel, -1, style=wx.RAISED_BORDER)
    self._tree.SetMinSize((640,400))
    pszr.Add(self._tree, 1, wx.EXPAND, 2)
    self._pdb_in = None
    szr.Layout()
    self.Fit()
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self.path_mgr = path_dialogs.manager()

  def LoadPDB (self, file_name) :
    from iotbx import file_reader
    f = file_reader.any_file(file_name,
      force_type="pdb",
      raise_sorry_if_errors=True)
    self._pdb_in = f
    hierarchy = f.file_object.construct_hierarchy()
    atoms = hierarchy.atoms()
    atoms.set_chemical_element_simple_if_necessary()
    self._tree.SetHierarchy(hierarchy)
    self._path_field.SetValue(f.file_name)
    self.Refresh()

  def OnOpen (self, event) :
    from iotbx import file_reader
    file_name = self.path_mgr.select_file(
      parent=self,
      message="Select PDB file",
      wildcard=file_reader.get_wildcard_strings(["pdb"]))
    self.LoadPDB(file_name)

  def OnSave (self, event) :
    save_header = self._header_box.GetValue()

  def OnClose (self, event) :
    if (self._tree.HaveUnsavedChanges()) :
      pass
    self.Destroy()

  def OnDestroy (self, event) :
    pass

if (__name__ == "__main__") :
  app = wxtbx.app.CCTBXApp(0)
  frame = PDBTreeFrame(None, -1, "PDB Editor")
  frame.Show()
  if (len(sys.argv) > 1) :
    frame.LoadPDB(sys.argv[1])
  app.MainLoop()
