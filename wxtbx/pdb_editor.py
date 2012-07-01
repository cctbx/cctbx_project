
import wxtbx.app
from wxtbx.phil_controls import simple_dialogs
from wx.lib.agw import customtreectrl
import wx
import sys

def format_residue_group (rg) :
  return "residue %s (resseq='%s', icode='%s')" % \
    (rg.resid().strip(), rg.resseq, rg.icode)

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
          self.ShowAtomMenu()
        elif (pdb_type == "atom_group") :
          self.ShowAtomGroupMenu()
        elif (pdb_type == "residue_group") :
          self.ShowResidueGroupMenu()
        elif (pdb_type == "chain") :
          self.ShowChainMenu()
        elif (pdb_type == "model") :
          self.ShowModelMenu()
        else :
          raise RuntimeError("Unrecognized object type '%s'" % pdb_type)

  def ShowAtomMenu (self) :
    labels_and_actions = [
      ("Delete atom", self.OnDeleteObject),
      ("Set name...", self.OnSetName),
      ("Set occupancy...", self.OnSetOccupancy),
      ("Set B-factor...", self.OnSetBfactor),
      ("Set element...", self.OnSetElement),
      ("Set charge...", self.OnSetCharge),
    ]
    self.ShowMenu(labels_and_actions)

  def ShowAtomGroupMenu (self) :
    pass

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

  #---------------------------------------------------------------------
  # EDITING ACTIONS
  def OnDeleteObject (self, event) :
    pass

  def OnSetName (self, event) :
    pass

  def OnSetOccupancy (self, event) :
    pass

  def OnSetBfactor (self, event) :
    pass

  def OnSetElement (self, event) :
    pass

  def OnSetCharge (self, event) :
    pass

  #---------------------------------------------------------------------
  # USER INPUT
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
    if (charge is not None) :
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

########################################################################
class PDBTreeFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
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
    self._tree = PDBTree(self.panel, -1, style=wx.RAISED_BORDER)
    self._tree.SetMinSize((640,400))
    pszr.Add(self._tree, 1, wx.EXPAND, 2)
    self._pdb_in = None
    szr.Layout()
    self.Fit()

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

if (__name__ == "__main__") :
  app = wxtbx.app.CCTBXApp(0)
  frame = PDBTreeFrame(None, -1, "PDB Editor")
  frame.Show()
  frame.LoadPDB(sys.argv[1])
  app.MainLoop()
