
from __future__ import absolute_import, division, print_function
import wxtbx.bitmaps
from wxtbx.phil_controls import simple_dialogs
from wxtbx.phil_controls import strctrl, intctrl
from wxtbx import symmetry_dialog
from wxtbx import path_dialogs
from wx.lib.embeddedimage import PyEmbeddedImage
from wx.lib.agw import customtreectrl
import wx
from libtbx.utils import Abort, Sorry, to_unicode, to_str
import string
import os

b_iso_icon = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAADnRFWHRTb2Z0d2FyZQBQeU1P"
    "TPa/er0AAAAYdEVYdFVSTABodHRwOi8vd3d3LnB5bW9sLm9yZ5iPN04AAAHrSURBVDiNxZLP"
    "TlNBFMa/mTNzL4XeIt6SCIKmRBKXpTSRujDxDdy4dmMXJr6BCzcufAX/4BPwFK4sxMREEzeG"
    "KpoiNo2tFHrnzp05LtqiWOKWb3Uy853fmZz5gPOW+Pdgfv5hOQiW6lovVaXMVQDJRNGmMc2d"
    "vb36q/8CFhYeV7S+/CwIrla1vgSiCEIoMDs412VjvjR2d+/UzgTE8b21XG7t+dTU9WoQLEOp"
    "GER5EGkADGuPYe0PHgzev2k2794c98lxQTRXVyquKnURRLMgilAoTGN1NUCpFGJmZhZaF0UY"
    "XrmxsrJVmQAAqEuZhxAhhNAgUigWBTY2gFoNiGMJohyILkhA3h83KQDI52+VpZyWgIAQAgCD"
    "mZEkQLsNCAFk2Z+ZQ+9fAADw3oA5hfcJmFM4Z9DpaGxvD729noNzQw9zilOAfv/1uzC85p3r"
    "Sef6kLIPIRSShJGmwWhACucO4VwX1u77M15w9MLa/QdEhdHXeRAZCKFH9wmc+4k0/eq9P355"
    "svxxEQTLLcCvC6EWh3vwYDbwfgDnfiHLOrD2Gw8GHxqt1qMnEwBjPn0Pw9Jb7w/XvTeLw30c"
    "wbkesqwNaz9zknxsHBw8PcnAqSCNFUW3y1JGdaK5qpRhBRBMVNg0prnT7W5NRPn89RuLodCO"
    "WgQHYQAAAABJRU5ErkJggg==")

b_aniso_icon = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAADnRFWHRTb2Z0d2FyZQBQeU1P"
    "TPa/er0AAAAYdEVYdFVSTABodHRwOi8vd3d3LnB5bW9sLm9yZ5iPN04AAAFnSURBVDiNzdLP"
    "ahNRHMXxz0xCtCOaYkVpFLLpSnHTCdLH6MJHsNZtt6661wew6COIO5/CjAbciBRBkSYtBNv0"
    "nzOxuS4yDbEVXepZ/S6X7zmH3738a0V/unxKeoGVeaIaCtr3ef5Xg8ekV9lokt7ELGo4xjbZ"
    "e1bWeAeVs/A6aZ03d2ks4Baux7F6FElCENE4JrzkNVSn4bVx5fYdNHEDV+p1cbNJv0+36+Jo"
    "JGbxlIlPh3ssJmwsYB5zuJwk4qUlWi2iyDAEOYZToZMGt2ldI51DHQkqMzPs7rK5Ke/1DELw"
    "Dd/pnDOIeDBbglUEnOzv0+nIi8JeCHbwFV9onzNAq1KCQxyhKAo/yrlfwp/J1qeecmIQ8/aI"
    "9BB7yMet5OV5Gx84OWB1evETgwFZlzQpoRpGOCjTd8gGrD4hmzaY/IMGW0NaOY3DEuriEz7S"
    "7vHo2Rm43N2vWubhJdLqOH0UePHqN+D/o5/JPnJ7ZjEAWAAAAABJRU5ErkJggg==")

def format_model(model):
  return "model '%s'" % model.id

def format_chain(chain):
  return "chain '%s'" % chain.id

def format_residue_group(rg):
  return "residue %s (resseq='%s', icode='%s') [%s]" % \
    (rg.resid().strip(), rg.resseq, rg.icode,
      ",".join([ ag.resname for ag in rg.atom_groups() ]))

def format_atom_group(ag):
  return "atom group %s (altloc='%s')" % (ag.resname, ag.altloc)

def format_atom(atom):
  if (atom.hetero):
    prefix = "hetatm"
  else :
    prefix = "atom"
  return "%s '%s' (x=%.3f, y=%.3f, z=%.3f, b=%.2f, occ=%.2f, elem='%s', charge='%s')" % (prefix, atom.name, atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.b,
      atom.occ, atom.element, atom.charge)
  #return "atom '%s'" % atom.name

def remove_objects_recursive(pdb_object):
  parent = pdb_object.parent()
  pdb_type = type(pdb_object).__name__
  if (pdb_type == 'atom'):
    parent.remove_atom(pdb_object)
    if (len(parent.atoms()) == 0):
      remove_objects_recursive(parent)
  elif (pdb_type == 'atom_group'):
    parent.remove_atom_group(pdb_object)
    if (len(parent.atom_groups()) == 0):
      remove_objects_recursive(parent)
  elif (pdb_type == 'residue_group'):
    parent.remove_residue_group(pdb_object)
    if (len(parent.residue_groups()) == 0):
      remove_objects_recursive(parent)
  elif (pdb_type == 'chain'):
    parent.remove_chain(pdb_object)
    if (len(parent.chains()) == 0):
      remove_objects_recursive(parent)
  else :
    if (len(parent.models()) == 1):
      raise Sorry("You must have at least one MODEL in the PDB file.")
    parent.remove_model(pdb_object)

def has_aniso_atoms(pdb_object):
  pdb_type = type(pdb_object).__name__
  if (pdb_type == 'atom'):
    return (pdb_object.uij != (-1,-1,-1,-1,-1,-1))
  else :
    for atom in pdb_object.atoms():
      if (atom.uij != (-1,-1,-1,-1,-1,-1)):
        return True
  return False

class PDBTree(customtreectrl.CustomTreeCtrl):
  max_states = 5 # maximum number of reverts possible
  def __init__(self, *args, **kwds):
    kwds = dict(kwds)
    kwds['agwStyle'] = wx.TR_HAS_VARIABLE_ROW_HEIGHT|wx.TR_HAS_BUTTONS| \
      wx.TR_TWIST_BUTTONS|wx.TR_HIDE_ROOT|wx.TR_MULTIPLE
    customtreectrl.CustomTreeCtrl.__init__(self, *args, **kwds)
    self.il = wx.ImageList(16,16)
    self.il.Add(b_iso_icon.GetBitmap())
    self.il.Add(b_aniso_icon.GetBitmap())
    self.SetImageList(self.il)
    self.Bind(wx.EVT_TREE_KEY_DOWN, self.OnChar)
    self.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnRightClick)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick)
    self.frame = self.GetTopLevelParent()
    # FIXME
    #self.Bind(wx.EVT_TREE_BEGIN_DRAG, self.OnStartDrag)
    #self.Bind(wx.EVT_TREE_END_DRAG, self.OnEndDrag)
    self.DeleteAllItems()
    self.path_mgr = path_dialogs.manager()
    self.dt = PDBTreeDropTarget(self)
    self.SetDropTarget(self.dt)

  def DeleteAllItems(self):
    customtreectrl.CustomTreeCtrl.DeleteAllItems(self)
    self.AddRoot("pdb_hierarchy")
    self._hierarchy_stack = []
    self._hierarchy_actions = []
    self._changes_made = False
    self._hierarchy = None
    self._hierarchy_start = None
    self._i_state = 0
    self._crystal_symmetry = None
    self._space_group = None
    self._unit_cell = None

  def PushState(self, action="edit", set_changes_flag=True):
    del self._hierarchy_stack[self._i_state+1:]
    del self._hierarchy_actions[self._i_state+1:]
    self._hierarchy_stack.append(self._hierarchy.deep_copy())
    self._hierarchy_actions.append(action)
    if (len(self._hierarchy_stack) > self.max_states):
      del self._hierarchy_stack[0]
      del self._hierarchy_actions[0]
    if (set_changes_flag):
      self._changes_made = True
    self._i_state = len(self._hierarchy_stack) - 1
    if wx.VERSION < (4,0):
      self.frame.EnableUndo(True)
    #print self._hierarchy_stack
    #print self._hierarchy_actions
    #print self._i_state

  def SetHierarchy(self, pdb_hierarchy):
    self.DeleteAllItems()
    self._hierarchy = pdb_hierarchy
    self._hierarchy_start = pdb_hierarchy
    self.PushState("starting model", set_changes_flag=False)
    self.PopulateTree(pdb_hierarchy)

  def SetCrystalSymmetry(self, crystal_symmetry):
    self._crystal_symmetry = crystal_symmetry
    if (crystal_symmetry is not None):
      self._unit_cell = crystal_symmetry.unit_cell()
      self._space_group = crystal_symmetry.space_group()
    else :
      self._unit_cell = self._space_group = None

  def SetState(self):
    #print "SetState(): n_states:", len(self._hierarchy_stack)
    #print "SetState(): i_state:", self._i_state
    #print "SetState(): action:", self._hierarchy_actions[self._i_state]
    customtreectrl.CustomTreeCtrl.DeleteAllItems(self)
    self.AddRoot("pdb_hierarchy")
    self._hierarchy = self._hierarchy_stack[self._i_state]
    self.PopulateTree(self._hierarchy)
    self.Refresh()

  def PopulateTree(self, pdb_hierarchy):
    root_node = self.GetRootItem()
    for model in pdb_hierarchy.models():
      self._InsertModelItem(root_node, model)

  def _InsertModelItem(self, root_node, model, index=None):
    model_node = None
    if (index is None):
      model_node = self.AppendItem(root_node, format_model(model), data=model)
    else :
      model_node = self.InsertItemByIndex(root_node, index,
        format_model(model), data=model)
    for chain in model.chains():
      self._InsertChainItem(model_node, chain)
    self.Expand(model_node)
    return model_node

  def _InsertChainItem(self, model_node, chain, index=None):
    chain_node = None
    if (index is None):
      chain_node = self.AppendItem(model_node, format_chain(chain), data=chain)
    else :
      chain_node = self.InsertItemByIndex(model_node, index,
        format_chain(chain), data=chain)
    for rg in chain.residue_groups():
      self._InsertResidueGroupItem(chain_node, rg)
    return chain_node

  def _InsertResidueGroupItem(self, chain_node, residue_group, index=None):
    rg_node = None
    if (index is None):
      rg_node = self.AppendItem(chain_node, format_residue_group(residue_group),
        data=residue_group)
    else :
      print("index =", index)
      rg_node = self.InsertItemByIndex(chain_node, index,
        format_residue_group(residue_group), data=residue_group)
    for ag in residue_group.atom_groups():
      self._InsertAtomGroupItem(rg_node, ag)
    return rg_node

  def _InsertAtomGroupItem(self, rg_node, atom_group, index=None):
    ag_node = None
    if (index is None):
      ag_node = self.AppendItem(rg_node, format_atom_group(atom_group),
        data=atom_group)
    else :
      ag_node = self.InsertItemByIndex(rg_node, index,
        format_atom_group(atom_group), data=atom_group)
    for atom in atom_group.atoms():
      self._InsertAtomItem(ag_node, atom)
    return ag_node

  def _InsertAtomItem(self, ag_node, atom, index=None):
    atom_node = None
    if (index is None):
      atom_node = self.AppendItem(ag_node, format_atom(atom), data=atom)
    else :
      atom_node = self.InsertItemByIndex(ag_node, index, format_atom(atom),
        data=atom)
    if (atom.uij == (-1,-1,-1,-1,-1,-1)):
      self.SetItemImage(atom_node, 0)
    else :
      self.SetItemImage(atom_node, 1)
    return atom_node

  def PropagateAtomChanges(self, node):
    child, cookie = self.GetFirstChild(node)
    if (child is None):
      pdb_object = self.GetItemPyData(node)
      if (type(pdb_object).__name__ == 'atom'):
        self.SetItemText(node, format_atom(pdb_object))
        if (pdb_object.uij == (-1,-1,-1,-1,-1,-1)):
          self.SetItemImage(node, 0)
        else :
          self.SetItemImage(node, 1)
      else :
        print("Can't modify object %s'" % type(pdb_object).__name__)
    else :
      while (child is not None):
        self.PropagateAtomChanges(child)
        child, cookie = self.GetNextChild(node, cookie)

  def FindItem(self, pdb_object):
    root_node = self.GetRootItem()
    return self._FindItem(root_node, pdb_object)

  def _FindItem(self, node, pdb_object):
    child, cookie = self.GetFirstChild(node)
    if (child is None):
      if (self.GetItemPyData(node) == pdb_object):
        return child
      return None
    while (child is not None):
      if (self.GetItemPyData(child) == pdb_object):
        return child
      item = self._FindItem(child, pdb_object)
      if (item is not None):
        return item
      child, cookie = self.GetNextChild(node, cookie)
    return None

  def _ApplyToAtoms(self, item, pdb_object, modify_action):
    assert (hasattr(modify_action, "__call__"))
    if (type(pdb_object).__name__ == 'atom'):
      modify_action(pdb_object)
      self.SetItemText(item, format_atom(pdb_object))
    else :
      for atom in pdb_object.atoms():
        modify_action(atom)
      self.PropagateAtomChanges(item)

  def GetSelectedObject(self, object_type=None):
    item = self.GetSelection()
    pdb_object = self.GetItemPyData(item)
    if (object_type is not None):
      assert (type(pdb_object).__name__ == object_type)
    return item, pdb_object

  def GetSelectionInfo(self):
    from scitbx.array_family import flex
    selection = flex.bool(self._hierarchy.atoms_size(), False)
    object_types = set([])
    for item in self.GetSelections():
      pdb_object = self.GetItemPyData(item)
      pdb_type = type(pdb_object).__name__
      object_types.add(pdb_type)
      if (pdb_type == 'atom'):
        selection[pdb_object.i_seq] = True
      else :
        for atom in pdb_object.atoms():
          selection[atom.i_seq] = True
    return list(object_types), selection.count(True)

  def FindInsertionIndex(self, pdb_object):
    pass

  def HaveUnsavedChanges(self):
    return self._changes_made

  def SaveChanges(self):
    self._changes_made = False

  def DeselectAll(self):
    items = self.GetSelections()
    for item in items :
      self.SelectItem(item, False)

  def Undo(self):
    if (self._i_state == 0):
      raise Sorry(("No more changes left to undo (only the previous %d "+
        "states will be remembered)") % self.max_states)
    else :
      print("Undoing action: %s" % self._hierarchy_actions[self._i_state])
      self._i_state -= 1
      self.SetState()
      if (self._i_state == 0):
        self.frame.EnableUndo(False)
      self.frame.EnableRedo(True)

  def Redo(self):
    if (self._i_state == len(self._hierarchy_stack) - 1):
      raise Sorry("No more changes to redo.")
    else :
      self._i_state += 1
      print("Restoring action: %s" % self._hierarchy_actions[self._i_state])
      self.SetState()
      self.frame.EnableUndo(True)
      if (self._i_state == len(self._hierarchy_stack) - 1):
        self.frame.EnableRedo(False)

  #---------------------------------------------------------------------
  # UI events
  def OnChar(self, event):
    evt2 = event.GetKeyEvent()
    key = evt2.GetKeyCode()
    shift_down = evt2.ShiftDown()
    all_items = self.GetSelections()
    last_item = None
    if (len(all_items) > 0):
      last_item = all_items[-1]
    if (last_item is None):
      first_item, cookie = self.GetFirstChild(self.GetRootItem())
      if (first_item is not None):
        self.SelectItem(first_item)
    elif (key == wx.WXK_SPACE) or (key == wx.WXK_TAB):
      for item in all_items :
        if (self.IsExpanded(item)):
          self.Collapse(item)
        else :
          self.Expand(item)
    elif (key == wx.WXK_LEFT):
      parent = self.GetItemParent(last_item)
      if (parent != self.GetRootItem()):
        if (not shift_down):
            self.DeselectAll()
        self.SelectItem(parent)
    elif (key == wx.WXK_RIGHT):
      first_child, cookie = self.GetFirstChild(last_item)
      if (first_child is not None):
        if (not shift_down):
          self.DeselectAll()
        self.SelectItem(first_child)
    elif (key == wx.WXK_DOWN):
      next_item = self.GetNextSibling(last_item)
      if (next_item is None):
        parent = self.GetItemParent(last_item)
        first_sibling, cookie = self.GetFirstChild(parent)
        # if the current node has no siblings, select the first child instead
        if (first_sibling == last_item):
          next_item, cookie = self.GetFirstChild(last_item)
        else :
          next_item = first_sibling
      if (next_item is not None):
        if (not shift_down):
          self.DeselectAll()
        self.SelectItem(next_item)
    elif (key == wx.WXK_UP):
      first_item = all_items[0]
      prev_item = self.GetPrevSibling(first_item)
      if (prev_item is None):
        parent = self.GetItemParent(first_item)
        prev_item = self.GetLastChild(parent)
      if (prev_item is not None):
        if (not shift_down):
          self.DeselectAll()
        self.SelectItem(prev_item)
    elif (key == wx.WXK_DELETE) or (key == wx.WXK_BACK):
      self.DeleteSelected()
    elif (key == wx.WXK_RETURN):
      self.ActionsForSelection(source_window=self)

  def OnRightClick(self, event):
    self.ActionsForSelection(source_window=self)

  # XXX is this useful or not?
  def OnDoubleClick(self, event):
    pass

  # TODO
  def OnStartDrag(self, event):
    items = self.GetSelections()
    object_types, n_atoms = self.GetSelectionInfo()
    if (len(object_types) == 1):
      event.Allow()
      data = wx.CustomDataObject(wx.CustomDataFormat('pdb_object'))
      data.SetData(object_types[0])
      drop_source = wx.DropSource(self)
      drop_source.SetData(data)
      result = drop_source.DoDragDrop(flags=wx.Drag_DefaultMove)
      print(result)
      self.Refresh()

  def OnDrop(self, x, y):
    print(x, y)
    return True

  def OnEndDrag(self, event):
    print(1)

  #---------------------------------------------------------------------
  # action menus
  def ActionsForSelection(self, source_window):
    all_sel = self.GetSelections()
    if (len(all_sel) == 0):
      raise Sorry("You must select an object in the PDB hierarchy to edit.")
    elif (len(all_sel) > 1):
      object_types, n_atoms = self.GetSelectionInfo()
      labels_and_actions = [
        ("Set occupancy...", self.OnSetOccupancy),
        ("Set B-factor...", self.OnSetBfactor),
        ("Set segment ID...", self.OnSetSegID),
        ("Convert to isotropic", self.OnSetIsotropic),
      ]
      if (len(object_types) == 1):
        labels_and_actions.extend([
          ("Delete object(s)...", self.OnDeleteObject),
          ("Apply rotation/translation...", self.OnMoveSites),
        ])
        if (self._crystal_symmetry is not None):
          labels_and_actions.extend([
            ("Apply symmetry operator...", self.OnApplySymop),
          ])
      if (object_types == ["residue_group"]):
        labels_and_actions.extend([
          ("Renumber residues...", self.OnRenumberResidues),
        ])
      labels_and_actions.extend([
          ("Toggle ATOM/HETATM...", self.OnSetAtomType),
          ("Reset element field...", self.OnResetElement),
          ("Set altloc...", self.OnSetAltlocOther),
      ])
      self.ShowMenu(labels_and_actions, source_window)
    else :
      sel = self.GetSelection()
      pdb_object = self.GetItemPyData(sel)
      if (pdb_object is not None):
        pdb_type = type(pdb_object).__name__
        if (pdb_type == "atom"):
          self.ShowAtomMenu(
            atom=pdb_object,
            source_window=source_window)
        elif (pdb_type == "atom_group"):
          self.ShowAtomGroupMenu(
            atom_group=pdb_object,
            source_window=source_window)
        elif (pdb_type == "residue_group"):
          self.ShowResidueGroupMenu(
            residue_group=pdb_object,
            source_window=source_window)
        elif (pdb_type == "chain"):
          self.ShowChainMenu(
            chain=pdb_object,
            source_window=source_window)
        elif (pdb_type == "model"):
          self.ShowModelMenu(
            model=pdb_object,
            source_window=source_window)
        else :
          raise RuntimeError("Unrecognized object type '%s'" % pdb_type)

  def ShowAtomMenu(self, atom, source_window):
    labels_and_actions = [
      ("Set name...", self.OnSetName),
      ("Set occupancy...", self.OnSetOccupancy),
      ("Set B-factor...", self.OnSetBfactor),
      ("Set element...", self.OnSetElement),
      ("Set charge...", self.OnSetCharge),
      ("Convert to isotropic", self.OnSetIsotropic),
      ("Delete atom", self.OnDeleteObject),
      ("Apply rotation/translation...", self.OnMoveSites),
      ("Toggle ATOM/HETATM...", self.OnSetAtomType),
      ("Reset element field...", self.OnResetElement),
    ]
    if (self._crystal_symmetry is not None):
      labels_and_actions.extend([
        ("Apply symmetry operator...", self.OnApplySymop),
      ])
    self.ShowMenu(labels_and_actions, source_window)

  def ShowAtomGroupMenu(self, atom_group, source_window):
    labels_and_actions = [
      ("Set altloc...", self.OnSetAltloc),
      ("Set occupancy...", self.OnSetOccupancy),
      ("Set B-factor...", self.OnSetBfactor),
      ("Set residue name...", self.OnSetResname),
      ("Convert to isotropic", self.OnSetIsotropic),
      ("Delete atom group", self.OnDeleteObject),
      ("Clone atom group", self.OnCloneAtoms),
      ("Apply rotation/translation...", self.OnMoveSites),
      ("Toggle ATOM/HETATM...", self.OnSetAtomType),
      ("Reset element field...", self.OnResetElement),
    ]
    if (atom_group.resname == "MET"):
      labels_and_actions.append(("Convert to SeMet", self.OnConvertMet))
    elif (atom_group.resname == "MSE"):
      labels_and_actions.append(("Convert to Met", self.OnConvertSeMet))
    labels_and_actions.append([
      ("Separate atoms...", self.OnSeparateAtoms),
    ])
    if (self._crystal_symmetry is not None):
      labels_and_actions.extend([
        ("Apply symmetry operator...", self.OnApplySymop),
      ])
    self.ShowMenu(labels_and_actions, source_window)

  def ShowResidueGroupMenu(self, residue_group, source_window):
    labels_and_actions = [
      ("Set residue number...", self.OnSetResseq),
      ("Set insertion code...", self.OnSetIcode),
      ("Set segment ID...", self.OnSetSegID),
      ("Set B-factor...", self.OnSetBfactor),
      ("Convert to isotropic", self.OnSetIsotropic),
      ("Delete residue", self.OnDeleteObject),
      ("Split residue", self.OnSplitResidue),
      ("Apply rotation/translation...", self.OnMoveSites),
      ("Insert residues after...", self.OnInsertAfter),
      ("Renumber residue...", self.OnRenumberResidues),
      ("Toggle ATOM/HETATM...", self.OnSetAtomType),
      ("Reset element field...", self.OnResetElement),
    ]
    if (len(residue_group.atom_groups()) == 1):
      labels_and_actions.append(("Set occupancy...", self.OnSetOccupancy))
    labels_and_actions.append([
      ("Separate atoms...", self.OnSeparateAtoms),
    ])
    if (self._crystal_symmetry is not None):
      labels_and_actions.extend([
        ("Apply symmetry operator...", self.OnApplySymop),
      ])
    self.ShowMenu(labels_and_actions, source_window)

  def ShowChainMenu(self, chain, source_window):
    labels_and_actions = [
      ("Set chain ID...", self.OnSetChainID),
      ("Set segment ID...", self.OnSetSegID),
      ("Set B-factor...", self.OnSetBfactor),
      ("Convert to isotropic", self.OnSetIsotropic),
      ("Delete chain", self.OnDeleteObject),
      ("Apply rotation/translation...", self.OnMoveSites),
      ("Toggle ATOM/HETATM...", self.OnSetAtomType),
      ("Reset element field...", self.OnResetElement),
      ("Add residues...", self.OnAddResidues),
      ("Renumber residues...", self.OnRenumberChain),
      ("Reset residue numbers...", self.OnResetNumbering),
    ]
    if (len(chain.conformers()) > 1):
      labels_and_actions.append(
        ("Delete alternate conformers", self.OnDeleteAltConfs))
    else :
      labels_and_actions.extend([
        ("Set altloc...", self.OnSetAltlocOther),
        ("Set occupancy...", self.OnSetOccupancy)])
    model = chain.parent()
    if (len(model.chains()) > 1):
      labels_and_actions.append(
        ("Merge with other chain...", self.OnMergeChain))
    if (self._crystal_symmetry is not None):
      labels_and_actions.extend([
        ("Apply symmetry operator...", self.OnApplySymop),
      ])
    self.ShowMenu(labels_and_actions, source_window)

  def ShowModelMenu(self, model, source_window):
    labels_and_actions = [
      ("Convert to isotropic", self.OnSetIsotropic),
      ("Delete model", self.OnDeleteObject),
      ("Apply rotation/translation...", self.OnMoveSites),
      ("Add chain(s)...", self.OnAddChain),
      ("Split model...", self.OnSplitModel),
      ("Toggle ATOM/HETATM...", self.OnSetAtomType),
      ("Reset element field...", self.OnResetElement),
    ]
    self.ShowMenu(labels_and_actions, source_window)

  def ShowMenu(self, items, source_window):
    def add_to_menu(m, label, action):
      item = m.Append(-1, label)
      source_window.Bind(wx.EVT_MENU, action, item)
    menu = wx.Menu()
    for entry in items :
      if isinstance(entry, list):
        advanced_menu = wx.Menu()
        menu.AppendMenu(-1, "Advanced options", advanced_menu)
        for label, action in entry :
          add_to_menu(advanced_menu, label, action)
      else :
        label, action = entry
        add_to_menu(menu, label, action)
    source_window.PopupMenu(menu)
    menu.Destroy()

  #---------------------------------------------------------------------
  # PROPERTY EDITING ACTIONS
  # atom
  def OnSetName(self, event):
    item, atom = self.GetSelectedObject('atom')
    new_name = self.GetNewName(atom.name)
    assert (new_name is not None) and (1 <= len(new_name) <= 4)
    if (new_name != atom.name):
      atom_group = atom.parent()
      for other_atom in atom_group.atoms():
        if (other_atom.name == new_name) and (atom != other_atom):
          confirm_action(("The atom group to which this atom belongs already "+
          "has another atom named \"%s\".  Are you sure you want to rename "+
          "the selected atom?") % new_name)
      atom.name = "%-4s" % new_name
      self.SetItemText(item, format_atom(atom))
      self.PushState("changed atom name ('%s' -> '%s')" % (atom.name,new_name))

  # atom
  def OnSetElement(self, event):
    item, atom = self.GetSelectedObject('atom')
    new_elem = self.GetNewElement(atom.element)
    assert (new_elem is None) or (len(new_elem) <= 2)
    if (new_elem != atom.element):
      if (new_elem is None) or (new_elem.isspace()):
        atom.element = '  '
      else :
        new_elem = new_elem.strip()
        atom.element = '%2s' % new_elem
      # TODO validate element symbol
      self.SetItemText(item, format_atom(atom))
      self.PushState("changed atom element to '%s'" % atom.element)

  # atom
  def OnSetCharge(self, event):
    item, atom = self.GetSelectedObject('atom')
    new_charge = self.GetNewCharge(atom.charge)
    assert (new_charge is None) or (-9 <= new_charge <= 9)
    if (new_charge in [0, None]):
      atom.charge = '  '
    elif (new_charge < 0):
      atom.charge = '%d-' % abs(new_charge)
    else :
      atom.charge = '%d+' % new_charge
    self.SetItemText(item, format_atom(atom))
    self.PushState("set atom charge to %s" % new_charge)

  # atom
  def OnSetXYZ(self, event):
    pass

  # atom_group
  def OnSetAltloc(self, event):
    """
    Change the altloc ID of an atom group.  This actually checks whether a
    split residue is having a blank altloc set, unlike the alternate function
    for arbitrary selections below.
    """
    item, atom_group = self.GetSelectedObject()
    new_altloc = self.GetNewAltloc(atom_group.altloc)
    assert (new_altloc is None) or (len(new_altloc) in [0,1])
    if (new_altloc != atom_group.altloc):
      if (new_altloc in [None, '']):
        rg = atom_group.parent()
        if (len(rg.atom_groups()) > 1):
          confirm_action("You have specified a blank altloc ID for this atom "+
            "group, but it is part of a residue containing multiple "+
            "conformations.  Are you sure this is what you want to do?")
        atom_group.altloc = ''
      else :
        atom_group.altloc = new_altloc
      self.SetItemText(item, format_atom_group(atom_group))
      self.PushState("set altloc to '%s'" % new_altloc)

  # atom_group
  def OnSetResname(self, event):
    item, atom_group = self.GetSelectedObject('atom_group')
    new_resname = self.GetNewResname(atom_group.resname)
    assert (new_resname is not None) and (len(new_resname) in [1,2,3])
    if (atom_group.resname != new_resname):
      atom_group.resname = new_resname
      self.SetItemText(item, format_atom_group(atom_group))
      rg_item = self.GetItemParent(item)
      self.SetItemText(rg_item, format_residue_group(atom_group.parent()))
      self.PushState("set resname to '%s" % new_resname)

  # atom_group (resname == MET)
  def OnConvertMet(self, event):
    item, atom_group = self.GetSelectedObject('atom_group')
    assert (atom_group.resname == "MET")
    atom_group.resname = "MSE"
    for atom in atom_group.atoms():
      if (atom.name == ' SD '):
        atom.name = ' SE '
        atom.element = 'Se'
        break
    self.SetItemText(item, format_atom_group(atom_group))
    rg_item = self.GetItemParent(item)
    self.SetItemText(rg_item, format_residue_group(atom_group.parent()))
    self.PropagateAtomChanges(item)
    self.PushState("converted MET to MSE")

  # atom_group (resname == MSE)
  def OnConvertSeMet(self, event):
    item, atom_group = self.GetSelectedObject('atom_group')
    assert (atom_group.resname == "MSE")
    atom_group.resname = "MET"
    for atom in atom_group.atoms():
      if (atom.name == ' SE '):
        atom.name = ' SD '
        atom.element = ' S'
        break
    self.SetItemText(item, format_atom_group(atom_group))
    rg_item = self.GetItemParent(item)
    self.SetItemText(rg_item, format_residue_group(atom_group.parent()))
    self.PropagateAtomChanges(item)
    self.PushState("converted MSE to MET")

  # residue_group
  def OnSetIcode(self, event):
    item, residue_group = self.GetSelectedObject('residue_group')
    new_icode = self.GetNewIcode(residue_group.icode)
    assert (new_icode is None) or (len(new_icode) == 1)
    if (new_icode != residue_group.icode):
      if (new_icode is None):
        residue_group.icode = ' '
      else :
        residue_group.icode = new_icode
      self.SetItemText(item, format_residue_group(residue_group))
      self.PushState("set insertion code to '%s'" % new_icode)

  def OnSetResseq(self, event):
    item, residue_group = self.GetSelectedObject('residue_group')
    new_resseq = self.GetNewResseq(residue_group.resseq_as_int())
    assert (new_resseq is not None)
    if (new_resseq != residue_group.resseq_as_int()):
      residue_group.resseq = "%4d" % new_resseq
      self.SetItemText(item, format_residue_group(residue_group))
      self.PushState("set residue number to '%s'" % new_resseq)

  def OnRenumberResidues(self, event):
    items = self.GetSelections()
    resseq_shift = self.GetResseqShift()
    if (resseq_shift is not None) and (resseq_shift != 0):
      for item in items :
        residue_group = self.GetItemPyData(item)
        assert (type(residue_group).__name__ == 'residue_group')
        resseq = residue_group.resseq_as_int()
        new_resseq = resseq + resseq_shift
        if (new_resseq > 9999) or (new_resseq < -999):
          raise NotImplementedError("Hybrid36 support not available.")
        else :
          residue_group.resseq = "%4d" % new_resseq
        self.SetItemText(item, format_residue_group(residue_group))
      self.PushState("renumber residues by '%d'" % resseq_shift)

  # chain
  def OnSetChainID(self, event):
    item, chain = self.GetSelectedObject('chain')
    new_id = self.GetNewChainID(chain.id)
    if (new_id != chain.id):
      if (new_id is None) or (new_id.isspace()):
        chain.id = ' '
      else :
        chain.id = "%2s" % new_id
      self.SetItemText(item, format_chain(chain))
      self.PushState("set chain ID to '%s'" % new_id)

  # chain
  def OnRenumberChain(self, event):
    item, chain = self.GetSelectedObject('chain')
    resseq_shift = self.GetResseqShift()
    if (resseq_shift is not None) and (resseq_shift != 0):
      child, cookie = self.GetFirstChild(item)
      while (child is not None):
        residue_group = self.GetItemPyData(child)
        assert (type(residue_group).__name__ == 'residue_group')
        resseq = residue_group.resseq_as_int()
        new_resseq = resseq + resseq_shift
        if (new_resseq > 9999) or (new_resseq < -999):
          raise NotImplementedError("Hybrid36 support not available.")
        else :
          residue_group.resseq = "%4d" % new_resseq
        self.SetItemText(child, format_residue_group(residue_group))
        child, cookie = self.GetNextChild(item, cookie)
      self.PushState("incremented residue numbers by %d" % resseq_shift)

  # chain
  def OnResetNumbering(self, event):
    item, chain = self.GetSelectedObject('chain')
    confirm = wx.MessageBox("This will renumber all residues in the chain "+
      "to start from 1, and reset the insertion code to blank.  Are you "+
      "sure you want to do this?", style=wx.YES_NO)
    if (confirm == wx.YES):
      i_res = 1
      child, cookie = self.GetFirstChild(item)
      while (child is not None):
        residue_group = self.GetItemPyData(child)
        assert (type(residue_group).__name__ == 'residue_group')
        residue_group.resseq = "%4d" % i_res
        i_res += 1
        self.SetItemText(child, format_residue_group(residue_group))
        child, cookie = self.GetNextChild(item, cookie)
      self.PushState("reset numbering for chain '%s'" % chain.id)

  # model
  def OnSetModelID(self, event):
    pass

  # all
  def OnSetOccupancy(self, event):
    items = self.GetSelections()
    new_occ = None
    for item in items :
      pdb_object = self.GetItemPyData(item)
      pdb_type = type(pdb_object).__name__
      occ = None
      if (pdb_type == 'atom'):
        occ = pdb_object.occ
      else :
        all_occ = pdb_object.atoms().extract_occ()
        if (all_occ.all_eq(all_occ[0])):
          occ = all_occ[0]
      if (new_occ is None):
        new_occ = self.GetNewOccupancy(occ)
      assert (0 <= new_occ <= 1.0)
      def apply_occ(atom) : atom.occ = new_occ
      self._ApplyToAtoms(item, pdb_object, apply_occ)

  # all
  def OnSetBfactor(self, event):
    items = self.GetSelections()
    new_b = None
    for item in items :
      pdb_object = self.GetItemPyData(item)
      pdb_type = type(pdb_object).__name__
      b_iso = None
      if (pdb_type == 'atom'):
        b_iso = pdb_object.b
      else :
        all_b = pdb_object.atoms().extract_b()
        if (all_b.all_eq(all_b[0])):
          b_iso = all_b[0]
      if (new_b is None):
        new_b = self.GetNewBiso(b_iso)
        assert (0 < new_b < 1000)
      def apply_b(atom) : atom.b = new_b
      self._ApplyToAtoms(item, pdb_object, apply_b)

  # all
  def OnSetIsotropic(self, event):
    items = self.GetSelections()
    for item in items :
      pdb_object = self.GetItemPyData(item)
      def set_isotropic(atom):
        atom.set_uij(new_uij=(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
      self._ApplyToAtoms(item, pdb_object, set_isotropic)

  # all
  def OnSetSegID(self, event):
    items = self.GetSelections()
    new_segid = None
    for item in items :
      pdb_object = self.GetItemPyData(item)
      pdb_type = type(pdb_object).__name__
      segid = None
      if (pdb_type == 'atom'):
        segid = pdb_object.segid
      else :
        from scitbx.array_family import flex
        segids = flex.std_string([ a.segid for a in pdb_object.atoms() ])
        if (segids.all_eq(segids[0])):
          segid = segids[0]
      if (segid is not None) and (segid.isspace()):
        segid = None
      if (new_segid is None):
        new_segid = self.GetNewSegID(segid)
      if (new_segid != segid):
        def apply_segid(atom) : atom.segid = new_segid
        self._ApplyToAtoms(item, pdb_object, apply_segid)

  # all
  def OnSetAltlocOther(self, event):
    """
    Set the altloc for an arbitrary collection of atoms, not just an atom
    group.
    """
    items = self.GetSelections()
    new_altloc = None
    for item in items :
      pdb_object = self.GetItemPyData(item)
      if (type(pdb_object).__name__ in ["residue_group", "chain"]):
        if (len(pdb_object.conformers()) > 1):
          raise Sorry("One or more residues in the selection contain "+
            "alternate conformers; you may only set the altloc identifier "+
            "for a single conformer at a time.")
      elif (type(pdb_object).__name__ == "atom"):
        atom_group = pdb_object.parent()
        ag_item = self.FindItem(atom_group)
        if (not ag_item in items):
          raise Sorry(("You have selected the atom '%s' to have its altloc "
            "identifier changed, but not the parent atom group.  Because the "+
            "altloc ID is an atom group property, you must set it for the "+
            "entire atom group, not individual atoms.") % pdb_object.id_str())
    new_altloc = self.GetNewAltloc("")
    modified = []
    for item in items :
      pdb_object = self.GetItemPyData(item)
      if (type(pdb_object) == 'atom_group'):
        pdb_object.altloc = new_altloc
        modified.append(pdb_object.memory_id())
        self.SetItemText(item, format_atom_group(pdb_object))
      else:
        atoms = []
        if (type(pdb_object) == 'atom'):
          atoms = [ pdb_object]
        else :
          atoms = pdb_object.atoms()
        for atom in atoms :
          atom_group = atom.parent()
          if (not atom_group.memory_id() in modified):
            atom_group.altloc = new_altloc
            modified.append(atom_group.memory_id())
            ag_item = self.FindItem(atom_group)
            self.SetItemText(ag_item, format_atom_group(atom_group))
    self.PushState("set altloc to %s" % new_altloc)

  # all
  def OnMoveSites(self, event):
    style = 0
    if (self._crystal_symmetry is not None):
      if (self._crystal_symmetry.unit_cell() is not None):
        style = simple_dialogs.RT_DIALOG_ENABLE_FRACTIONAL
    dlg = simple_dialogs.RTDialog(
      parent=self,
      title="Rotation/translation operator",
      wxtbxStyle=style)
    rt = None
    fractional = False
    if (dlg.ShowModal() == wx.ID_OK):
      rt = dlg.GetMatrix()
      fractional = dlg.IsFractional()
    wx.CallAfter(dlg.Destroy)
    if (rt is None):
      raise Abort()
    self.ApplyTransformation(rt, fractional=fractional)
    self.PushState("moved sites")

  # all
  def OnApplySymop(self, event):
    if (self._space_group is None):
      raise Sorry("No space group defined.")
    space_group = self._crystal_symmetry.space_group()
    space_group_info = space_group.info()
    dlg = simple_dialogs.SymopChoiceDialog(
      parent=self,
      title="Apply symmetry operator",
      label="Valid symmetry operators",
      caption="Please select a symmetry operator from the list of choices "+
        "supported for space group '%s'" % str(space_group_info))
    dlg.SetSpaceGroup(space_group)
    if (dlg.ShowModal() == wx.ID_OK):
      smx = dlg.GetValue()
    wx.CallAfter(dlg.Destroy)
    if (smx is None):
      raise Abort()
    from scitbx import matrix
    r = smx.r().as_double()
    t = smx.t().as_double()
    rt = matrix.rt((r,t))
    self.ApplyTransformation(rt, fractional=True)
    self.PushState("applied symop '%s'" % str(smx))

  def ApplyTransformation(self, rt, fractional=False):
    items = self.GetSelections()
    unit_cell = None
    if (fractional):
      unit_cell = self._crystal_symmetry.unit_cell()
      assert (unit_cell is not None)
    # FIXME
    for item in items :
      pdb_object = self.GetItemPyData(item)
      if has_aniso_atoms(pdb_object):
        confirm_action("This action will invalidate any ANISOU records "+
          "present in the selected atoms, which will automatically be "+
          "converted to isotropic.  Are you sure you want to continue?")
        break
    for item in items :
      pdb_object = self.GetItemPyData(item)
      pdb_type = type(pdb_object).__name__
      if (pdb_type == 'atom'):
        from scitbx.array_family import flex
        sites = flex.vec3_double([pdb_object.xyz])
        if (fractional):
          sites = unit_cell.fractionalize(sites_cart=sites)
        sites = rt.r.elems * sites + rt.t.elems
        if (fractional):
          sites = unit_cell.orthogonalize(sites_frac=sites)
        pdb_object.xyz = sites[0]
        if (pdb_object.uij != (-1,-1,-1,-1,-1,-1)):
          pdb_object.uij = (-1,-1,-1,-1,-1,-1)
        self.SetItemText(item, format_atom(pdb_object))
      else :
        atoms = pdb_object.atoms()
        sites = atoms.extract_xyz()
        if (fractional):
          sites = unit_cell.fractionalize(sites_cart=sites)
        sites = rt.r.elems * sites + rt.t.elems
        if (fractional):
          sites = unit_cell.orthogonalize(sites_frac=sites)
        atoms.set_xyz(sites)
        for atom in atoms :
          if (atom.uij != (-1,-1,-1,-1,-1,-1)):
            atom.uij = (-1,-1,-1,-1,-1,-1)
        self.PropagateAtomChanges(item)

  # all
  def OnResetElement(self, event):
    items = self.GetSelections()
    n_changed = n_unassigned = 0
    for item in items :
      pdb_object = self.GetItemPyData(item)
      pdb_type = type(pdb_object).__name__
      if (pdb_type == 'atom'):
        pdb_object.element = ''
        pdb_object.set_chemical_element_simple_if_necessary()
        if (pdb_object.element == ''):
          n_unassigned += 1
        else :
          n_changed += 1
      else :
        for atom in pdb_object.atoms():
          atom.element = ''
          atom.set_chemical_element_simple_if_necessary()
          if (atom.element == ''):
            n_unassigned += 1
          else :
            n_changed += 1
      self.PropagateAtomChanges(item)
    self.PushState("reset chemical element field")
    if (n_unassigned > 0):
      wx.MessageBox(("WARNING: %d atoms now have blank element fields because "+
        "the element could not be automatically determined from the atom "+
        "name.  You may need to set these manually!") % n_unassigned)
    elif (n_changed > 0):
      wx.MessageBox("Chemical element reset for %d atoms." % n_changed)

  # all
  def OnSetAtomType(self, event):
    items = self.GetSelections()
    n_hetatm = n_atom = 0
    for item in items :
      pdb_object = self.GetItemPyData(item)
      pdb_type = type(pdb_object).__name__
      if (pdb_type == 'atom'):
        if (pdb_object.hetero):
          n_hetatm += 1
        else :
          n_atom += 1
      else :
        for atom in pdb_object.atoms():
          if (atom.hetero):
            n_hetatm += 1
          else :
            n_atom += 1
    atom_type = "ATOM"
    if (n_atom == 0) and (n_hetatm > 0):
      atom_type = "HETATM"
    new_type = self.GetAtomType(atom_type)
    for item in items :
      pdb_object = self.GetItemPyData(item)
      pdb_type = type(pdb_object).__name__
      if (pdb_type == 'atom'):
        if (new_type == "HETATM"):
          pdb_object.hetero = True
        else :
          pdb_object.hetero = False
      else :
        for atom in pdb_object.atoms():
          if (new_type == "HETATM"):
            atom.hetero = True
          else :
            atom.hetero = False
      self.PropagateAtomChanges(item)
    self.PushState("set record types to %s" % new_type)

  #---------------------------------------------------------------------
  # HIERARCHY EDITING
  def OnDeleteObject(self, event):
    self.DeleteSelected()

  def DeleteSelected(self):
    object_types, n_atoms = self.GetSelectionInfo()
    if (len(object_types) > 1):
      raise Sorry(("Multiple object types selected (%s) - you may only delete "+
        "one type of object at a time.") % (" ".join(object_types)))
    confirm_action("Are you sure you want to delete the selected %d atom(s)?"
      % n_atoms)
    for item in self.GetSelections():
      pdb_object = self.GetItemPyData(item)
      self._DeleteObject(item, pdb_object)
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("deleted %d atoms" % n_atoms)

  def _DeleteObject(self, item, pdb_object):
    remove_objects_recursive(pdb_object)
    self._RemoveItem(item)

  def _RemoveItem(self, item):
    parent_item = self.GetItemParent(item)
    self.DeleteChildren(item)
    self.Delete(item)
    while (not self.HasChildren(parent_item)):
      del_item = parent_item
      parent_item = self.GetItemParent(del_item)
      self.Delete(del_item)

  # chain
  def OnDeleteAltConfs(self, event):
    item, chain = self.GetSelectedObject('chain')
    n_alt_atoms = 0
    n_alt_residues = 0
    for residue_group in chain.residue_groups():
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1):
        n_alt_residues += 1
        for ag in atom_groups[1:] :
          n_alt_atoms += len(ag.atoms())
    if (n_alt_atoms == 0):
      raise Sorry("No alternate conformations found in this chain.")
    confirm_action(("There are %d residues with alternate conformations; "+
      "removing these will delete %d atoms from the model.  Are you sure "+
      "you want to continue?") % (n_alt_residues, n_alt_atoms))
    # TODO more control over what happens to remaining atom_groups
    child, cookie = self.GetFirstChild(item)
    while (child is not None):
      residue_group = self.GetItemPyData(child)
      assert (type(residue_group).__name__ == 'residue_group')
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1):
        first_group = atom_groups[0]
        first_group.altloc = ''
        for atom in first_group.atoms():
          atom.occ = 1.0
        child2, cookie2 = self.GetFirstChild(child)
        self.SetItemText(child2, format_atom_group(first_group))
        for atom_group in atom_groups[1:] :
          residue_group.remove_atom_group(atom_group)
          item2 = self.FindItem(atom_group)
          assert (item2 is not None)
          self.DeleteChildren(item2)
          self.Delete(item2)
        self.SetItemText(child, format_residue_group(residue_group))
      child, cookie = self.GetNextChild(item, cookie)
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("deleted %d atoms with alternate conformations" %
      n_alt_atoms)

  # atom_group
  def OnCloneAtoms(self, event):
    item, target_atom_group = self.GetSelectedObject('atom_group')
    residue_group = target_atom_group.parent()
    atom_groups = residue_group.atom_groups()
    assert (len(atom_groups) > 0)
    start_occ = 1/(len(atom_groups) + 1)
    new_occ = self.GetNewOccupancy(start_occ, new=True)
    self._AddAtomGroup(
      item=self.GetItemParent(item),
      residue_group=residue_group,
      new_group=target_atom_group.detached_copy(),
      new_occ=new_occ)
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("cloned atom group")

  # atom_group, residue_group
  def OnSeparateAtoms(self, event):
    item, pdb_object = self.GetSelectedObject()
    assert (type(pdb_object).__name__ in ["atom_group", "residue_group"])
    confirm_action("This action will divide all atoms in the selected "+
      "residue into separate residues, and renumber them starting "+
      "from 1.  Are you sure you want to do this?")
    if (type(pdb_object).__name__ == "atom_group"):
      self.SeparateAtoms(item, pdb_object)
    else :
      for atom_group in pdb_object.atom_groups():
        ag_item = self.FindItem(atom_group)
        self.SeparateAtoms(ag_item, atom_group)

  def SeparateAtoms(self, item, atom_group):
    import iotbx.pdb.hierarchy
    rg_item = self.GetItemParent(item)
    chain_item = self.GetItemParent(rg_item)
    assert (not None in [rg_item, chain_item])
    rg = atom_group.parent()
    chain = rg.parent()
    i = 1
    for atom in atom_group.atoms():
      atom_new = atom.detached_copy()
      atom_group.remove_atom(atom)
      atom_item = self.FindItem(atom)
      self.Delete(atom_item)
      ag_new = iotbx.pdb.hierarchy.atom_group(resname=atom_group.resname,
        altloc=atom_group.altloc)
      ag_new.append_atom(atom_new)
      rg_new = iotbx.pdb.hierarchy.residue_group(resseq="%d" % i)
      rg_new.append_atom_group(ag_new)
      chain.append_residue_group(rg_new)
      self._InsertResidueGroupItem(chain_item, rg_new)
      i += 1
    self.Delete(item)
    if (len(rg.atoms()) == 0):
      chain.remove_residue_group(rg)
      self.Delete(rg_item)
    self.PushState("split residue")

  # residue_group
  def OnSplitResidue(self, event):
    item, residue_group = self.GetSelectedObject('residue_group')
    atom_groups = residue_group.atom_groups()
    assert (len(atom_groups) > 0)
    start_occ = 1/(len(atom_groups) + 1)
    new_occ = self.GetNewOccupancy(start_occ, new=True)
    self._AddAtomGroup(
      item=item,
      residue_group=residue_group,
      new_group=atom_groups[0].detached_copy(),
      new_occ=new_occ)
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("split residue")

  def _AddAtomGroup(self, item, residue_group, new_group, new_occ):
    atom_groups = residue_group.atom_groups()
    for atom in new_group.atoms():
      atom.occ = new_occ
    for atom_group in atom_groups :
      for atom in atom_group.atoms():
        atom.occ = max(0, atom.occ - new_occ/len(atom_groups))
    if (len(atom_groups) == 1):
      atom_groups[0].altloc = 'A'
      new_group.altloc = 'B'
    else :
      new_altloc = None
      for char in string.uppercase :
        for atom_group in atom_groups :
          if (atom_group.altloc == char):
            break
        else :
          new_altloc = char
      new_group.altloc = char
    residue_group.append_atom_group(new_group)
    self._InsertAtomGroupItem(item, new_group)
    self.PropagateAtomChanges(item)
    child, cookie = self.GetFirstChild(item)
    while (child is not None):
      atom_group = self.GetItemPyData(child)
      assert (type(atom_group).__name__ == 'atom_group')
      self.SetItemText(child, format_atom_group(atom_group))
      child, cookie = self.GetNextChild(item, cookie)
    self._hierarchy.atoms().reset_i_seq()

  # residue_group
  def OnInsertAfter(self, event):
    item, residue_group = self.GetSelectedObject('residue_group')
    new_residues = self.GetResiduesFromFile()
    chain = residue_group.parent()
    index = None
    for k, other_rg in enumerate(chain.residue_groups()):
      if (other_rg == residue_group):
        index = k+1
        break
    assert (index is not None)
    self._AddResidues(
      chain=chain,
      chain_node=self.GetItemParent(item),
      residues=new_residues,
      index=index)
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("inserted %d residues" % len(new_residues))

  # chain
  def OnAddResidues(self, event):
    item, chain = self.GetSelectedObject('chain')
    new_residues = self.GetResiduesFromFile()
    index = self.GetInsertionIndex(chain)
    self._AddResidues(
      chain=chain,
      chain_node=item,
      residues=new_residues,
      index=index)
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("added %d residues" % len(new_residues))

  # chain
  def OnMergeChain(self, event):
    item, chain = self.GetSelectedObject('chain')
    model = chain.parent()
    assert (len(model.chains()) > 1)
    target_chain = self.GetChainForMerging(model.chains(), chain)
    assert (chain != target_chain)
    target_item = self.FindItem(target_chain)
    assert (target_item is not None)
    index = self.GetInsertionIndex(target_chain)
    merge_residues = [ rg.detached_copy() for rg in chain.residue_groups() ]
    self._DeleteObject(item, chain)
    self._AddResidues(
      chain=target_chain,
      chain_node=target_item,
      residues=merge_residues,
      index=index)
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("merged chain '%s' with chain '%s'" % (chain.id,
      target_chain.id))

  # model
  def OnAddChain(self, event):
    item, model = self.GetSelectedObject('model')
    new_chains = self.GetChainsFromFile()
    for chain in new_chains :
      model.append_chain(chain)
      self._InsertChainItem(item, chain)
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("added %d chain(s) to model" % len(new_chains))
    self.Refresh()

  # model
  def OnSplitModel(self, event):
    item, model = self.GetSelectedObject('model')
    n_models = len(self._hierarchy.models())
    new_model = model.detached_copy()
    new_model.id = str(n_models + 1)
    self._hierarchy.append_model(new_model)
    root_node = self.GetRootItem()
    self._InsertModelItem(root_node, new_model)
    child, cookie = self.GetFirstChild(root_node)
    i_model = 1
    while (child is not None):
      other_model = self.GetItemPyData(child)
      assert (type(other_model).__name__ == 'model')
      other_model.id = str(i_model)
      self.SetItemText(child, format_model(other_model))
      child, cookie = self.GetNextChild(root_node, cookie)
      i_model += 1
    self._hierarchy.atoms().reset_i_seq()
    self.PushState("split model '%s'" % model.id)

  def _AddResidues(self, chain, chain_node, residues, index=None):
    for rg in residues :
      if (index is None):
        chain.append_residue_group(rg)
      else :
        chain.insert_residue_group(index, rg)
      self._InsertResidueGroupItem(chain_node, rg, index=index)
      if (index is not None):
        index += 1

  #---------------------------------------------------------------------
  # USER INPUT
  def GetNewName(self, name=None):
    dlg = simple_dialogs.StringDialog(
      parent=self,
      title="Set atom name",
      label="New name",
      caption="Please specify the atom name.  This is four characters in "+
        "length, but spaces will be added to the end if necessary.  Note that "+
        "leading spaces are significant, since they determine the column "+
        "alignment and the identity of the atom.  (For instance, 'CA  ' and "+
        "' CA ' have very different meanings.)",
      value=name)
    dlg.SetMinLength(1)
    dlg.SetMaxLength(4)
    dlg.SetOptional(False)
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewOccupancy(self, occ=None, new=False):
    desc_str = "selected"
    if (new):
      desc_str = "new"
    dlg = simple_dialogs.FloatDialog(
      parent=self,
      title="Set new occupancy",
      label="New occupancy",
      caption=("Please specify the occupancy for the %s atom(s); this "+
        "value represents the fraction of unit cells in which the atom(s) "+
        "is/are present, and must be a value between 0 (no contribution to "+
        "F_calc) and 1.0; note that atoms with multiple conformers should "+
        "always have a sum of occupancies of 1.0.  The precision will be "+
        "truncated to two digits after the decimal point.") % desc_str,
      value=occ)
    dlg.SetMin(0.0)
    dlg.SetMax(1.0)
    dlg.SetOptional(False)
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewBiso(self, b=None):
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
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewCharge(self, charge=None):
    if (charge is not None) and (charge.isspace()):
      charge = None
    elif (charge is not None):
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
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewElement(self, elem=None):
    if (elem is not None) and (elem.isspace()):
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
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewResseq(self, resseq):
    dlg = simple_dialogs.IntegerDialog(
      parent=self,
      title="Set residue number",
      label="Residue number",
      caption="The residue number can be any value, but the official PDB "+
        "format limits it to a range from -999 to 9999.  If you specify a "+
        "value outside of this range, it will be converted to Hybrid36 "+
        "encoding, which is recognized by Phenix, Coot, and CCP4, but not "+
        "by the PDB.",
      value=resseq)
    dlg.SetOptional(False)
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewResname(self, resname):
    dlg = simple_dialogs.StringDialog(
      parent=self,
      title="Set residue name",
      label="Residue code",
      caption="The residue name can be one to three characters, but is "+
        "nearly always three except for elemental ions.  Only alphanumeric "+
        "characters are allowed.",
      value=resname)
    dlg.SetMinLength(1)
    dlg.SetMaxLength(3)
    dlg.SetOptional(False)
    result = simple_dialogs.get_phil_value_from_dialog(dlg)
    if (not result.isalnum()):
      raise Sorry("Residue code must be alphanumeric (value: '%s')" % result)
    return result

  def GetNewIcode(self, icode):
    if (icode.isspace()):
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
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewAltloc(self, altloc_id):
    if (altloc_id.isspace()):
      altloc_id = None
    dlg = simple_dialogs.StringDialog(
      parent=self,
      title="Set altloc ID",
      label="New altloc",
      caption="Please specify an altloc (i.e. conformation) identifier.  By "+
        "default this is left blank, unless the selected residue group "+
        "either has multiple conformations, or is present at partial "+
        "occupancy and interacts with one conformation of an adjacent "+
        "group.  Non-blank altloc IDs must be a single character, usually "+
        "uppercase letters.",
      value=altloc_id)
    #dlg.SetMinLength(1)
    dlg.SetMaxLength(1)
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewChainID(self, chain_id):
    if (chain_id.isspace()):
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
      value=chain_id)
    #dlg.SetMinLength(1)
    dlg.SetMaxLength(2)
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewSegID(self, segid):
    if (segid is not None) and (segid.isspace()):
      segid = None
    dlg = simple_dialogs.StringDialog(
      parent=self,
      title="Set segment ID",
      label="New ID",
      caption="The segment ID (segid) is an optional field for disambiguating "+
        "between chains with identical IDs, or for otherwise flagging part of "+
        "the model.  It may be up to four characters in length (with spaces "+
        "significant).  Note that while most programs in Phenix should "+
        "preserve the segid and use it for atom selections, it is no longer "+
        "accepted as part of the official PDB format.",
      value=segid)
    dlg.SetMaxLength(4)
    dlg.SetOptional(True)
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetResseqShift(self):
    dlg = simple_dialogs.IntegerDialog(
      parent=self,
      title="Shift residue numbers",
      label="Increment by",
      caption="You may shift the residue numbers by any amount, even if this "+
        "will make them negative; note however that the official PDB format "+
        "is limited to residue numbers between -999 and 999.  (PHENIX and "+
        "other programs can handle a much larger range using Hybrid36 "+
        "encoding, but this is not supported by the PDB.)",
      value=None)
    dlg.SetOptional(False)
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def _GetPDBFile(self):
    from iotbx import file_reader
    file_name = self.path_mgr.select_file(
      parent=self,
      message="Select PDB file to insert",
      wildcard=file_reader.get_wildcard_strings(["pdb"]))
    pdb_in = file_reader.any_file(to_str(file_name),
      force_type="pdb",
      raise_sorry_if_errors=True)
    hierarchy = pdb_in.file_object.hierarchy
    if (len(hierarchy.models()) > 1):
      raise Sorry("Multi-MODEL PDB files not supported for this action.")
    return hierarchy

  def GetChainForMerging(self, chains, merge_chain):
    assert (len(chains) > 1)
    dlg = SelectChainDialog(
      parent=self,
      title="Select chain to merge with",
      message=("You have selected to merge the chain '%s' (%d residues) with "+
        "another chain in the model; please select a target chain.") %
        (merge_chain.id, len(merge_chain.residue_groups())),
      chains=chains,
      exclude_chain=merge_chain)
    target_chain = None
    if (dlg.ShowModal() == wx.ID_OK):
      target_chain = dlg.GetChain()
    wx.CallAfter(dlg.Destroy)
    if (target_chain is None) : raise Abort()
    return target_chain

  def GetChainsFromFile(self):
    hierarchy = self._GetPDBFile()
    chains = hierarchy.models()[0].chains()
    if (len(chains) == 1):
      return [ chains[0].detached_copy() ]
    else :
      dlg = SelectChainDialog(
        parent=self,
        title="Select chains to add",
        message="Multiple chains were found in this PDB file; you may insert "+
          "them all at once, or one at a time.",
        chains=chains,
        allow_all_chains=True)
      add_chains = None
      if (dlg.ShowModal() == wx.ID_OK):
        add_chains = dlg.GetChains()
      wx.CallAfter(dlg.Destroy)
      if (add_chains is None) : raise Abort()
      return [ c.detached_copy() for c in add_chains ]

  def GetResiduesFromFile(self):
    hierarchy = self._GetPDBFile()
    chains = hierarchy.models()[0].chains()
    if (len(chains) == 1):
      return [ rg.detached_copy() for rg in chains[0].residue_groups() ]
    else :
      dlg = SelectChainDialog(
        parent=self,
        title="Select residues to add",
        message="Multiple chains were found in this PDB file; please select "+
          "the one containing the residues you wish to insert into the "+
          "selected chain.",
        chains=chains,
        allow_all_chains=False,
        show_residue_range=True)
      add_chains = resseq_start = resseq_end = None
      if (dlg.ShowModal() == wx.ID_OK):
        add_chains = dlg.GetChains()
        resseq_start, resseq_end = dlg.GetResseqRange()
      wx.CallAfter(dlg.Destroy)
      if (add_chains is None) : raise Abort()
      return [ rg.detached_copy() for rg in add_chains[0].residue_groups() ]

  def GetInsertionIndex(self, chain):
    dlg = AddResiduesDialog(
      parent=self,
      title="Insert residues")
    insert_type = resid = None
    if (dlg.ShowModal() == wx.ID_OK):
      insert_type = dlg.GetInsertionType()
      resid = dlg.GetResID()
    wx.CallAfter(dlg.Destroy)
    if (insert_type is None):
      raise Abort()
    elif (insert_type == ADD_RESIDUES_START):
      return 0
    elif (insert_type == ADD_RESIDUES_END):
      return None
    else :
      if (resid is None):
        raise Sorry("You need to specify a residue ID to insert after.")
      for k, residue_group in enumerate(chain.residue_groups()):
        if (residue_group.resid().strip() == resid.strip()):
          return k + 1
      raise Sorry("The residue ID %s was not found in the target chain." %
        resid.strip())

  def GetAtomType(self, atom_type):
    dlg = simple_dialogs.ChoiceDialog(
      parent=self,
      title="Set atom record type",
      label="Record type",
      value=None,
      caption="Please select the label for the selected atom record(s).  "+
        "This will not affect the refinement or visualization of these "+
        "atoms, but HETATM records may be handled differently by some "+
        "programs (including Phaser).")
    if (atom_type in ["ATOM", None]):
      dlg.SetChoices(["*ATOM", "HETATM"])
    else :
      dlg.SetChoices(["ATOM", "*HETATM"])
    new_type = None
    if (dlg.ShowModal() == wx.ID_OK):
      new_type = dlg.GetPhilValue()
    wx.CallAfter(dlg.Destroy)
    if (new_type is None):
      raise Abort()
    return new_type

########################################################################
# AUXILARY GUI CLASSES
class PDBTreeDropTarget(wx.DropTarget):
  def __init__(self, tree):
    wx.DropTarget.__init__(self)
    self.tree = tree
    #self.df = wx.CustomDataFormat("pdb_object")
    #self.cdo = wx.CustomDataObject(self.df)
    #self.SetDataObject(self.cdo)

  def OnDrop(self, x, y):
    return self.tree.OnDrop(x,y)

  def OnEnter(self, x, y, d):
    return d

  def OnDragOver(self, x, y, d):
    print(x, y, d)
    return d

  def OnData(self, x, y, d):
    return d

########################################################################
# INPUT DIALOGS
ADD_RESIDUES_START = 0
ADD_RESIDUES_END = 1
ADD_RESIDUES_AFTER = 2
class AddResiduesDialog(wx.Dialog):
  def __init__(self, *args, **kwds):
    super(AddResiduesDialog, self).__init__(*args, **kwds)
    style = self.GetWindowStyle()
    style |= wx.WS_EX_VALIDATE_RECURSIVELY|wx.RAISED_BORDER|wx.CAPTION
    self.SetWindowStyle(style)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    caption_txt = wx.StaticText(self, -1, "Please specify where to insert "+
      "the new residue.  The residue ID is the combination of residue number "+
      "and insertion code, but you only need to specify the residue number "+
      "if the insertion code is blank or irrelevant.")
    caption_txt.Wrap(480)
    szr.Add(caption_txt, 0, wx.ALL, 5)
    box = wx.BoxSizer(wx.HORIZONTAL)
    txt1 = wx.StaticText(self, -1, "Insert residues:")
    box.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._insert_choice = wx.Choice(self, -1,
      choices=["at the start of the chain",
               "at the end of the chain",
               "after residue"])
    self._insert_choice.SetSelection(ADD_RESIDUES_END)
    self.Bind(wx.EVT_CHOICE, self.OnChangeInsertion, self._insert_choice)
    box.Add(self._insert_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt2 = wx.StaticText(self, -1, "residue ID:")
    box.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._resid_txt = strctrl.StrCtrl(
      parent=self,
      value=None,
      size=(80,-1))
    self._resid_txt.SetOptional(True)
    self._resid_txt.Enable(False)
    self._resid_txt.SetMaxLength(5)
    self._resid_txt.SetMinLength(1)
    box.Add(self._resid_txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr.Add(box)
    cancel_btn = wx.Button(self, wx.ID_CANCEL)
    ok_btn = wx.Button(self, wx.ID_OK)
    ok_btn.SetDefault()
    szr4 = wx.StdDialogButtonSizer()
    szr4.Add(cancel_btn)
    szr4.Add(ok_btn, 0, wx.LEFT, 5)
    szr.Add(szr4, 0, wx.ALL|wx.ALIGN_RIGHT, 5)
    szr.Layout()
    self.Fit()
    self.Centre(wx.BOTH)

  def OnChangeInsertion(self, event):
    insert_type = self._insert_choice.GetSelection()
    if (insert_type == ADD_RESIDUES_AFTER):
      self._resid_txt.SetOptional(False)
      self._resid_txt.Enable()
    else :
      self._resid_txt.SetOptional(True)
      self._resid_txt.Enable(False)

  def GetInsertionType(self):
    return self._insert_choice.GetSelection()

  def GetResID(self):
    return self._resid_txt.GetPhilValue()

class SelectChainDialog(wx.Dialog):
  def __init__(self, parent, title, message, chains, allow_all_chains=False,
      show_residue_range=False, exclude_chain=None):
    self._chains = []
    for chain in chains :
      if (chain != exclude_chain):
        self._chains.append(chain)
    self.allow_all_chains = allow_all_chains
    super(SelectChainDialog, self).__init__(parent=parent, title=title,
      style=wx.WS_EX_VALIDATE_RECURSIVELY|wx.RAISED_BORDER|wx.CAPTION)
    style = self.GetWindowStyle()
    style |= wx.WS_EX_VALIDATE_RECURSIVELY|wx.RAISED_BORDER|wx.CAPTION
    self.SetWindowStyle(style)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    caption_txt = wx.StaticText(self, -1, message)
    caption_txt.Wrap(480)
    szr.Add(caption_txt, 0, wx.ALL, 5)
    box = wx.FlexGridSizer(rows=2, cols=2)
    txt1 = wx.StaticText(self, -1, "Select chain:")
    box.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    choices = []
    if (allow_all_chains):
      choices.append("All")
    for chain in self._chains :
      n_rg = len(chain.residue_groups())
      choices.append("chain '%s' (%d residues)" % (chain.id, n_rg))
    self._chain_choice = wx.Choice(self, -1, choices=choices)
    box.Add(self._chain_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    if (show_residue_range):
      txt2 = wx.StaticText(self, -1, "Residue range (optional):")
      box.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      box2 = wx.BoxSizer(wx.HORIZONTAL)
      box.Add(box2)
      self._resseq_start = intctrl.IntCtrl(
        parent=self,
        name="Starting resseq",
        value=None)
      self._resseq_start.SetOptional(True)
      self._resseq_end = intctrl.IntCtrl(
        parent=self,
        name="Ending resseq",
        value=None)
      self._resseq_end.SetOptional(False)
      txt3 = wx.StaticText(self, -1, "to")
      box2.Add(self._resseq_start, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      box2.Add(txt3, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      box2.Add(self._resseq_end, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    else :
      self._resseq_start = self._resseq_end = None
    szr.Add(box)
    cancel_btn = wx.Button(self, wx.ID_CANCEL)
    ok_btn = wx.Button(self, wx.ID_OK)
    ok_btn.SetDefault()
    szr4 = wx.StdDialogButtonSizer()
    szr4.Add(cancel_btn)
    szr4.Add(ok_btn, 0, wx.LEFT, 5)
    szr.Add(szr4, 0, wx.ALL|wx.ALIGN_RIGHT, 5)
    szr.Layout()
    self.Fit()
    self.Centre(wx.BOTH)

  def GetChains(self):
    chain_sel = self._chain_choice.GetSelection()
    if (self.allow_all_chains):
      if (chain_sel == 0):
        return self._chains
      else :
        chain_sel -= 1
    return [ self._chains[chain_sel] ]

  def GetChain(self):
    assert (not self.allow_all_chains)
    return self.GetChains()[0]

  def GetResseqRange(self):
    assert (not None in [self._resseq_start, self._resseq_end])
    return self._resseq_start.GetPhilValue(), self._resseq_end.GetPhilValue()

########################################################################
class PDBTreeFrame(wx.Frame):
  def __init__(self, *args, **kwds):
    wx.Frame.__init__(self, *args, **kwds)
    self.statusbar = self.CreateStatusBar()
    self.statusbar.SetStatusText("Right-click any item for a list of editing actions")
    # panel setup
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
    self._header_box.SetValue(False)
    szr3.Add(self._header_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._serial_box = wx.CheckBox(self.panel, -1,
      "Reset atom serial numbers")
    self._serial_box.SetValue(True)
    szr3.Add(self._serial_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    pszr.Add(szr3, 0, wx.EXPAND)
    self._tree = PDBTree(self.panel, -1, style=wx.RAISED_BORDER)
    self._tree.SetMinSize((640,400))
    pszr.Add(self._tree, 1, wx.EXPAND, 2)
    # toolbar setup
    self.toolbar = self.CreateToolBar(style=wx.TB_TEXT)
    if wx.VERSION < (4,0):
      self.toolbar.AddTool = self.toolbar.AddLabelTool
    bmp = wxtbx.bitmaps.fetch_custom_icon_bitmap("phenix.pdbtools")
    btn = self.toolbar.AddTool(-1, "Load file", bmp,
      shortHelp="Load file", kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnOpen, btn)
    bmp = wxtbx.bitmaps.fetch_icon_bitmap("actions", "save_all")
    btn = self.toolbar.AddTool(-1, "Save file", bmp,
      shortHelp="Save file", kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnSave, btn)
    self.toolbar.AddSeparator()
    bmp = wxtbx.bitmaps.fetch_custom_icon_bitmap("tools")
    btn = self.toolbar.AddTool(-1, "Edit...", bmp, shortHelp="Edit...",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnEditModel, btn)
    bmp = wxtbx.bitmaps.fetch_custom_icon_bitmap("symmetry")
    btn = self.toolbar.AddTool(-1, "Symmetry", bmp, shortHelp="Symmetry",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnEditSymmetry, btn)
    bmp = wxtbx.bitmaps.fetch_icon_bitmap("actions", "editdelete")
    btn = self.toolbar.AddTool(-1, "Delete", bmp, shortHelp="Delete",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self._tree.OnDeleteObject, btn)
    self.toolbar.Realize()
    # menu setup
    menu_bar = wx.MenuBar()
    file_menu = wx.Menu()
    menu_bar.Append(file_menu, "File")
    item = file_menu.Append(-1, "&Open file\tCtrl-O")
    self.Bind(wx.EVT_MENU, self.OnOpen, item)
    item = file_menu.Append(-1, "&Quit\tCtrl-Q")
    self.Bind(wx.EVT_MENU, self.OnClose, item)
    edit_menu = wx.Menu()
    #menu_bar.Append(edit_menu, "Edit")
    self.undo_item = edit_menu.Append(-1, "&Undo\tCtrl-Z")
    self.redo_item = edit_menu.Append(-1, "&Redo\tCtrl-Y")
    self.undo_item.Enable(False)
    self.redo_item.Enable(False)
    self.Bind(wx.EVT_MENU, self.OnUndo, self.undo_item)
    self.Bind(wx.EVT_MENU, self.OnRedo, self.redo_item)
    self.SetMenuBar(menu_bar)
    #
    self._pdb_in = None
    self._crystal_symmetry = None
    szr.Layout()
    self.Fit()
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy, self)
    self.path_mgr = self._tree.path_mgr
    self._callback = None
    self._callback_program = None

  def LoadPDB(self, file_name):
    from iotbx import file_reader
    f = file_reader.any_file(file_name,
      force_type="pdb",
      raise_sorry_if_errors=True)
    self._pdb_in = f
    self._crystal_symmetry = f.crystal_symmetry()
    hierarchy = f.file_object.hierarchy
    self._hierarchy = hierarchy
    atoms = hierarchy.atoms()
    atoms.reset_i_seq()
    atoms.set_chemical_element_simple_if_necessary()
    self._tree.SetHierarchy(hierarchy)
    self._tree.SetCrystalSymmetry(self._crystal_symmetry)
    self._path_field.SetValue(os.path.abspath(to_unicode(f.file_name)))
    self._tree.SetFocus()
    self.Refresh()

  def OnOpen(self, event):
    from iotbx import file_reader
    file_name = self.path_mgr.select_file(
      parent=self,
      message="Select PDB file",
      wildcard=file_reader.get_wildcard_strings(["pdb"]))
    self.LoadPDB(to_str(file_name))

  def OnSave(self, event):
    if (self._pdb_in is None):
      return
    self.Save()

  def Save(self):
    from iotbx import file_reader
    input_file = os.path.abspath(self._pdb_in.file_name)
    output_file = os.path.splitext(input_file)[0] + "_edited.pdb"
    file_name = self.path_mgr.select_file(
      parent=self,
      message="Save PDB file",
      style=wx.FD_SAVE,
      wildcard=file_reader.get_wildcard_strings(["pdb"]),
      current_file=output_file)
    save_header = self._header_box.GetValue()
    f = open(file_name, "w")
    for method in ["title_section",] :
      section_lines = getattr(self._pdb_in.file_object.input, method)()
      for line in section_lines :
        f.write(line + "\n")
    if (save_header):
      wx.MessageBox("Warning: any atom selections present in the PDB header "+
        "may be invalidated by model modifications!  Please use caution "+
        "when using the modified PDB file in any downstream applications.")
      for method in [
          "remark_section",
          "heterogen_section",
          "secondary_structure_section",] :
        section_lines = getattr(self._pdb_in.file_object.input, method)()
        for line in section_lines :
          f.write(line + "\n")
    reset_serial = self._serial_box.GetValue()
    if (reset_serial):
      self._hierarchy.atoms().reset_serial()
    pdb_out = self._hierarchy.as_pdb_string(
      crystal_symmetry=self._crystal_symmetry)
    for line in pdb_out.splitlines():
      if (not line.startswith("BREAK")):
        f.write("%s\n" % line)
    f.write("END")
    f.close()
    self._tree.SaveChanges()
    if (self._callback is not None):
      confirm_action(("Modified structure saved to %s.  Do you want "+
        "to update the file path in %s?") % (to_str(file_name),
                                             self._callback_program))
      self._callback(old_file_name=self._pdb_in.file_name,
                     new_file_name=to_str(file_name))
    else :
      wx.MessageBox("Modified structure saved to %s." % file_name)

  def OnEditSymmetry(self, event):
    dlg = symmetry_dialog.SymmetryDialog(parent=self,
      title="Edit model symmetry (CRYST1 record)",
      caption="Please enter the symmetry information for the PDB file; this "+
        "is optional or ignored in some cases (such as molecular replacement "+
        "search models), but required for many crystallography programs.")
    dlg.SetSymmetry(self._crystal_symmetry)
    if (dlg.ShowModal() == wx.ID_OK):
      old_symm = self._crystal_symmetry
      symm = dlg.GetSymmetry(allow_incomplete=True)
      #if (symm.space_group() is None) and (old_symm.space_group() is not None):
      self._crystal_symmetry = symm
    wx.CallAfter(dlg.Destroy)

  def OnEditModel(self, event):
    self._tree.ActionsForSelection(source_window=event.GetEventObject())

  def OnClose(self, event):
    if (self._pdb_in is not None) and (self._tree.HaveUnsavedChanges()):
      confirm = wx.MessageBox("You have unsaved changes - do you want to "+
        "save the modified PDB file now?", style=wx.YES_NO)
      if (confirm == wx.YES):
        self.Save()
    self.Destroy()

  def OnUndo(self, event):
    self._tree.Undo()

  def OnRedo(self, event):
    self._tree.Redo()

  def EnableRedo(self, enable=True):
    self.redo_item.Enable(enable)

  def EnableUndo(self, enable=True):
    self.undo_item.Enable(enable)

  def OnDestroy(self, event):
    pass

  def SetCallback(self, callback, program_name):
    self._callback = callback
    self._callback_program = program_name

def confirm_action(msg):
  confirm = wx.MessageBox(
    message=msg,
    style=wx.YES_NO|wx.ICON_QUESTION)
  if (confirm == wx.NO):
    raise Abort()
  return True
