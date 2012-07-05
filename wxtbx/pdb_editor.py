
import wxtbx.app
import wxtbx.bitmaps
from wxtbx.phil_controls import simple_dialogs
from wxtbx import symmetry_dialog
from wxtbx import path_dialogs
from libtbx.utils import Abort, Sorry
from wx.lib.agw import customtreectrl
import wx
import string
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
    self._hierarchy = None

  def SetHierarchy (self, pdb_hierarchy) :
    self.DeleteAllItems()
    self._hierarchy = pdb_hierarchy
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
    self.ActionsForSelection(source_window=self)

  def OnDoubleClick (self, event) :
    pass
  #  self.ActionsForSelection()

  def ActionsForSelection (self, source_window) :
    if (self.flag_multiple_selections) :
      pass
    else :
      sel = self.GetSelection()
      pdb_object = self._node_lookup.get(sel, None)
      if (pdb_object is not None) :
        pdb_type = type(pdb_object).__name__
        if (pdb_type == "atom") :
          self.ShowAtomMenu(
            atom=pdb_object,
            source_window=source_window)
        elif (pdb_type == "atom_group") :
          self.ShowAtomGroupMenu(
            atom_group=pdb_object,
            source_window=source_window)
        elif (pdb_type == "residue_group") :
          self.ShowResidueGroupMenu(
            residue_group=pdb_object,
            source_window=source_window)
        elif (pdb_type == "chain") :
          self.ShowChainMenu(
            chain=pdb_object,
            source_window=source_window)
        elif (pdb_type == "model") :
          self.ShowModelMenu(
            model=pdb_object,
            source_window=source_window)
        else :
          raise RuntimeError("Unrecognized object type '%s'" % pdb_type)

  def ShowAtomMenu (self, atom, source_window) :
    labels_and_actions = [
      ("Set name...", self.OnSetName),
      ("Set occupancy...", self.OnSetOccupancy),
      ("Set B-factor...", self.OnSetBfactor),
      ("Set element...", self.OnSetElement),
      ("Set charge...", self.OnSetCharge),
      ("Delete atom", self.OnDeleteObject),
      ("Apply rotation/translation...", self.OnMoveSites),
    ]
    self.ShowMenu(labels_and_actions, source_window)

  def ShowAtomGroupMenu (self, atom_group, source_window) :
    labels_and_actions = [
      ("Set altloc...", self.OnSetAltloc),
      ("Set occupancy...", self.OnSetOccupancy),
      ("Set B-factor...", self.OnSetBfactor),
      ("Set residue name...", self.OnSetResname),
      ("Delete atom group", self.OnDeleteObject),
      ("Clone atom group", self.OnCloneAtoms),
      ("Apply rotation/translation...", self.OnMoveSites),
    ]
    if (atom_group.resname == "MET") :
      labels_and_actions.append(("Convert to SeMet", self.OnConvertMet))
    elif (atom_group.resname == "MSE") :
      labels_and_actions.append(("Convert to Met", self.OnConvertSeMet))
    self.ShowMenu(labels_and_actions, source_window)

  def ShowResidueGroupMenu (self, residue_group, source_window) :
    labels_and_actions = [
      ("Set residue number...", self.OnSetResseq),
      ("Set insertion code...", self.OnSetIcode),
      ("Delete residue", self.OnDeleteObject),
      ("Split residue", self.OnSplitResidue),
      ("Apply rotation/translation...", self.OnMoveSites),
    ]
    self.ShowMenu(labels_and_actions, source_window)

  def ShowChainMenu (self, residue_group, source_window) :
    labels_and_actions = [
      ("Set chain ID...", self.OnSetChainID),
      ("Set segment ID...", self.OnSetSegID),
      ("Delete chain", self.OnDeleteObject),
      ("Delete alternate conformers", self.OnDeleteAltConfs),
      ("Apply rotation/translation...", self.OnMoveSites),
    ]
    self.ShowMenu(labels_and_actions, source_window)

  def ShowModelMenu (self, residue_group, source_window) :
    labels_and_actions = [
      ("Set model ID...", self.OnSetChainID),
      ("Delete model", self.OnDeleteObject),
      ("Apply rotation/translation...", self.OnMoveSites),
    ]
    self.ShowMenu(labels_and_actions, source_window)

  def ShowMenu (self, items, source_window) :
    menu = wx.Menu()
    for label, action in items :
      item = menu.Append(-1, label)
      source_window.Bind(wx.EVT_MENU, action, item)
    source_window.PopupMenu(menu)
    menu.Destroy()

  def GetSelectedObject (self, object_type=None) :
    item = self.GetSelection()
    pdb_object = self._node_lookup.get(item, None)
    if (object_type is not None) :
      assert (type(pdb_object).__name__ == object_type)
    return item, pdb_object

  #---------------------------------------------------------------------
  # PROPERTY EDITING ACTIONS
  # atom
  def OnSetName (self, event) :
    item, atom = self.GetSelectedObject('atom')
    new_name = self.GetNewName(atom.name)
    assert (new_name is not None) and (1 <= len(new_name) <= 4)
    if (new_name != atom.name) :
      atom_group = atom.parent()
      for other_atom in atom_group.atoms() :
        if (other_atom.name == new_name) and (atom is not other_atom) :
          confirm_action(("The atom group to which this atom belongs already "+
          "has another atom named \"%s\".  Are you sure you want to rename "+
          "the selected atom?") % new_name)
      self._changes_made = True
      atom.name = "%-4s" % new_name
      self.SetItemText(item, format_atom(atom))

  # atom
  def OnSetElement (self, event) :
    item, atom = self.GetSelectedObject('atom')
    new_elem = self.GetNewElement(atom.element)
    assert (new_elem is None) or (len(new_elem) <= 2)
    if (new_elem != atom.element) :
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

  # atom
  def OnSetXYZ (self, event) :
    pass

  # atom_group
  def OnSetAltloc (self, event) :
    item, atom_group = self.GetSelectedObject('atom_group')
    new_altloc = self.GetNewAltloc(atom_group.altloc)
    assert (new_altloc is None) or (len(new_altloc) in [0,1])
    if (new_altloc != atom_group.altloc) :
      if (new_altloc in [None, '']) :
        rg = atom_group.parent()
        if (len(rg.atom_groups()) > 1) :
          confirm_action("You have specified a blank altloc ID for this atom "+
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
    if (atom_group.resname != new_resname) :
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
    item, residue_group = self.GetSelectedObject('residue_group')
    new_resseq = self.GetNewResseq(residue_group.resseq_as_int())
    assert (new_resseq is not None)
    if (new_resseq != residue_group.resseq_as_int()) :
      self._changes_made = True
      if (new_resseq > 9999) or (new_resseq < -999) :
        raise NotImplementedError("Hybrid36 support not available.")
      else :
        residue_group.resseq = "%4d" % new_resseq
      self.SetItemText(item, format_residue_group(residue_group))

  # residue_group
  def OnSetIcode (self, event) :
    item, residue_group = self.GetSelectedObject('residue_group')
    new_icode = self.GetNewIcode(residue_group.icode)
    assert (new_icode is None) or (len(new_icode) == 1)
    if (new_icode != residue_group.icode) :
      self._changes_made = True
      if (new_icode is None) :
        residue_group.icode = ' '
      else :
        residue_group.icode = new_icode
      self.SetItemText(item, format_residue_group(residue_group))

  # chain
  def OnSetChainID (self, event) :
    item, chain = self.GetSelectedObject('chain')
    new_id = self.GetNewChainID(chain.id)
    if (new_id != chain.id) :
      sefl._changes_made = True
      if (new_id is None) or (new_id.isspace()) :
        chain.id = ' '
      else :
        chain.id = "%2s" % new_id

  # chain
  def OnRenumber (self, event) :
    item, chain = self.GetSelectedObject('chain')
    resseq_shift = self.GetResseqShift()
    if (resseq_shift is not None) and (resseq_shift != 0) :
      self._changes_made = True
      child, cookie = self.GetFirstChild(item)
      while (child is not None) :
        residue_group = self._node_lookup.get(child, None)
        assert (type(residue_group).__name__ == 'residue_group')
        resseq = residue_group.resseq_as_int()
        new_resseq = resseq + resseq_shift
        if (new_resseq > 9999) or (new_resseq < -999) :
          raise NotImplementedError("Hybrid36 support not available.")
        else :
          residue_group.resseq = "%4d" % new_resseq
        self.SetItemText(child, format_residue_group(residue_group))
        child, cookie = self.GetNextChild(node, cookie)

  # chain
  def OnDeleteAltConfs (self, event) :
    item, chain = self.GetSelectedObject('chain')
    n_alt_atoms = 0
    n_alt_residues = 0
    for residue_group in chain.residue_groups() :
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1) :
        n_alt_residues += 1
        for ag in atom_groups[1:] :
          n_alt_atoms += len(ag.atoms())
    if (len(n_alt_atoms) == 0) :
      raise Sorry("No alternate conformations found in this chain.")
    confirm_action(("There are %d residues with alternate conformations; "+
      "removing these will delete %d atoms from the model.  Are you sure "+
      "you want to continue?") % (n_alt_residues, n_alt_atoms))
    # TODO more control over what happens to remaining atom_groups
    child, cookie = self.GetFirstChild(item)
    while (child is not None) :
      residue_group = self._node_lookup.get(child, None)
      assert (type(residue_group).__name__ == 'residue_group')
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1) :
        first_group = atom_groups[0]
        first_group.altloc = ''
        for atom in first_group.atoms() :
          atom.occ = 1.0
        child2, cookie2 = self.GetFirstChild(child)
        self.SetItemText(child2, format_atom_group(first_group))
        for atom_group in atom_groups[1:] :
          residue_group.remove_atom_group(atom_group)
          for item2, pdb_object in self._node_lookup.iteritems() :
            if (pdb_object is atom_group) :
              self.DeleteChildren(item2)
              self.DeleteItem(item2)
        self.SetItemText(child, format_residue_group(residue_group))
      child, cookie = self.GetNextchild(item, cookie)

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

  # all
  def OnSetSegID (self, event) :
    if (self.flag_multiple_selections) :
      assert 0
    else :
      item = self.GetSelection()
      pdb_object = self._node_lookup.get(item, None)
      pdb_type = type(pdb_object).__name__
      segid = None
      if (pdb_type == 'atom') :
        segid = pdb_object.segid
      else :
        from scitbx.array_family import flex
        segids = flex.std_string([ a.segid for a in pdb_object.atoms() ])
        if (segids.all_eq(segids[0])) :
          segid = segids[0]
      if (segid is not None) and (segid.isspace()) :
        segid = None
      new_segid = self.GetNewSegID(segid)
      if (new_segid != segid) :
        self._changes_made = True
        if (pdb_type == 'atom') :
          pdb_object.segid = new_segid
          self.SetItemText(item, format_atom(pdb_object))
        else :
          for atom in pdb_object.atoms() :
            atom.segid = new_segid
          self.PropagateAtomChanges(item)

  # all
  def OnMoveSites (self, event) :
    if (self.flag_multiple_selections) :
      assert 0
    else :
      item = self.GetSelection()
      pdb_object = self._node_lookup.get(item, None)
      pdb_type = type(pdb_object).__name__
      rt = simple_dialogs.get_rt_matrix(self)
      if (pdb_type == 'atom') :
        from scitbx.array_family import flex
        v = flex.vec3_double([pdb_object.xyz])
        v = rt.r.elems * v + rt.t.elems
        pdb_object.xyz = v[0]
        self.SetItemText(item, format_atom(pdb_object))
      else :
        atoms = pdb_object.atoms()
        sites = atoms.extract_xyz()
        sites = rt.r.elems * v + rt.t.elems
        atoms.set_xyz(sites)
        self.PropagateAtomChanges(item)

  #---------------------------------------------------------------------
  # HIERARCHY EDITING
  def OnDeleteObject (self, event) :
    if (self.flag_multiple_selections) :
      assert 0
    else :
      item = self.GetSelection()
      pdb_object = self._node_lookup.get(item, None)
      pdb_type = type(pdb_object).__name__
      n_atoms = None
      if (pdb_type == 'atom') :
        n_atoms = 1
      else :
        n_atoms = len(pdb_object.atoms())
      confirm_action("Are you sure you want to delete the selected %d atom(s)?"
        % n_atoms)
      parent = pdb_object.parent()
      clean_up = True
      if (pdb_type == 'atom') :
        parent.remove_atom(pdb_object)
      elif (pdb_type == 'atom_group') :
        parent.remove_atom_group(pdb_object)
      elif (pdb_type == 'residue_group') :
        parent.remove_residue_group(pdb_object)
      elif (pdb_type == 'chain') :
        parent.remove_chain(pdb_object)
      else :
        if (len(self._hierarchy.models()) == 1) :
          raise Sorry("You must have at least one MODEL in the PDB file.")
        parent.remove_model(pdb_object)
      parent_item = self.GetItemParent(item)
      self.DeleteChildren(item)
      self.DeleteItem(item)
      while (not self.HasChildren(parent_item)) :
        del_item = parent_item
        parent_itme = self.GetItemParent(del_item)
        self.DeleteItem(del_item)

  # model
  def OnAddChain (self, event) :
    pass

  # chain
  def OnAddResidues (self, event) :
    pass

  # residue_group
  def OnSplitResidue (self, event) :
    print 1
    item, residue_group = self.GetSelectedObject('residue_group')
    atom_groups = residue_group.atom_groups()
    self._changes_made = True
    assert (len(atom_groups) > 0)
    start_occ = 1/(len(atom_groups) + 1)
    new_occ = self.GetNewOccupancy(start_occ, new=True)
    new_group = atom_groups[0].detached_copy()
    for atom in new_group.atoms() :
      atom.occ = new_occ
    for atom_group in atom_groups :
      for atom in atom_group.atoms() :
        atom.occ = max(0, atom.occ - new_occ/len(atom_groups))
    if (len(atom_groups) == 1) :
      atom_groups[0].altloc = 'A'
      new_group.altloc = 'B'
    else :
      new_altloc = None
      for char in string.uppercase :
        for atom_group in atom_groups :
          if (atom_group.altloc == char) :
            break
        else :
          new_altloc = char
      new_group.altloc = char
    residue_group.append_atom_group(new_group)
    ag_item = self.AppendItem(item, format_atom_group(new_group))
    self._node_lookup[ag_item] = new_group
    for atom in new_group.atoms() :
      new_item = self.AppendItem(ag_item, format_atom(atom))
      self._node_lookup[new_item] = atom
    child, cookie = self.GetFirstChild(item)
    while (child is not None) :
      atom_group = self._node_lookup.get(child, None)
      assert (type(atom_group).__name__ == 'atom_group')
      self.SetItemText(child, format_atom_group(atom_group))
      child, cookie = self.GetNextChild(item, cookie)

  # atom_group
  def OnCloneAtoms (self, event) :
    pass

  #---------------------------------------------------------------------
  # USER INPUT
  def GetNewName (self, name=None) :
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

  def GetNewOccupancy (self, occ=None, new=False) :
    desc_str = "selected"
    if (new) :
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
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewCharge (self, charge=None) :
    if (charge is not None) and (charge.isspace()) :
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
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewElement (self, elem=None) :
    if (elem is not None) and (elem.isspace()) :
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

  def GetNewResseq (self, resseq) :
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

  def GetNewIcode (self, icode) :
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
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewChainID (self, chain_id) :
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
      value=chain_id)
    #dlg.SetMinLength(1)
    dlg.SetMaxLength(2)
    return simple_dialogs.get_phil_value_from_dialog(dlg)

  def GetNewSegID (self, segid) :
    if (segid.isspace()) :
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

ADD_RESIDUES_START = 0
ADD_RESIDUES_END = 1
ADD_RESIDUES_AFTER = 2
class AddResiduesDialog (wx.Dialog) :
  def __init__ (self, *args, **kwds) :
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
    szr.Add(caption_txt, 0, wx.ALL, 5)
    box = wx.BoxSizer(wx.HORIZONTAL)
    txt1 = wx.StaticText(self, -1, "Insert residues:")
    box.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._insert_choice = wx.Choice(self, -1,
      choices=["at the end of the chain",
               "at the start of the chain",
               "after residue"])
    self._insert_choice.SetSelection(0)
    self.Bind(wx.EVT_CHOICE, self.OnChangeInsertion, self._insert_choice)
    box.Add(self._insert_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    txt2 = wx.StaticText(self, -1, "residue ID:")
    box.Add(txt2, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self._resid_txt = phil_controls.strctrl(
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

  def OnChangeInsertion (self, event) :
    insert_type = self._insert_choice.GetSelection()
    if (insert_type == ADD_RESIDUES_AFTER) :
      self._resid_txt.SetOptional(False)
      self._resid_txt.Enable()
    else :
      self._resid_txt.SetOptional(True)
      self._resid_txt.Enable(False)

  def GetInsertionType (self) :
    return self._insert_choice.GetSelection()

  def GetResID (self) :
    return self._resid_txt.GetPhilValue()

########################################################################
class PDBTreeFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.statusbar = self.CreateStatusBar()
    self.statusbar.SetStatusText("Right-click any item for a list of editing actions")
    # toolbar setup
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    bmp = wxtbx.bitmaps.fetch_custom_icon_bitmap("phenix.pdbtools")
    btn = self.toolbar.AddLabelTool(-1, "Load file", bmp,
      shortHelp="Load file", kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnOpen, btn)
    bmp = wxtbx.bitmaps.fetch_icon_bitmap("actions", "save_all")
    btn = self.toolbar.AddLabelTool(-1, "Save file", bmp,
      shortHelp="Save file", kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnSave, btn)
    bmp = wxtbx.bitmaps.fetch_custom_icon_bitmap("tools")
    btn = self.toolbar.AddLabelTool(-1, "Edit...", bmp, shortHelp="Edit...",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnEditModel, btn)
    bmp = wxtbx.bitmaps.fetch_custom_icon_bitmap("symmetry")
    btn = self.toolbar.AddLabelTool(-1, "Symmetry", bmp, shortHelp="Symmetry",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnEditSymmetry, btn)
    self.toolbar.Realize()
    #
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
    self._crystal_symmetry = None
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
    self._crystal_symmetry = f.file_object.crystal_symmetry()
    hierarchy = f.file_object.construct_hierarchy()
    self._hierarchy = hierarchy
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
    if (self._pdb_in is None) :
      return
    file_name = self.path_mgr.select_file(
      parent=self,
      message="Save PDB file",
      style=wx.SAVE,
      wildcard=file_reader.get_wildcard_strings(["pdb"]),
      current_file=self._pdb_in.file_name)
    save_header = self._header_box.GetValue()
    f = open(file_name, "w")
    if (save_header) :
      for method in [
          "title_section",
          "remark_section",
          "heterogen_section",
          "secondary_structure_section",] :
        section_lines = getattr(self._pdb_in.file_object, method)()
        for line in section_lines :
          f.write(line + "\n")
    f.write(self._hierarchy.as_pdb_string(
      crystal_symmetry=self._crystal_symmetry))
    f.close()
    wx.MessageBox("Modified structure saved to %s." % file_name)
    self._tree.SaveChanges()

  def OnEditSymmetry (self, event) :
    dlg = symmetry_dialog.SymmetryDialog(parent=self,
      title="Edit model symmetry (CRYST1 record)",
      caption="Please enter the symmetry information for the PDB file; this "+
        "is optional or ignored in some cases (such as molecular replacement "+
        "search models), but required for many crystallography programs.")
    dlg.SetSymmetry(self._crystal_symmetry)
    if (dlg.ShowModal() == wx.ID_OK) :
      old_symm = self._crystal_symmetry
      symm = dlg.GetSymmetry(allow_incomplete=True)
      #if (symm.space_group() is None) and (old_symm.space_group() is not None) :
      self._crystal_symmetry = symm
    wx.CallAfter(dlg.Destroy)

  def OnEditModel (self, event) :
    self._tree.ActionsForSelection(source_window=event.GetEventObject())

  def OnClose (self, event) :
    if (self._tree.HaveUnsavedChanges()) :
      pass
    self.Destroy()

  def OnDestroy (self, event) :
    pass

def confirm_action (msg) :
  confirm = wx.MessageBox(
    message=msg,
    style=wx.YES_NO)
  if (confirm == wx.NO) :
    raise Abort()
  return True

if (__name__ == "__main__") :
  app = wxtbx.app.CCTBXApp(0)
  frame = PDBTreeFrame(None, -1, "PDB Editor")
  frame.Show()
  if (len(sys.argv) > 1) :
    frame.LoadPDB(sys.argv[1])
  app.MainLoop()
