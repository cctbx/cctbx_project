import boost.python
ext = boost.python.import_ext("iotbx_pdb_hierarchy_v2_ext")
from iotbx_pdb_hierarchy_v2_ext import *

class _root(boost.python.injector, ext.root):

  def only_model(self):
    assert self.models_size() == 1
    return self.models()[0]

  def only_chain(self):
    return self.only_model().only_chain()

  def only_residue_group(self):
    return self.only_chain().only_residue_group()

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_atom(self):
    return self.only_atom_group().only_atom()

class _model(boost.python.injector, ext.model):

  def only_chain(self):
    assert self.chains_size() == 1
    return self.chains()[0]

  def only_residue_group(self):
    return self.only_chain().only_residue_group()

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_atom(self):
    return self.only_atom_group().only_atom()

class _chain(boost.python.injector, ext.chain):

  def only_residue_group(self):
    assert self.residue_groups_size() == 1
    return self.residue_groups()[0]

  def only_atom_group(self):
    return self.only_residue_group().only_atom_group()

  def only_atom(self):
    return self.only_atom_group().only_atom()

class _residue_group(boost.python.injector, ext.residue_group):

  def only_atom_group(self):
    assert self.atom_groups_size() == 1
    return self.atom_groups()[0]

  def only_atom(self):
    return self.only_atom_group().only_atom()

  def move_blank_altloc_atom_groups_to_front(self):
    blank_altloc_char = ' '
    n_blank_altloc_atom_groups = 0
    for i_atom_group,atom_group in enumerate(self.atom_groups()):
      if (atom_group.altloc in ["", blank_altloc_char]):
        if (i_atom_group != n_blank_altloc_atom_groups):
          self.remove_atom_group(
            i=i_atom_group)
          self.insert_atom_group(
            i=n_blank_altloc_atom_groups, atom_group=atom_group)
        n_blank_altloc_atom_groups += 1
    return n_blank_altloc_atom_groups

  def process_blank_altloc(self):
    blank_altloc_char = ' '
    n_blank_altloc_atom_groups = self.move_blank_altloc_atom_groups_to_front()
    if (n_blank_altloc_atom_groups == 0):
      return (0, 0)
    atom_groups = self.atom_groups()
    blank_name_set = set()
    for atom_group in atom_groups[:n_blank_altloc_atom_groups]:
      atom_group.altloc = ""
      for atom in atom_group.atoms():
        blank_name_set.add((atom_group.resname, atom.name))
    blank_but_alt = set()
    for atom_group in atom_groups[n_blank_altloc_atom_groups:]:
      for atom in atom_group.atoms():
        if ((atom_group.resname, atom.name) in blank_name_set):
          blank_but_alt.add((atom_group.resname, atom.name))
    n_blank_but_alt_atom_groups = 0
    if (len(blank_but_alt) != 0):
      n_atom_groups_removed = 0
      for i_atom_group in xrange(n_blank_altloc_atom_groups):
        atom_group = atom_groups[i_atom_group]
        new_atom_group = None
        n_atoms_removed = 0
        for i_atom,atom in enumerate(atom_group.atoms()):
          if ((atom_group.resname, atom.name) in blank_but_alt):
            atom_group.remove_atom(i=i_atom-n_atoms_removed)
            n_atoms_removed += 1
            if (new_atom_group is None):
              new_atom_group = ext.atom_group(
                altloc=blank_altloc_char, resname=atom_group.resname)
              self.insert_atom_group(
                i=n_blank_altloc_atom_groups+n_blank_but_alt_atom_groups,
                atom_group=new_atom_group)
              n_blank_but_alt_atom_groups += 1
            new_atom_group.append_atom(atom)
        if (atom_group.atoms_size() == 0):
          self.remove_atom_group(i=i_atom_group-n_atom_groups_removed)
          n_atom_groups_removed += 1
          n_blank_altloc_atom_groups -= 1
    return (n_blank_altloc_atom_groups, n_blank_but_alt_atom_groups)

class _atom_group(boost.python.injector, ext.atom_group):

  def only_atom(self):
    assert self.atoms_size() == 1
    return self.atoms()[0]
