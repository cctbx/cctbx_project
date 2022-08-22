from __future__ import absolute_import, division, print_function
import iotbx
from string import ascii_letters

from mmtbx.ligands.ready_set_basics import construct_xyz

# class smart_add_atoms(list):
#   def __init__(self): pass

#   def append(self, item):
#     for chain1 in item:
#       remove = []
#       for atom1 in chain1.atoms():
#         for chain_list in self:
#           for chain2 in chain_list:
#             for atom2 in chain2.atoms():
#               if atom1.quote()==atom2.quote():
#                 remove.append(atom1)
#       if remove:
#         for atom in remove:
#           remove_atom_from_chain(chain1, atom)
#     list.append(self, item)

def _new_atom(name, element, xyz, occ, b, hetero, segid=' '*4):
  # altloc???
  atom = iotbx.pdb.hierarchy.atom()
  atom.name = name
  atom.element = "H"
  atom.xyz = xyz
  atom.occ = occ
  atom.b = b
  atom.hetero = hetero
  atom.segid = segid
  return atom

def _add_atom_to_chain(atom, ag, icode=None):
  rg = _add_atom_to_residue_group(atom, ag, icode=icode)
  chain = ag.parent().parent()
  tc = iotbx.pdb.hierarchy.chain()
  tc.id = chain.id
  tc.append_residue_group(rg)
  return tc

def _add_atom_to_residue_group(atom, ag, icode=None):
  tag = iotbx.pdb.hierarchy.atom_group()
  tag.resname = ag.resname
  tag.append_atom(atom)
  rg = iotbx.pdb.hierarchy.residue_group()
  rg.resseq = ag.parent().resseq
  if icode is not None: rg.icode=icode
  rg.append_atom_group(tag)
  for i, c in enumerate(ascii_letters):
    if c==ag.parent().parent().id:
      break
  atom.tmp = i
  return rg

def new_atom_with_inheritance(name, element, xyz, parent=None):
  occ=1
  b=20
  hetero=False
  if parent:
    occ=parent.occ
    b=parent.b
    hetero=parent.hetero
  return _new_atom(name, element, xyz, occ, b, hetero)

def add_hydrogens_to_atom_group_using_bad(ag,
                                          atom_name,
                                          atom_element,
                                          bond_atom,
                                          angle_atom,
                                          dihedral_atom,
                                          bond_length,
                                          angle,
                                          dihedral,
                                          append_to_end_of_model=False,
                                          ):
  rc = []
  if ag.get_atom(atom_name.strip()): return []
  if type(bond_atom)==type(''):
    ba = ag.get_atom(bond_atom.strip())
    if ba is None: return []
  else: ba = bond_atom
  if type(angle_atom)==type(''):
    aa = ag.get_atom(angle_atom.strip())
    if aa is None: return []
  else: aa = angle_atom
  if type(dihedral_atom)==type(''):
    da = ag.get_atom(dihedral_atom.strip())
    if da is None: return []
  else: da = dihedral_atom
  ro2 = construct_xyz(ba, bond_length,
                      aa, angle,
                      da, dihedral,
                      period=1,
                     )
  atom = _new_atom(atom_name, atom_element, ro2[0], ba.occ, ba.b, ba.hetero)
  if append_to_end_of_model:
    chain = _add_atom_to_chain(atom, ag)
    rc.append(chain)
  else:
    ag.append_atom(atom)
  return rc
