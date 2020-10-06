from __future__ import absolute_import, division, print_function
import iotbx
from mmtbx.ligands.ready_set_basics import construct_xyz
import six
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

def generate_atom_group_atom_names(rg, names, return_Nones=False):
  '''
  Generate all alt. loc. groups of names
  '''
  atom_groups = rg.atom_groups()
  atom_altlocs = {}
  for ag in atom_groups:
    for atom in ag.atoms():
      atom_altlocs.setdefault(atom.parent().altloc, [])
      atom_altlocs[atom.parent().altloc].append(atom)
  if len(atom_altlocs)>1 and '' in atom_altlocs:
    for key in atom_altlocs:
      if key=='': continue
      for atom in atom_altlocs['']:
        atom_altlocs[key].append(atom)
    del atom_altlocs['']
  for key, value in six.iteritems(atom_altlocs):
    atoms=[]
    for name in names:
      for atom in value:
        if atom.name.strip()==name.strip():
          atoms.append(atom)
          break
      else:
        if return_Nones:
          atoms.append(None)
        else:
          print('not all atoms found. missing %s from %s' % (name, names))
          break
    if len(atoms)!=len(names):
      yield None, None
    else:
      yield atoms[0].parent(), atoms

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
