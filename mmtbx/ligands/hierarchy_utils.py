from __future__ import absolute_import, division, print_function
import iotbx
from string import ascii_letters

from mmtbx.ligands.ready_set_basics import construct_xyz

def _new_atom(name, element, xyz, occ, b, hetero, segid=' '*4):
  # altloc???
  atom = iotbx.pdb.hierarchy.atom()
  atom.name = name
  atom.element = element
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
    chain = _add_atom_to_chain(atom, ag, icode=ba.parent().parent().icode)
    rc.append(chain)
  else:
    ag.append_atom(atom)
  return rc

def attempt_to_squash_alt_loc(hierarchy):
  indices = hierarchy.altloc_indices()
  altlocs = [_f for _f in indices if _f]
  if len(altlocs)==0: return hierarchy
  elif len(altlocs)>1: return None
  squash_hierarchy = hierarchy.deep_copy()
  for rg in squash_hierarchy.residue_groups():
    if len(rg.atom_groups())==1: continue
    ags = rg.atom_groups()
    detached_ag = ags[1].detached_copy()
    for atom in detached_ag.atoms():
      ags[0].append_atom(atom.detached_copy())
    rg.remove_atom_group(ags[1])
  return squash_hierarchy

def get_bonds_as_dict(geometry_restraints_manager, include_non_zero_origin_id=True):
  bonds={}
  for bond in geometry_restraints_manager.get_all_bond_proxies():
    if not hasattr(bond, 'get_proxies_with_origin_id'): continue
    for p in bond.get_proxies_with_origin_id():
      tmp=bonds.setdefault(p.i_seqs[0], [])
      tmp.append(p.i_seqs[1])
      tmp=bonds.setdefault(p.i_seqs[1], [])
      tmp.append(p.i_seqs[0])
    if include_non_zero_origin_id:
      for p in bond.get_proxies_without_origin_id(0):
        tmp=bonds.setdefault(p.i_seqs[0], [])
        tmp.append(p.i_seqs[1])
        tmp=bonds.setdefault(p.i_seqs[1], [])
        tmp.append(p.i_seqs[0])
  return bonds

def simple_valence_check(ph, geometry_restraints_manager):
  from mmtbx.ligands.chemistry import get_valences
  bonds = get_bonds_as_dict(geometry_restraints_manager.geometry)
  for atom in ph.atoms():
    if atom.element_is_hydrogen(): continue
    if atom.parent().resname in ['HOH']: continue
    number_of_bonds = len(bonds.get(atom.i_seq, None))
    # if number_of_bonds is None: continue
    # print(atom.quote(), number_of_bonds, get_valences(atom.element, atom.charge_as_int()))
    v = get_valences(atom.element, charge=atom.charge_as_int())
    if number_of_bonds not in v:
      print(atom.quote(), number_of_bonds, v)

def main(filename):
  from iotbx import pdb
  pdb_inp=pdb.input(filename)
  ph=pdb_inp.construct_hierarchy()
  print('is_hierarchy_altloc_consistent')
  print(ph.is_hierarchy_altloc_consistent())

if __name__ == '__main__':
  import sys
  main(*tuple(sys.argv[1:]))
