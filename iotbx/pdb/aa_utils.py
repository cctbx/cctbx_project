"""Utilities for working with atom names"""
from __future__ import division, print_function
from iotbx.pdb import modified_aa_names
from iotbx.pdb.amino_acid_codes import three_letter_given_one_letter
from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter

from mmtbx.chemical_components import get_atom_names, get_bond_pairs

def get_aa_parent(code):
  one = modified_aa_names.lookup.get(code.upper(), False)
  if not one: return None
  return three_letter_given_one_letter.get(one, None)

def get_aa_children(code):
  one = one_letter_given_three_letter.get(code, None)
  if one is None: return None
  if one not in modified_aa_names.lookup.values(): return None
  rc=[]
  for key, item in modified_aa_names.lookup.items():
    if item==one: rc.append(key)
  rc.sort()
  return rc

def standard_polymerise(code):
  names = get_atom_names(code, heavy_atom_only=True)
  names = set(names)
  return set(['N','CA','C','O']).issubset(names)

def compare_atom_names(child, parent):
  c_names = set(get_atom_names(child))
  p_names = get_atom_names(parent, heavy_atom_only=True)
  p_names.remove('OXT')
  return set(p_names).issubset(c_names)

def _is_specific_number_heavy_bonds(bonds, name, number, verbose=False):
  n_bonds = 0
  for bond in bonds:
    if verbose: print(bond)
    if 'OXT' in bond: continue
    if name in bond:
      n_bonds+=1
    if verbose: print(n_bonds)
  return n_bonds!=number

def is_n_terminal(bonds):
  return _is_specific_number_heavy_bonds(bonds, 'N', 1)

def is_c_terminal(bonds):
  return _is_specific_number_heavy_bonds(bonds, 'C', 2)

def is_ca_mod(bonds, n_bonds=3):
  return _is_specific_number_heavy_bonds(bonds, 'CA', n_bonds)

def is_standard_bonding(bonds, include_cb=True):
  standard = [('CA', 'N'), ('C', 'CA'), ('C', 'O')]
  if include_cb:
    standard.append(('CA', 'CB'))
  return set(standard).issubset(set(bonds))

def is_not_standard_main_chain(code):
  assert standard_polymerise(code)
  bonds = get_bond_pairs(code, heavy_atom_only=True, use_tuple=True)
  if is_n_terminal(bonds):
    return 'N terminal'
  elif is_c_terminal(bonds):
    return 'C terminal'
  elif is_ca_mod(bonds):
    return 'CA mod'
  elif not is_standard_bonding(bonds):
    return 'bonding mod'
  return False

def get_aa_type(code):
  parent = get_aa_parent(code)
  if not standard_polymerise(code):
    return 'non-polymer'
  exact_subset = True
  if not compare_atom_names(code, parent):
    exact_subset = False
  standard_main_chain = True
  if is_not_standard_main_chain(code):
    standard_main_chain = False
    return is_not_standard_main_chain(code)
  return 'ok'

def get_useable_sorted_on_size(code):
  rc = get_aa_children(code)
  tmp = []
  overall = {}
  for code in rc:
    aa_type = get_aa_type(code)
    overall.setdefault(aa_type, [])
    overall[aa_type].append(code)
    if aa_type in ['ok']:
      tmp.append(code)
  def myFunc(e): return len(get_atom_names(e, heavy_atom_only=True))
  tmp.sort(key=myFunc)
  return tmp

def tst_types(code):
  rc = get_aa_children(code)
  overall = {}
  for code in rc:
    aa_type = get_aa_type(code)
    overall.setdefault(aa_type, [])
    overall[aa_type].append(code)
  print(overall)
  for key, item in overall.items():
    print(key, item)

if __name__ == '__main__':
  import os, sys
  if len(sys.argv)>1:
    aa_list = sys.argv[1:]
  else:
    aa_list = ['CYS', 'DAL', 'NWM', 'GLY']
  for resname in aa_list:
    rc = get_aa_parent(resname)
    print('  Parent %s : %s' % (resname, rc))
    rc = get_aa_children(resname)
    print('  Children %s : %s' % (resname, rc))
    print(get_useable_sorted_on_size(resname))
    tst_types(resname)
    continue
    for code in rc:
      print(code,get_aa_type(code))
    #   cmd = 'phenix.reel --chem %s' % code
    #   os.system(cmd)
