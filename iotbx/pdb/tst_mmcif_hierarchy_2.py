"""Test mmCIF labels"""
from __future__ import absolute_import, division, print_function

import libtbx.load_env

import iotbx.pdb
import os.path

import libtbx.load_env

def exercise_cif_label_asym_id_simple():
  # 1 chain, no water
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/p9.pdb", test=os.path.isfile)
  h = iotbx.pdb.input(pdb_fname).construct_hierarchy()
  n_atoms = h.atoms_size()
  for i in range(n_atoms):
    assert h.get_label_asym_id_iseq(i) == 'A'
  for rg in h.only_model().only_chain().residue_groups():
    asym_id = h.get_label_asym_id(rg)
    if asym_id != 'A':
      print ('Error:', rg.id_str(), asym_id)
      assert 0
  for i in range(10):
    assert h.get_label_asym_id_iseq(10) == 'A' # SIC! getting for the same atom!

def exercise_cif_label_asym_id_one_chain_ligand_water():
  '''1 chain, ligand, water
  waters are the same label_asym_id if in the same chain
  each ligand gets its own label_asym_id
  '''
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/one_chain_ligand_water.pdb",
    test=os.path.isfile)

  h = iotbx.pdb.input(pdb_fname).construct_hierarchy()
  assert h.only_model().chains_size() == 1
  # 83 protein atoms
  for i in range(84):
    assert h.get_label_asym_id_iseq(i) == 'A'
  for i in range(84, 96):
    assert h.get_label_asym_id_iseq(i) == 'B'
  assert h.get_label_asym_id_iseq(96) == 'C'
  assert h.get_label_asym_id_iseq(97) == 'C'

  # cb = h.as_cif_block()
  # cb.show()


def exercise_cif_label_asym_id_two_chains_ligands_waters():
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/two_chains_ligand_water.pdb",
    test=os.path.isfile)
  h = iotbx.pdb.input(pdb_fname).construct_hierarchy()
  # print ('n chains:', h.only_model().chains_size())
  # for a in h.atoms():
  #   print (a.i_seq, a.id_str())
  # 83 protein atoms
  for i in range(26):
    assert h.get_label_asym_id_iseq(i) == 'A'
  # CA
  assert h.get_label_asym_id_iseq(26) == 'B'
  # BEN - 1 residue in AC
  for i in range(27, 45):
    assert h.get_label_asym_id_iseq(i) == 'C', i
  # SO4 - first
  for i in range(45, 50):
    assert h.get_label_asym_id_iseq(i) == 'D'
  # SO4 - second
  for i in range(50, 55):
    assert h.get_label_asym_id_iseq(i) == 'E'
  # GOL
  for i in range(55, 61):
    assert h.get_label_asym_id_iseq(i) == 'F'
  # water
  for i in range(61, 65):
    assert h.get_label_asym_id_iseq(i) == 'G'

  # cb = h.as_cif_block()
  # cb.show()

def exercise_cif_label_asym_id_two_chains():
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/one_chain_ligand_water.pdb",
    test=os.path.isfile)
  with open(pdb_fname, 'r') as f:
    full_model = f.readlines()
  second_chain = (''.join(full_model[5:])).replace(' A ', ' B ')
  full_model = ''.join(full_model) + second_chain
  h = iotbx.pdb.input(lines=full_model, source_info=None).construct_hierarchy()
  # print ('n chains:', h.only_model().chains_size())
  # for a in h.atoms():
  #   print (a.i_seq, a.id_str())
  for i in range(84):
    assert h.get_label_asym_id_iseq(i) == 'A'
  for i in range(84, 96):
    assert h.get_label_asym_id_iseq(i) == 'B'
  assert h.get_label_asym_id_iseq(96) == 'C'
  assert h.get_label_asym_id_iseq(97) == 'C'

  for i in range(98, 182):
    assert h.get_label_asym_id_iseq(i) == 'D'
  for i in range(182, 194):
    assert h.get_label_asym_id_iseq(i) == 'E'
  assert h.get_label_asym_id_iseq(194) == 'F'
  assert h.get_label_asym_id_iseq(195) == 'F'
  # cb = h.as_cif_block()
  # cb.show()

def exercise_cif_label_asym_id_ligand_before_protein():
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/ligand_before_prot.pdb",
    test=os.path.isfile)
  h = iotbx.pdb.input(pdb_fname).construct_hierarchy()
  # for i in range(h.atoms_size()):
  #   print (i, h.get_label_asym_id_iseq(i))
  h.get_label_asym_id_iseq(0) == 'A'
  h.get_label_asym_id_iseq(1) == 'B'
  h.get_label_asym_id_iseq(2) == 'C'
  for i in range(3, 23):
    assert h.get_label_asym_id_iseq(i) == 'D'

if __name__ == '__main__':
  exercise_cif_label_asym_id_simple()
  exercise_cif_label_asym_id_one_chain_ligand_water()
  exercise_cif_label_asym_id_two_chains_ligands_waters()
  exercise_cif_label_asym_id_two_chains()
  exercise_cif_label_asym_id_ligand_before_protein()
  print ('OK')
