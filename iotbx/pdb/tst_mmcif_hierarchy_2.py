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

pdb_str_unknown_residue_in_chain = """\
ATOM      1  N   GLY C   4     -21.642  55.906  20.798  1.00 67.42           N
ATOM      2  CA  GLY C   4     -20.299  56.304  20.417  1.00 57.03           C
ATOM      3  C   GLY C   4     -19.307  55.155  20.368  1.00 49.28           C
ATOM      4  O   GLY C   4     -18.458  55.105  19.483  1.00 48.20           O
HETATM    5  O3  HP3 C   5     -14.220  56.995  18.819  1.00 44.52           O
HETATM    6  O4  HP3 C   5     -12.986  54.899  19.271  1.00 37.25           O
HETATM    7  P   HP3 C   5     -13.657  56.096  19.903  1.00 43.88           P
HETATM    8  O1  HP3 C   5     -12.687  56.861  20.781  1.00 39.50           O
HETATM    9  CE2 HP3 C   5     -15.039  55.553  20.898  1.00 47.32           C
HETATM   10  NE1 HP3 C   5     -15.776  56.425  21.622  1.00 39.88           N
HETATM   11  ND1 HP3 C   5     -16.726  55.591  22.236  1.00 50.21           N
HETATM   12  CD2 HP3 C   5     -15.404  54.221  20.962  1.00 45.89           C
HETATM   13  NG  HP3 C   5     -16.446  54.278  21.799  1.00 47.37           N
HETATM   14  CB  HP3 C   5     -17.339  53.258  22.329  1.00 38.02           C
HETATM   15  CA  HP3 C   5     -18.519  53.059  21.378  1.00 48.49           C
HETATM   16  N   HP3 C   5     -19.414  54.234  21.320  1.00 45.36           N
HETATM   17  C   HP3 C   5     -19.284  51.855  21.843  1.00 53.75           C
HETATM   18  O   HP3 C   5     -20.297  51.977  22.516  1.00 43.83           O
ATOM     19  N   ALA C   6     -18.767  50.683  21.488  1.00 52.53           N
ATOM     20  CA  ALA C   6     -19.317  49.413  21.937  1.00 51.73           C
ATOM     21  C   ALA C   6     -18.996  49.182  23.422  1.00 50.17           C
ATOM     22  O   ALA C   6     -18.163  49.879  23.996  1.00 44.08           O
ATOM     23  CB  ALA C   6     -18.768  48.277  21.078  1.00 42.73           C
ATOM     24  N   GLY C   7     -19.666  48.212  24.044  1.00 53.16           N
ATOM     25  CA  GLY C   7     -19.394  47.871  25.429  1.00 57.47           C
ATOM     26  C   GLY C   7     -19.954  48.819  26.469  1.00 55.41           C
ATOM     27  O   GLY C   7     -19.599  48.696  27.649  1.00 49.98           O
ATOM     28  N   ALA C   8     -20.814  49.758  26.077  1.00 54.83           N
ATOM     29  CA  ALA C   8     -21.403  50.715  27.012  1.00 55.88           C
ATOM     30  C   ALA C   8     -22.659  50.143  27.663  1.00 61.52           C
ATOM     31  O   ALA C   8     -22.687  49.901  28.870  1.00 59.53           O
ATOM     32  CB  ALA C   8     -21.723  52.026  26.303  1.00 56.85           C
"""

pdb_str_dna = """\
ATOM      1  P    DA A   1       0.000   0.000   0.000  1.00 20.00           P
ATOM      2  O3'  DA A   1       1.000   0.000   0.000  1.00 20.00           O
ATOM      3  P    DC A   2       2.500   0.000   0.000  1.00 20.00           P
ATOM      4  O3'  DC A   2       3.500   0.000   0.000  1.00 20.00           O
ATOM      5  P   5MC A   3       5.000   0.000   0.000  1.00 20.00           P
"""

def get_label_ids(h):
  rgs = h.only_chain().residue_groups()
  asym_ids = [h.get_label_asym_id(rg) for rg in rgs]
  seq_ids = [h.get_label_seq_id(rg.atom_groups()[0]) for rg in rgs]
  return asym_ids, seq_ids

def exercise_cif_label_asym_id_unknown_residue_in_chain():
  '''A residue with a non-CCD name (HP3, from a user restraints file) that is
  peptide-bonded in the middle of a chain is part of the polymer: whole chain
  keeps one label_asym_id and the residue gets a label_seq_id.
  (wwPDB report: depositor-named residues mid-chain split label_asym_id)'''
  h = iotbx.pdb.input(
    lines=pdb_str_unknown_residue_in_chain, source_info=None).construct_hierarchy()
  asym_ids, seq_ids = get_label_ids(h)
  assert asym_ids == ['A']*5, asym_ids
  assert seq_ids == ['1', '2', '3', '4', '5'], seq_ids
  # round trip through mmCIF must keep a single chain
  h2 = iotbx.pdb.input(
    lines=h.as_mmcif_string(), source_info=None).construct_hierarchy()
  assert h2.only_model().chains_size() == 1, h2.only_model().chains_size()

def exercise_cif_label_asym_id_unknown_residue_not_linked():
  '''The same unknown residue moved away from the chain is a ligand: it gets
  its own label_asym_id and no label_seq_id.'''
  h = iotbx.pdb.input(
    lines=pdb_str_unknown_residue_in_chain, source_info=None).construct_hierarchy()
  for atom in h.only_chain().residue_groups()[1].atoms():
    x, y, z = atom.xyz
    atom.set_xyz(new_xyz=(x+20, y, z))
  asym_ids, seq_ids = get_label_ids(h)
  assert asym_ids == ['A', 'B', 'C', 'C', 'C'], asym_ids
  assert seq_ids == ['1', '.', '3', '4', '5'], seq_ids

def exercise_cif_label_asym_id_d_amino_acids():
  '''D-amino acids (residue class d_amino_acid) are polymer residues'''
  h = iotbx.pdb.input(
    lines=pdb_str_unknown_residue_in_chain, source_info=None).construct_hierarchy()
  for ag in h.atom_groups():
    ag.resname = 'DAL'
  asym_ids, seq_ids = get_label_ids(h)
  assert asym_ids == ['A']*5, asym_ids
  assert seq_ids == ['1', '2', '3', '4', '5'], seq_ids

def exercise_cif_label_seq_id_nucleic_acid():
  '''Nucleic acid residues get label_seq_id like amino acids'''
  h = iotbx.pdb.input(lines=pdb_str_dna, source_info=None).construct_hierarchy()
  asym_ids, seq_ids = get_label_ids(h)
  assert asym_ids == ['A']*3, asym_ids
  assert seq_ids == ['1', '2', '3'], seq_ids

if __name__ == '__main__':
  exercise_cif_label_asym_id_simple()
  exercise_cif_label_asym_id_one_chain_ligand_water()
  exercise_cif_label_asym_id_two_chains_ligands_waters()
  exercise_cif_label_asym_id_two_chains()
  exercise_cif_label_asym_id_ligand_before_protein()
  exercise_cif_label_asym_id_unknown_residue_in_chain()
  exercise_cif_label_asym_id_unknown_residue_not_linked()
  exercise_cif_label_asym_id_d_amino_acids()
  exercise_cif_label_seq_id_nucleic_acid()
  print ('OK')
