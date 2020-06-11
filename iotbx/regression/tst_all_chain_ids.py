from __future__ import absolute_import, division, print_function
from iotbx.regression.ncs.tst_mtrix_biomt_cmdl import pdb_str_0
import mmtbx.model
import iotbx.pdb

def test_1():
  """
  Chains A,B,C are present in the model. The rest of them are being generated
  by expansion. The problem is that iotbx/pdb/utils.py:all_chain_ids()
  generate chain ids with leading white space for single character ids.
  This results in inconsistent selection behavior demonstrated in this test.

  Many instruments quote chain ID when generate selection strings such as:
  iotbx/pdb/secondary_structure.py
  nqh_flip in mmtbx/validation/clashscore.py

  maybe others

  Therefore they don't work properly when the model is expanded using MTRIX
  and/or BIOMT matrices internally.

  To reproduce nqh_flip problem:
  phenix.fetch_pdb 1qgt
  phenix.molprobity 1qgt.pdb
  """
  # Loading and expanding the model
  inp = iotbx.pdb.input(lines=pdb_str_0, source_info=None)
  model = mmtbx.model.manager(model_input=inp)
  model.expand_with_BIOMT_records()
  h = model.get_hierarchy()
  asc = model.get_atom_selection_cache()
  chain_ids = [chain.id for chain in h.only_model().chains()]
  print (chain_ids)
  # Note leading whitespaces starting with chain D:
  # assert chain_ids == ['A', 'B', 'C', ' D', ' E', ' F', ' G', ' H', ' I', ' J', ' K', ' L', ' M', ' N', ' O', ' P', ' Q', ' R', ' S', ' T', ' U', ' V', ' W', ' X', ' Y', ' Z', ' 0']
  assert chain_ids == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '0']

  # Set of test cases
  A_selections = ["chain A", "chain 'A'", "chain ' A'"]
  D_selections = ["chain D", "chain 'D'", "chain ' D'"]

  n_atoms_A = []
  for a_sel in A_selections:
    chA = asc.selection(a_sel)
    n_atoms = h.select(chA).atoms_size()
    print("chA:", n_atoms)
    n_atoms_A.append(n_atoms)
  assert n_atoms_A == [44, 44, 0]

  n_atoms_D = []
  for d_sel in D_selections:
    chD = asc.selection(d_sel)
    n_atoms = h.select(chD).atoms_size()
    print("chD:", n_atoms)
    n_atoms_D.append(n_atoms)
  assert n_atoms_D == [44, 44, 0]
  print ("Is intenal expansion consistent?", n_atoms_A == n_atoms_D)
  assert n_atoms_A == n_atoms_D

  # BUT!!!! if we dump the expanded model and re-load it, everyting is consistent.
  expanded_lines = model.model_as_pdb()
  exp_h = iotbx.pdb.input(lines=expanded_lines, source_info=None).construct_hierarchy()

  exp_asc = exp_h.atom_selection_cache()
  # literally copy-paste from above
  n_atoms_A = []
  for a_sel in A_selections:
    chA = exp_asc.selection(a_sel)
    n_atoms = exp_h.select(chA).atoms_size()
    print("exp chA:", n_atoms)
    n_atoms_A.append(n_atoms)
  assert n_atoms_A == [44, 44, 0]

  n_atoms_D = []
  for d_sel in D_selections:
    chD = exp_asc.selection(d_sel)
    n_atoms = exp_h.select(chD).atoms_size()
    print("exp chD:", n_atoms)
    n_atoms_D.append(n_atoms)
  assert n_atoms_D == [44, 44, 0]
  print ("Is reading from pdb consistent?", n_atoms_A == n_atoms_D)
  assert n_atoms_A == n_atoms_D

if (__name__ == "__main__"):
  test_1()
