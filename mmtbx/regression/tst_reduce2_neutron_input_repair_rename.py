"""The second half of tst_reduce2_neutron_input_repair (cases b to e: names
permuted, one isotope per conformer, H and D apart), in its own file to keep
each test under a minute.
"""
from __future__ import absolute_import, division, print_function
from mmtbx.regression.tst_reduce2_neutron_input_repair import (
  test_rename,
  test_rename_refused,
  test_relabel,
  test_conformer,
  test_conformer_carbon,
  test_conformer_refused,
  test_different_parents_refused)

if __name__ == "__main__":
  test_rename()
  test_rename_refused()
  test_relabel()
  test_conformer()
  test_conformer_carbon()
  test_conformer_refused()
  test_different_parents_refused()
  print("OK")
