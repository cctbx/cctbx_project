"""The second half of tst_reduce2_protonation (supplied H, neutron options,
qmready charge), in its own file to keep each test under a minute.
"""
from __future__ import absolute_import, division, print_function
from mmtbx.regression.tst_reduce2_protonation import (
  test_supplied_h,
  test_supplied_half_pair,
  test_supplied_h_only_conformers,
  test_supplied_pair_relabelled,
  test_acid_out_of_plane_line,
  test_supplied_h_per_conformer,
  test_neutron_options,
  test_qmready_charge)

if __name__ == "__main__":
  test_supplied_h()
  test_supplied_half_pair()
  test_supplied_h_only_conformers()
  test_supplied_pair_relabelled()
  test_acid_out_of_plane_line()
  test_supplied_h_per_conformer()
  test_neutron_options()
  test_qmready_charge()
  print("OK")
