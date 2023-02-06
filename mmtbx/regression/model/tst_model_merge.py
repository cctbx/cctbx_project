from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb
from libtbx.utils import format_cpu_times
from iotbx.regression.ncs.tst_ncs import pdb_str_5
from libtbx.utils import null_out


def exercise_merge():
  inp = iotbx.pdb.input(lines=pdb_str_5, source_info=None)
  ab = mmtbx.model.manager(model_input=inp)
  ab.set_log(null_out())
  a = ab.apply_selection_string('chain A')
  b = ab.apply_selection_string('chain B')
  assert a.overall_counts().n_residues == 3
  a.merge_other_model(b)
  assert a.overall_counts().n_residues == 5

  try:
    ab.merge_other_model(b)
    raise AssertionError("Failed to fail in duplicate merge")
  except Exception as e:
    pass # ok

def run():
  exercise_merge()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
