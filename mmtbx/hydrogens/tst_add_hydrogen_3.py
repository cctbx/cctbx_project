from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from libtbx.utils import null_out
#from libtbx.test_utils import approx_equal

# ------------------------------------------------------------------------------

def run():
  test_000(pdb_str = pdb_str_000)

# ------------------------------------------------------------------------------

def test_000(pdb_str):
  '''
    Make sure reduce does not crash for single_atom_residue models
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  # initial model (has no H atoms)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  number_h_expected = model_initial.get_hd_selection().count(True)
  assert(number_h_expected == 0)
  # place H atoms
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model_initial)
  reduce_add_h_obj.run()
  # We don't expect H atoms to be placed
  # (not enough restraints for single atom residues)
  model_h_added = reduce_add_h_obj.get_model()
  number_h_placed = model_h_added.get_hd_selection().count(True)
  assert(number_h_placed == 0)

# ------------------------------------------------------------------------------

pdb_str_000 = """
REMARK Make sure reduce does not crash for single_atom_residue models
CRYST1   22.029   33.502   24.035  90.00  90.00  90.00 P 1
ATOM      6  P     A A  10     -62.272  56.445  13.820  1.00 15.00           P
ATOM      7  P     G A  11     -63.673  51.410  11.026  1.00 15.00           P
ATOM      8  P     U A  12     -62.888  45.926   9.711  1.00 15.00           P
ATOM      9  P     U A  13     -60.326  41.305  11.244  1.00 15.00           P
ATOM     10  P     U A  14     -57.909  36.481  13.207  1.00 15.00           P
ATOM     11  P     G A  15     -62.106  32.943  15.800  1.00 15.00           P
ATOM     12  P     A A  16     -65.446  37.240  15.291  1.00 15.00           P
ATOM     13  P     U A  17     -66.286  42.354  18.232  1.00 15.00           P
ATOM     14  P     C A  18     -64.629  46.517  21.258  1.00 15.00           P
ATOM     15  P     A A  19     -60.460  50.019  23.746  1.00 15.00           P
ATOM     16  P     U A  20     -54.257  51.133  23.481  1.00 15.00           P
"""

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
