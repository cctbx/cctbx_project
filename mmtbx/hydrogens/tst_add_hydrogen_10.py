from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen

def run_00():
  """
  Exercise adding H to water
  """
  pdb_str = """
CRYST1   25.876   39.544   32.276  90.00 109.64  90.00 P 1 21 1
HETATM  846  O   HOH A 511       4.272   0.767 -10.242  1.00  7.02           O
  """
  pdb_inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  o = reduce_hydrogen.place_hydrogens(model = model, exclude_water = False)
  o.run()
  cntr = 0
  expected = [
    [5.122, 0.767, -10.242],
    [4.068, 1.592, -10.242]
  ]
  for a in o.get_model().get_hierarchy().atoms():
    if a.element_is_hydrogen():
      cntr+=1
      xyz = [round(_,3) for _ in a.xyz]
      assert xyz in expected
  assert cntr==2

def run_01():
  """
  Exercise adding H to water: make sure workaround does not remove them
  """
  pdb_str = """
CRYST1   19.465   21.432   29.523  90.00  90.00  90.00 P 21 21 21    4
HETATM  151  O   HOH A1006       9.937  14.244   1.856  0.50  8.38           O
HETATM  165  O   HOH A1106       9.290  13.738   1.763  0.50 18.99           O
HETATM  157  O   HOH A1012      -0.833  19.856   2.677  0.50 12.16           O
HETATM  166  O   HOH A1112      -0.886  20.218   1.931  0.50  9.08           O
END
  """
  pdb_inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  o = reduce_hydrogen.place_hydrogens(model = model, exclude_water = False)
  o.run()
  assert o.get_model().get_xray_structure().hd_selection().count(True)==8

if (__name__ == "__main__"):
  t0 = time.time()
  run_00()
  run_01()
  print("OK. Time: %8.3f"%(time.time()-t0))
