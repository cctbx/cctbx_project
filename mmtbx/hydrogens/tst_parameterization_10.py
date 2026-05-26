from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen

def run():
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

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
