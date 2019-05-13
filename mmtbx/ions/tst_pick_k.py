
from __future__ import division
from __future__ import print_function
import os
from libtbx import easy_run
import time

def exercise():
  # FIXME
  print("Temporarily disabled, skipping")
  return
  from mmtbx.regression.make_fake_anomalous_data import generate_potassium_inputs
  base = "tst_pick_k"
  mtz_file, pdb_file = generate_potassium_inputs(file_base=base, anonymize=True)
  time.sleep(2)
  args = [pdb_file, mtz_file, "nproc=1", "elements=K,MG", "use_phaser=False"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_k = 0
  for line in result.stdout_lines :
    if ("Probable cation: K+1" in line):
      n_k += 1
  assert n_k == 3
  os.remove(pdb_file)
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")
  print("OK")

if (__name__ == "__main__"):
  exercise()
