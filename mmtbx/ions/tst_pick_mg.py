
from __future__ import division
from __future__ import print_function
import os
from libtbx import easy_run
import time

def exercise():
  from mmtbx.regression.make_fake_anomalous_data import generate_magnessium_inputs
  base = "tst_pick_mg"
  mtz_file, pdb_file = generate_magnessium_inputs(file_base=base, anonymize=True)
  time.sleep(2)
  args = [ "\"%s\"" % pdb_file, "\"%s\"" % mtz_file, "nproc=1",
           "use_phaser=False", "elements=MG" ]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_mg = 0
  for line in result.stdout_lines :
    if ("Probable cation: MG+2" in line):
      n_mg += 1
  if n_mg != 2:
    print("\n".join(result.stdout_lines))
    raise RuntimeError("Expected 2 MG+2, found %d" % n_mg)
  os.remove(pdb_file)
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + ".pdb")
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")
  print("OK")

if __name__ == "__main__":
  exercise()
