
from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run
import time

def exercise():
  from mmtbx.regression.make_fake_anomalous_data import generate_zinc_inputs
  base = "tst_pick_approx_zn"
  mtz_file, pdb_file = generate_zinc_inputs(file_base=base, anonymize = True)
  time.sleep(2)
  args = ["\"%s\"" % pdb_file, "\"%s\"" % mtz_file, "wavelength=1.54",
          "nproc=1", "elements=CA,ZN", "use_phaser=False"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_zn = 0
  for line in result.stdout_lines:
    if "Probable cation: ZN+2" in line:
      n_zn += 1
  if n_zn != 1:
    print("\n".join(result.stdout_lines))
    raise RuntimeError("Expected 1 ZN+2, found %d" % n_zn)
  os.remove(pdb_file)
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + ".pdb")
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")
  print("OK")

if (__name__ == "__main__"):
  exercise()
