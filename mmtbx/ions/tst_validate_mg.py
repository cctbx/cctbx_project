
from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run
import time

def exercise():
  from mmtbx.regression.make_fake_anomalous_data import generate_magnessium_inputs
  base = "tst_validate_mg"
  mtz_file, pdb_file = generate_magnessium_inputs(file_base=base, anonymize=False)
  time.sleep(2)
  args = ["\"%s\"" % pdb_file, "\"%s\"" % mtz_file, "nproc=1"]
  result = easy_run.fully_buffered("mmtbx.validate_ions %s" % " ".join(args)
    ).raise_if_errors()
  n_mg, n_bad = 0, 0
  for line in result.stdout_lines :
    if "| MG" in line:
      n_mg += 1
    if "!!!" in line:
      n_bad += 1
  assert n_mg == 2 and n_bad == 0
  for ext in [".pdb", ".mtz", "_fmodel.eff"]:
    os.remove(base + ext)
  print("OK")

if (__name__ == "__main__"):
  exercise()
