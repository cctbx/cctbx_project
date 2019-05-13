
from __future__ import division
from __future__ import print_function
import os
from libtbx import easy_run
import time

def exercise():
  from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
  base = "tst_validate_ca"
  mtz_file, pdb_file = generate_calcium_inputs(file_base=base, anonymize=False)
  time.sleep(2)
  args = ["\"%s\"" % pdb_file, "\"%s\"" % mtz_file, "wavelength=1.12",
          "nproc=1"]
  result = easy_run.fully_buffered("mmtbx.validate_ions %s" % " ".join(args)
    ).raise_if_errors()
  n_ca, n_bad = 0, 0
  for line in result.stdout_lines:
    if "| CA" in line:
      n_ca += 1
    if "!!!" in line:
      n_bad += 1
  assert n_ca == 1 and n_bad == 0
  for ext in [".pdb", ".mtz", "_fmodel.eff"]:
    os.remove(base + ext)
  print("OK")

if (__name__ == "__main__"):
  print("WARNING: TEST TOO SLOW. MAKE IT RUN UNDER 300s AND ENABLE BACK.")
  if 0: #XXX FIXME disabled
    t0 = time.time()
    exercise()
    print("Time: %6.2f"%(time.time()-t0))
    print("OK")
