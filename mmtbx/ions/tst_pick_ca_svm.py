
from __future__ import division
from libtbx import easy_run
import time

def exercise () :
  try :
    import svm
    import svmutil
  except ImportError :
    print "libsvm not available, skipping this test"
    return
  from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
  mtz_file, pdb_file = generate_calcium_inputs(
    file_base = "tst_ions_pick_ca_svm", anonymize = True)
  args = ["\"%s\"" % pdb_file, "\"%s\"" % mtz_file, "wavelength=1.1",
    "nproc=1", "use_phaser=False", "use_svm=True", "elements=CA"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_ca = 0
  for line in result.stdout_lines:
    if ("Final choice: CA" in line) :
      n_ca += 1
  if (n_ca != 1) :
    print "\n".join(result.stdout_lines)
    raise RuntimeError("Expected 1 Ca2+, found %d" % n_ca)
  print "OK"

if (__name__ == "__main__") :
  exercise()
