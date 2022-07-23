
from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run

def exercise():
  try :
    from libsvm import svm      # import dependency
    from libsvm import svmutil  # import dependency
  except ImportError :
    print("libsvm not available, skipping this test")
    return
  from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
  base = "tst_pick_ca_svm"
  mtz_file, pdb_file = generate_calcium_inputs(file_base=base, anonymize=True)
  args = ["\"%s\"" % pdb_file, "\"%s\"" % mtz_file, "wavelength=1.1",
    "nproc=1", "use_phaser=False", "use_svm=True", "elements=CA,ZN"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  os.remove(pdb_file)
  os.remove(mtz_file)
  # "zn_frag_hoh.pdb" => "zn_frag_fmodel.eff"
  os.remove(os.path.splitext(pdb_file)[0][:-4] + ".pdb")
  os.remove(os.path.splitext(pdb_file)[0][:-4] + "_fmodel.eff")
  n_ca = 0
  for line in result.stdout_lines:
    if ("Final choice: CA" in line):
      n_ca += 1
  if (n_ca != 1):
    print("\n".join(result.stdout_lines))
    raise RuntimeError("Expected 1 Ca2+, found %d" % n_ca)
  print("OK")

if (__name__ == "__main__"):
  exercise()
