
from __future__ import division
from libtbx import easy_run

def exercise () :
  from mmtbx.regression import make_fake_anomalous_data
  mtz_file, pdb_file = make_fake_anomalous_data.generate_potassium_inputs()
  args = [pdb_file, mtz_file, "nproc=1", "elements=K,MG", "use_phaser=False"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_k = 0
  for line in result.stdout_lines :
    if ("Probable cation: K+1" in line) :
      n_k += 1
  assert n_k == 3
  print "OK"

if (__name__ == "__main__") :
  exercise()
