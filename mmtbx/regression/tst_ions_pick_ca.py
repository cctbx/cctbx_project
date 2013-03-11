
from __future__ import division
from libtbx import easy_run

def exercise () :
  from mmtbx.regression import make_fake_anomalous_data
  mtz_file, pdb_file = make_fake_anomalous_data.generate_calcium_inputs()
  args = [pdb_file, mtz_file, "wavelength=1.12", "nproc=1", "use_phaser=False"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_ca = 0
  for line in result.stdout_lines:
    #print line
    if "Probable cation: CA+2" in line:
      n_ca += 1

  assert n_ca == 1
  print "OK"

if (__name__ == "__main__") :
  exercise()
