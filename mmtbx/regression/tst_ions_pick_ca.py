
from __future__ import division
from libtbx import easy_run

def exercise () :
  from mmtbx.regression import make_fake_anomalous_data
  mtz_file, pdb_file = make_fake_anomalous_data.generate_calcium_inputs()
  args = ["\"%s\"" % pdb_file, "\"%s\"" % mtz_file, "wavelength=1.1",
    "nproc=1", "use_phaser=False", "fpp_ratio_max=1.2"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_ca = 0
  for line in result.stdout_lines:
    #print line
    if "Probable cation: CA+2" in line:
      n_ca += 1
  if (n_ca != 1) :
    print "\n".join(result.stdout_lines)
    raise RuntimeError("Expected 1 Ca2+, found %d" % n_ca)
  print "OK"

if (__name__ == "__main__") :
  exercise()
