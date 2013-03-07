
from __future__ import division
from libtbx import easy_run

def exercise () :
  from mmtbx.regression import make_fake_anomalous_data
  mtz_file, pdb_file = make_fake_anomalous_data.generate_calcium_inputs()
  args = ["ca_frag.pdb", "ca_frag.mtz", "wavelength=1.12", "nproc=1"]
  result = easy_run.fully_buffered("mmtbx.validate_ions %s" % " ".join(args)
    ).raise_if_errors()
  n_ca, n_bad = 0, 0
  for line in result.stdout_lines:
    print line
    if "| CA" in line:
      n_ca += 1
    if "!!!" in line:
      n_bad += 1

  assert n_ca == 1 and n_bad == 0
  print "OK"

if (__name__ == "__main__") :
  exercise()
