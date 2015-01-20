
from __future__ import division
from mmtbx.regression import tst_build_alt_confs
from scitbx.array_family import flex
from libtbx import easy_run
import libtbx.load_env # import dependency
import os

def exercise () :
  pdb_file = "tmp_ringer.pdb"
  mtz_file = "tmp_ringer.mtz"
  open(pdb_file, "w").write(tst_build_alt_confs.pdb_raw)
  args = [
    "phenix.fmodel",
    pdb_file,
    "high_resolution=2.0",
    "type=real",
    "r_free_flags_fraction=0.1",
    "random_seed=12345",
    "label=F",
    "output.file_name=%s" % mtz_file,
  ]
  easy_run.fully_buffered(args).raise_if_errors()
  result = easy_run.fully_buffered(
    "phenix.maps \"%s\" \"%s\" output.prefix=tmp_ringer" %
    (pdb_file, mtz_file)).raise_if_errors()
  assert (result.return_code == 0)
  result = easy_run.fully_buffered(
    "mmtbx.ringer \"%s\" tmp_ringer_map_coeffs.mtz" % pdb_file).raise_if_errors()
  _lines1 = open("tmp_ringer_ringer.csv").read().splitlines()
  lines1 = []
  for line in _lines1 :
    if ("2mFo-DFc" in line) :
      lines1.append(line)
  os.remove("tmp_ringer_ringer.csv")
  assert (result.return_code == 0)
  # Now with ccp4 map as input
  result2 = easy_run.fully_buffered(
    "phenix.mtz2map \"%s\" tmp_ringer_map_coeffs.mtz" % pdb_file)
  assert (result2.return_code == 0)
  result3 = easy_run.fully_buffered(
    "mmtbx.ringer \"%s\" tmp_ringer_map_coeffs_2mFo-DFc.ccp4" % pdb_file)
  assert (result3.return_code == 0)
  lines2 = open("tmp_ringer_ringer.csv").read().splitlines()
  assert len(lines1) == len(lines2)
  for line1, line2 in zip(lines1, lines2) :
    fields1 = line1.split(",")
    fields2 = line2.split(",")
    rho1 = flex.double([ float(x) for x in fields1[4:] ])
    rho2 = flex.double([ float(x) for x in fields2[4:] ])
    cc = flex.linear_correlation(x=rho1, y=rho2).coefficient()
    assert (cc >= 0.99), cc

if (__name__ == "__main__") :
  exercise()
  print "OK"
