from __future__ import division
from libtbx.easy_run import call
from smtbx.regression.test_data import fnames
import os

ciffile=fnames.thpp_cif
insfile=fnames.thpp_ins
hklfile=fnames.thpp_hkl
outfile=fnames.thpp_out

cl1 = "smtbx.anom_refine {} {} F -e13000 -t --outfile={} --overwrite".format(
    ciffile, hklfile, outfile)
cl2 = "smtbx.anom_refine {} {} F -e13000 -t --outfile={} --overwrite".format(
    insfile, hklfile, outfile)

out1 = "13000.0 0.152 -0.053 0.744 1.303"
out2 = "13000.0 -0.032 0.001 0.006 0.008"

def run():
  for cl,out in zip((cl1,cl2), (out1,out2)):
    print(cl)
    assert call(cl)==0, "smtbx.anom_refine refinement failed"
    with open(outfile) as f: result = f.readlines()[0].strip()
    assert out == result, "{} | {}".format(result, out)
  os.remove(outfile)


run()
