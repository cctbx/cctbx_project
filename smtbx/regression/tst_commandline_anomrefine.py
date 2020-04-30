from __future__ import division
from libtbx.easy_run import call
from libtbx.test_utils import approx_equal
from smtbx.regression.test_data import fnames
import os

ciffile=fnames.thpp_cif
insfile=fnames.thpp_ins
hklfile=fnames.thpp_hkl
outfile=fnames.thpp_out

cl1 = "smtbx.anom_refine {} {} F -e13000 -t -o{} -O -d0 -s0 -c25".format(
    ciffile, hklfile, outfile)
cl2 = "smtbx.anom_refine {} {} F -e13000 -T -o{} -O -d0 -s0 -c25".format(
    insfile, hklfile, outfile)

out1 = [13000, .152, -.053, .743, 1.303]
out2 = [13000, -.032, .001, .005, .008, .015, .016, 3.674, 4.279]

def run():
  for cl,out in zip((cl1,cl2), (out1,out2)):
    print(cl)
    assert call(cl)==0, "smtbx.anom_refine refinement failed"
    with open(outfile) as f: r_lines = [l for l in f.readlines()]
    results = []
    for line in r_lines: results.extend([float(r) for r in line.split()])
    for i in range(len(results)):
      assert \
          approx_equal(out[i], results[i], eps=.0015), \
          "Wrong result for '{}'".format(cl)
  os.remove(outfile)


run()
