from cctbx import maptbx
from cctbx.maptbx import bcr
import sys
from libtbx.test_utils import approx_equal

def run():
  o = maptbx.atom_curves(scattering_type="C")
  r = o.image(d_min = 2, b_iso = 0, radius_max=5, radius_step=0.03)
  #
  if 0:
    for dist, val in zip(r.radii, r.image_values):
      print "dist, val:", dist, val
    sys.stdout.flush()
  #
  bpeak,cpeak,rpeak,npeak,curve,curres = bcr.get_BCR(
    dens=r.image_values, dist=r.radii,  mxp=5, epsc=0.001, kpres=0)
  #
  if 0:
    print "bpeak:",bpeak
    print "cpeak:",cpeak
    print "rpeak:",rpeak
  #
  assert approx_equal(bpeak, [47.407781705050496, 23.21986190711327, 15.081095199013678, 9.663325738278937, 3.9912812852140402, 0.019352186565330765])
  assert approx_equal(cpeak, [14.737557556181374, -13.257929632841835, 7.709460378947266, -5.514409032143366, 3.127463490356433, 0.21703597852557188])
  assert approx_equal(rpeak, [0.0, 1.7264594573760674, 2.8646162024878516, 3.906036940851505, 4.8535995904354445, 4.980925062238841])

if (__name__ == "__main__"):
  run()
