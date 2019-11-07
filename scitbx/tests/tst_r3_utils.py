from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal

def exercise_clash_detector_simple():
  from scitbx.array_family import flex
  sites_cart = flex.vec3_double([
    (1,2.0,3),
    (1,3.3,3)])
  from scitbx.r3_utils import clash_detector_simple
  cd = clash_detector_simple(n_sites=2, threshold=1.2)
  assert approx_equal(cd.threshold_sq, 1.2**2)
  assert not cd.has_clash(sites_cart=sites_cart)
  sites_cart[1] = (1,3.1,3)
  assert cd.has_clash(sites_cart=sites_cart)
  cd.add_exclusion(i=0, j=1)
  assert not cd.has_clash(sites_cart=sites_cart)

def run(args):
  assert len(args) == 0
  exercise_clash_detector_simple()
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
