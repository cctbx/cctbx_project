from scitbx.math.superpose import kabsch_rotation, least_squares_fit
from scitbx.math import euler_angles_as_matrix
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import random

def random_rotation():
  return euler_angles_as_matrix([random.uniform(0,360) for i in xrange(3)])

def exercise():
  reference = flex.vec3_double(flex.random_double(10*3)*10-5)
  other = reference.deep_copy()
  assert approx_equal(kabsch_rotation(reference, other), [1,0,0,0,1,0,0,0,1])
  # pure rotations
  for n_sites in [1,2,3,7,10,30]:
    reference = flex.vec3_double(flex.random_double(n_sites*3)*10-5)
    other = reference.deep_copy()
    for i_trial in xrange(10):
      c = random_rotation()
      other_c = tuple(c) * other
      kc = kabsch_rotation(reference, other_c)
      assert approx_equal(kc.elems * other_c, reference)
  # rotations + local shifts
  for n_sites in [1,2,3,7,10,30]:
    reference = flex.vec3_double(flex.random_double(n_sites*3)*10-5)
    other = reference.deep_copy()
    for i_trial in xrange(10):
      other_s = other + flex.vec3_double(flex.random_double(n_sites*3)*0.5)
      ks = kabsch_rotation(reference, other_s)
      other_ks = ks.elems * other_s
      rms_ks = reference.rms_difference(other_ks)
      c = random_rotation()
      other_cs = tuple(c) * other_s
      kcs = kabsch_rotation(reference, other_cs)
      other_kcs = tuple(kcs) * other_cs
      rms_kcs = reference.rms_difference(other_kcs)
      assert approx_equal(rms_kcs, rms_ks)
  # degenerate situations
  reference = flex.vec3_double([])
  kc = kabsch_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  reference = flex.vec3_double([(0,0,0)])
  kc = kabsch_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  reference = flex.vec3_double([(-1,0,0), (1,0,0)])
  kc = kabsch_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  other = flex.vec3_double([(0,-1,0), (0,1,0)])
  kc = kabsch_rotation(reference, other)
  assert approx_equal(kc.elems * other, reference)
  reference = flex.vec3_double([(-2,0,0), (-1,0,0), (1,0,0), (2,0,0)])
  kc = kabsch_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  other = flex.vec3_double([(0,-2,0), (0,-1,0), (0,1,0), (0,2,0)])
  kc = kabsch_rotation(reference, other)
  assert approx_equal(kc.elems * other, reference)
  reference = flex.vec3_double([(-1,0,0), (1,0,0), (0,-1,0), (0,1,0)])
  kc = kabsch_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  other = flex.vec3_double([(0,0,-1), (0,0,1), (0,-1,0), (0,1,0)])
  kc = kabsch_rotation(reference, other)
  assert approx_equal(kc.elems * other, reference)
  # global shifts
  for n_sites in [1,2,3,7,10,30]:
    reference = flex.vec3_double(flex.random_double(n_sites*3)*10-5)
    other = reference + list(flex.random_double(3)*100-50)
    for i_trial in xrange(10):
      s = least_squares_fit(reference, other)
      assert approx_equal(reference, s.other_sites_best_fit())
      c = random_rotation()
      s = least_squares_fit(reference, tuple(c)*other)
      assert approx_equal(reference, s.other_sites_best_fit())
  print "OK"

if (__name__ == "__main__"):
  exercise()
