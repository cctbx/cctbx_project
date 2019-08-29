from __future__ import absolute_import, division, print_function
from scitbx.math.superpose import kabsch_rotation, kearsley_rotation, least_squares_fit
from scitbx.math import euler_angles_as_matrix
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import random
from six.moves import range

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def random_rotation():
  return euler_angles_as_matrix([random.uniform(0,360) for i in range(3)])

def exercise_rotation():
  reference = flex.vec3_double(flex.random_double(10*3)*10-5)
  other = reference.deep_copy()
  assert approx_equal(kabsch_rotation(reference, other), [1,0,0,0,1,0,0,0,1])
  assert approx_equal(kearsley_rotation(reference, other), [1,0,0,0,1,0,0,0,1])
  # pure rotations
  for n_sites in [1,2,3,7,10,30]:
    reference = flex.vec3_double(flex.random_double(n_sites*3)*10-5)
    other = reference.deep_copy()
    for i_trial in range(10):
      c = random_rotation()
      other_c = tuple(c) * other
      kc = kabsch_rotation(reference, other_c)
      assert approx_equal(kc.elems * other_c, reference)
      kc = kearsley_rotation(reference, other_c)
      assert approx_equal(kc.elems * other_c, reference)
  for n_sites in [1,2,3,7,10,30]:
    reference = flex.vec3_double(flex.random_double(n_sites*3)*10-5)
    other = reference.deep_copy()
    for i_trial in range(10):
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

      other_s = other + flex.vec3_double(flex.random_double(n_sites*3)*0.5)
      ks = kearsley_rotation(reference, other_s)
      other_ks = ks.elems * other_s
      rms_ks = reference.rms_difference(other_ks)
      c = random_rotation()
      other_cs = tuple(c) * other_s
      kcs = kearsley_rotation(reference, other_cs)
      other_kcs = tuple(kcs) * other_cs
      rms_kcs = reference.rms_difference(other_kcs)
      assert approx_equal(rms_kcs, rms_ks)

  #degenerate systems
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


  reference = flex.vec3_double([])
  kc = kearsley_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  reference = flex.vec3_double([(0,0,0)])
  kc = kearsley_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  reference = flex.vec3_double([(-1,0,0), (1,0,0)])
  kc = kearsley_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  other = flex.vec3_double([(0,-1,0), (0,1,0)])
  kc = kearsley_rotation(reference, other)
  assert approx_equal(kc.elems * other, reference)
  reference = flex.vec3_double([(-2,0,0), (-1,0,0), (1,0,0), (2,0,0)])
  kc = kearsley_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  other = flex.vec3_double([(0,-2,0), (0,-1,0), (0,1,0), (0,2,0)])
  kc = kearsley_rotation(reference, other)
  assert approx_equal(kc.elems * other, reference)
  reference = flex.vec3_double([(-1,0,0), (1,0,0), (0,-1,0), (0,1,0)])
  kc = kearsley_rotation(reference, reference)
  assert approx_equal(kc.elems * reference, reference)
  other = flex.vec3_double([(0,0,-1), (0,0,1), (0,-1,0), (0,1,0)])
  kc = kearsley_rotation(reference, other)
  assert approx_equal(kc.elems * other, reference)

def exercise(method):
  assert method in ["kearsley", "kabsch"]
  # global shifts
  for n_sites in [1,3,7,10,30]:
    reference = flex.vec3_double(flex.random_double(n_sites*3)*10-5)
    other = reference + list(flex.random_double(3)*100-50)
    for i_trial in range(10):
      s = least_squares_fit(reference, other, method)
      assert approx_equal(reference, s.other_sites_best_fit())
      c = random_rotation()
      s = least_squares_fit(reference, tuple(c)*other, method)
      if method == "kearsley": # Kabsch fails in special cases
        assert approx_equal(s.r.determinant(), 1)
      assert approx_equal(reference, s.other_sites_best_fit())
      assert approx_equal(s.rt().r, s.r)
      assert approx_equal(s.rt().t, s.t)
      assert approx_equal(reference, s.rt() * s.other_sites)

if (__name__ == "__main__"):
  exercise_rotation()
  exercise("kabsch")
  exercise("kearsley")
  print("OK")
