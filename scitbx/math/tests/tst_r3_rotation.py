from __future__ import absolute_import, division, print_function
import scitbx.math
from scitbx.math import r3_rotation_vector_to_vector as vector_to_vector
from scitbx.math import r3_rotation_vector_to_001 as vector_to_001
from scitbx.math import r3_rotation_vector_to_010 as vector_to_010
from scitbx.math import r3_rotation_vector_to_100 as vector_to_100
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.utils import format_cpu_times
from libtbx.test_utils import Exception_expected, approx_equal
import math
import time
import sys
from six.moves import range

def exercise_axis_and_angle(
      axis_range=2,
      angle_max_division=12,
      angle_min_power=-30):
  from_matrix = scitbx.math.r3_rotation_axis_and_angle_from_matrix(
    r=[1,0,0,0,1,0,0,0,1])
  assert approx_equal(from_matrix.axis, [1/3**0.5]*3)
  assert approx_equal(from_matrix.angle(deg=True), 0)
  assert approx_equal(from_matrix.angle(deg=False), 0)
  from_matrix = scitbx.math.r3_rotation_axis_and_angle_from_matrix(r=[0]*9)
  assert approx_equal(from_matrix.axis, [0,0,1])
  assert approx_equal(from_matrix.angle(deg=True), 90)
  assert approx_equal(from_matrix.angle(deg=False), math.pi/2)
  #
  angles = []
  for d in range(1,angle_max_division+1):
    angles.append(360/d)
    angles.append(-360/d)
  for p in range(-angle_min_power+1):
    angles.append(10**(-p))
    angles.append(-10**(-p))
  hex_orth = matrix.sqr([
    8.7903631196301042, -4.3951815598150503, 0,
    0, 7.6126777700894994, 0,
    0, 0, 14.943617303371177])
  rot_axes = []
  rot_angles = []
  rot_matrices = []
  for u in range(-axis_range, axis_range+1):
    for v in range(-axis_range, axis_range+1):
      for w in range(axis_range+1):
        for axis in [(u,v,w), (hex_orth*matrix.col((u,v,w))).elems]:
          for angle in angles:
            try:
              r = scitbx.math.r3_rotation_axis_and_angle_as_matrix(
                axis=axis, angle=angle, deg=True)
              rot_axes.append(axis)
              rot_angles.append(angle)
              rot_matrices.append(r)
            except RuntimeError:
              assert axis == (0,0,0)
              try:
                scitbx.math.r3_rotation_axis_and_angle_as_unit_quaternion(
                  axis=axis, angle=angle, deg=True)
              except RuntimeError: pass
              else: raise Exception_expected
            else:
              q = scitbx.math.r3_rotation_axis_and_angle_as_unit_quaternion(
                axis=axis, angle=angle, deg=True)
              assert approx_equal(abs(matrix.col(q)), 1)
              rq = scitbx.math.r3_rotation_unit_quaternion_as_matrix(q=q)
              assert approx_equal(rq, r)
              from_matrix = scitbx.math.r3_rotation_axis_and_angle_from_matrix(
                r=r)
              rr = from_matrix.as_matrix()
              assert approx_equal(rr, r)
              qq = from_matrix.as_unit_quaternion()
              assert approx_equal(abs(matrix.col(qq)), 1)
              rq = scitbx.math.r3_rotation_unit_quaternion_as_matrix(q=qq)
              assert approx_equal(rq, r)
              qq = scitbx.math.r3_rotation_matrix_as_unit_quaternion(r=r)
              assert approx_equal(abs(matrix.col(qq)), 1)
              rq = scitbx.math.r3_rotation_unit_quaternion_as_matrix(q=qq)
              assert approx_equal(rq, r)
              rm = matrix.sqr(r)
              assert rm.is_r3_rotation_matrix()
              qm = rm.r3_rotation_matrix_as_unit_quaternion()
              assert approx_equal(abs(qm), 1)
              assert approx_equal(qm, qq)
              rqmm = qm.unit_quaternion_as_r3_rotation_matrix()
              assert approx_equal(rqmm, r)
              for deg in [False, True]:
                rr = scitbx.math.r3_rotation_axis_and_angle_as_matrix(
                  axis=from_matrix.axis,
                  angle=from_matrix.angle(deg=deg),
                  deg=deg,
                  min_axis_length=1-1.e-5)
                assert approx_equal(rr, r)
                qq = scitbx.math.r3_rotation_axis_and_angle_as_unit_quaternion(
                  axis=from_matrix.axis,
                  angle=from_matrix.angle(deg=deg),
                  deg=deg,
                  min_axis_length=1-1.e-5)
                qq = from_matrix.as_unit_quaternion()
                assert approx_equal(abs(matrix.col(qq)), 1)
                rq = scitbx.math.r3_rotation_unit_quaternion_as_matrix(q=qq)
                assert approx_equal(rq, r)
                for conv in [
                  scitbx.math.r3_rotation_axis_and_angle_as_matrix,
                  scitbx.math.r3_rotation_axis_and_angle_as_unit_quaternion]:
                  try:
                    conv(
                      axis=from_matrix.axis,
                      angle=from_matrix.angle(deg=deg),
                      deg=deg,
                      min_axis_length=1+1.e-5)
                  except RuntimeError: pass
                  else: raise Exception_expected
  rot_matrices_vectorized = scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      rot_axes, rot_angles, deg=True)
  for a,b in zip(rot_matrices, rot_matrices_vectorized):
    assert a==b
  #
  for i_trial in range(100):
    r = flex.random_double_r3_rotation_matrix()
    from_matrix = scitbx.math.r3_rotation_axis_and_angle_from_matrix(r=r)
    rr = from_matrix.as_matrix()
    assert approx_equal(rr, r)
    assert approx_equal(math.cos(from_matrix.angle()), (r[0]+r[4]+r[8]-1)/2)
    q = matrix.col(from_matrix.as_unit_quaternion())
    assert approx_equal(abs(q), 1)
    rq = scitbx.math.r3_rotation_unit_quaternion_as_matrix(q=q)
    assert approx_equal(rq, r)
    rq = matrix.col(q).unit_quaternion_as_r3_rotation_matrix()
    assert approx_equal(rq, r)

def unit_quaternion_matrix_timings(n_trials=50, n_repeats=500):
  cpp_um = scitbx.math.r3_rotation_unit_quaternion_as_matrix
  cpp_mu = scitbx.math.r3_rotation_matrix_as_unit_quaternion
  times = [0,0,0,0]
  for i_trial in range(n_trials):
    qm = matrix.col.random(n=4, a=-1, b=1).normalize()
    rm = qm.unit_quaternion_as_r3_rotation_matrix()
    q = qm.elems
    r = cpp_um(q=q)
    t0 = time.time()
    for i_repeat in range(n_repeats):
      qm.unit_quaternion_as_r3_rotation_matrix()
    times[0] += time.time()-t0
    t0 = time.time()
    for i_repeat in range(n_repeats):
      rm.r3_rotation_matrix_as_unit_quaternion()
    times[1] += time.time()-t0
    t0 = time.time()
    for i_repeat in range(n_repeats):
      cpp_um(q=q)
    times[2] += time.time()-t0
    t0 = time.time()
    for i_repeat in range(n_repeats):
      cpp_mu(r=r)
    times[3] += time.time()-t0
  print("times unit quaternion <-> matrix")
  print("     py: %.2f %.2f s" % (times[0], times[1]))
  print("    c++: %.2f %.2f s" % (times[2], times[3]))
  if (times[2] != 0 and times[3] != 0):
    print("  ratio: %.2f %.2f" % (times[0]/times[2], times[1]/times[3]))
  sys.stdout.flush()

def check_vector_to_vector(g, t):
  assert approx_equal(abs(matrix.col(g)), 1)
  assert approx_equal(abs(matrix.col(t)), 1)
  r = matrix.sqr(vector_to_vector(given_unit_vector=g, target_unit_vector=t))
  assert approx_equal(r * matrix.col(g), t)
  assert approx_equal(r.determinant(), 1)

def check_vector_to_001(g):
  assert approx_equal(abs(matrix.col(g)), 1)
  r = matrix.sqr(vector_to_001(given_unit_vector=g))
  assert approx_equal(r * matrix.col(g), (0,0,1))
  assert approx_equal(r.determinant(), 1)
  p = matrix.col(g).vector_to_001_rotation()
  assert approx_equal(p * matrix.col(g), (0,0,1))
  assert approx_equal(p.determinant(), 1)
  assert approx_equal(p, r)

def check_vector_to_010(g):
  assert approx_equal(abs(matrix.col(g)), 1)
  r = matrix.sqr(vector_to_010(given_unit_vector=g))
  assert approx_equal(r * matrix.col(g), (0,1,0))
  assert approx_equal(r.determinant(), 1)

def check_vector_to_100(g):
  assert approx_equal(abs(matrix.col(g)), 1)
  r = matrix.sqr(vector_to_100(given_unit_vector=g))
  assert approx_equal(r * matrix.col(g), (1,0,0))
  assert approx_equal(r.determinant(), 1)

def exercise_vector_to_vector(angle_exponent_step=10, n_trials=10):
  principal_vectors = [matrix.col(v) for v in ((1,0,0), (0,1,0), (0,0,1))]
  for g0 in principal_vectors:
    for t0 in principal_vectors:
      for g in [g0, -g0]:
        for t in [t0, -t0]:
          check_vector_to_vector(g=g, t=t)
  if (sys.platform.startswith("osf")):
    max_exp = 300
  else:
    max_exp = 340
  for ig,g0 in enumerate(principal_vectors):
    for it,t0 in enumerate(principal_vectors):
      if (ig == it): continue
      axis = g0.cross(t0)
      for e in range(0, max_exp, angle_exponent_step):
        angle = 10**(-e)
        for angle in [angle, -angle]:
          r = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
            axis=axis,
            angle=angle))
          for g in [g0, -g0]:
            for t in [t0, -t0]:
              check_vector_to_vector(g, r*t)
              check_vector_to_vector(r*g, t)
              check_vector_to_vector(r*g, r*t)
  for i_trial in range(n_trials):
    g = matrix.col(flex.random_double_point_on_sphere())
    check_vector_to_vector(g, g)
    check_vector_to_vector(g, -g)
    t = matrix.col(flex.random_double_point_on_sphere())
    check_vector_to_vector(g, t)
    for e in range(0, max_exp, angle_exponent_step):
      angle = 10**(-e)
      for angle in [angle, -angle]:
        rt = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          axis=g,
          angle=angle)) * t
        check_vector_to_vector(rt, t)
        check_vector_to_vector(rt, -t)
  #
  check_vector_to_001((0,0,1))
  check_vector_to_001((0,0,-1))
  for e in range(0, max_exp, angle_exponent_step):
    angle = 10**(-e)
    for angle in [angle, -angle]:
      rg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
        axis=(1,1,0),
        angle=angle)) * matrix.col((0,0,1))
      check_vector_to_001(rg)
      check_vector_to_001(-rg)
  for i_trial in range(n_trials):
    g = matrix.col(flex.random_double_point_on_sphere())
    check_vector_to_001(g)
  #
  check_vector_to_010((0,1,0))
  check_vector_to_010((0,-1,0))
  for e in range(0, max_exp, angle_exponent_step):
    angle = 10**(-e)
    for angle in [angle, -angle]:
      rg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
        axis=(1,0,1),
        angle=angle)) * matrix.col((0,1,0))
      check_vector_to_010(rg)
      check_vector_to_010(-rg)
  for i_trial in range(n_trials):
    g = matrix.col(flex.random_double_point_on_sphere())
    check_vector_to_010(g)
  #
  check_vector_to_100((1,0,0))
  check_vector_to_100((-1,0,0))
  for e in range(0, max_exp, angle_exponent_step):
    angle = 10**(-e)
    for angle in [angle, -angle]:
      rg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
        axis=(0,1,1),
        angle=angle)) * matrix.col((1,0,0))
      check_vector_to_100(rg)
      check_vector_to_100(-rg)
  for i_trial in range(n_trials):
    g = matrix.col(flex.random_double_point_on_sphere())
    check_vector_to_100(g)
  #
  rvv = vector_to_vector((0,0,-1), (0,0,1))
  rv1 = vector_to_001((0,0,-1))
  assert approx_equal(rv1, rvv)
  rvv = vector_to_vector((0,-1,0), (0,1,0))
  rv1 = vector_to_010((0,-1,0))
  assert approx_equal(rv1, rvv)
  rvv = vector_to_vector((-1,0,0), (1,0,0))
  rv1 = vector_to_100((-1,0,0))
  assert approx_equal(rv1, rvv)

def exercise():
  exercise_axis_and_angle()
  unit_quaternion_matrix_timings()
  exercise_vector_to_vector()
  print(format_cpu_times())

if (__name__ == "__main__"):
  exercise()
