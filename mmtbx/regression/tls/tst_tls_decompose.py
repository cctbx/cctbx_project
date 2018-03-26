from __future__ import division
import libtbx.utils
from libtbx.test_utils import approx_equal
import scitbx.random
from scitbx import matrix
from scitbx.math import basic_statistics
from scitbx.array_family import flex
import itertools
import math, time, sys

from libtbx.utils import null_out, multi_out

from mmtbx.tls import analysis as tls_analysis
from mmtbx.tls.decompose import decompose_tls_matrices

d2r = math.pi/180

def dot(a,b):
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def generate_random_tls_motions(g):
  """Given a random number generator, return tls motions"""

  # Generate random input
  tx, ty, tz = sorted(abs(g(3))*2.00 + 0.10, reverse=True)      # on order of 1
  dx, dy, dz = sorted(abs(g(3))*0.01 + 0.01, reverse=True) # on order of
  sx, sy, sz = g(3)*5.0                             # on order of 5, don't need sorting, can be negative
  # Generate three rotation axes
  l_x = flex.vec3_double([tuple(g(3))])
  l_y = flex.vec3_double([tuple(g(3))])
  l_z = l_x.cross(l_y)
  l_x = l_y.cross(l_z)
  l_x = (l_x*(1.0/l_x.norm()))
  l_y = (l_y*(1.0/l_y.norm()))
  l_z = (l_z*(1.0/l_z.norm()))
  assert approx_equal(l_x.dot(l_y)[0], 0.0)
  assert approx_equal(l_y.dot(l_z)[0], 0.0)
  assert approx_equal(l_x.dot(l_z)[0], 0.0)
  assert approx_equal(l_x.dot(l_x)[0], 1.0)
  assert approx_equal(l_y.dot(l_y)[0], 1.0)
  assert approx_equal(l_z.dot(l_z)[0], 1.0)
  l_x = l_x[0]
  l_y = l_y[0]
  l_z = l_z[0]
  # Generate three points
  w_lx = tuple(g(3)*10.0)
  w_ly = tuple(g(3)*10.0)
  w_lz = tuple(g(3)*10.0)
  # Generate three vibration axes
  v_x = flex.vec3_double([tuple(g(3))])
  v_y = flex.vec3_double([tuple(g(3))])
  v_z = v_x.cross(v_y)
  v_x = v_y.cross(v_z)
  v_x = (v_x*(1.0/v_x.norm()))
  v_y = (v_y*(1.0/v_y.norm()))
  v_z = (v_z*(1.0/v_z.norm()))
  assert approx_equal(v_x.dot(v_y)[0], 0.0)
  assert approx_equal(v_y.dot(v_z)[0], 0.0)
  assert approx_equal(v_x.dot(v_z)[0], 0.0)
  assert approx_equal(v_x.dot(v_x)[0], 1.0)
  assert approx_equal(v_y.dot(v_y)[0], 1.0)
  assert approx_equal(v_z.dot(v_z)[0], 1.0)
  v_x = v_x[0]
  v_y = v_y[0]
  v_z = v_z[0]

  return dx, dy, dz, l_x, l_y, l_z, sx, sy, sz, tx, ty, tz, v_x, v_y, v_z, w_lx, w_ly, w_lz

def run_and_compare_dtm_and_tao(dx, dy, dz, l_x, l_y, l_z, sx, sy, sz, tx, ty, tz, v_x, v_y, v_z, w_M_lx, w_M_ly, w_M_lz):

  o = tls_analysis.tls_from_motions(
        dx=dx, dy=dy, dz=dz,
        l_x=l_x, l_y=l_y, l_z=l_z,
        sx=sx, sy=sy, sz=sz,
        tx=tx, ty=ty, tz=tz,
        v_x=v_x, v_y=v_y, v_z=v_z,
        w_M_lx=w_M_lx, w_M_ly=w_M_ly, w_M_lz=w_M_lz)

  # Add a bit of noise and convert to degrees
  T = tuple(flex.double(o.T_M.as_sym_mat3()))
  L = tuple(flex.double(o.L_M.as_sym_mat3()))
  S = tuple(flex.double(o.S_M))

  try:
    dtm = decompose_tls_matrices(T=tuple(T),
                                 L=tuple(L),
                                 S=tuple(S),
                                 l_and_s_in_degrees=False,
                                 verbose=False)
  except Exception as e:
    print e
    # Should never raise error
    raise

  try:
    tao = tls_analysis.run(T=matrix.sym(sym_mat3=T),
                           L=matrix.sym(sym_mat3=L),
                           S=matrix.sqr(S),
                           log=null_out())
  except Exception as e:
    # Raises error when invalid
    print 'TAO error:', e
    print 'DTM error:', dtm.error()
    # Return failure status
    return (dtm.is_valid(), False)

  # One succeeded one failed
  if not dtm.is_valid():
    print 'TAO error:', 'none'
    print 'DTM error:', dtm.error()
    return (dtm.is_valid(), True)

  # both suceeded
  check_equal_dtm_and_tao_results(dtm=dtm, tao=tao)

  return (True, True)

def check_equal_dtm_and_tao_results(dtm, tao):

  # tao orders axes by eigenvalues in INCREASING order.
  # dtm orders axes by eigenvalues in DECREASING order.
  # Therefore eigenvalues for x and z swap. As this is
  # a non-cyclic permutation the direction of some axes
  # changes as both systems are forced to be right-handed.

  # extract reult object
  tar = tao.result

  # Swap x and z
  assert approx_equal(dtm.l_amplitudes, (tar.dz, tar.dy, tar.dx))
  assert approx_equal(dtm.s_amplitudes, (tar.sz, tar.sy, tar.sx))
  assert approx_equal(dtm.v_amplitudes, (tar.tz, tar.ty, tar.tx))

  # Always check the vibration axes
  checklist = [(dtm.v_axis_directions[0], tar.v_z_M),
               (dtm.v_axis_directions[1], tar.v_y_M),
               (dtm.v_axis_directions[2], tar.v_x_M)]

  # Special case where all eigenvalues are the same
  if not (approx_equal(dtm.l_amplitudes[0], dtm.l_amplitudes[1], out=null_out()) and
          approx_equal(dtm.l_amplitudes[1], dtm.l_amplitudes[2], out=null_out()) and
          approx_equal(dtm.l_amplitudes[2], dtm.l_amplitudes[0], out=null_out())):
      checklist.extend([(dtm.l_axis_directions[0], tar.l_z),
                        (dtm.l_axis_directions[1], tar.l_y),
                        (dtm.l_axis_directions[2], tar.l_x)])

  # Check that each vector is the same (direction may be inverted)
  for a,b in checklist:
    d = dot(a,b) # calcuate normalised dot product
    d = d/abs(d) # i.e. +/-1
    assert approx_equal(abs(d), 1.0)
    assert approx_equal(flex.double(a)*d, b)

  # Check that intersection points of rotation axes are the same
  assert approx_equal(dtm.l_axis_intersections, (tar.w_M_lz, tar.w_M_ly, tar.w_M_lx))

def exercise_decompose_tls_matrices():

  T = (9.598824505, 15.650977723, 8.860020860, -0.005192226, -1.725114669, -0.869540262)
  L = (1.23001438881,4.61188308667,2.18102205791,0.435192109608,0.517157269537,-1.55699937912)
  S = (-0.725668016697,-2.01548281778,0.431134160979,3.6329322896,-0.252099816856,-0.771964293756,-0.968895922456,1.70374138375,0.977767833552)

  dtm_deg_input = decompose_tls_matrices(T=T, L=L, S=S, l_and_s_in_degrees=True, verbose=False)
  dtm_rad_input = decompose_tls_matrices(T=T, L=tuple([v*d2r*d2r for v in L]), S=tuple([v*d2r for v in S]),
                                         l_and_s_in_degrees=False, verbose=False)

  for dtm in [dtm_deg_input,dtm_rad_input]:

    assert dtm.is_valid()
    assert dtm.error() == 'none'

    assert approx_equal(dtm.v_amplitudes, (3.06921018016, 2.94143684759, 1.41201837191))
    assert approx_equal(dtm.l_amplitudes, (0.0404764451411, 0.024588606101, 0.0141767022316))
    assert approx_equal(dtm.s_amplitudes, (-2.8463226589406196, 13.419545390835827, -4.109522433370094))

    assert approx_equal(dtm.l_axis_directions, [(0.04063819985118014,  0.9008939255667365, -0.4321327013659231),
                                                (0.6515256136265412,   0.30400592381522407, 0.6950501946433874),
                                                (0.7575373994077218,  -0.3097911121420477, -0.5746012141793462)])

    assert approx_equal(dtm.v_axis_directions, [(0.5480708503205247,  -0.74962057386549,    0.3710624452386817),
                                                (0.40867651787831666, -0.1470759298532635, -0.9007508948608665),
                                                (0.7297957568825613,   0.6453198169089266,  0.22574429592093403)])

    assert approx_equal(dtm.l_axis_intersections, [(-13.015849649869374, 35.879232351401456, -60.666685853674494),
                                                   (-65.65635765706187, -50.022061840205446, -84.73832882566137),
                                                   (-110.51255057376397, 51.16255655992905, -197.9357742701504)])

  # Check that same as output from "mmtbx.tls.analysis.run"
  # TAO needs input in radians
  tao = tls_analysis.run(T=matrix.sym(sym_mat3=T),
                         L=matrix.sym(sym_mat3=L)*d2r*d2r,
                         S=matrix.sqr(S)*d2r,
                         log=null_out())

  check_equal_dtm_and_tao_results(dtm=dtm_rad_input, tao=tao)
  check_equal_dtm_and_tao_results(dtm=dtm_deg_input, tao=tao)

def exercise_tls_analysis_and_decompose_similar_results():
  """Test random TLS matrices to see that different methods give same results"""

  from scitbx.random import set_random_seed, variate, normal_distribution
  set_random_seed(100)
  g = variate(normal_distribution(0,1))

  #for noise in (0.001, 0.1):
  n_success = 0
  n_run = 0
  for i in xrange(24):
    n_run += 1

    dx, dy, dz, l_x, l_y, l_z, sx, sy, sz, tx, ty, tz, \
            v_x, v_y, v_z, w_M_lx, w_M_ly, w_M_lz = generate_random_tls_motions(g=g)

    ii = i//3

    # Different tests (three of each)
    if   ii==0:
        pass
    elif ii==1:
        dx = 1e-9
    elif ii==2:
        dy = 1e-9
    elif ii==3:
        dz = 1e-9
    elif ii==4:
        dx = dy = 1e-9
    elif ii==5:
        dy = dz = 1e-9
    elif ii==6:
        dz = dx = 1e-9
    elif ii==7:
        dx = dy = dz = 1e-9

    suc1, suc2 = run_and_compare_dtm_and_tao(dx, dy, dz, l_x, l_y, l_z, sx, sy, sz, tx, ty, tz, v_x, v_y, v_z, w_M_lx, w_M_ly, w_M_lz)

    assert suc1 is suc2, 'DTM {} while TAO {}'.format(*[['failed', 'succeeded'][int(v)] for v in [suc1,suc2]])

def run():
  libtbx.utils.show_times_at_exit()
  exercise_decompose_tls_matrices
  exercise_tls_analysis_and_decompose_similar_results()

if __name__ == '__main__':
  run()
