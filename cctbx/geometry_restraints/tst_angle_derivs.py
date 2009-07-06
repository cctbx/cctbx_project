from cctbx.geometry_restraints import angle, angle_delta_deg
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import math
import sys

def derivs_fd(a, order, eps=1.e-6):
  result = flex.vec3_double()
  sites0 = a.sites
  sites = [list(site) for site in sites0]
  for i_site in xrange(3):
    ds = []
    for i_dim in xrange(3):
      samples = []
      for signed_eps in [eps, -eps]:
        sites[i_site][i_dim] = sites0[i_site][i_dim] + signed_eps
        a_eps = angle(sites=sites, angle_ideal=a.angle_ideal, weight=a.weight)
        if (order == 1):
          samples.append(a_eps.residual())
        elif (order == 2):
          samples.append(a_eps.gradients()[i_site][i_dim])
        else:
          raise RuntimeError
      sites[i_site][i_dim] = sites0[i_site][i_dim]
      ds.append((samples[0]-samples[1])/(2*eps))
    result.append(tuple(ds))
  return result

def compare_derivs(out, ana, fin, expect_failure, eps=1.e-6):
  s = 1 / max(1, flex.max(flex.abs(ana.as_double())))
  print >> out, "  ana:", list(ana*s)
  print >> out, "  fin:", list(fin*s)
  if (not expect_failure):
    assert approx_equal(ana*s, fin*s, eps=eps)
  else:
    approx_equal(ana*s, fin*s, eps=eps)
  print >> out

def check_derivs(out, a, expect_failure=False):
  print >> out, "sites:", a.sites
  print >> out
  gc = a.grads_and_curvs()
  print >> out, "grads:", a.sites
  g_fd = derivs_fd(a=a, order=1)
  compare_derivs(out=out, ana=gc[:3], fin=g_fd, expect_failure=expect_failure)
  print >> out, "curvs:", a.sites
  c_fd = derivs_fd(a=a, order=2)
  compare_derivs(out=out, ana=gc[3:], fin=c_fd, expect_failure=expect_failure)
  print >> out

def run(args):
  assert args in [[], ["--verbose"]]
  if (len(args) == 0):
    out = null_out()
  else:
    out = sys.stdout
  #
  mt = flex.mersenne_twister(seed=0)
  #
  for i_trial in xrange(10):
    l0 = mt.random_double() + 0.5
    l1 = mt.random_double() + 0.5
    l2 = mt.random_double() + 0.5
    angle_model = mt.random_double() * 178 + 1 \
                + 180 * (mt.random_size_t() % 3 - 1)
    v = matrix.col(mt.random_double_point_on_sphere())
    axis = v.ortho()
    site1 = v * l1
    site0 = site1 + v * l0
    r = axis.axis_and_angle_as_r3_rotation_matrix(angle=angle_model, deg=True)
    site2 = site1 + (r * v) * l2
    a = angle(
      sites=[site0, site1, site2],
      angle_ideal=mt.random_double() * 720 - 360,
      weight=mt.random_double() * 10 + 0.1)
    assert approx_equal(min(
      abs(angle_delta_deg(angle_1=a.angle_model, angle_2= angle_model)),
      abs(angle_delta_deg(angle_1=a.angle_model, angle_2=-angle_model))), 0)
    check_derivs(out=out, a=a)
  #
  for site2 in [(0,2.3,0), (0,0,2.5)]:
    perm = flex.size_t([0,1,2])
    while True:
      a = angle(
        sites=tuple(flex.vec3_double(
          [(1.2,0,0), (0,0,0), site2]).select(perm)),
        angle_ideal=mt.random_double() * 720 - 360,
        weight=mt.random_double() * 10 + 0.1)
      check_derivs(out=out, a=a)
      if (not perm.next_permutation()):
        break
  #
  for site0 in [(1,0,0),(0,1,0),(0,0,1),(1,1,1)]:
    perm = flex.size_t([0,1,2])
    while True:
      a = angle(
        sites=tuple(flex.vec3_double(
          [site0, (0,0,0), -matrix.col(site0)]).select(perm)),
        angle_ideal=180,
        weight=1.3)
      check_derivs(out=out, a=a, expect_failure=True)
      if (not perm.next_permutation()):
        break
  #
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
