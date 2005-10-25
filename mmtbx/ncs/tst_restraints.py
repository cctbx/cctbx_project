from mmtbx import ncs
import mmtbx.ncs.restraints
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times

def exercise_pair_registry_basic():
  registry = ncs.restraints.pair_registry(30, 3)
  selection_pairs = registry.selection_pairs()
  assert zip(*selection_pairs[0]) == []
  assert zip(*selection_pairs[1]) == []
  assert registry.enter(i_seq=0, j_seq=10, j_ncs=2) == 1
  assert registry.enter(i_seq=10, j_seq=0, j_ncs=2) == 0
  assert registry.enter(i_seq=0, j_seq=10, j_ncs=1) == -2
  assert registry.enter(i_seq=21, j_seq=1, j_ncs=1) == 1
  assert registry.enter(i_seq=0, j_seq=20, j_ncs=1) == 1
  selection_pairs = registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert zip(*selection_pairs[0]) == [(0, 20), (1, 21)]
  assert zip(*selection_pairs[1]) == [(0, 10)]

def adp_iso_residual_sum(weight, average_power, u_isos):
  n = u_isos.size()
  u_ave = flex.sum(u_isos) / n
  u_diff = u_isos - u_ave
  return weight * flex.sum_sq(u_diff) / u_ave**average_power

def adp_iso_analytical_gradients(weight, average_power, u_isos):
  """
um = (u1 + u2)/2
r=w*((u1-um)^2+(u2-um)^2)/um^mp
D[r,u1]

um = (u1 + u2 + u3)/3
r=w*((u1-um)^2+(u2-um)^2+(u3-um)^2)/um^mp
D[r,u1]

um = (u1 + u2 + u3 + u4)/4
r=w*((u1-um)^2+(u2-um)^2+(u3-um)^2+(u4-um)^2)/um^mp
D[r,u1]
"""
  n = u_isos.size()
  u_ave = flex.sum(u_isos) / n
  u_diff = u_isos - u_ave
  n_ave_pow = n*u_ave**average_power
  return weight * (  2/n_ave_pow*(u_diff*n-flex.sum(u_diff))
                   - average_power/(n_ave_pow*u_ave)*flex.sum_sq(u_diff))

def adp_iso_finite_difference_gradients(
      weight,
      average_power,
      u_isos,
      eps=1.e-6):
  result = flex.double()
  for i_u_iso in xrange(u_isos.size()):
    rs = []
    for signed_eps in [eps, -eps]:
      u_isos_eps = u_isos.deep_copy()
      u_isos_eps[i_u_iso] += signed_eps
      r = adp_iso_residual_sum(
        weight=weight, average_power=average_power, u_isos=u_isos_eps)
      rs.append(r)
    result.append((rs[0]-rs[1])/(2*eps))
  return result

def exercise_adp_iso_analytical():
  mersenne_twister = flex.mersenne_twister(seed=0)
  for weight in [1.234, 2.134]:
    for average_power in [0.345, 0.589]:
      for size in xrange(2,20):
        for i in xrange(10):
          u_isos = mersenne_twister.random_double(size=size) + 1.e-3
          a = adp_iso_analytical_gradients(
            weight=weight, average_power=average_power, u_isos=u_isos)
          f = adp_iso_finite_difference_gradients(
            weight=weight, average_power=average_power, u_isos=u_isos)
          assert approx_equal(a, f)

def exercise_pair_registry_adp_iso():
  mersenne_twister = flex.mersenne_twister(seed=0)
  for n_seq in xrange(2,20):
    registry = ncs.restraints.pair_registry(n_seq=n_seq, n_ncs=2)
    for j_seq in xrange(1,n_seq):
      assert registry.enter(i_seq=0, j_seq=j_seq, j_ncs=1) == 1
    selection_pairs = registry.selection_pairs()
    assert zip(*selection_pairs[0]) == zip([0]*(n_seq-1), xrange(1,n_seq))
    weight = 2.134
    average_power = 0.589
    u_isos = mersenne_twister.random_double(size=n_seq) + 1.e-3
    gradients_in = mersenne_twister.random_double(size=n_seq)
    gradients = gradients_in.deep_copy()
    registry_residual_sum = registry.adp_iso_residual_sum(
      weight=2.134,
      average_power=0.589,
      u_isos=u_isos,
      u_average_min=1.e-6,
      gradients=gradients)
    gradients -= gradients_in
    assert approx_equal(
      registry_residual_sum,
      adp_iso_residual_sum(
        weight=weight, average_power=average_power, u_isos=u_isos))
    assert approx_equal(
      gradients,
      adp_iso_analytical_gradients(
        weight=weight, average_power=average_power, u_isos=u_isos))

if (__name__ == "__main__"):
  exercise_pair_registry_basic()
  exercise_adp_iso_analytical()
  exercise_pair_registry_adp_iso()
  print format_cpu_times()
