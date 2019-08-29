from __future__ import absolute_import, division, print_function
from mmtbx import ncs
import mmtbx.ncs.cartesian_restraints
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal
from libtbx.utils import format_cpu_times
from six.moves import zip
from six.moves import range

def exercise_pair_registry_basic():
  registry = ncs.cartesian_restraints.pair_registry(30, 3)
  assert registry.n_seq() == 30
  assert registry.number_of_additional_isolated_sites == 0
  selection_pairs = registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == []
  assert list(zip(*selection_pairs[1])) == []
  registry = ncs.cartesian_restraints.pair_registry(30, 3)
  assert registry.enter(i_seq=0, j_seq=20, j_ncs=1) == (0, 0)
  assert registry.enter(i_seq=0, j_seq=20, j_ncs=1) == (0, 1)
  assert registry.enter(i_seq=20, j_seq=0, j_ncs=1) == (2, 1)
  assert registry.enter(i_seq=20, j_seq=1, j_ncs=1) == (2, 1)
  assert registry.enter(i_seq=21, j_seq=1, j_ncs=1) == (0, 0)
  assert registry.enter(i_seq=21, j_seq=2, j_ncs=1) == (3, 1)
  assert registry.enter(i_seq=0, j_seq=20, j_ncs=2) == (1, 1)
  assert registry.enter(i_seq=0, j_seq=10, j_ncs=2) == (0, 0)
  assert registry.enter(i_seq=0, j_seq=11, j_ncs=2) == (3, 10)
  assert registry.enter(i_seq=0, j_seq=10, j_ncs=2) == (0, 1)
  selection_pairs = registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(0, 20), (21, 1)]
  assert list(zip(*selection_pairs[1])) == [(0, 10)]
  selection = flex.bool(30, False)
  sel_registry = registry.proxy_select(iselection=selection.iselection())
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == []
  assert list(zip(*selection_pairs[1])) == []
  selection = flex.bool(30, True)
  sel_registry = registry
  for i in range(10):
    sel_registry = sel_registry.proxy_select(iselection=selection.iselection())
    selection_pairs = sel_registry.selection_pairs()
    assert len(selection_pairs) == 2
    assert list(zip(*selection_pairs[0])) == [(0, 20), (21, 1)]
    assert list(zip(*selection_pairs[1])) == [(0, 10)]
  selection[0] = False
  sel_registry = registry.proxy_select(iselection=selection.iselection())
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(20, 0)]
  assert list(zip(*selection_pairs[1])) == []
  selection[0] = True
  selection[1] = False
  sel_registry = registry.proxy_select(iselection=selection.iselection())
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(0, 19)]
  assert list(zip(*selection_pairs[1])) == [(0, 9)]
  selection[1] = True
  selection[2] = False
  selection[4] = False
  sel_registry = registry.proxy_select(iselection=selection.iselection())
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(0, 18), (19, 1)]
  assert list(zip(*selection_pairs[1])) == [(0, 8)]
  sel_registry = sel_registry.proxy_select(
    iselection=flex.size_t([0,1,8,18,19,23,27]))
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(0, 3), (4, 1)]
  assert list(zip(*selection_pairs[1])) == [(0, 2)]
  assert sel_registry.enter(i_seq=0, j_seq=3, j_ncs=1) == (0, 1)
  assert sel_registry.enter(i_seq=0, j_seq=2, j_ncs=1) == (1, 2)
  assert sel_registry.enter(i_seq=3, j_seq=0, j_ncs=1) == (2, 1)
  assert sel_registry.enter(i_seq=0, j_seq=4, j_ncs=1) == (3, 3)
  assert sel_registry.enter(i_seq=6, j_seq=5, j_ncs=1) == (0, 0)
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(0, 3), (4, 1), (6, 5)]
  assert list(zip(*selection_pairs[1])) == [(0, 2)]
  sel_registry = sel_registry.proxy_select(
    iselection=flex.size_t([0,2,4,1]))
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(2, 3)]
  assert list(zip(*selection_pairs[1])) == [(0, 1)]
  registry.register_additional_isolated_sites(number=10)
  assert registry.number_of_additional_isolated_sites == 10
  registry.register_additional_isolated_sites(number=3)
  assert registry.number_of_additional_isolated_sites == 13
  selection.resize(43, False)
  sel_registry = registry.proxy_select(iselection=selection.iselection())
  assert sel_registry.number_of_additional_isolated_sites == 0
  selection[35] = True
  sel_registry = registry.proxy_select(iselection=selection.iselection())
  assert sel_registry.number_of_additional_isolated_sites == 1
  selection[38] = True
  selection[42] = True
  sel_registry = registry.proxy_select(iselection=selection.iselection())
  assert sel_registry.number_of_additional_isolated_sites == 3
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(0, 18), (19, 1)]
  assert list(zip(*selection_pairs[1])) == [(0, 8)]
  #
  iselection = flex.size_t([
    # random permutation with omissions first 30
    8, 10, 6, 12, 20, 15, 16, 22, 11, 23, 14, 27,
    3, 13, 5, 28, 17, 0, 19, 26, 9, 4, 1, 29, 18,
    # random permutation with omissions other 13
    32, 37, 34, 41, 38, 42, 33, 35, 30])
  sel_registry = registry.proxy_select(iselection=iselection)
  assert sel_registry.number_of_additional_isolated_sites == 9
  selection_pairs = sel_registry.selection_pairs()
  assert len(selection_pairs) == 2
  assert list(zip(*selection_pairs[0])) == [(17, 4)]
  assert list(zip(*selection_pairs[1])) == [(17, 1)]
  #
  iselection = flex.size_t([0,1,30])
  sel_registry = registry.proxy_select(iselection=iselection)
  assert sel_registry.number_of_additional_isolated_sites == 1
  iselection = flex.size_t([0,30,1])
  try: registry.proxy_select(iselection=iselection)
  except RuntimeError as e:
    assert str(e).endswith(
      "): MMTBX_ASSERT("
      "result_i_seq < result_n_seq || iselection[result_i_seq] >= n_seq)"
      " failure.")
  else:
    raise Exception_expected

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
  for i_u_iso in range(u_isos.size()):
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
      for size in range(2,20):
        for i in range(10):
          u_isos = mersenne_twister.random_double(size=size) + 1.e-3
          a = adp_iso_analytical_gradients(
            weight=weight, average_power=average_power, u_isos=u_isos)
          f = adp_iso_finite_difference_gradients(
            weight=weight, average_power=average_power, u_isos=u_isos)
          assert approx_equal(a, f)

def exercise_pair_registry_adp_iso():
  mersenne_twister = flex.mersenne_twister(seed=0)
  for n_seq in range(2,20):
    registry = ncs.cartesian_restraints.pair_registry(n_seq=n_seq, n_ncs=n_seq)
    for j_seq in range(1,n_seq):
      assert registry.enter(i_seq=0, j_seq=j_seq, j_ncs=j_seq) == (0, 0)
    selection_pairs = registry.selection_pairs()
    for j in range(1,n_seq):
      assert list(zip(*selection_pairs[j-1])) == [(0,j)]
    weight = 2.134
    average_power = 0.589
    u_isos = mersenne_twister.random_double(size=n_seq) + 1.e-3
    gradients_in = mersenne_twister.random_double(size=n_seq)
    gradients = gradients_in.deep_copy()
    registry_residual_sum = registry.adp_iso_residual_sum(
      weight=weight,
      average_power=average_power,
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

def exercise():
  exercise_pair_registry_basic()
  exercise_adp_iso_analytical()
  exercise_pair_registry_adp_iso()
  print(format_cpu_times())

if (__name__ == "__main__"):
  exercise()
