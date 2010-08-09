from __future__ import division
from cctbx.array_family import flex
from cctbx import geometry_restraints
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.crystal import direct_space_asu
from scitbx import matrix
from scitbx import stl
from libtbx.test_utils import \
  Exception_expected, approx_equal, not_approx_equal, eps_eq, show_diff
from libtbx.utils import null_out
from cStringIO import StringIO
import math
import sys

def finite_difference_gradients(restraint_type, sites_cart, proxy, unit_cell=None, eps=1.e-8):
  def residual(restraint_type, sites_cart, proxy, unit_cell):
    if unit_cell is None:
      return restraint_type(
        sites_cart=sites_cart,
        proxy=proxy).residual()
    else:
      return restraint_type(
        unit_cell=unit_cell,
        sites_cart=sites_cart,
        proxy=proxy).residual()
  result = []
  for i in xrange(len(sites_cart)):
    result_i = []
    for j in xrange(3):
      h = [0,0,0]
      h[j] = eps
      h = matrix.col(h)
      sites_cart[i]=matrix.col(sites_cart[i]) + h
      qp = residual(restraint_type,sites_cart,proxy,unit_cell)
      sites_cart[i]=matrix.col(sites_cart[i]) - 2*h
      qm = residual(restraint_type,sites_cart,proxy,unit_cell)
      dq = (qp-qm)/2
      result_i.append(dq/(eps))
    result.append(result_i)
  return result

def exercise_bond_similarity():
  # test without symmetry operations
  i_seqs=((0,2),
          (1,3),
          (4,5))
  weights=(1,2,3)
  p = geometry_restraints.bond_similarity_proxy(
    i_seqs=i_seqs,
    weights=weights)
  assert tuple(p.i_seqs) == i_seqs
  assert approx_equal(p.weights, weights)
  assert p.sym_ops == None
  #
  expected_deltas = \
    (-0.033333333333333, 0.066666666666666, -0.033333333333333)
  expected_rms_deltas = math.sqrt(
    sum([delta * delta for delta in expected_deltas])
    /len(expected_deltas))
  expected_residual = sum([weights[i] * expected_deltas[i]
                          * expected_deltas[i]
                          for i in range(3)])\
                    / sum([w for w in weights])
  expected_gradients = (
    ((0,0,0.011111111111), (0,0,-0.011111111111)),
    ((0,-0.044444444444,0), (0,0.044444444444,0)),
    ((0.033333333333,0,0), (-0.033333333333,0,0)))
  sites_array=[
    ((1,2,3),(1,2,4.5)),((2,4,6),(2,5.6,6)),((4,14,19),(5.5,14,19))]
  b = geometry_restraints.bond_similarity(
    sites_array=sites_array,
    weights=weights)
  assert approx_equal(b.sites_array, sites_array)
  assert approx_equal(b.weights, weights)
  assert approx_equal(b.mean_distance(), 1.533333333333333)
  assert approx_equal(b.deltas(), expected_deltas)
  assert approx_equal(b.rms_deltas(), expected_rms_deltas)
  assert approx_equal(b.residual(), expected_residual)
  assert approx_equal(b.gradients(), expected_gradients)
  #
  sites_cart = flex.vec3_double(
    [(1,2,3),(2,4,6),(1,2,4.5),(2,5.6,6),(4,14,19),(5.5,14,19)])
  b = geometry_restraints.bond_similarity(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(b.sites_array, sites_array)
  assert approx_equal(b.weights, weights)
  assert approx_equal(b.mean_distance(), 1.533333333333333)
  assert approx_equal(b.deltas(), expected_deltas)
  assert approx_equal(b.rms_deltas(), expected_rms_deltas)
  assert approx_equal(b.residual(), expected_residual)
  assert approx_equal(b.gradients(), expected_gradients)
  # test with symmetry operations
  unit_mx = sgtbx.rt_mx()
  i_seqs=((0,2),
          (1,3),
          (4,5))
  sym_ops=(unit_mx,
           unit_mx,
           sgtbx.rt_mx('1+x,y,z'),)
  weights=(1,2,3)
  p = geometry_restraints.bond_similarity_proxy(
    i_seqs=i_seqs,
    sym_ops=sym_ops,
    weights=weights)
  assert tuple(p.i_seqs) == i_seqs
  assert tuple(p.sym_ops) == sym_ops
  assert approx_equal(p.weights, weights)
  #
  expected_deltas = \
    (-0.033333333333333, 0.066666666666666, -0.033333333333333)
  expected_rms_deltas = math.sqrt(
    sum([delta * delta for delta in expected_deltas])
    /len(expected_deltas))
  expected_residual = sum([weights[i] * expected_deltas[i]
                          * expected_deltas[i]
                          for i in range(3)])\
                    / sum([w for w in weights])
  expected_gradients = (
    ((0,0,0.011111111111), (0,0,-0.011111111111)),
    ((0,-0.044444444444,0), (0,0.044444444444,0)),
    ((0.033333333333,0,0), (-0.033333333333,0,0)))
  sites_array=[
    ((1,2,3),(1,2,4.5)),((2,4,6),(2,5.6,6)),((14,24,29),(15.5,24,29))]
  b = geometry_restraints.bond_similarity(
    sites_array=sites_array,
    weights=weights)
  assert approx_equal(b.sites_array, sites_array)
  assert approx_equal(b.weights, weights)
  assert approx_equal(b.mean_distance(), 1.533333333333333)
  assert approx_equal(b.deltas(), expected_deltas)
  assert approx_equal(b.rms_deltas(), expected_rms_deltas)
  assert approx_equal(b.residual(), expected_residual)
  assert approx_equal(b.gradients(), expected_gradients)
  #
  unit_cell = uctbx.unit_cell([15,25,30,90,90,90])
  sites_cart = flex.vec3_double(
    [(1,2,3),(2,4,6),(1,2,4.5),(2,5.6,6),(14,24,29),(0.5,24,29)])
  b = geometry_restraints.bond_similarity(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(b.sites_array, sites_array)
  assert approx_equal(b.weights, weights)
  assert approx_equal(b.mean_distance(), 1.533333333333333)
  assert approx_equal(b.deltas(), expected_deltas)
  assert approx_equal(b.rms_deltas(), expected_rms_deltas)
  assert approx_equal(b.residual(), expected_residual)
  assert approx_equal(b.gradients(), expected_gradients)
  #
  proxies = geometry_restraints.shared_bond_similarity_proxy([p,p])
  assert eps_eq(geometry_restraints.bond_similarity_residuals(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [expected_residual]*2)
  assert eps_eq(geometry_restraints.bond_similarity_deltas_rms(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [expected_rms_deltas]*2)
  residual_sum = geometry_restraints.bond_similarity_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert eps_eq(residual_sum, 2*expected_residual)
  gradient_array = flex.vec3_double(sites_cart.size(), (0,0,0))
  residual_sum = geometry_restraints.bond_similarity_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=gradient_array)
  assert eps_eq(residual_sum, 2*expected_residual)
  fd_grads = finite_difference_gradients(
    restraint_type=geometry_restraints.bond_similarity,
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p)
  for g,e in zip(gradient_array, fd_grads):
    assert approx_equal(g, matrix.col(e)*2)
  # check proxies with and without sym_ops are happy side-by-side
  p_sym = geometry_restraints.bond_similarity_proxy(
    i_seqs=i_seqs,
    sym_ops=sym_ops,
    weights=weights)
  assert p_sym.sym_ops == sym_ops
  restraint_sym = geometry_restraints.bond_similarity(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p_sym)
  p_no_sym = geometry_restraints.bond_similarity_proxy(
    i_seqs=i_seqs,
    weights=weights)
  assert p_no_sym.sym_ops == None
  restraint_no_sym = geometry_restraints.bond_similarity(
    sites_cart=sites_cart,
    proxy=p_no_sym)
  proxies = geometry_restraints.shared_bond_similarity_proxy([p_sym,p_no_sym])
  assert approx_equal(geometry_restraints.bond_similarity_deltas_rms(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [restraint_sym.rms_deltas(),restraint_no_sym.rms_deltas()])
  assert approx_equal(geometry_restraints.bond_similarity_residuals(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [restraint_sym.residual(),restraint_no_sym.residual()])
  residual_sum = geometry_restraints.bond_similarity_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, restraint_sym.residual() + restraint_no_sym.residual())

def exercise_bond():
  p = geometry_restraints.bond_params(
    distance_ideal=3.5,
    weight=1)
  assert approx_equal(p.distance_ideal, 3.5)
  assert approx_equal(p.weight, 1)
  assert approx_equal(p.slack, 0)
  p.distance_ideal = 35
  assert approx_equal(p.distance_ideal, 35)
  p.distance_ideal = 3.5
  assert approx_equal(p.distance_ideal, 3.5)
  p.weight = 10
  assert approx_equal(p.weight, 10)
  p.weight = 1
  assert approx_equal(p.weight, 1)
  p.slack = 3
  assert approx_equal(p.slack, 3)
  p.slack = 0
  assert approx_equal(p.slack, 0)
  #
  c = p.scale_weight(factor=2)
  assert c.distance_ideal == p.distance_ideal
  assert approx_equal(c.weight, 2)
  assert c.slack == p.slack
  #
  t = geometry_restraints.bond_params_table()
  assert t.size() == 0
  d = geometry_restraints.bond_params_dict()
  assert len(d) == 0
  p = geometry_restraints.bond_params(distance_ideal=3, weight=2)
  d[10] = p
  assert approx_equal(d[10].distance_ideal, 3)
  t.append(d)
  t.append(d)
  assert approx_equal(t[1][10].distance_ideal, 3)
  t[0][13] = p
  assert approx_equal(t[0][13].distance_ideal, 3)
  t[0][13].distance_ideal = 5
  assert approx_equal(t[0][13].distance_ideal, 5)
  assert approx_equal(t[1][10].distance_ideal, 3)
  t[1][1] = geometry_restraints.bond_params(distance_ideal=4, weight=5)
  while (t.size() < 14):
    t.append(geometry_restraints.bond_params_dict())
  s = t.proxy_select(iselection=flex.size_t([1]))
  assert approx_equal(s[0][0].distance_ideal, 4)
  assert approx_equal(s[0][0].weight, 5)
  #
  rest = t.proxy_remove(selection=flex.bool([True]*14))
  assert [p.size() for p in rest] == [0]*14
  rest = t.proxy_remove(selection=flex.bool([False]*14))
  assert [p.size() for p in rest] == [2,2]+[0]*12
  rest = t.proxy_remove(selection=flex.bool([False,True]+[False]*12))
  assert [p.size() for p in rest] == [2,1]+[0]*12
  rest = t.proxy_remove(
    selection=flex.bool([True,False]+[False]*8+[True]+[False]*3))
  assert [p.size() for p in rest] == [1,2]+[0]*12
  rest = t.proxy_remove(
    selection=flex.bool([True,True]+[False]*8+[True]+[False]*3))
  assert [p.size() for p in rest] == [1]+[0]*13
  #
  p = geometry_restraints.bond_params(distance_ideal=2.8, weight=2)
  assert t[3].size() == 0
  t.update(i_seq=4, j_seq=3, params=p)
  assert t[3].size() == 1
  assert approx_equal(t[3][4].distance_ideal, 2.8)
  p = geometry_restraints.bond_params(distance_ideal=3.8, weight=3)
  t.update(i_seq=3, j_seq=5, params=p)
  assert t[3].size() == 2
  assert approx_equal(t[3][5].distance_ideal, 3.8)
  p = geometry_restraints.bond_params(distance_ideal=1.8, weight=4)
  t.update(i_seq=3, j_seq=4, params=p)
  assert t[3].size() == 2
  assert approx_equal(t[3][4].distance_ideal, 1.8)
  #
  assert geometry_restraints.bond_params_table().mean_residual(1) == 0
  assert t.mean_residual(bond_stretch_factor=0) == 0
  assert approx_equal(t.mean_residual(0.5), 9.26166666667)
  assert approx_equal(t.mean_residual(1.0), 37.0466666667)
  assert approx_equal(t.mean_residual(2.0), 148.186666667)
  #
  p = geometry_restraints.bond_simple_proxy(
    i_seqs=[1,0],
    distance_ideal=3.5,
    weight=1)
  assert p.i_seqs == (1,0)
  assert approx_equal(p.distance_ideal, 3.5)
  assert approx_equal(p.weight, 1)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1)
  assert approx_equal(p.distance_ideal, 3.5)
  assert approx_equal(p.weight, 1)
  p.distance_ideal = 35
  assert approx_equal(p.distance_ideal, 35)
  p.distance_ideal = 3.5
  assert approx_equal(p.distance_ideal, 3.5)
  p.weight = 10
  assert approx_equal(p.weight, 10)
  p.weight = 1
  assert approx_equal(p.weight, 1)
  b = geometry_restraints.bond(
    sites=[(1,2,3),(2,4,6)],
    distance_ideal=3.5,
    weight=1)
  assert approx_equal(b.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(b.distance_ideal, 3.5)
  assert approx_equal(b.weight, 1)
  assert approx_equal(b.distance_model**2, 14)
  assert approx_equal(b.delta, -0.241657386774)
  assert approx_equal(b.residual(), 0.0583982925824)
  assert approx_equal(b.gradients(),
    ((-0.12917130661302928, -0.25834261322605856, -0.38751391983908784),
     ( 0.12917130661302928,  0.25834261322605856,  0.38751391983908784)))
  b = geometry_restraints.bond(
    sites=[(1,2,3),(1,2,3)],
    distance_ideal=3.5,
    weight=1)
  assert approx_equal(b.distance_model, 0)
  assert approx_equal(b.gradients(), [(0,0,0), (0,0,0)])
  sites_cart = flex.vec3_double([(1,2,3),(2,4,6)])
  b = geometry_restraints.bond(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(b.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(b.distance_ideal, 3.5)
  assert approx_equal(b.weight, 1)
  assert approx_equal(b.distance_model**2, 14)
  proxies = geometry_restraints.shared_bond_simple_proxy([p,p])
  for proxy in proxies:
    assert approx_equal(proxy.weight, 1)
    proxy.weight = 12
  for proxy in proxies:
    assert approx_equal(proxy.weight, 12)
    proxy.weight = 1
  tab = geometry_restraints.extract_bond_params(
    n_seq=2, bond_simple_proxies=proxies)
  assert tab[0].keys() == [1]
  assert approx_equal(tab[0].values()[0].distance_ideal, 3.5)
  assert approx_equal(tab[0].values()[0].weight, 1)
  assert tab[1].keys() == []
  assert approx_equal(geometry_restraints.bond_distances_model(
    sites_cart=sites_cart,
    proxies=proxies), [14**0.5]*2)
  assert approx_equal(geometry_restraints.bond_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [-0.241657386774]*2)
  assert approx_equal(geometry_restraints.bond_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [0.0583982925824]*2)
  residual_sum = geometry_restraints.bond_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*0.0583982925824)
  #
  for sym_op in (None, sgtbx.rt_mx("-y,z,x")):
    if sym_op is None:
      p = geometry_restraints.bond_simple_proxy(
        i_seqs=[1,0],
        distance_ideal=3.5,
        weight=1)
    else:
      for proxy_t in (geometry_restraints.bond_simple_proxy,
                      geometry_restraints.bond_sym_proxy):
        p = proxy_t(
          i_seqs=[1,0],
          rt_mx_ji=sgtbx.rt_mx("-y,z,x"),
          distance_ideal=3.5,
          weight=1)
      assert p.rt_mx_ji == sym_op
    assert p.i_seqs == (1,0)
    assert approx_equal(p.distance_ideal, 3.5)
    assert approx_equal(p.weight, 1)
  unit_cell = uctbx.unit_cell([15,25,30,90,90,90])
  sites_cart = flex.vec3_double([[1,2,3],[2,3,4]])
  b = geometry_restraints.bond(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(b.sites, [(2,3,4),(-1.2,2.5,2)])
  assert approx_equal(b.distance_ideal, 3.5)
  assert approx_equal(b.weight, 1)
  assert approx_equal(b.distance_model**2, 14.49)
  #
  sites_cart = flex.vec3_double([[1,2,3],[2,3,4]])
  asu_mappings = direct_space_asu.non_crystallographic_asu_mappings(
    sites_cart=sites_cart)
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=5)
  p = geometry_restraints.bond_asu_proxy(
    pair=pair_generator.next(),
    distance_ideal=2,
    weight=10)
  p = geometry_restraints.bond_asu_proxy(pair=p, params=p)
  assert pair_generator.at_end()
  assert p.i_seq == 0
  assert p.j_seq == 1
  assert p.j_sym == 0
  assert approx_equal(p.distance_ideal, 2)
  assert approx_equal(p.weight, 10)
  p.distance_ideal = 35
  assert approx_equal(p.distance_ideal, 35)
  p.distance_ideal = 2
  assert approx_equal(p.distance_ideal, 2)
  p.weight = 1
  assert approx_equal(p.weight, 1)
  p.weight = 10
  assert approx_equal(p.weight, 10)
  assert p.as_simple_proxy().i_seqs == (0,1)
  assert approx_equal(p.as_simple_proxy().distance_ideal, 2)
  assert approx_equal(p.as_simple_proxy().weight, 10)
  sym_proxies = geometry_restraints.shared_bond_asu_proxy([p,p])
  for proxy in sym_proxies:
    assert approx_equal(proxy.distance_ideal, 2)
    proxy.distance_ideal = -4
  for proxy in sym_proxies:
    assert approx_equal(proxy.distance_ideal, -4)
    proxy.distance_ideal = 2
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  sorted_asu_proxies.process(proxies=sym_proxies)
  assert approx_equal(
    geometry_restraints.bond_distances_model(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies),
    [3**.5]*2)
  assert approx_equal(
    geometry_restraints.bond_deltas(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies),
    [2-3**.5]*2)
  assert approx_equal(
    geometry_restraints.bond_residuals(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies),
    [10*(2-3**.5)**2]*2)
  assert approx_equal(
    geometry_restraints.bond_residual_sum(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies,
      gradient_array=None),
    (10*(2-3**.5)**2)*2)
  gradient_array = flex.vec3_double(2, [0,0,0])
  assert approx_equal(
    geometry_restraints.bond_residual_sum(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies,
      gradient_array=gradient_array),
    (10*(2-3**.5)**2)*2)
  assert approx_equal(gradient_array,
    [[ 6.1880215351700611]*3,
     [-6.1880215351700611]*3])
  for disable_cache in [False, True]:
    gradient_array = flex.vec3_double(2, [0,0,0])
    assert approx_equal(
      geometry_restraints.bond_residual_sum(
        sites_cart=sites_cart,
        sorted_asu_proxies=sorted_asu_proxies,
        gradient_array=gradient_array,
        disable_cache=disable_cache),
      (10*(2-3**.5)**2)*2)
    assert approx_equal(gradient_array,
      [[ 6.1880215351700611]*3,
       [-6.1880215351700611]*3])
  #
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  sorted_asu_proxies.push_back(proxy=sym_proxies[0])
  assert sorted_asu_proxies.simple.size() == 0
  assert sorted_asu_proxies.asu.size() == 1
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  sorted_asu_proxies.process(proxy=proxies[0])
  sorted_asu_proxies.process(proxy=sym_proxies[0])
  assert sorted_asu_proxies.simple.size() == 2
  assert sorted_asu_proxies.asu.size() == 0
  assert sorted_asu_proxies.n_total() == 2
  residual_0 = geometry_restraints.bond(
    sites_cart=sites_cart,
    proxy=proxies[0]).residual()
  residual_1 = geometry_restraints.bond(
    sites_cart=sites_cart,
    asu_mappings=asu_mappings,
    proxy=sym_proxies[0]).residual()
  assert approx_equal(residual_1, 10*(2-3**.5)**2)
  gradient_array = flex.vec3_double(2, [0,0,0])
  assert approx_equal(geometry_restraints.bond_residual_sum(
    sites_cart=sites_cart,
    sorted_asu_proxies=sorted_asu_proxies,
    gradient_array=gradient_array), residual_0+residual_1)
  assert approx_equal(gradient_array,
    [(5.1354626519124107, 5.1354626519124107, 5.1354626519124107),
     (-5.1354626519124107, -5.1354626519124107, -5.1354626519124107)])
  sorted_asu_proxies.process(sorted_asu_proxies.simple.deep_copy())
  assert sorted_asu_proxies.simple.size() == 4
  sorted_asu_proxies.process(sorted_asu_proxies.asu.deep_copy())
  assert sorted_asu_proxies.asu.size() == 0
  sorted_asu_proxies.push_back(sorted_asu_proxies.asu.deep_copy())
  assert sorted_asu_proxies.asu.size() == 0
  #
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  assert pair_asu_table.table()[0].keys() == []
  geometry_restraints.add_pairs(
    pair_asu_table=pair_asu_table,
    bond_simple_proxies=proxies)
  assert pair_asu_table.table()[0].keys() == [1]
  #
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    pair_asu_table=pair_asu_table)
  assert sorted_asu_proxies.simple.size() == 0
  assert sorted_asu_proxies.asu.size() == 0
  pair_asu_table.add_all_pairs(distance_cutoff=2)
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    pair_asu_table=pair_asu_table)
  assert sorted_asu_proxies.simple.size() == 1
  assert sorted_asu_proxies.asu.size() == 0
  #
  def sign(x):
    if (x < 0): return -1
    return 1
  mt = flex.mersenne_twister(seed=0)
  for slack in [0, 1/3., 2/3., 1]:
    for ish in xrange(9):
      sh = ish / 2.
      site1 = matrix.col((1,2,3)) \
            + sh * matrix.col(mt.random_double_point_on_sphere())
      b = geometry_restraints.bond(
        sites=[(1,2,3),site1],
        distance_ideal=2,
        weight=1,
        slack=slack)
      assert approx_equal(b.distance_model, sh)
      assert approx_equal(
        b.delta_slack,
        sign(b.delta) * max(0, (abs(b.delta) - b.slack)))
      #
      for i in xrange(3):
        rs = []
        eps = 1.e-6
        for signed_eps in [eps, -eps]:
          site0 = [1,2,3]
          site0[i] += signed_eps
          be = geometry_restraints.bond(
            sites=[site0,site1],
            distance_ideal=2,
            weight=1,
            slack=slack)
          rs.append(be.residual())
        g_fin = (rs[0]-rs[1])/(2*eps)
        g_ana = b.gradients()[0][i]
        assert approx_equal(g_ana, g_fin)

  #
  sites_cart = flex.vec3_double([(1,2,3),(2,3,4)])
  gradients_inp = flex.vec3_double([(10,20,30),(20,30,40)])
  site_symmetry_table_indices = flex.size_t([0,1])
  home_sites_cart = flex.vec3_double([(1.1,1.9,3.05),(2.1,3.2,4.3)])
  iselection = flex.size_t([0,1])
  def get(weight, slack, no_grads=False, no_site_sym=False):
    gradients = None
    if (not no_grads):
      gradients = gradients_inp.deep_copy()
    site_sym = None
    if (not no_site_sym):
      site_sym = site_symmetry_table_indices
    residual_sum = geometry_restraints \
      .home_restraints_summation_skip_special_positions(
        sites_cart=sites_cart,
        gradients=gradients,
        site_symmetry_table_indices=site_sym,
        home_sites_cart=home_sites_cart,
        iselection=iselection,
        weight=weight,
        slack=slack)
    if (not no_grads):
      return [residual_sum, gradients]
    return residual_sum
  assert approx_equal(get(1, 0),
    [0.0225, [(9.8,20.2,29.9), (20,30,40)]])
  assert approx_equal(get(1, 0, True), 0.0225)
  assert approx_equal(get(1, 0, True, True), 0.1625)
  assert approx_equal(get(1, 0, False, True),
    [0.1625, [(9.8,20.2,29.9), (19.8,29.6,39.4)]])
  assert approx_equal(get(3, 0),
    [0.0675, [(9.4,20.6,29.7), (20,30,40)]])
  assert approx_equal(get(1, 3),
    [0.0, [(10,20,30), (20,30,40)]])
  site_symmetry_table_indices = flex.size_t([1,0])
  assert approx_equal(get(1, 0),
    [0.14, [(10,20,30), (19.8,29.6,39.4)]])
  site_symmetry_table_indices = flex.size_t([0,0])
  iselection = flex.size_t([1])
  assert approx_equal(get(1, 0),
    [0.14, [(10,20,30), (19.8,29.6,39.4)]])
  iselection = flex.size_t([])
  assert approx_equal(get(1, 0),
    [0.0, [(10,20,30), (20,30,40)]])
  iselection = flex.size_t([0,1])

class py_nonbonded_cos(object): # prototype

  def __init__(self, max_residual, exponent=1):
    self.max_residual = max_residual
    self.exponent = exponent

  def residual_and_gradients(self, proxy, sites_cart, gradient_array=None):
    vdw_distance = proxy.vdw_distance
    i, j = proxy.i_seqs
    diff_vec = matrix.col(sites_cart[i]) - matrix.col(sites_cart[j])
    d = abs(diff_vec)
    if (d >= vdw_distance): return 0
    x = d / vdw_distance
    m = self.max_residual
    e = self.exponent
    pi = math.pi
    pix = pi * x
    cpixpo = math.cos(pix) + 1
    result = m * (cpixpo/2)**e
    if (gradient_array is not None and d != 0):
      """
      r=m ((Cos[Pi x]+1)/2)^e
      g=D[r,x]
      """
      drdx = -(e*m*pi*cpixpo**(e-1)*math.sin(pix)) / 2**e
      drdd = drdx / vdw_distance
      gradient_0 = diff_vec * drdd / d
      gradient_array[i] = matrix.col(gradient_array[i]) + gradient_0
      gradient_array[j] = matrix.col(gradient_array[j]) - gradient_0
    return result

def nb_cos_finite_difference_gradients(nbf, proxy, sites_cart, eps=1.e-6):
  result = []
  for i in proxy.i_seqs:
    sc = list(sites_cart[i])
    for c in xrange(3):
      scc0 = sc[c]
      rs = []
      for signed_eps in [eps, -eps]:
        sc[c] = scc0 + signed_eps
        sites_cart[i] = sc
        rs.append(nbf.residual_and_gradients(
          proxy=proxy,
          sites_cart=sites_cart))
      sc[c] = scc0
      result.append((rs[0]-rs[1])/(2*eps))
    sites_cart[i] = sc
  return flex.vec3_double(flex.double(result))

def exercise_nonbonded_cos(verbose=0):
  if (verbose): log = sys.stdout
  else: log = null_out()
  for exponent in [1, 2, 1.5]:
    nbf = py_nonbonded_cos(max_residual=13, exponent=exponent)
    sites_cart = flex.vec3_double([(1,2,3), (0.3,2.4,3.5)])
    proxy = geometry_restraints.nonbonded_simple_proxy(
      i_seqs=(0,1), vdw_distance=2)
    gradient_array = flex.vec3_double(2, (0,0,0))
    r = nbf.residual_and_gradients(
      proxy=proxy,
      sites_cart=sites_cart,
      gradient_array=gradient_array)
    fd_grads = nb_cos_finite_difference_gradients(
      nbf=nbf, proxy=proxy, sites_cart=sites_cart)
    print >> log, list(gradient_array)
    print >> log, list(fd_grads)
    assert approx_equal(gradient_array, fd_grads)
    nc = geometry_restraints.nonbonded_cos(
      sites=list(sites_cart),
      vdw_distance=proxy.vdw_distance,
      function=geometry_restraints.cos_repulsion_function(
        max_residual=nbf.max_residual, exponent=nbf.exponent))
    assert approx_equal(nc.residual(), r)
    print >> log, nc.gradients()
    assert approx_equal(nc.gradients(), gradient_array)
    print >> log
    sc0 = matrix.col(sites_cart[0])
    v01 = matrix.col(sites_cart[1]) - sc0
    v01 *= 1/abs(v01)
    for i in xrange(21,-1,-1):
      for eps in [0, 1.e-3, -1.e-3]:
        d = i/10 + eps
        sites_cart[1] = sc0 + d * v01
        gradient_array = flex.vec3_double(2, (0,0,0))
        r = nbf.residual_and_gradients(
          proxy=proxy,
          sites_cart=sites_cart,
          gradient_array=gradient_array)
        print >> log, "d, r:", d, r
        if (d == 2):
          assert abs(r) <= 1.e-8
        else:
          fd_grads = nb_cos_finite_difference_gradients(
            nbf=nbf, proxy=proxy, sites_cart=sites_cart)
          print >> log, list(gradient_array)
          print >> log, list(fd_grads)
          assert approx_equal(gradient_array, fd_grads)
          nc = geometry_restraints.nonbonded_cos(
            sites=list(sites_cart),
            vdw_distance=proxy.vdw_distance,
            function=geometry_restraints.cos_repulsion_function(
              max_residual=nbf.max_residual, exponent=nbf.exponent))
          assert approx_equal(nc.residual(), r)
          print >> log, nc.gradients()
          assert approx_equal(nc.gradients(), gradient_array)
          print >> log

def exercise_nonbonded():
  params = geometry_restraints.nonbonded_params()
  assert params.distance_table.size() == 0
  assert params.radius_table.size() == 0
  assert approx_equal(params.factor_1_4_interactions, 2/3.)
  assert approx_equal(params.const_shrink_1_4_interactions, 0)
  assert approx_equal(params.default_distance, 0)
  assert approx_equal(params.minimum_distance, 0)
  params = geometry_restraints.nonbonded_params(
    factor_1_4_interactions=.5,
    const_shrink_1_4_interactions=.1,
    default_distance=1,
    minimum_distance=2)
  assert approx_equal(params.factor_1_4_interactions, .5)
  assert approx_equal(params.const_shrink_1_4_interactions, .1)
  assert approx_equal(params.default_distance, 1)
  assert approx_equal(params.minimum_distance, 2)
  params.factor_1_4_interactions = .4
  assert approx_equal(params.factor_1_4_interactions, .4)
  params.const_shrink_1_4_interactions = .2
  assert approx_equal(params.const_shrink_1_4_interactions, .2)
  params.default_distance = .3
  assert approx_equal(params.default_distance, .3)
  params.minimum_distance = .6
  assert approx_equal(params.minimum_distance, .6)
  #
  p = geometry_restraints.nonbonded_simple_proxy(
    i_seqs=[0,1],
    vdw_distance=5)
  assert p.i_seqs == (0,1)
  assert approx_equal(p.vdw_distance, 5)
  r = geometry_restraints.nonbonded_prolsq(
    sites=[(1,2,3),(2,4,6)],
    vdw_distance=5,
    function=geometry_restraints.prolsq_repulsion_function())
  assert approx_equal(r.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(r.vdw_distance, 5)
  assert approx_equal(r.function.c_rep, 16)
  assert approx_equal(r.delta, 3.74165738677)
  assert approx_equal(r.residual(), 40.1158130612)
  assert approx_equal(r.gradients(),
    [(34.081026602378813, 68.162053204757626, 102.24307980713644),
     (-34.081026602378813, -68.162053204757626, -102.24307980713644)])
  sites_cart = flex.vec3_double([(1,2,3),(2,4,6)])
  r = geometry_restraints.nonbonded_prolsq(
    sites_cart=sites_cart,
    proxy=p,
    function=geometry_restraints.prolsq_repulsion_function())
  assert approx_equal(r.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(r.vdw_distance, 5)
  assert approx_equal(r.function.c_rep, 16)
  assert approx_equal(r.delta, 3.74165738677)
  proxies = geometry_restraints.shared_nonbonded_simple_proxy([p,p])
  for proxy in proxies:
    assert approx_equal(proxy.vdw_distance, 5)
  assert approx_equal(geometry_restraints.nonbonded_deltas(
      sites_cart=sites_cart,
      proxies=proxies),
    [3.74165738677]*2)
  assert approx_equal(geometry_restraints.nonbonded_residuals(
      sites_cart=sites_cart,
      proxies=proxies,
      function=geometry_restraints.prolsq_repulsion_function()),
    [40.1158130612]*2)
  assert approx_equal(geometry_restraints.nonbonded_residuals(
      sites_cart=sites_cart,
      proxies=proxies,
      function=geometry_restraints.inverse_power_repulsion_function(10)),
    [1.3363062095621219]*2)
  assert approx_equal(geometry_restraints.nonbonded_residuals(
      sites_cart=sites_cart,
      proxies=proxies,
      function=geometry_restraints.cos_repulsion_function(max_residual=13)),
    [1.9279613709216095]*2)
  assert approx_equal(geometry_restraints.nonbonded_residuals(
      sites_cart=sites_cart,
      proxies=proxies,
      function=geometry_restraints.gaussian_repulsion_function(
        max_residual=12, norm_height_at_vdw_distance=0.2)),
    [4.8725695136639997]*2)
  residual_sum = geometry_restraints.nonbonded_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None,
    function=geometry_restraints.prolsq_repulsion_function())
  assert approx_equal(residual_sum, 2*40.1158130612)
  residual_sum = geometry_restraints.nonbonded_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None,
    function=geometry_restraints.inverse_power_repulsion_function(10))
  assert approx_equal(residual_sum, 2*1.3363062095621219)
  residual_sum = geometry_restraints.nonbonded_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None,
    function=geometry_restraints.cos_repulsion_function(max_residual=13))
  assert approx_equal(residual_sum, 2*1.9279613709216095)
  residual_sum = geometry_restraints.nonbonded_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None,
    function=geometry_restraints.gaussian_repulsion_function(
      max_residual=12, norm_height_at_vdw_distance=0.2))
  assert approx_equal(residual_sum, 2*4.8725695136639997)
  #
  sites_cart = flex.vec3_double([[1,2,3],[2,3,4]])
  asu_mappings = direct_space_asu.non_crystallographic_asu_mappings(
    sites_cart=sites_cart)
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=5)
  p = geometry_restraints.nonbonded_asu_proxy(
    pair=pair_generator.next(),
    vdw_distance=2)
  assert pair_generator.at_end()
  assert p.i_seq == 0
  assert p.j_seq == 1
  assert p.j_sym == 0
  assert approx_equal(p.vdw_distance, 2)
  p.vdw_distance = 3
  assert approx_equal(p.vdw_distance, 3)
  p.vdw_distance = 2
  sym_proxies = geometry_restraints.shared_nonbonded_asu_proxy([p,p])
  for proxy in sym_proxies:
    assert approx_equal(proxy.vdw_distance, 2)
    proxy.vdw_distance = 3
  for proxy in sym_proxies:
    assert approx_equal(proxy.vdw_distance, 3)
    proxy.vdw_distance = 2
  f = geometry_restraints.prolsq_repulsion_function(
    c_rep=1, k_rep=4, irexp=2, rexp=3)
  assert approx_equal(f.c_rep, 1)
  assert approx_equal(f.k_rep, 4)
  assert approx_equal(f.irexp, 2)
  assert approx_equal(f.rexp, 3)
  assert approx_equal(f.residual(vdw_distance=3, delta=2.9), 2492774.43588)
  assert approx_equal(f.residual(vdw_distance=3, delta=3), 2460375.0)
  r = geometry_restraints.nonbonded_prolsq(
    sites=list(sites_cart),
    vdw_distance=p.vdw_distance,
    function=f)
  assert approx_equal(r.function.c_rep, 1)
  assert approx_equal(r.diff_vec, [-1,-1,-1])
  assert approx_equal(r.delta**2, 3)
  assert approx_equal(r.residual(), 226981)
  assert approx_equal(r.gradients(),
    [(22326.0, 22326.0, 22326.0), (-22326.0, -22326.0, -22326.0)])
  f = geometry_restraints.prolsq_repulsion_function()
  assert approx_equal(f.residual(vdw_distance=3, delta=2.9), 0.0016)
  assert approx_equal(f.residual(vdw_distance=3, delta=3), 0)
  r = geometry_restraints.nonbonded_prolsq(
    sites=list(sites_cart),
    vdw_distance=p.vdw_distance,
    function=f)
  assert approx_equal(r.function.c_rep, 16)
  assert approx_equal(r.function.k_rep, 1)
  assert approx_equal(r.function.irexp, 1)
  assert approx_equal(r.function.rexp, 4)
  assert approx_equal(r.diff_vec, [-1,-1,-1])
  assert approx_equal(r.delta**2, 3)
  assert approx_equal(r.residual(), 0.0824764182859)
  assert approx_equal(r.gradients(),
    [(0.71084793153727288, 0.71084793153727288, 0.71084793153727288),
     (-0.71084793153727288, -0.71084793153727288, -0.71084793153727288)])
  #
  assert approx_equal(f.residual(vdw_distance=3, delta=2.9), 0.0016)
  assert approx_equal(f.residual(vdw_distance=3, delta=3), 0)
  for irexp in [1,2,3,4,5]:
    f = geometry_restraints.inverse_power_repulsion_function(
      nonbonded_distance_cutoff=100,
      k_rep=4,
      irexp=irexp)
    assert approx_equal(f.k_rep, 4)
    assert approx_equal(f.irexp, irexp)
    if (irexp == 1):
      assert approx_equal(f.residual(vdw_distance=3, delta=2.9), 4.13793103448)
      assert approx_equal(f.residual(vdw_distance=3, delta=3), 4)
    r = geometry_restraints.nonbonded_inverse_power(
      sites=list(sites_cart),
      vdw_distance=p.vdw_distance,
      function=f)
    assert approx_equal(r.diff_vec, [-1,-1,-1])
    assert approx_equal(r.delta**2, 3)
    assert approx_equal(r.residual(), f.k_rep*p.vdw_distance/r.delta**f.irexp)
    g = -f.irexp*f.k_rep*p.vdw_distance/(r.delta**(f.irexp+1))/r.delta
    assert approx_equal(r.gradients(), [(-g,-g,-g),(g,g,g)])
  for nonbonded_distance_cutoff in [1,2]:
    f = geometry_restraints.inverse_power_repulsion_function(
      nonbonded_distance_cutoff=nonbonded_distance_cutoff,
      k_rep=4,
      irexp=2)
    r = geometry_restraints.nonbonded_inverse_power(
      sites=list(sites_cart),
      vdw_distance=p.vdw_distance,
      function=f)
    if (nonbonded_distance_cutoff == 1):
      assert approx_equal(r.residual(), 0)
      assert approx_equal(r.gradients(), [(0,0,0),(0,0,0)])
    else:
      assert not_approx_equal(r.residual(), 0)
      assert not_approx_equal(r.gradients(), [(0,0,0),(0,0,0)])
  #
  for exponent in [1,2,3]:
    f = geometry_restraints.cos_repulsion_function(
      max_residual=13,
      exponent=exponent)
    assert approx_equal(f.max_residual, 13)
    assert approx_equal(f.exponent, exponent)
    if (exponent == 1):
      assert approx_equal(
        f.residual(vdw_distance=3, delta=2.9), 0.0356076801062)
      assert approx_equal(f.residual(vdw_distance=3, delta=3), 0)
    r = geometry_restraints.nonbonded_cos(
      sites=list(sites_cart),
      vdw_distance=p.vdw_distance,
      function=f)
    assert approx_equal(r.diff_vec, [-1,-1,-1])
    assert approx_equal(r.delta**2, 3)
    pynbc = py_nonbonded_cos(max_residual=f.max_residual, exponent=f.exponent)
    gradient_array = flex.vec3_double(sites_cart.size(), (0,0,0))
    pr = pynbc.residual_and_gradients(
      proxy=geometry_restraints.nonbonded_simple_proxy(
        i_seqs=(0,1), vdw_distance=r.vdw_distance),
      sites_cart=sites_cart,
      gradient_array=gradient_array)
    assert approx_equal(r.residual(), pr)
    assert approx_equal(r.gradients(), gradient_array)
  #
  expected_residuals = iter([
    1.39552236695,
    2.66705891104,
    3.89565157972])
  expected_gradients = iter([
    2.45678379644,
    2.88800522497,
    2.92825483126])
  for norm_height_at_vdw_distance in [0.1, 0.2, 0.3]:
    f = geometry_restraints.gaussian_repulsion_function(
      max_residual=12,
      norm_height_at_vdw_distance=norm_height_at_vdw_distance)
    assert approx_equal(f.max_residual, 12)
    assert approx_equal(
      f.norm_height_at_vdw_distance(), norm_height_at_vdw_distance)
    assert approx_equal(
      f.residual(vdw_distance=3, delta=3), 12*norm_height_at_vdw_distance)
    assert approx_equal(
      f.residual(vdw_distance=3, delta=2.9), expected_residuals.next())
    r = geometry_restraints.nonbonded_gaussian(
      sites=list(sites_cart),
      vdw_distance=p.vdw_distance,
      function=f)
    assert approx_equal(r.diff_vec, [-1,-1,-1])
    assert approx_equal(r.delta**2, 3)
    e = expected_gradients.next()
    assert approx_equal(r.gradients(), [[e]*3,[-e]*3])
  #
  sorted_asu_proxies = geometry_restraints.nonbonded_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  sorted_asu_proxies.process(proxies=sym_proxies)
  assert approx_equal(
    flex.pow2(geometry_restraints.nonbonded_deltas(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies)),
    [3]*2)
  assert approx_equal(
    geometry_restraints.nonbonded_residuals(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies,
      function=geometry_restraints.prolsq_repulsion_function()),
    [0.0824764182859]*2)
  assert approx_equal(
    geometry_restraints.nonbonded_residual_sum(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies,
      gradient_array=None,
      function=geometry_restraints.prolsq_repulsion_function()),
    0.0824764182859*2)
  for disable_cache in [False, True]:
    gradient_array = flex.vec3_double(2, [0,0,0])
    assert approx_equal(
      geometry_restraints.nonbonded_residual_sum(
        sites_cart=sites_cart,
        sorted_asu_proxies=sorted_asu_proxies,
        gradient_array=gradient_array,
        function=geometry_restraints.prolsq_repulsion_function(),
        disable_cache=disable_cache),
      0.0824764182859*2)
    assert approx_equal(gradient_array,
      [(1.4216958630745458, 1.4216958630745458, 1.4216958630745458),
       (-1.4216958630745458, -1.4216958630745458, -1.4216958630745458)])
  #
  sorted_asu_proxies = geometry_restraints.nonbonded_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  sorted_asu_proxies.process(proxy=proxies[0])
  sorted_asu_proxies.process(proxy=sym_proxies[0], sym_excl_flag=False)
  assert sorted_asu_proxies.simple.size() == 2
  assert sorted_asu_proxies.asu.size() == 0
  assert sorted_asu_proxies.n_total() == 2
  residual_0 = geometry_restraints.nonbonded_prolsq(
    sites_cart=sites_cart,
    proxy=proxies[0],
    function=geometry_restraints.prolsq_repulsion_function()).residual()
  residual_1 = geometry_restraints.nonbonded_prolsq(
    sites_cart=sites_cart,
    asu_mappings=asu_mappings,
    proxy=sym_proxies[0],
    function=geometry_restraints.prolsq_repulsion_function()).residual()
  gradient_array = flex.vec3_double(2, [0,0,0])
  assert approx_equal(geometry_restraints.nonbonded_residual_sum(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies,
      gradient_array=gradient_array,
      function=geometry_restraints.prolsq_repulsion_function()),
    residual_0+residual_1)
  assert approx_equal(gradient_array,
    [(1290.2817767146657, 1290.2817767146657, 1290.2817767146657),
     (-1290.2817767146657, -1290.2817767146657, -1290.2817767146657)])
  sorted_asu_proxies.process(sorted_asu_proxies.simple.deep_copy())
  assert sorted_asu_proxies.simple.size() == 4
  sorted_asu_proxies.process(sorted_asu_proxies.asu.deep_copy())
  assert sorted_asu_proxies.asu.size() == 0
  #
  f = geometry_restraints.prolsq_repulsion_function()
  assert approx_equal((f.c_rep, f.k_rep, f.irexp, f.rexp), (16,1,1,4))
  c = f.customized_copy()
  assert approx_equal((c.c_rep, c.k_rep, c.irexp, c.rexp), (16,1,1,4))
  c = f.customized_copy(c_rep=8, k_rep=2, irexp=3, rexp=5)
  assert approx_equal((c.c_rep, c.k_rep, c.irexp, c.rexp), (8,2,3,5))

def exercise_angle():
  # test without symmetry operations
  p = geometry_restraints.angle_proxy(
    i_seqs=[2,1,0],
    angle_ideal=95,
    weight=1)
  def check(p):
    assert p.i_seqs == (2,1,0)
    assert p.sym_ops is None
    assert approx_equal(p.angle_ideal, 95)
    assert approx_equal(p.weight, 1)
    assert p.slack == 0.0
  check(p)
  p = geometry_restraints.angle_proxy(
    i_seqs=[2,1,0],
    sym_ops=None,
    angle_ideal=95,
    weight=1)
  check(p)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2)
  c = geometry_restraints.angle_proxy(
    i_seqs=[3,4,5],
    proxy=p)
  assert c.i_seqs == (3,4,5)
  assert c.sym_ops is None
  assert approx_equal(c.angle_ideal, 95)
  assert approx_equal(c.weight, 1)
  c = p.scale_weight(factor=3.14)
  assert c.i_seqs == (0,1,2)
  assert c.sym_ops is None
  assert approx_equal(c.angle_ideal, 95)
  assert approx_equal(c.weight, 3.14)
  a = geometry_restraints.angle(
    sites=[(1,0,0),(0,0,0),(0,1,0)],
    angle_ideal=95,
    weight=1)
  assert approx_equal(a.sites, [(1,0,0),(0,0,0),(0,1,0)])
  assert approx_equal(a.angle_ideal, 95)
  assert approx_equal(a.weight, 1)
  assert a.have_angle_model
  assert approx_equal(a.angle_model, 90)
  assert approx_equal(a.delta, 5)
  assert approx_equal(a.residual(), 25)
  assert approx_equal(a.gradients(epsilon=1.e-100),
    ((0.0, 572.95779513082323, 0.0),
     (-572.95779513082323, -572.95779513082323, 0.0),
     (572.95779513082323, 0.0, 0.0)))
  assert a.slack == 0.0
  sites_cart = flex.vec3_double([(1,0,0),(0,0,0),(0,1,0)])
  a = geometry_restraints.angle(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(a.sites, [(1,0,0),(0,0,0),(0,1,0)])
  assert approx_equal(a.angle_ideal, 95)
  assert approx_equal(a.weight, 1)
  assert a.have_angle_model
  assert approx_equal(a.angle_model, 90)
  proxies = geometry_restraints.shared_angle_proxy([p,p])
  for proxy in proxies:
    assert approx_equal(proxy.weight, 1)
  assert approx_equal(geometry_restraints.angle_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [5]*2)
  assert approx_equal(geometry_restraints.angle_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [25]*2)
  residual_sum = geometry_restraints.angle_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*25)
  # test with symmetry operations
  sym_ops = (sgtbx.rt_mx('-1+x,+y,+z'),sgtbx.rt_mx(),sgtbx.rt_mx())
  p = geometry_restraints.angle_proxy(
    i_seqs=[2,1,0],
    sym_ops=sym_ops,
    angle_ideal=95,
    weight=1)
  assert p.i_seqs == (2,1,0)
  assert p.sym_ops == sym_ops
  c = geometry_restraints.angle_proxy(
    i_seqs=[3,4,5],
    proxy=p)
  assert c.i_seqs == (3,4,5)
  assert c.sym_ops == sym_ops
  assert approx_equal(c.angle_ideal, 95)
  assert approx_equal(c.weight, 1)
  c = p.scale_weight(factor=5.82)
  assert c.sym_ops == sym_ops
  assert approx_equal(c.weight, 5.82)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2)
  assert p.sym_ops == (sgtbx.rt_mx(),sgtbx.rt_mx(),sgtbx.rt_mx('-1+x,+y,+z'))
  assert approx_equal(p.angle_ideal, 95)
  assert approx_equal(p.weight, 1)
  unit_cell = uctbx.unit_cell([15,25,30,90,90,90])
  sites_cart = flex.vec3_double([(1,0,0),(0,1,0),(14,0,0)])
  a = geometry_restraints.angle(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(a.sites, [(1,0,0),(0,1,0),(-1,0,0)])
  assert approx_equal(a.angle_ideal, 95)
  assert approx_equal(a.weight, 1)
  assert a.have_angle_model
  assert approx_equal(a.angle_model, 90)
  assert approx_equal(a.delta, 5)
  assert approx_equal(a.residual(), a.weight*a.delta**2)
  proxies = geometry_restraints.shared_angle_proxy()
  for i in range(10): proxies.append(p)
  assert approx_equal(geometry_restraints.angle_deltas(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [5]*10)
  assert approx_equal(geometry_restraints.angle_residuals(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [25]*10)
  gradient_array = flex.double([0]*10)
  residual_sum = geometry_restraints.angle_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 10*25)
  # check proxies with and without sym_ops are happy side-by-side
  p_sym = geometry_restraints.angle_proxy(
    i_seqs=[2,1,0],
    sym_ops=sym_ops,
    angle_ideal=95,
    weight=1)
  assert p_sym.sym_ops == sym_ops
  angle_sym = geometry_restraints.angle(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p_sym)
  p_no_sym = geometry_restraints.angle_proxy(
    i_seqs=[2,1,0],
    angle_ideal=95,
    weight=1)
  assert p_no_sym.sym_ops == None
  angle_no_sym = geometry_restraints.angle(
    sites_cart=sites_cart,
    proxy=p_no_sym)
  proxies = geometry_restraints.shared_angle_proxy([p_sym,p_no_sym])
  assert approx_equal(geometry_restraints.angle_deltas(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [angle_sym.delta,angle_no_sym.delta])
  assert approx_equal(geometry_restraints.angle_residuals(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [angle_sym.residual(),angle_no_sym.residual()])
  residual_sum = geometry_restraints.angle_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, angle_sym.residual() + angle_no_sym.residual())
  #
  proxies = geometry_restraints.shared_angle_proxy([
    geometry_restraints.angle_proxy([0,1,2], 1, 2),
    geometry_restraints.angle_proxy([1,2,3], 2, 3),
    geometry_restraints.angle_proxy([2,3,0], 3, 4),
    geometry_restraints.angle_proxy([3,1,2], 4, 5)])
  selected = proxies.proxy_select(n_seq=4, iselection=flex.size_t([0,2]))
  assert selected.size() == 0
  selected = proxies.proxy_select(n_seq=4, iselection=flex.size_t([0,1,2]))
  assert selected.size() == 1
  assert selected[0].i_seqs == (0,1,2)
  selected = proxies.proxy_select(n_seq=4, iselection=flex.size_t([0,2,3]))
  assert selected.size() == 1
  assert selected[0].i_seqs == (1,2,0)
  selected = proxies.proxy_select(n_seq=4, iselection=flex.size_t([1,2,3]))
  assert selected.size() == 2
  assert selected[0].i_seqs == (0,1,2)
  assert selected[1].i_seqs == (2,0,1)
  assert approx_equal(selected[0].angle_ideal, 2)
  assert approx_equal(selected[1].angle_ideal, 4)
  assert approx_equal(selected[0].weight, 3)
  assert approx_equal(selected[1].weight, 5)
  selected = proxies.proxy_select(n_seq=4, iselection=flex.size_t([0,1,2,3]))
  assert selected.size() == 4
  #
  rest = proxies.proxy_remove(selection=flex.bool([True,True,True,True]))
  assert rest.size() == 0
  rest = proxies.proxy_remove(selection=flex.bool([False,True,True,True]))
  assert rest.size() == 2
  assert rest[0].i_seqs == (0,1,2)
  assert rest[1].i_seqs == (2,3,0)
  rest = proxies.proxy_remove(selection=flex.bool([True,True,True,False]))
  assert rest.size() == 3
  #
  unit_cell = uctbx.unit_cell([15,11.5,16.25,90,99.5,90])
  sites_cart = flex.vec3_double(
    [(12.87,0.10,9.04),(12.54,0.44,7.73),(13.47,0.34,6.71)])
  rt_mx = sgtbx.rt_mx('2-X,-Y,1-Z')
  u_mx = sgtbx.rt_mx()
  i_seqs = [2,0,1]
  a_sym_ops = [u_mx,rt_mx,rt_mx]
  b_sym_ops = [rt_mx,u_mx,u_mx]
  a_proxy = geometry_restraints.angle_proxy(
    i_seqs=i_seqs,
    sym_ops=a_sym_ops,
    angle_ideal=120,
    weight=1)
  b_proxy = geometry_restraints.angle_proxy(
    i_seqs=i_seqs,
    sym_ops=b_sym_ops,
    angle_ideal=120,
    weight=1)
  a_angle = geometry_restraints.angle(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=a_proxy)
  b_angle = geometry_restraints.angle(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=b_proxy)
  assert approx_equal(a_angle.delta, b_angle.delta)
  assert approx_equal(a_angle.residual(), b_angle.residual())
  assert not_approx_equal(a_angle.gradients(), b_angle.gradients())
  assert approx_equal(a_angle.grads_and_curvs(), [
    (82.47780595428577, 9.9176398449691465, -275.30239023150972),
    (-287.7152475786869, 62.177141285587624, 345.71504772803047),
    (205.23744162440113, -72.094781130556768, -70.412657496520751),
    (548.64354849364315, -141.3499058192476, 4777.1475618734567),
    (5313.4607140709404, -424.64322327702632, 7689.6945986452274),
    (2641.6450110035271, 222.22748382126554, 413.72502597622014)])
  a_gradient_array = flex.vec3_double(sites_cart.size())
  a_residual_sum = geometry_restraints.angle_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=geometry_restraints.shared_angle_proxy([a_proxy]),
    gradient_array=a_gradient_array)
  b_gradient_array = flex.vec3_double(sites_cart.size())
  b_residual_sum = geometry_restraints.angle_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=geometry_restraints.shared_angle_proxy([b_proxy]),
    gradient_array=b_gradient_array)
  for a,b in zip(a_gradient_array, b_gradient_array):
    assert approx_equal(a, b)
  fd_grads = finite_difference_gradients(
    restraint_type=geometry_restraints.angle,
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=a_proxy)
  for g,e in zip(a_gradient_array, fd_grads):
    assert approx_equal(g, e)
  #exercise slack parameter
  a = geometry_restraints.angle(
    sites=[(1,0,0),(0,0,0),(0,1,0)],
    angle_ideal=95,
    weight=1,
    slack=0)
  assert approx_equal(a.angle_model, 90)
  assert approx_equal(a.delta, 5)
  assert approx_equal(a.delta_slack, 5)
  assert approx_equal(a.residual(), 25)
  a = geometry_restraints.angle(
    sites=[(1,0,0),(0,0,0),(0,1,0)],
    angle_ideal=95,
    weight=1,
    slack=3)
  assert approx_equal(a.delta, 5)
  assert approx_equal(a.delta_slack, 2)
  a = geometry_restraints.angle(
    sites=[(1,0,0),(0,0,0),(0,1,0)],
    angle_ideal=95,
    weight=1,
    slack=5)
  assert approx_equal(a.delta, 5)
  assert approx_equal(a.delta_slack, 0)
  a = geometry_restraints.angle(
    sites=[(1,0,0),(0,0,0),(0,1,0)],
    angle_ideal=95,
    weight=1,
    slack=7)
  assert approx_equal(a.delta, 5)
  assert approx_equal(a.delta_slack, 0)
  a = geometry_restraints.angle(
    sites=[(1,0,0),(0,0,0),(0,1,0)],
    angle_ideal=85,
    weight=1,
    slack=3)
  assert approx_equal(a.delta, -5)
  assert approx_equal(a.delta_slack, -2)
  a = geometry_restraints.angle(
    sites=[(1,0,0),(0,0,0),(0,1,0)],
    angle_ideal=85,
    weight=1,
    slack=0)
  assert approx_equal(a.delta, -5)
  assert approx_equal(a.delta_slack, -5)

def exercise_dihedral():
  u_mx = sgtbx.rt_mx() # unit matrix
  sym_ops = (u_mx, sgtbx.rt_mx('1+X,+Y,+Z'), u_mx, sgtbx.rt_mx('+X,-1+Y,2+Z'))
  p = geometry_restraints.dihedral_proxy(
    i_seqs=[3,2,1,0],
    sym_ops=sym_ops,
    angle_ideal=-40,
    weight=1,
    periodicity=2,
    alt_angle_ideals=[180])
  assert p.i_seqs == (3,2,1,0)
  assert p.sym_ops == sym_ops
  assert p.alt_angle_ideals == (180,)
  c = geometry_restraints.dihedral_proxy(
    i_seqs=[6,8,5,3],
    proxy=p)
  assert c.i_seqs == (6,8,5,3)
  assert c.sym_ops == sym_ops
  assert approx_equal(c.angle_ideal, -40)
  assert approx_equal(c.weight, 1)
  assert c.periodicity == 2
  assert c.alt_angle_ideals == (180,)
  c = p.scale_weight(factor=0.37)
  assert c.i_seqs == (3,2,1,0)
  assert c.sym_ops == sym_ops
  assert approx_equal(c.weight, 0.37)
  assert c.alt_angle_ideals == (180,)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2,3)
  assert p.sym_ops == (
    sgtbx.rt_mx('+X,-1+Y,2+Z'), u_mx, sgtbx.rt_mx('1+X,+Y,+Z'), u_mx)
  assert approx_equal(p.angle_ideal, -40)
  assert approx_equal(p.weight, 1)
  assert p.periodicity == 2
  assert p.alt_angle_ideals == (180,)
  p.angle_ideal = 50
  p.weight = 2
  p.periodicity = 3
  p.alt_angle_ideals = None
  assert approx_equal(p.angle_ideal, 50)
  assert approx_equal(p.weight, 2)
  assert p.periodicity == 3
  assert p.alt_angle_ideals is None
  p.angle_ideal = -40
  p.weight = 1
  p.periodicity = 2
  p.alt_angle_ideals = (25,-25)
  assert approx_equal(p.angle_ideal, -40)
  assert approx_equal(p.weight, 1)
  assert p.periodicity == 2
  assert p.alt_angle_ideals == (25,-25,)
  #
  u_mx = sgtbx.rt_mx() # unit matrix
  sym_ops = (u_mx, u_mx, sgtbx.rt_mx('+X,1+Y,+Z'), sgtbx.rt_mx('+X,1+Y,+Z'))
  unit_cell = uctbx.unit_cell([15,25,30,90,90,90])
  sites_cart = flex.vec3_double([(2,24,0),(1,24,0),(1,1,0),(0,1,0)])
  for periodicity in [0, 1, -1]:
    p = geometry_restraints.dihedral_proxy(
      i_seqs=[0,1,2,3],
      sym_ops=sym_ops,
      angle_ideal=175,
      weight=1,
      periodicity=periodicity)
    d = geometry_restraints.dihedral(
      unit_cell=unit_cell,
      sites_cart=sites_cart,
      proxy=p)
    assert approx_equal(d.sites, [(2,24,0),(1,24,0),(1,26,0),(0,26,0)])
    assert approx_equal(d.angle_ideal, 175)
    assert approx_equal(d.weight, 1)
    assert d.have_angle_model
    assert approx_equal(d.angle_model, 180)
    assert approx_equal(d.delta, -5)
    if (periodicity <= 0):
      assert approx_equal(d.residual(), 25)
      assert approx_equal(d.gradients(epsilon=1.e-100),
        ((0, 0, 572.95779513082323),
         (0, 0, -572.95779513082323),
         (0, 0, -572.95779513082323),
         (0, 0, 572.95779513082323)))
    else:
      assert approx_equal(d.residual(), 36.5308983192)
      assert approx_equal(d.gradients(epsilon=1.e-100),
        ((0.0, 0.0, 836.69513037751835),
         (0.0, 0.0, -836.69513037751835),
         (0.0, 0.0, -836.69513037751835),
         (0.0, 0.0, 836.69513037751835)))
  #
  p = geometry_restraints.dihedral_proxy(
    i_seqs=[0,1,2,3],
    sym_ops=sym_ops,
    angle_ideal=175,
    weight=1,
    periodicity=0)
  proxies = geometry_restraints.shared_dihedral_proxy([p,p])
  assert proxies.count_harmonic() == 2
  assert approx_equal(geometry_restraints.dihedral_deltas(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [-5]*2)
  assert approx_equal(geometry_restraints.dihedral_residuals(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [25]*2)
  residual_sum = geometry_restraints.dihedral_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*25)
  # check proxies with and without sym_ops are happy side-by-side
  p_sym = geometry_restraints.dihedral_proxy(
    i_seqs=[0,1,2,3],
    sym_ops=sym_ops,
    angle_ideal=175,
    weight=1)
  assert p_sym.sym_ops == sym_ops
  dihedral_sym = geometry_restraints.dihedral(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p_sym)
  p_no_sym = geometry_restraints.dihedral_proxy(
    i_seqs=[0,1,2,3],
    angle_ideal=175,
    weight=1)
  assert p_no_sym.sym_ops == None
  dihedral_no_sym = geometry_restraints.dihedral(
    sites_cart=sites_cart,
    proxy=p_no_sym)
  proxies = geometry_restraints.shared_dihedral_proxy([p_sym,p_no_sym])
  assert approx_equal(geometry_restraints.dihedral_deltas(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [dihedral_sym.delta,dihedral_no_sym.delta])
  assert approx_equal(geometry_restraints.dihedral_residuals(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [dihedral_sym.residual(),dihedral_no_sym.residual()])
  residual_sum = geometry_restraints.dihedral_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(
    residual_sum,
    dihedral_sym.residual() + dihedral_no_sym.residual())
  #
  unit_cell = uctbx.unit_cell([15,11.5,16.25,90,99.5,90])
  sites_cart = flex.vec3_double(
    [(12.87,0.10,9.04),(12.54,0.44,7.73),(13.47,0.34,6.71)])
  rt_mx = sgtbx.rt_mx('2-X,-Y,1-Z')
  u_mx = sgtbx.rt_mx()
  i_seqs = [2,0,1,2]
  a_sym_ops = [u_mx,u_mx,u_mx,rt_mx]
  b_sym_ops = [u_mx,rt_mx,rt_mx,rt_mx]
  a_proxy = geometry_restraints.dihedral_proxy(
    i_seqs=i_seqs,
    sym_ops=a_sym_ops,
    angle_ideal=0,
    weight=1)
  b_proxy = geometry_restraints.dihedral_proxy(
    i_seqs=i_seqs,
    sym_ops=b_sym_ops,
    angle_ideal=0,
    weight=1)
  a_dihedral = geometry_restraints.dihedral(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=a_proxy)
  b_dihedral = geometry_restraints.dihedral(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=b_proxy)
  assert approx_equal(a_dihedral.delta, b_dihedral.delta)
  assert approx_equal(a_dihedral.residual(), b_dihedral.residual())
  assert not_approx_equal(a_dihedral.gradients(), b_dihedral.gradients())
  a_gradient_array = flex.vec3_double(sites_cart.size())
  a_residual_sum = geometry_restraints.dihedral_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=geometry_restraints.shared_dihedral_proxy([a_proxy]),
    gradient_array=a_gradient_array)
  b_gradient_array = flex.vec3_double(sites_cart.size())
  b_residual_sum = geometry_restraints.dihedral_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=geometry_restraints.shared_dihedral_proxy([b_proxy]),
    gradient_array=b_gradient_array)
  for a,b in zip(a_gradient_array, b_gradient_array):
    assert approx_equal(a, b)
  fd_grads = finite_difference_gradients(
    restraint_type=geometry_restraints.dihedral,
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=a_proxy)
  for g,e in zip(a_gradient_array, fd_grads):
    assert approx_equal(g, e)
  #
  p = geometry_restraints.dihedral_proxy(
    i_seqs=[3,2,1,0],
    angle_ideal=-40,
    weight=1,
    periodicity=-2,
    alt_angle_ideals=None)
  assert p.i_seqs == (3,2,1,0)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2,3)
  assert approx_equal(p.angle_ideal, -40)
  assert approx_equal(p.weight, 1)
  assert p.periodicity == -2
  assert p.alt_angle_ideals is None
  p.angle_ideal = 50
  p.weight = 2
  p.periodicity = 3
  assert approx_equal(p.angle_ideal, 50)
  assert approx_equal(p.weight, 2)
  assert p.periodicity == 3
  p.angle_ideal = -40
  p.weight = 1
  p.periodicity = -2
  assert approx_equal(p.angle_ideal, -40)
  assert approx_equal(p.weight, 1)
  assert p.periodicity == -2
  d = geometry_restraints.dihedral(
    sites=[(1,0,0),(0,0,0),(0,1,0),(1,0,1)],
    angle_ideal=-40,
    weight=1)
  assert approx_equal(d.sites, [(1,0,0),(0,0,0),(0,1,0),(1,0,1)])
  assert approx_equal(d.angle_ideal, -40)
  assert approx_equal(d.weight, 1)
  assert d.have_angle_model
  assert approx_equal(d.angle_model, -45)
  assert approx_equal(d.delta, 5)
  assert approx_equal(d.residual(), 25)
  assert approx_equal(d.gradients(epsilon=1.e-100),
    ((0, 0, -572.95779513082323),
     (286.47889756541161, 0, 286.47889756541161),
     (0, 0, 0),
     (-286.47889756541161, 0, 286.47889756541161)))
  sites_cart = flex.vec3_double([(1,0,0),(0,0,0),(0,1,0),(-1,0,-1)])
  d = geometry_restraints.dihedral(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(d.sites, [(1,0,0),(0,0,0),(0,1,0),(-1,0,-1)])
  assert approx_equal(d.angle_ideal, -40)
  assert approx_equal(d.weight, 1)
  assert d.periodicity == -2
  assert d.have_angle_model
  assert approx_equal(d.angle_model, 135)
  assert approx_equal(d.delta, 5)
  proxies = geometry_restraints.shared_dihedral_proxy([p,p])
  for proxy in proxies:
    assert proxy.periodicity == -2
  proxies[1].periodicity = 3
  assert proxies[1].periodicity == 3
  assert proxies[0].periodicity == -2
  proxies[1].periodicity = -2
  assert proxies[1].periodicity == -2
  assert approx_equal(geometry_restraints.dihedral_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [5]*2)
  assert approx_equal(geometry_restraints.dihedral_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [25]*2)
  residual_sum = geometry_restraints.dihedral_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*25)
  #
  sites_cart = flex.vec3_double((
    (44.14, -3.376, 8.756),
    (43.598, -2.045, 8.726),
    (42.178, -2.036, 9.302),
    (41.818, -0.984, 10.006)))
  r_orig = geometry_restraints.dihedral(
    sites=list(sites_cart), angle_ideal=0, weight=1)
  perm = flex.size_t(xrange(4))
  n_perms = 0
  n_equiv = 0
  n_equiv_direct = 0
  while 1:
    sites_perm = sites_cart.select(perm)
    r = geometry_restraints.dihedral(
      sites=list(sites_perm), angle_ideal=0, weight=1)
    if (   abs(r.angle_model - r_orig.angle_model) < 1.e-6
        or abs(r.angle_model + r_orig.angle_model) < 1.e-6):
      n_equiv += 1
      p = geometry_restraints.dihedral_proxy(
        i_seqs=list(perm),
        angle_ideal=r.angle_model,
        weight=1)
      rp = geometry_restraints.dihedral(
        sites_cart=sites_cart,
        proxy=p)
      assert approx_equal(rp.angle_model, r.angle_model)
      if (abs(p.angle_ideal - r_orig.angle_model) < 1.e-6):
        n_equiv_direct += 1
      p_sorted = p.sort_i_seqs()
      assert p_sorted.i_seqs == (0,1,2,3)
      assert approx_equal(p_sorted.angle_ideal, r_orig.angle_model)
    #
    p = geometry_restraints.dihedral_proxy(
      i_seqs=list(perm),
      angle_ideal=12,
      weight=1,
      alt_angle_ideals=[30,-40])
    p_sorted = p.sort_i_seqs()
    if (p_sorted.angle_ideal > 0):
      assert approx_equal(p_sorted.alt_angle_ideals, [30,-40])
    else:
      assert approx_equal(p_sorted.alt_angle_ideals, [-30,40])
    #
    n_perms += 1
    if (not perm.next_permutation()):
      break
  assert n_perms == 24
  assert n_equiv == 4
  assert n_equiv_direct == 2
  #
  proxies = geometry_restraints.shared_dihedral_proxy([
    geometry_restraints.dihedral_proxy([0,1,2,3], 1, 2, 3),
    geometry_restraints.dihedral_proxy([1,2,3,4], 2, 3, 4),
    geometry_restraints.dihedral_proxy([2,3,0,4], 3, 4, 5),
    geometry_restraints.dihedral_proxy([3,1,2,4], 4, 5, 6)])
  selected = proxies.proxy_select(n_seq=5, iselection=flex.size_t([0,2,4]))
  assert selected.size() == 0
  selected = proxies.proxy_select(n_seq=5, iselection=flex.size_t([1,2,3,4]))
  assert selected.size() == 2
  assert selected[0].i_seqs == (0,1,2,3)
  assert selected[1].i_seqs == (2,0,1,3)
  assert approx_equal(selected[0].angle_ideal, 2)
  assert approx_equal(selected[1].angle_ideal, 4)
  assert approx_equal(selected[0].weight, 3)
  assert approx_equal(selected[1].weight, 5)
  assert selected[0].periodicity == 4
  assert selected[1].periodicity == 6
  #
  rest = proxies.proxy_remove(selection=flex.bool([True,True,True,True,True]))
  assert rest.size() == 0
  rest = proxies.proxy_remove(selection=flex.bool([False,True,True,True,True]))
  assert rest.size() == 2
  assert rest[0].i_seqs == (0,1,2,3)
  assert rest[1].i_seqs == (2,3,0,4)
  rest = proxies.proxy_remove(selection=flex.bool([True,True,True,True,False]))
  assert rest.size() == 3
  #
  def get_d(angle_ideal, angle_model, periodicity, alt_angle_ideals=None):
    a = angle_model * math.pi / 180
    c, s = math.cos(a), math.sin(a)
    d = geometry_restraints.dihedral(
      sites=[(1,0,-1),(0,0,-1),(0,0,0),(c,s,0)],
      angle_ideal=angle_ideal,
      weight=1/15**2,
      periodicity=periodicity,
      alt_angle_ideals=alt_angle_ideals)
    v = math.fmod(d.angle_model-angle_model, 360)
    if (v < 0): v += 360
    if (v > 360-1e-8): v -= 360
    assert approx_equal(v, 0)
    return d
  #
  for periodicity in xrange(1,6):
    f = open("plot_geo_restr_dihedral_periodicity_%d.xy" % periodicity, "w")
    for signed_periodicity in [periodicity, -periodicity]:
      for angle_model in xrange(0, 720+1, 1):
        d = get_d(
          angle_ideal=70,
          angle_model=angle_model,
          periodicity=signed_periodicity)
        print >> f, angle_model, d.residual()
      print >> f, "&"
    f.close()
  #
  intersection_angle = 120
  for angle_ideal in xrange(0, 720+5, 5):
    for periodicity in xrange(1,6):
      for signed_periodicity in [periodicity, -periodicity]:
        residuals = []
        for offset in [0, intersection_angle, -intersection_angle]:
          d = get_d(
            angle_ideal=angle_ideal,
            angle_model=angle_ideal + offset / periodicity,
            periodicity=signed_periodicity)
          residuals.append(d.residual())
        assert approx_equal(residuals[0], 0)
        assert approx_equal(residuals[1], residuals[2])
        assert approx_equal(residuals[1], d.weight * d.delta**2)
      #
      for offset in [intersection_angle, -intersection_angle]:
        for offset2 in [30, -30]:
          residuals = []
          for signed_periodicity in [periodicity, -periodicity]:
            d = get_d(
              angle_ideal=angle_ideal,
              angle_model=angle_ideal + (offset + offset2) / periodicity,
              periodicity=signed_periodicity)
            residuals.append(d.residual())
          if ((offset > 0) == (offset2 < 0)):
            assert residuals[0] > residuals[1]
          else:
            assert residuals[0] < residuals[1]

  #test alt_angle_ideals
  d = get_d(angle_ideal=180, angle_model=170, periodicity=2)
  assert approx_equal(d.delta, 10.)
  assert d.alt_angle_ideals is None
  d = get_d(angle_ideal=180, angle_model=170, periodicity=2,
    alt_angle_ideals=(25,))
  assert approx_equal(d.delta, 10.)
  d = get_d(angle_ideal=180, angle_model=10, periodicity=2,
    alt_angle_ideals=(25,))
  assert approx_equal(d.delta, -10.)
  d = get_d(angle_ideal=180, angle_model=10, periodicity=2,
    alt_angle_ideals=(15,))
  assert approx_equal(d.delta, 5.)
  d = get_d(angle_ideal=180, angle_model=10, periodicity=2,
    alt_angle_ideals=(15,345))
  assert approx_equal(d.delta, 5.)
  d = get_d(angle_ideal=180, angle_model=-10, periodicity=2,
    alt_angle_ideals=(15,345))
  assert approx_equal(d.delta, -5.)
  d = get_d(angle_ideal=30, angle_model=28, periodicity=1)
  assert approx_equal(d.delta, 2.)
  d = get_d(angle_ideal=30, angle_model=-30, periodicity=1)
  assert approx_equal(d.delta, 60.)
  d = get_d(angle_ideal=30, angle_model=28, periodicity=1,
    alt_angle_ideals=(-30,))
  assert approx_equal(d.delta, 2.)
  d = get_d(angle_ideal=30, angle_model=-28, periodicity=1,
    alt_angle_ideals=(-30,))
  assert approx_equal(d.delta, -2.)

def exercise_chirality():
  p = geometry_restraints.chirality_proxy(
    i_seqs=[0,2,3,1],
    volume_ideal=4,
    both_signs=False,
    weight=1)
  assert p.i_seqs == (0,2,3,1)
  assert approx_equal(p.volume_ideal, 4)
  assert not p.both_signs
  assert approx_equal(p.weight, 1)
  c = geometry_restraints.chirality_proxy(
    i_seqs=[9,0,4,6],
    proxy=p)
  assert c.i_seqs == (9,0,4,6)
  assert approx_equal(c.volume_ideal, 4)
  assert not c.both_signs
  assert approx_equal(c.weight, 1)
  c = p.scale_weight(factor=9.32)
  assert c.i_seqs == (0,2,3,1)
  assert approx_equal(c.weight, 9.32)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2,3)
  assert approx_equal(p.volume_ideal, 4)
  p = geometry_restraints.chirality_proxy(
    i_seqs=[0,2,1,3],
    volume_ideal=-4,
    both_signs=False,
    weight=1)
  assert approx_equal(p.volume_ideal, -4)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2,3)
  assert approx_equal(p.volume_ideal, 4)
  c = geometry_restraints.chirality(
    sites=[(1,0,0),(0,0,0),(0,1,0),(1,0,1)],
    volume_ideal=4,
    both_signs=False,
    weight=1)
  assert approx_equal(c.sites, [(1,0,0),(0,0,0),(0,1,0),(1,0,1)])
  assert approx_equal(c.volume_ideal, 4)
  assert approx_equal(c.weight, 1)
  assert not c.both_signs
  assert approx_equal(c.volume_model, -1)
  assert approx_equal(c.delta_sign, -1)
  assert approx_equal(c.delta, 5)
  assert approx_equal(c.residual(), 25)
  assert approx_equal(c.gradients(),
    ((10, 0, -10),
     (-10, -10, 0),
     (-0, 10, -0),
     (-0, -0, 10)))
  sites_cart = flex.vec3_double([(1,0,0),(0,0,0),(0,1,0),(-1,0,-1)])
  c = geometry_restraints.chirality(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(c.sites, [(1,0,0),(0,0,0),(0,1,0),(-1,0,-1)])
  assert approx_equal(c.volume_ideal, 4)
  assert approx_equal(c.weight, 1)
  assert not c.both_signs
  assert approx_equal(c.volume_model, 1)
  assert approx_equal(c.delta_sign, -1)
  assert approx_equal(c.delta, 3)
  proxies = geometry_restraints.shared_chirality_proxy([p,p])
  for proxy in proxies:
    assert proxy.volume_ideal == 4
  assert approx_equal(geometry_restraints.chirality_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [3]*2)
  assert approx_equal(geometry_restraints.chirality_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [9]*2)
  residual_sum = geometry_restraints.chirality_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*9)
  #
  proxies = geometry_restraints.shared_chirality_proxy([
    geometry_restraints.chirality_proxy([0,1,2,3], 1, False, 2),
    geometry_restraints.chirality_proxy([1,2,3,4], 2, True, 3),
    geometry_restraints.chirality_proxy([2,3,0,4], 3, True, 4),
    geometry_restraints.chirality_proxy([3,1,2,4], 4, False, 5)])
  selected = proxies.proxy_select(n_seq=5, iselection=flex.size_t([0,2,4]))
  assert selected.size() == 0
  selected = proxies.proxy_select(n_seq=5, iselection=flex.size_t([1,2,3,4]))
  assert selected.size() == 2
  assert selected[0].i_seqs == (0,1,2,3)
  assert selected[1].i_seqs == (2,0,1,3)
  assert approx_equal(selected[0].volume_ideal, 2)
  assert approx_equal(selected[1].volume_ideal, 4)
  assert selected[0].both_signs
  assert not selected[1].both_signs
  assert approx_equal(selected[0].weight, 3)
  assert approx_equal(selected[1].weight, 5)
  #
  rest = proxies.proxy_remove(selection=flex.bool([True,True,True,True,True]))
  assert rest.size() == 0
  rest = proxies.proxy_remove(selection=flex.bool([False,True,True,True,True]))
  assert rest.size() == 2
  assert rest[0].i_seqs == (0,1,2,3)
  assert rest[1].i_seqs == (2,3,0,4)
  rest = proxies.proxy_remove(selection=flex.bool([True,True,True,True,False]))
  assert rest.size() == 3

def exercise_planarity():
  weights = flex.double([1, 2, 3, 4])
  u_mx = sgtbx.rt_mx()
  sym_ops = (u_mx, sgtbx.rt_mx('1+x,y,z'), u_mx, sgtbx.rt_mx('1+x,1+y,z'))
  p = geometry_restraints.planarity_proxy(
    i_seqs=flex.size_t([3,1,0,2]),
    sym_ops=sym_ops,
    weights=weights)
  assert tuple(p.i_seqs) == (3,1,0,2)
  assert tuple(p.sym_ops) == sym_ops
  assert approx_equal(p.weights, weights)
  c = geometry_restraints.planarity_proxy(
    i_seqs=flex.size_t([8,6,3,1]),
    proxy=p)
  assert tuple(c.i_seqs) == (8,6,3,1)
  assert c.sym_ops == sym_ops
  c = p.scale_weights(factor=4.94)
  assert c.i_seqs.id() == p.i_seqs.id()
  assert c.weights.id() != p.weights.id()
  assert approx_equal(c.weights, p.weights*4.94)
  p = p.sort_i_seqs()
  assert tuple(p.i_seqs) == (0,1,2,3)
  assert tuple(p.weights) == (3,2,4,1)
  assert tuple(p.sym_ops) == (u_mx, sgtbx.rt_mx('1+x,y,z'), sgtbx.rt_mx('1+x,1+y,z'), u_mx)
  #
  unit_cell = uctbx.unit_cell([15,25,30,90,90,90])
  sites_cart = flex.vec3_double([(1,24,1.1),(1,1,1),(14,1,1),(14,24,0.9)])
  expected_residual = 0.04
  expected_gradients = [(0,0,0.4), (0,0,0), (0,0,0), (0,0,-0.4)]
  p = geometry_restraints.planarity_proxy(
    i_seqs=(0,1,2,3),
    sym_ops=[u_mx, sgtbx.rt_mx('x,y,z'), sgtbx.rt_mx('x,1+y,z'), sgtbx.rt_mx('1-x,y,z')],
    weights=(2,2,2,2))
  planarity = geometry_restraints.planarity(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(planarity.deltas(), (0.1, 0, 0, -0.1))
  assert approx_equal(planarity.residual(), 0.04)
  assert approx_equal(planarity.gradients(), expected_gradients)
  #
  proxies = geometry_restraints.shared_planarity_proxy([p,p])
  for proxy in proxies:
    assert tuple(proxy.i_seqs) == (0,1,2,3)
  assert eps_eq(geometry_restraints.planarity_deltas_rms(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [0.070710678118654821]*2)
  assert eps_eq(geometry_restraints.planarity_residuals(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [expected_residual]*2)
  residual_sum = geometry_restraints.planarity_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert eps_eq(residual_sum, 2*expected_residual)
  gradient_array = flex.vec3_double(proxy.i_seqs.size(), (0,0,0))
  residual_sum = geometry_restraints.planarity_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=gradient_array)
  assert eps_eq(residual_sum, 2*0.04)
  for g,e in zip(gradient_array, expected_gradients):
    assert eps_eq(g, matrix.col(e)*2)
  #
  sites_cart = flex.vec3_double([
    (-6.9010753374697966, 1.3017288659588333, -1.4469233441387523),
    (-4.947324488687852, -1.0193474269570115, 0.16296067326855093),
    (-6.9598378855214706, -0.66835111494675281, -1.7153810358296142),
    (-4.846552160625774, 0.96315156534510438, 0.51500258491293438)])
  weights = flex.double([1, 2, 3, 4])
  expected_residual = 0.000428526964094
  expected_gradients = [
    (0.019669677238598002, 0.0024733761183690019, -0.020428665017957027),
    (0.020015197633649708, 0.0025168238009701843, -0.020787517902107686),
    (-0.020270795833584803, -0.002548964159754365, 0.0210529787910856),
    (-0.019414079038662237, -0.0024412357595847371, 0.020163204128978415)]
  p = geometry_restraints.planarity_proxy(
    i_seqs=flex.size_t([0,1,2,3]),
    weights=weights)
  assert tuple(p.i_seqs) == (0,1,2,3)
  assert approx_equal(p.weights, weights)
  perm = flex.size_t([3,1,0,2])
  p = geometry_restraints.planarity_proxy(
    i_seqs=flex.size_t([0,1,2,3]).select(perm),
    weights=weights.select(perm))
  assert tuple(p.i_seqs) == (3,1,0,2)
  assert not_approx_equal(p.weights, weights)
  p = p.sort_i_seqs()
  assert tuple(p.i_seqs) == (0,1,2,3)
  assert approx_equal(p.weights, weights)
  for i_constructor in xrange(2):
    if (i_constructor == 0):
      l = geometry_restraints.planarity(
        sites=sites_cart,
        weights=weights)
    else:
      l = geometry_restraints.planarity(
        sites_cart=sites_cart,
        proxy=p)
    assert approx_equal(l.sites, sites_cart)
    assert approx_equal(l.weights, weights)
    assert eps_eq(l.deltas(), (0.014233272168667327, 0.007241647943016986,
                               -0.0048894168534149443, -0.0035120793736139956))
    assert eps_eq(l.rms_deltas(), 0.00853329655308)
    assert eps_eq(l.residual(), expected_residual)
    assert eps_eq(l.gradients(), expected_gradients)
    assert eps_eq(l.normal(),
      (0.69097523765119184, 0.086887122267422581, -0.71763768639680903))
    assert eps_eq(l.residual(), l.lambda_min())
    assert eps_eq(l.center_of_mass(),
      (-5.7061446613913009, 0.11105869285849694, -0.42071347654387559))
    assert eps_eq(l.center_of_mass(), l.sites.mean_weighted(weights=l.weights))
    assert eps_eq(l.residual_tensor(),
      (10.250312599815825, 8.7000194514224525, 10.265208176541265,
       2.7229147081229312, 10.19874296603952, 3.6750425846794936))
    assert eps_eq(l.eigensystem().values(),
      [21.998140770294835, 7.2169709305206142, 0.00042852696409348911])
  proxies = geometry_restraints.shared_planarity_proxy([p,p])
  for proxy in proxies:
    assert tuple(proxy.i_seqs) == (0,1,2,3)
  assert eps_eq(geometry_restraints.planarity_deltas_rms(
    sites_cart=sites_cart,
    proxies=proxies), [0.0085332965530764398]*2)
  assert eps_eq(geometry_restraints.planarity_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [expected_residual]*2)
  residual_sum = geometry_restraints.planarity_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert eps_eq(residual_sum, 2*expected_residual)
  gradient_array = flex.vec3_double(proxy.i_seqs.size(), (0,0,0))
  residual_sum = geometry_restraints.planarity_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=gradient_array)
  assert eps_eq(residual_sum, 2*expected_residual)
  for g,e in zip(gradient_array, expected_gradients):
    assert eps_eq(g, matrix.col(e)*2)
  #
  def make_proxy(i_seqs, weights):
    return geometry_restraints.planarity_proxy(
      flex.size_t(i_seqs),
      flex.double(weights))
  proxies = geometry_restraints.shared_planarity_proxy([
    make_proxy([0,1,2,3], [2,3,4,5]),
    make_proxy([1,2,3,4], [3,4,5,6]),
    make_proxy([2,3,0,4], [4,5,6,7]),
    make_proxy([3,1,2,4], [5,6,7,8])])
  selected = proxies.proxy_select(n_seq=5, iselection=flex.size_t([0,2,4]))
  assert selected.size() == 0
  selected = proxies.proxy_select(n_seq=5, iselection=flex.size_t([1,2,3,4]))
  assert selected.size() == 2
  assert list(selected[0].i_seqs) == [0,1,2,3]
  assert list(selected[1].i_seqs) == [2,0,1,3]
  assert approx_equal(selected[0].weights, [3,4,5,6])
  assert approx_equal(selected[1].weights, [5,6,7,8])
  #
  for i_remove in range(10,15):
    sel = flex.size_t(range(10,i_remove)+range(i_remove+1,15))
    pp = geometry_restraints.planarity_proxy(
      i_seqs=flex.size_t([10, 11, 12, 13, 14]),
      weights=flex.double(range(5))+13)
    pps = geometry_restraints.shared_planarity_proxy()
    pps.append(pp)
    selected = pps.proxy_select(20, sel)
    assert list(selected[0].i_seqs) == [0,1,2,3]
    assert approx_equal(selected[0].weights,
      pp.weights[:i_remove-10].concatenate(pp.weights[i_remove+1-10:]))
  #
  rest = proxies.proxy_remove(selection=flex.bool([True,True,True,True,True]))
  assert rest.size() == 0
  rest = proxies.proxy_remove(selection=flex.bool([False,True,True,True,True]))
  assert rest.size() == 2
  assert list(rest[0].i_seqs) == [0,1,2,3]
  assert list(rest[1].i_seqs) == [2,3,0,4]
  rest = proxies.proxy_remove(selection=flex.bool([True,True,True,True,False]))
  assert rest.size() == 3
  #
  unit_cell = uctbx.unit_cell([15,11.5,16.25,90,99.5,90])
  sites_cart = flex.vec3_double(
    [(12.87,0.10,9.04),(12.54,0.44,7.73),(13.47,0.34,6.71)])
  rt_mx = sgtbx.rt_mx('2-X,-Y,1-Z')
  u_mx = sgtbx.rt_mx()
  i_seqs = flex.size_t([0,1,2,0,1,2])
  sym_ops = (u_mx,u_mx,u_mx,rt_mx,rt_mx,rt_mx)
  weights = flex.double([1]*6)
  p = geometry_restraints.planarity_proxy(
    i_seqs=i_seqs,
    sym_ops=sym_ops,
    weights=weights)
  planarity = geometry_restraints.planarity(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p)
  gradient_array = flex.vec3_double(sites_cart.size())
  residual_sum = geometry_restraints.planarity_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=geometry_restraints.shared_planarity_proxy([p]),
    gradient_array=gradient_array)
  fd_grads = finite_difference_gradients(
    restraint_type=geometry_restraints.planarity,
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p)
  for g,e in zip(gradient_array, fd_grads):
    assert approx_equal(g, e)
  # check proxies with and without sym_ops are happy side-by-side
  p_sym = geometry_restraints.planarity_proxy(
    i_seqs=i_seqs,
    sym_ops=sym_ops,
    weights=weights)
  assert p_sym.sym_ops == sym_ops
  restraint_sym = geometry_restraints.planarity(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxy=p_sym)
  p_no_sym = geometry_restraints.planarity_proxy(
    i_seqs=i_seqs,
    weights=weights)
  assert p_no_sym.sym_ops == None
  restraint_no_sym = geometry_restraints.planarity(
    sites_cart=sites_cart,
    proxy=p_no_sym)
  proxies = geometry_restraints.shared_planarity_proxy([p_sym,p_no_sym])
  assert approx_equal(geometry_restraints.planarity_deltas_rms(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [restraint_sym.rms_deltas(),restraint_no_sym.rms_deltas()])
  assert approx_equal(geometry_restraints.planarity_residuals(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies), [restraint_sym.residual(),restraint_no_sym.residual()])
  residual_sum = geometry_restraints.planarity_residual_sum(
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, restraint_sym.residual() + restraint_no_sym.residual())

def exercise_proxy_show():
  if sys.platform.startswith("win") and sys.version_info[:2] < (2,6):
    # This appears to be a windows-specific bug with string formatting
    # for python versions prior to 2.6, where the exponent is printed
    # with 3 digits rather than 2.
    print "Skipping exercise_proxy_show()"
    return
  # zeolite AHT
  crystal_symmetry = crystal.symmetry(
    unit_cell=(15.794, 9.206, 8.589, 90, 90, 90),
    space_group_symbol="C m c m")
  sites_cart_cry = crystal_symmetry.unit_cell().orthogonalization_matrix() \
    * flex.vec3_double([(0.1681, 0.6646, 0.4372), (0.0000, 0.6644, 0.5629)])
  asu_mappings = crystal_symmetry.asu_mappings(buffer_thickness=3.0)
  asu_mappings.process_sites_cart(original_sites=sites_cart_cry)
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  pair_asu_table.add_all_pairs(distance_cutoff=2.9)
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  sio = StringIO()
  sorted_asu_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart_cry,
    f=sio)
  assert not show_diff(sio.getvalue(), """\
Bond restraints: 0
""")
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    pair_asu_table=pair_asu_table)
  mt = flex.mersenne_twister(seed=5)
  for proxy in sorted_asu_proxies.asu:
    proxy.distance_ideal = 2.9 + (mt.random_double()-0.5)*0.2
    proxy.weight = 1+mt.random_double()*100
    if (mt.random_double() > 0.5):
      proxy.slack = mt.random_double()*0.1
  sio = StringIO()
  sorted_asu_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart_cry,
    site_labels=["Si1", "Si2"],
    f=sio)
  assert not show_diff(sio.getvalue(), """\
Bond restraints: 3
Sorted by residual:
bond Si1
     Si2
  ideal  model  slack  delta    sigma   weight residual sym.op.
  2.979  2.866  0.004  0.112 1.71e-01 3.42e+01 4.01e-01 x,y,z
bond Si2
     Si1
  ideal  model  slack  delta    sigma   weight residual sym.op.
  2.822  2.866  0.042 -0.045 1.29e-01 6.05e+01 4.35e-04 x,y,z
bond Si2
     Si1
  ideal  model  delta    sigma   weight residual sym.op.
  2.867  2.866  0.001 1.26e-01 6.33e+01 6.17e-05 -x,y,z
""")
  sio = StringIO()
  sorted_asu_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart_cry,
    f=sio,
    prefix="&",
    max_items=2)
  assert not show_diff(sio.getvalue(), """\
&Bond restraints: 3
&Sorted by residual:
&bond 0
&     1
&  ideal  model  slack  delta    sigma   weight residual sym.op.
&  2.979  2.866  0.004  0.112 1.71e-01 3.42e+01 4.01e-01 x,y,z
&bond 1
&     0
&  ideal  model  slack  delta    sigma   weight residual sym.op.
&  2.822  2.866  0.042 -0.045 1.29e-01 6.05e+01 4.35e-04 x,y,z
&... (remaining 1 not shown)
""")
  sio = StringIO()
  sorted_asu_proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart_cry,
    f=sio,
    prefix="*",
    max_items=0)
  assert not show_diff(sio.getvalue(), """\
*Bond restraints: 3
""")
  #
  for unit_cell in [None, uctbx.unit_cell([15,11.5,16.25,90,99.5,90])]:
    simple_proxies = geometry_restraints.shared_bond_simple_proxy([
      geometry_restraints.bond_simple_proxy((0,1), 2.979, 10),
      geometry_restraints.bond_simple_proxy(
        (1,0), sgtbx.rt_mx("x,y,z"), 2.822, 250)
    ])
    sio = StringIO()
    simple_proxies.show_sorted(
      "residual", sites_cart_cry, unit_cell=unit_cell, f=sio)
    assert not show_diff(sio.getvalue(), """\
Bond restraints: 2
Sorted by residual:
bond 1
     0
  ideal  model  delta    sigma   weight residual sym.op.
  2.822  2.866 -0.044 6.32e-02 2.50e+02 4.86e-01 x,y,z
bond 0
     1
  ideal  model  delta    sigma   weight residual
  2.979  2.866  0.113 3.16e-01 1.00e+01 1.27e-01
""")
  mt = flex.mersenne_twister(seed=73)
  sites_cart = flex.vec3_double(mt.random_double(size=18))
  site_labels = ["a", "ba", "c", "dada", "e", "f"]
  #
  proxies = geometry_restraints.shared_angle_proxy()
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=flex.vec3_double(),
    f=sio)
  assert not show_diff(sio.getvalue(), """\
Bond angle restraints: 0
""")
  proxies = geometry_restraints.shared_angle_proxy([
    geometry_restraints.angle_proxy(
      i_seqs=[2,1,0],
      angle_ideal=59,
      weight=2),
    geometry_restraints.angle_proxy(
      i_seqs=[3,0,1],
      angle_ideal=99,
      weight=8)])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=sio,
    prefix="+")
  assert not show_diff(sio.getvalue(), """\
+Bond angle restraints: 2
+Sorted by residual:
+angle 3
+      0
+      1
+    ideal   model   delta    sigma   weight residual
+    99.00   99.72   -0.72 3.54e-01 8.00e+00 4.19e+00
+angle 2
+      1
+      0
+    ideal   model   delta    sigma   weight residual
+    59.00   58.06    0.94 7.07e-01 2.00e+00 1.76e+00
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="delta",
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=sio,
    prefix="@",
    max_items=1)
  assert not show_diff(sio.getvalue(), """\
@Bond angle restraints: 2
@Sorted by delta:
@angle c
@      ba
@      a
@    ideal   model   delta    sigma   weight residual
@    59.00   58.06    0.94 7.07e-01 2.00e+00 1.76e+00
@... (remaining 1 not shown)
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=sio,
    max_items=0)
  assert not show_diff(sio.getvalue(), """\
Bond angle restraints: 2
""")
  #
  proxies = geometry_restraints.shared_dihedral_proxy()
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=flex.vec3_double(),
    f=sio)
  assert not show_diff(sio.getvalue(), """\
Dihedral angle restraints: 0
""")
  proxies = geometry_restraints.shared_dihedral_proxy([
    geometry_restraints.dihedral_proxy(
      i_seqs=[0,1,3,4],
      angle_ideal=59,
      weight=2,
      periodicity=-1),
    geometry_restraints.dihedral_proxy(
      i_seqs=[3,2,0,5],
      angle_ideal=99,
      weight=8,
      periodicity=-1)])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=sio,
    prefix="-")
  assert not show_diff(sio.getvalue(), """\
-Dihedral angle restraints: 2
-  sinusoidal: 0
-    harmonic: 2
-Sorted by residual:
-dihedral 3
-         2
-         0
-         5
-    ideal   model   delta  harmonic     sigma   weight residual
-    99.00   16.67   82.33    -1      3.54e-01 8.00e+00 5.42e+04
-dihedral 0
-         1
-         3
-         4
-    ideal   model   delta  harmonic     sigma   weight residual
-    59.00 -159.79 -141.21    -1      7.07e-01 2.00e+00 3.99e+04
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="delta",
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=sio,
    prefix="^",
    max_items=1)
  assert not show_diff(sio.getvalue(), """\
^Dihedral angle restraints: 2
^  sinusoidal: 0
^    harmonic: 2
^Sorted by delta:
^dihedral a
^         ba
^         dada
^         e
^    ideal   model   delta  harmonic     sigma   weight residual
^    59.00 -159.79 -141.21    -1      7.07e-01 2.00e+00 3.99e+04
^... (remaining 1 not shown)
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=sio,
    max_items=0)
  assert not show_diff(sio.getvalue(), """\
Dihedral angle restraints: 2
  sinusoidal: 0
    harmonic: 2
""")
  #
  proxies = geometry_restraints.shared_chirality_proxy()
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=flex.vec3_double(),
    f=sio)
  assert not show_diff(sio.getvalue(), """\
Chirality restraints: 0
""")
  proxies = geometry_restraints.shared_chirality_proxy([
    geometry_restraints.chirality_proxy(
      i_seqs=[0,1,3,4],
      volume_ideal=0.09,
      both_signs=False,
      weight=2),
    geometry_restraints.chirality_proxy(
      i_seqs=[3,2,0,5],
      volume_ideal=0.16,
      both_signs=True,
      weight=8)])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=sio,
    prefix="$")
  assert not show_diff(sio.getvalue(), """\
$Chirality restraints: 2
$Sorted by residual:
$chirality 3
$          2
$          0
$          5
$  both_signs  ideal   model   delta    sigma   weight residual
$    True       0.16    0.05    0.11 3.54e-01 8.00e+00 9.34e-02
$chirality 0
$          1
$          3
$          4
$  both_signs  ideal   model   delta    sigma   weight residual
$    False      0.09   -0.04    0.13 7.07e-01 2.00e+00 3.33e-02
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="delta",
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=sio,
    prefix="*",
    max_items=1)
  assert not show_diff(sio.getvalue(), """\
*Chirality restraints: 2
*Sorted by delta:
*chirality a
*          ba
*          dada
*          e
*  both_signs  ideal   model   delta    sigma   weight residual
*    False      0.09   -0.04    0.13 7.07e-01 2.00e+00 3.33e-02
*... (remaining 1 not shown)
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=sio,
    max_items=0)
  assert not show_diff(sio.getvalue(), """\
Chirality restraints: 2
""")
  #
  proxies = geometry_restraints.shared_planarity_proxy()
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=flex.vec3_double(),
    f=sio)
  assert not show_diff(sio.getvalue(), """\
Planarity restraints: 0
""")
  proxies = geometry_restraints.shared_planarity_proxy([
    geometry_restraints.planarity_proxy(
      i_seqs=[0,2,4,1],
      weights=[0.31,0.2,0.31,0.4]),
    geometry_restraints.planarity_proxy(
      i_seqs=[0,2,3,4,5],
      weights=[0.01,0.11,0.21,0.31,0.41])])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=sio,
    prefix=":")
  assert not show_diff(sio.getvalue(), """\
:Planarity restraints: 2
:Sorted by residual:
:           delta    sigma   weight rms_deltas residual
:plane 0    0.004 1.80e+00 3.10e-01   1.46e-01 2.52e-02
:      2   -0.196 2.24e+00 2.00e-01
:      4   -0.115 1.80e+00 3.10e-01
:      1    0.184 1.58e+00 4.00e-01
:           delta    sigma   weight rms_deltas residual
:plane 0   -0.332 1.00e+01 1.00e-02   1.78e-01 9.86e-03
:      2    0.152 3.02e+00 1.10e-01
:      3   -0.143 2.18e+00 2.10e-01
:      4   -0.030 1.80e+00 3.10e-01
:      5    0.063 1.56e+00 4.10e-01
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="rms_deltas",
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=sio,
    prefix="<",
    max_items=1)
  assert not show_diff(sio.getvalue(), """\
<Planarity restraints: 2
<Sorted by rms_deltas:
<              delta    sigma   weight rms_deltas residual
<plane a      -0.332 1.00e+01 1.00e-02   1.78e-01 9.86e-03
<      c       0.152 3.02e+00 1.10e-01
<      dada   -0.143 2.18e+00 2.10e-01
<      e      -0.030 1.80e+00 3.10e-01
<      f       0.063 1.56e+00 4.10e-01
<... (remaining 1 not shown)
""")
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    sites_cart=sites_cart,
    f=sio,
    max_items=0)
  assert not show_diff(sio.getvalue(), """\
Planarity restraints: 2
""")
  #
  unit_cell = uctbx.unit_cell([15,11.5,16.25,90,99.5,90])
  sites_cart = flex.vec3_double(
    [(12.87,0.10,9.04),(12.54,0.44,7.73),(13.47,0.34,6.71),(1,2,3)])
  rt_mx = sgtbx.rt_mx('2-X,-Y,1-Z')
  u_mx = sgtbx.rt_mx()
  p = geometry_restraints.angle_proxy(
    i_seqs=[2,0,1],
    sym_ops=[u_mx,rt_mx,rt_mx],
    angle_ideal=120,
    weight=1)
  proxies = geometry_restraints.shared_angle_proxy([p,p])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    f=sio,
    prefix='!')
  assert not show_diff(sio.getvalue(), """\
!Bond angle restraints: 2
!Sorted by residual:
!angle 2
!      0  -x+2,-y,-z+1
!      1  -x+2,-y,-z+1
!    ideal   model   delta    sigma   weight residual
!   120.00  122.78   -2.78 1.00e+00 1.00e+00 7.73e+00
!angle 2
!      0  -x+2,-y,-z+1
!      1  -x+2,-y,-z+1
!    ideal   model   delta    sigma   weight residual
!   120.00  122.78   -2.78 1.00e+00 1.00e+00 7.73e+00
""")
  # test proxies with and without proxy.sym_ops side by side
  p2 = geometry_restraints.angle_proxy(
    i_seqs=[2,1,0],
    angle_ideal=59,
    weight=2)
  proxies = geometry_restraints.shared_angle_proxy([p,p2])
  sio = StringIO()
  proxies.show_sorted(
    by_value="delta",
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    f=sio,
    prefix='~')
  assert not show_diff(sio.getvalue(), """\
~Bond angle restraints: 2
~Sorted by delta:
~angle 2
~      1
~      0
~    ideal   model   delta    sigma   weight residual
~    59.00  121.08  -62.08 7.07e-01 2.00e+00 7.71e+03
~angle 2
~      0  -x+2,-y,-z+1
~      1  -x+2,-y,-z+1
~    ideal   model   delta    sigma   weight residual
~   120.00  122.78   -2.78 1.00e+00 1.00e+00 7.73e+00
""")
  #
  p = geometry_restraints.dihedral_proxy(
    i_seqs=[2,0,1,2],
    sym_ops=[u_mx,u_mx,u_mx,rt_mx],
    angle_ideal=0,
    weight=1)
  proxies = geometry_restraints.shared_dihedral_proxy([p,p])
  sio = StringIO()
  proxies.show_sorted(
    by_value="delta",
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    f=sio,
    prefix='%')
  assert not show_diff(sio.getvalue(), """\
%Dihedral angle restraints: 2
%  sinusoidal: 0
%    harmonic: 2
%Sorted by delta:
%dihedral 2
%         0
%         1
%         2  -x+2,-y,-z+1
%    ideal   model   delta  harmonic     sigma   weight residual
%     0.00    5.16   -5.16     0      1.00e+00 1.00e+00 2.67e+01
%dihedral 2
%         0
%         1
%         2  -x+2,-y,-z+1
%    ideal   model   delta  harmonic     sigma   weight residual
%     0.00    5.16   -5.16     0      1.00e+00 1.00e+00 2.67e+01
""")
  #
  p = geometry_restraints.planarity_proxy(
    i_seqs=flex.size_t([0,1,2,0,1,2]),
    sym_ops=[u_mx,u_mx,u_mx,rt_mx,rt_mx,rt_mx],
    weights=flex.double([1]*6))
  proxies = geometry_restraints.shared_planarity_proxy([p,p])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    f=sio,
    prefix=">")
  assert not show_diff(sio.getvalue(), """\
>Planarity restraints: 2
>Sorted by residual:
>           delta    sigma   weight rms_deltas residual sym.op.
>plane 0    0.017 1.00e+00 1.00e+00   1.60e-02 1.53e-03
>      1   -0.015 1.00e+00 1.00e+00
>      2    0.016 1.00e+00 1.00e+00
>      0   -0.017 1.00e+00 1.00e+00                     -x+2,-y,-z+1
>      1    0.015 1.00e+00 1.00e+00                     -x+2,-y,-z+1
>      2   -0.016 1.00e+00 1.00e+00                     -x+2,-y,-z+1
>           delta    sigma   weight rms_deltas residual sym.op.
>plane 0    0.017 1.00e+00 1.00e+00   1.60e-02 1.53e-03
>      1   -0.015 1.00e+00 1.00e+00
>      2    0.016 1.00e+00 1.00e+00
>      0   -0.017 1.00e+00 1.00e+00                     -x+2,-y,-z+1
>      1    0.015 1.00e+00 1.00e+00                     -x+2,-y,-z+1
>      2   -0.016 1.00e+00 1.00e+00                     -x+2,-y,-z+1
""")
  # test proxies with and without proxy.sym_ops side by side
  p2 = geometry_restraints.planarity_proxy(
    i_seqs=[0,2,3,1],
    weights=[0.31,0.2,0.31,0.4])
  proxies = geometry_restraints.shared_planarity_proxy([p,p2])
  sio = StringIO()
  proxies.show_sorted(
    by_value="rms_deltas",
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    f=sio,
    prefix="#")
  assert not show_diff(sio.getvalue(), """\
#Planarity restraints: 2
#Sorted by rms_deltas:
#           delta    sigma   weight rms_deltas residual sym.op.
#plane 0   -0.053 1.80e+00 3.10e-01   6.02e-02 4.53e-03
#      2   -0.071 2.24e+00 2.00e-01
#      3   -0.005 1.80e+00 3.10e-01
#      1    0.081 1.58e+00 4.00e-01
#           delta    sigma   weight rms_deltas residual sym.op.
#plane 0    0.017 1.00e+00 1.00e+00   1.60e-02 1.53e-03
#      1   -0.015 1.00e+00 1.00e+00
#      2    0.016 1.00e+00 1.00e+00
#      0   -0.017 1.00e+00 1.00e+00                     -x+2,-y,-z+1
#      1    0.015 1.00e+00 1.00e+00                     -x+2,-y,-z+1
#      2   -0.016 1.00e+00 1.00e+00                     -x+2,-y,-z+1
""")
  #
  unit_cell = uctbx.unit_cell([15,25,30,90,90,90])
  sites_cart = flex.vec3_double(
    [(1,2,3),(2,4,6),(1,2,4.5),(2,5.6,6),(14,24,29),(0.5,24,29)])
  p = geometry_restraints.bond_similarity_proxy(
    i_seqs=[(0,2),(1,3),(4,5)],
    sym_ops=[u_mx,u_mx,sgtbx.rt_mx('1+x,y,z')],
    weights=(1,2,3))
  proxies = geometry_restraints.shared_bond_similarity_proxy([p,p])
  sio = StringIO()
  proxies.show_sorted(
    by_value="residual",
    unit_cell=unit_cell,
    sites_cart=sites_cart,
    f=sio,
    prefix=">")
  assert not show_diff(sio.getvalue(), """\
>Bond similarity restraints: 2
>Sorted by residual:
>            delta    sigma   weight rms_deltas residual sym.op.
>bond 0-2   -0.033 1.00e+00 1.00e+00   4.71e-02 2.22e-03
>     1-3    0.067 7.07e-01 2.00e+00
>     4-5   -0.033 5.77e-01 3.00e+00                     x+1,y,z
>            delta    sigma   weight rms_deltas residual sym.op.
>bond 0-2   -0.033 1.00e+00 1.00e+00   4.71e-02 2.22e-03
>     1-3    0.067 7.07e-01 2.00e+00
>     4-5   -0.033 5.77e-01 3.00e+00                     x+1,y,z
""")
  sio = StringIO()
  p = geometry_restraints.bond_similarity_proxy(
    i_seqs=[(0,2),(1,3),(4,5)],
    weights=(1,2,3))
  proxies = geometry_restraints.shared_bond_similarity_proxy([p,p])
  proxies.show_sorted(
    by_value="rms_deltas",
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=sio,
    prefix=">")
  assert not show_diff(sio.getvalue(), """\
>Bond similarity restraints: 2
>Sorted by rms_deltas:
>                delta    sigma   weight rms_deltas residual
>bond a-c       -6.033 1.00e+00 1.00e+00   5.98e+00 3.56e+01
>     ba-dada   -5.933 7.07e-01 2.00e+00
>     e-f        5.967 5.77e-01 3.00e+00
>                delta    sigma   weight rms_deltas residual
>bond a-c       -6.033 1.00e+00 1.00e+00   5.98e+00 3.56e+01
>     ba-dada   -5.933 7.07e-01 2.00e+00
>     e-f        5.967 5.77e-01 3.00e+00
""")

def exercise():
  #exercise_bond_similarity()
  exercise_bond()
  exercise_nonbonded()
  exercise_nonbonded_cos()
  exercise_angle()
  exercise_dihedral()
  exercise_chirality()
  exercise_planarity()
  exercise_proxy_show()
  print "OK"

if (__name__ == "__main__"):
  exercise()
