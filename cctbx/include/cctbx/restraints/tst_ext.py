from cctbx import restraints
from cctbx import crystal
from cctbx import sgtbx
from cctbx.crystal import direct_space_asu
from cctbx.array_family import flex
from scitbx import matrix
from scitbx import stl
from libtbx.test_utils import approx_equal, eps_eq
from libtbx.itertbx import count

def exercise_bond():
  p = restraints.bond_params(
    distance_ideal=3.5,
    weight=1)
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
  p = restraints.bond_simple_proxy(
    i_seqs=[0,1],
    distance_ideal=3.5,
    weight=1)
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
  b = restraints.bond(
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
  sites_cart = flex.vec3_double([(1,2,3),(2,4,6)])
  b = restraints.bond(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(b.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(b.distance_ideal, 3.5)
  assert approx_equal(b.weight, 1)
  assert approx_equal(b.distance_model**2, 14)
  proxies = restraints.shared_bond_simple_proxy([p,p])
  for proxy in proxies:
    assert approx_equal(proxy.weight, 1)
    proxy.weight = 12
  for proxy in proxies:
    assert approx_equal(proxy.weight, 12)
    proxy.weight = 1
  assert approx_equal(restraints.bond_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [-0.241657386774]*2)
  assert approx_equal(restraints.bond_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [0.0583982925824]*2)
  residual_sum = restraints.bond_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*0.0583982925824)
  #
  sites_cart = flex.vec3_double([[1,2,3],[2,3,4]])
  asu_mappings = direct_space_asu.non_crystallographic_asu_mappings(
    sites_cart=sites_cart)
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=5)
  p = restraints.bond_asu_proxy(
    pair=pair_generator.next(),
    distance_ideal=2,
    weight=10)
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
  sym_proxies = restraints.shared_bond_asu_proxy([p,p])
  for proxy in sym_proxies:
    assert approx_equal(proxy.distance_ideal, 2)
    proxy.distance_ideal = -4
  for proxy in sym_proxies:
    assert approx_equal(proxy.distance_ideal, -4)
    proxy.distance_ideal = 2
  assert approx_equal(
    restraints.bond_deltas(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies),
    [2-3**.5]*2)
  assert approx_equal(
    restraints.bond_residuals(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies),
    [10*(2-3**.5)**2]*2)
  assert approx_equal(
    restraints.bond_residual_sum(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies,
      gradient_array=None),
    (10*(2-3**.5)**2)*2)
  gradient_array = flex.vec3_double(2, [0,0,0])
  assert approx_equal(
    restraints.bond_residual_sum(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies,
      gradient_array=gradient_array),
    (10*(2-3**.5)**2)*2)
  assert approx_equal(gradient_array,
    [[ 6.1880215351700611]*3,
     [-6.1880215351700611]*3])
  for disable_cache in [00000, 0001]:
    gradient_array = flex.vec3_double(2, [0,0,0])
    assert approx_equal(
      restraints.bond_residual_sum(
        sites_cart=sites_cart,
        asu_mappings=asu_mappings,
        proxies=sym_proxies,
        gradient_array=gradient_array,
        disable_cache=disable_cache),
      (10*(2-3**.5)**2)*2)
    assert approx_equal(gradient_array,
      [[ 6.1880215351700611]*3,
       [-6.1880215351700611]*3])
  #
  sorted_asu_proxies = restraints.bond_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  assert sorted_asu_proxies.asu_mappings().is_locked()
  sorted_asu_proxies.push_back(proxy=sym_proxies[0])
  assert sorted_asu_proxies.simple.size() == 0
  assert sorted_asu_proxies.sym.size() == 1
  sorted_asu_proxies = restraints.bond_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  assert not sorted_asu_proxies.process(proxy=proxies[0])
  assert not sorted_asu_proxies.process(proxy=sym_proxies[0])
  assert sorted_asu_proxies.simple.size() == 2
  assert sorted_asu_proxies.sym.size() == 0
  assert sorted_asu_proxies.n_total() == 2
  residual_0 = restraints.bond(
    sites_cart=sites_cart,
    proxy=proxies[0]).residual()
  residual_1 = restraints.bond(
    sites_cart=sites_cart,
    asu_mappings=asu_mappings,
    proxy=sym_proxies[0]).residual()
  assert approx_equal(residual_1, 10*(2-3**.5)**2)
  gradient_array = flex.vec3_double(2, [0,0,0])
  assert approx_equal(restraints.bond_residual_sum(
    sites_cart=sites_cart,
    sorted_asu_proxies=sorted_asu_proxies,
    gradient_array=gradient_array), residual_0+residual_1)
  assert approx_equal(gradient_array,
    [(5.1354626519124107, 5.1354626519124107, 5.1354626519124107),
     (-5.1354626519124107, -5.1354626519124107, -5.1354626519124107)])

def exercise_bond_tables():
  assert "map_indexing_suite_bond_params_dict_entry" in restraints.__dict__
  assert "map_indexing_suite_bond_sym_dict_entry" in restraints.__dict__
  assert "map_indexing_suite_bond_asu_dict_entry" in restraints.__dict__
  t = restraints.bond_params_table()
  assert t.size() == 0
  d = restraints.bond_params_dict()
  assert len(d) == 0
  p = restraints.bond_params(distance_ideal=3, weight=2)
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
  #
  d = restraints.bond_sym_dict()
  assert len(d) == 0
  sym_ops = sgtbx.space_group("P 41").all_ops()
  for i,j_sym in enumerate([10,18,13]):
    d[j_sym] = restraints.bond_sym_ops(sym_ops[:i])
    assert len(d) == i+1
    assert len(d[j_sym]) == i
    assert [str(s) for s in sym_ops[:i]] == [str(s) for s in d[j_sym]]
    d[j_sym] = sym_ops[:i]
    assert [str(s) for s in sym_ops[:i]] == [str(s) for s in d[j_sym]]
  assert [item.key() for item in d] == [10,13,18]
  assert d[13].size() == 2
  d[13].append(sym_ops[-1])
  assert d[13].size() == 3
  del d[13][0]
  assert d[13].size() == 2
  d[13].clear()
  assert d[13].size() == 0
  t = restraints.bond_sym_table()
  t.append(d)
  assert t.size() == 1
  assert len(t[0][10]) == 0
  t.append(d)
  assert t.size() == 2
  assert len(t[1][18]) == 1
  t = restraints.bond_sym_table(3)
  for d in t:
    assert len(d) == 0
  t[1][10] = sym_ops[:2]
  assert len(t[1]) == 1
  assert len(t[1][10]) == 2
  #
  t = restraints.bond_asu_table(3)
  for d in t:
    assert len(d) == 0
  t[1][10] = restraints.bond_asu_dict()
  assert t[1][10].size() == 0
  t[1][10].append(restraints.bond_asu_j_sym_groups())
  assert t[1][10].size() == 1
  assert t[1][10][0].size() == 0
  t[1][10][0].insert(3)
  assert t[1][10][0].size() == 1
  t[1][10].append(restraints.bond_asu_j_sym_groups())
  assert t[1][10][1].size() == 0
  t[1][10][1].insert([4,5,4])
  assert t[1][10][1].size() == 2

def exercise_repulsion():
  p = restraints.repulsion_simple_proxy(
    i_seqs=[0,1],
    vdw_radius=5)
  assert p.i_seqs == (0,1)
  assert approx_equal(p.vdw_radius, 5)
  r = restraints.repulsion(
    sites=[(1,2,3),(2,4,6)],
    vdw_radius=5)
  assert approx_equal(r.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(r.vdw_radius, 5)
  assert approx_equal(r.function.c_rep, 16)
  assert approx_equal(r.delta, 3.74165738677)
  assert approx_equal(r.residual(), 40.1158130612)
  assert approx_equal(r.gradients(),
    [(34.081026602378813, 68.162053204757626, 102.24307980713644),
     (-34.081026602378813, -68.162053204757626, -102.24307980713644)])
  sites_cart = flex.vec3_double([(1,2,3),(2,4,6)])
  r = restraints.repulsion(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(r.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(r.vdw_radius, 5)
  assert approx_equal(r.function.c_rep, 16)
  assert approx_equal(r.delta, 3.74165738677)
  proxies = restraints.shared_repulsion_simple_proxy([p,p])
  for proxy in proxies:
    assert approx_equal(proxy.vdw_radius, 5)
  assert approx_equal(restraints.repulsion_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [3.74165738677]*2)
  assert approx_equal(restraints.repulsion_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [40.1158130612]*2)
  residual_sum = restraints.repulsion_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*40.1158130612)
  #
  sites_cart = flex.vec3_double([[1,2,3],[2,3,4]])
  asu_mappings = direct_space_asu.non_crystallographic_asu_mappings(
    sites_cart=sites_cart)
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=5)
  p = restraints.repulsion_asu_proxy(
    pair=pair_generator.next(),
    vdw_radius=2)
  assert pair_generator.at_end()
  assert p.i_seq == 0
  assert p.j_seq == 1
  assert p.j_sym == 0
  assert approx_equal(p.vdw_radius, 2)
  p.vdw_radius = 3
  assert approx_equal(p.vdw_radius, 3)
  p.vdw_radius = 2
  sym_proxies = restraints.shared_repulsion_asu_proxy([p,p])
  for proxy in sym_proxies:
    assert approx_equal(proxy.vdw_radius, 2)
    proxy.vdw_radius = 3
  for proxy in sym_proxies:
    assert approx_equal(proxy.vdw_radius, 3)
    proxy.vdw_radius = 2
  f = restraints.repulsion_function(
    c_rep=1, k_rep=4, irexp=2, rexp=3)
  assert approx_equal(f.c_rep, 1)
  assert approx_equal(f.k_rep, 4)
  assert approx_equal(f.irexp, 2)
  assert approx_equal(f.rexp, 3)
  r = restraints.repulsion(
    sites=list(sites_cart),
    vdw_radius=p.vdw_radius,
    function=f)
  assert approx_equal(r.function.c_rep, 1)
  assert approx_equal(r.diff_vec, [-1,-1,-1])
  assert approx_equal(r.delta**2, 3)
  assert approx_equal(r.residual(), 226981)
  assert approx_equal(r.gradients(),
    [(22326.0, 22326.0, 22326.0), (-22326.0, -22326.0, -22326.0)])
  r = restraints.repulsion(
    sites=list(sites_cart),
    vdw_radius=p.vdw_radius,
    function=restraints.repulsion_function())
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
  assert approx_equal(
    flex.pow2(restraints.repulsion_deltas(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies)),
    [3]*2)
  assert approx_equal(
    flex.pow2(restraints.repulsion_deltas(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies,
      function=restraints.repulsion_function())),
    [3]*2)
  assert approx_equal(
    restraints.repulsion_residuals(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies),
    [0.0824764182859]*2)
  assert approx_equal(
    restraints.repulsion_residuals(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies,
      function=restraints.repulsion_function()),
    [0.0824764182859]*2)
  assert approx_equal(
    restraints.repulsion_residual_sum(
      sites_cart=sites_cart,
      asu_mappings=asu_mappings,
      proxies=sym_proxies,
      gradient_array=None),
    0.0824764182859*2)
  for disable_cache in [00000, 0001]:
    gradient_array = flex.vec3_double(2, [0,0,0])
    assert approx_equal(
      restraints.repulsion_residual_sum(
        sites_cart=sites_cart,
        asu_mappings=asu_mappings,
        proxies=sym_proxies,
        gradient_array=gradient_array,
        function=restraints.repulsion_function(),
        disable_cache=disable_cache),
      0.0824764182859*2)
    assert approx_equal(gradient_array,
      [(1.4216958630745458, 1.4216958630745458, 1.4216958630745458),
       (-1.4216958630745458, -1.4216958630745458, -1.4216958630745458)])
  #
  sorted_asu_proxies = restraints.repulsion_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  assert sorted_asu_proxies.asu_mappings().is_locked()
  assert not sorted_asu_proxies.process(proxy=proxies[0])
  assert not sorted_asu_proxies.process(proxy=sym_proxies[0])
  assert sorted_asu_proxies.simple.size() == 2
  assert sorted_asu_proxies.sym.size() == 0
  assert sorted_asu_proxies.n_total() == 2
  residual_0 = restraints.repulsion(
    sites_cart=sites_cart,
    proxy=proxies[0]).residual()
  residual_1 = restraints.repulsion(
    sites_cart=sites_cart,
    asu_mappings=asu_mappings,
    proxy=sym_proxies[0]).residual()
  gradient_array = flex.vec3_double(2, [0,0,0])
  assert approx_equal(restraints.repulsion_residual_sum(
    sites_cart=sites_cart,
    sorted_asu_proxies=sorted_asu_proxies,
    gradient_array=gradient_array), residual_0+residual_1)
  assert approx_equal(gradient_array,
    [(1290.2817767146657, 1290.2817767146657, 1290.2817767146657),
     (-1290.2817767146657, -1290.2817767146657, -1290.2817767146657)])

def exercise_angle():
  p = restraints.angle_proxy(
    i_seqs=[0,1,2],
    angle_ideal=95,
    weight=1)
  assert p.i_seqs == (0,1,2)
  assert approx_equal(p.angle_ideal, 95)
  assert approx_equal(p.weight, 1)
  a = restraints.angle(
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
  sites_cart = flex.vec3_double([(1,0,0),(0,0,0),(0,1,0)])
  a = restraints.angle(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(a.sites, [(1,0,0),(0,0,0),(0,1,0)])
  assert approx_equal(a.angle_ideal, 95)
  assert approx_equal(a.weight, 1)
  assert a.have_angle_model
  assert approx_equal(a.angle_model, 90)
  proxies = restraints.shared_angle_proxy([p,p])
  for proxy in proxies:
    assert approx_equal(proxy.weight, 1)
  assert approx_equal(restraints.angle_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [5]*2)
  assert approx_equal(restraints.angle_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [25]*2)
  residual_sum = restraints.angle_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*25)

def exercise_dihedral():
  p = restraints.dihedral_proxy(
    i_seqs=[0,1,2,3],
    angle_ideal=-40,
    weight=1,
    periodicity=2)
  assert p.i_seqs == (0,1,2,3)
  assert approx_equal(p.angle_ideal, -40)
  assert approx_equal(p.weight, 1)
  assert p.periodicity == 2
  d = restraints.dihedral(
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
  d = restraints.dihedral(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(d.sites, [(1,0,0),(0,0,0),(0,1,0),(-1,0,-1)])
  assert approx_equal(d.angle_ideal, -40)
  assert approx_equal(d.weight, 1)
  assert d.periodicity == 2
  assert d.have_angle_model
  assert approx_equal(d.angle_model, 135)
  assert approx_equal(d.delta, 5)
  proxies = restraints.shared_dihedral_proxy([p,p])
  for proxy in proxies:
    assert proxy.periodicity == 2
  assert approx_equal(restraints.dihedral_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [5]*2)
  assert approx_equal(restraints.dihedral_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [25]*2)
  residual_sum = restraints.dihedral_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*25)

def exercise_chirality():
  p = restraints.chirality_proxy(
    i_seqs=[0,1,2,3],
    volume_ideal=4,
    both_signs=00000,
    weight=1)
  assert p.i_seqs == (0,1,2,3)
  assert approx_equal(p.volume_ideal, 4)
  assert not p.both_signs
  assert approx_equal(p.weight, 1)
  c = restraints.chirality(
    sites=[(1,0,0),(0,0,0),(0,1,0),(1,0,1)],
    volume_ideal=4,
    both_signs=00000,
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
  c = restraints.chirality(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(c.sites, [(1,0,0),(0,0,0),(0,1,0),(-1,0,-1)])
  assert approx_equal(c.volume_ideal, 4)
  assert approx_equal(c.weight, 1)
  assert not c.both_signs
  assert approx_equal(c.volume_model, 1)
  assert approx_equal(c.delta_sign, -1)
  assert approx_equal(c.delta, 3)
  proxies = restraints.shared_chirality_proxy([p,p])
  for proxy in proxies:
    assert proxy.volume_ideal == 4
  assert approx_equal(restraints.chirality_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [3]*2)
  assert approx_equal(restraints.chirality_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [9]*2)
  residual_sum = restraints.chirality_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert approx_equal(residual_sum, 2*9)

def exercise_planarity():
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
  p = restraints.planarity_proxy(
    i_seqs=flex.size_t([0,1,2,3]),
    weights=weights)
  assert tuple(p.i_seqs) == (0,1,2,3)
  assert approx_equal(p.weights, weights)
  for i_constructor in xrange(2):
    if (i_constructor == 0):
      l = restraints.planarity(
        sites=sites_cart,
        weights=weights)
    else:
      l = restraints.planarity(
        sites_cart=sites_cart,
        proxy=p)
    assert approx_equal(l.sites, sites_cart)
    assert approx_equal(l.weights, weights)
    assert eps_eq(l.deltas(), (0.014233272168667327, 0.007241647943016986,
                               -0.0048894168534149443, -0.0035120793736139956))
    assert eps_eq(l.residual(), expected_residual)
    assert eps_eq(l.gradients(), expected_gradients)
    assert eps_eq(l.normal(),
      (0.69097523765119184, 0.086887122267422581, -0.71763768639680903))
    assert eps_eq(l.residual(), l.lambda_min())
    assert eps_eq(l.center_of_mass(),
      (-5.7061446613913009, 0.11105869285849694, -0.42071347654387559))
    assert eps_eq(l.residual_tensor(),
      (10.250312599815825, 8.7000194514224525, 10.265208176541265,
       2.7229147081229312, 10.19874296603952, 3.6750425846794936))
    assert eps_eq(l.eigensystem().values(),
      [21.998140770294835, 7.2169709305206142, 0.00042852696409348911])
  proxies = restraints.shared_planarity_proxy([p,p])
  for proxy in proxies:
    assert tuple(proxy.i_seqs) == (0,1,2,3)
  assert eps_eq(restraints.planarity_deltas_rms(
    sites_cart=sites_cart,
    proxies=proxies), [0.0085332965530764398]*2)
  assert eps_eq(restraints.planarity_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [expected_residual]*2)
  residual_sum = restraints.planarity_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=None)
  assert eps_eq(residual_sum, 2*expected_residual)
  gradient_array = flex.vec3_double(proxy.i_seqs.size(), (0,0,0))
  residual_sum = restraints.planarity_residual_sum(
    sites_cart=sites_cart,
    proxies=proxies,
    gradient_array=gradient_array)
  assert eps_eq(residual_sum, 2*expected_residual)
  for g,e in zip(gradient_array, expected_gradients):
    assert eps_eq(g, matrix.col(e)*2)

def exercise_bonded_interactions():
  bond_proxies = restraints.shared_bond_simple_proxy()
  for pairs in [(0,1),(1,2),(2,3),(2,4),(3,5),(3,9),(5,6),(5,7),
                (4,8),(4,9),(8,10),(9,10),(9,11),(10,11)]:
    bond_proxies.append(restraints.bond_simple_proxy(
      i_seqs=pairs,
      distance_ideal=1,
      weight=1))
  bond_sets = restraints.bond_sets(
    n_sites=12,
    proxies=bond_proxies)
  expected_1_2_all = [
    (1,),
    (0, 2),
    (1, 3, 4),
    (2, 5, 9),
    (2, 8, 9),
    (3, 6, 7),
    (5,),
    (5,),
    (4, 10),
    (3, 4, 10, 11),
    (8, 9, 11),
    (9, 10)]
  assert [tuple(s) for s in bond_sets] == expected_1_2_all
  expected_1_2_elim = [
    (1,),
    (2,),
    (3, 4),
    (5, 9),
    (8, 9),
    (6, 7),
    (),
    (),
    (10,),
    (10, 11),
    (11,),
    ()]
  expected_1_3_all = [
    (2,),
    (3, 4),
    (0, 5, 8, 9),
    (1, 4, 6, 7, 10, 11),
    (1, 3, 10, 11),
    (2, 9),
    (3, 7),
    (3, 6),
    (2, 9, 11),
    (2, 5, 8),
    (3, 4),
    (3, 4, 8)]
  expected_1_3_elim = [
    (2,),
    (3, 4),
    (5, 8, 9),
    (4, 6, 7, 10, 11),
    (10, 11),
    (9,),
    (7,),
    (),
    (9, 11),
    (),
    (),
    ()]
  expected_1_4_all = [
    (3, 4),
    (5, 8, 9),
    (6, 7, 10, 11),
    (0, 8),
    (0, 5),
    (1, 4, 10, 11),
    (2, 9),
    (2, 9),
    (1, 3),
    (1, 6, 7),
    (2, 5),
    (2, 5)]
  expected_1_4_elim = [
    (3, 4),
    (5, 8, 9),
    (6, 7, 10, 11),
    (8,),
    (5,),
    (10, 11),
    (9,),
    (9,),
    (),
    (),
    (),
    ()]
  for i_seq_0 in xrange(bond_sets.size()):
    bonded_interactions = restraints.bonded_interactions(
      bond_sets=bond_sets,
      i_seq_0=i_seq_0)
    assert bonded_interactions.i_seq_0() == i_seq_0
    assert not bonded_interactions.eliminate_redundant_interactions()
    assert tuple(bonded_interactions.interactions_1_2()) \
        == expected_1_2_all[i_seq_0]
    assert tuple(bonded_interactions.interactions_1_3()) \
        == expected_1_3_all[i_seq_0]
    assert tuple(bonded_interactions.interactions_1_4()) \
        == expected_1_4_all[i_seq_0]
    assert bonded_interactions.interaction_type_of(j_seq=20) == 0
    assert bonded_interactions.interaction_type_of(
      j_seq=expected_1_2_all[i_seq_0][0]) == 2
    assert bonded_interactions.interaction_type_of(
      j_seq=expected_1_3_all[i_seq_0][0]) == 3
    assert bonded_interactions.interaction_type_of(
      j_seq=expected_1_4_all[i_seq_0][0]) == 4
    all_i_seqs = {}
    for i_seq in bond_sets[i_seq_0]:
      all_i_seqs[i_seq] = 0
    for i_seq in expected_1_3_all[i_seq_0]:
      assert all_i_seqs.get(i_seq, None) is None
      all_i_seqs[i_seq] = 0
    for i_seq in expected_1_4_all[i_seq_0]:
      assert all_i_seqs.get(i_seq, None) is None
      all_i_seqs[i_seq] = 0
    bonded_interactions = restraints.bonded_interactions(
      bond_sets=bond_sets,
      i_seq_0=i_seq_0,
      eliminate_redundant_interactions=0001)
    assert bonded_interactions.i_seq_0() == i_seq_0
    assert bonded_interactions.eliminate_redundant_interactions()
    assert tuple(bonded_interactions.interactions_1_2()) \
        == expected_1_2_elim[i_seq_0]
    assert tuple(bonded_interactions.interactions_1_3()) \
        == expected_1_3_elim[i_seq_0]
    assert tuple(bonded_interactions.interactions_1_4()) \
        == expected_1_4_elim[i_seq_0]

def exercise():
  exercise_bond()
  exercise_repulsion()
  exercise_angle()
  exercise_dihedral()
  exercise_chirality()
  exercise_planarity()
  exercise_bond_tables()
  exercise_bonded_interactions()
  print "OK"

if (__name__ == "__main__"):
  exercise()
