from cctbx import geometry_restraints
from cctbx import crystal
from cctbx import sgtbx
from cctbx.crystal import direct_space_asu
from cctbx.array_family import flex
from scitbx import matrix
from scitbx import stl
from libtbx.test_utils import approx_equal, eps_eq
from libtbx.itertbx import count

def exercise_bond():
  p = geometry_restraints.bond_params(
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
  #
  p = geometry_restraints.bond_simple_proxy(
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
  for disable_cache in [00000, 0001]:
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
  assert sorted_asu_proxies.asu_mappings().is_locked()
  sorted_asu_proxies.push_back(proxy=sym_proxies[0])
  assert sorted_asu_proxies.simple.size() == 0
  assert sorted_asu_proxies.asu.size() == 1
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  assert not sorted_asu_proxies.process(proxy=proxies[0])
  assert not sorted_asu_proxies.process(proxy=sym_proxies[0])
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
  r = geometry_restraints.nonbonded(
    sites=[(1,2,3),(2,4,6)],
    vdw_distance=5)
  assert approx_equal(r.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(r.vdw_distance, 5)
  assert approx_equal(r.function.c_rep, 16)
  assert approx_equal(r.delta, 3.74165738677)
  assert approx_equal(r.residual(), 40.1158130612)
  assert approx_equal(r.gradients(),
    [(34.081026602378813, 68.162053204757626, 102.24307980713644),
     (-34.081026602378813, -68.162053204757626, -102.24307980713644)])
  sites_cart = flex.vec3_double([(1,2,3),(2,4,6)])
  r = geometry_restraints.nonbonded(
    sites_cart=sites_cart,
    proxy=p)
  assert approx_equal(r.sites, [(1,2,3),(2,4,6)])
  assert approx_equal(r.vdw_distance, 5)
  assert approx_equal(r.function.c_rep, 16)
  assert approx_equal(r.delta, 3.74165738677)
  proxies = geometry_restraints.shared_nonbonded_simple_proxy([p,p])
  for proxy in proxies:
    assert approx_equal(proxy.vdw_distance, 5)
  assert approx_equal(geometry_restraints.nonbonded_deltas(
    sites_cart=sites_cart,
    proxies=proxies), [3.74165738677]*2)
  assert approx_equal(geometry_restraints.nonbonded_residuals(
    sites_cart=sites_cart,
    proxies=proxies), [40.1158130612]*2)
  residual_sum = geometry_restraints.nonbonded_residual_sum(
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
  f = geometry_restraints.repulsion_function(
    c_rep=1, k_rep=4, irexp=2, rexp=3)
  assert approx_equal(f.c_rep, 1)
  assert approx_equal(f.k_rep, 4)
  assert approx_equal(f.irexp, 2)
  assert approx_equal(f.rexp, 3)
  r = geometry_restraints.nonbonded(
    sites=list(sites_cart),
    vdw_distance=p.vdw_distance,
    function=f)
  assert approx_equal(r.function.c_rep, 1)
  assert approx_equal(r.diff_vec, [-1,-1,-1])
  assert approx_equal(r.delta**2, 3)
  assert approx_equal(r.residual(), 226981)
  assert approx_equal(r.gradients(),
    [(22326.0, 22326.0, 22326.0), (-22326.0, -22326.0, -22326.0)])
  r = geometry_restraints.nonbonded(
    sites=list(sites_cart),
    vdw_distance=p.vdw_distance,
    function=geometry_restraints.repulsion_function())
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
  sorted_asu_proxies = geometry_restraints.nonbonded_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  sorted_asu_proxies.process(proxies=sym_proxies)
  assert approx_equal(
    flex.pow2(geometry_restraints.nonbonded_deltas(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies)),
    [3]*2)
  assert approx_equal(
    flex.pow2(geometry_restraints.nonbonded_deltas(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies,
      function=geometry_restraints.repulsion_function())),
    [3]*2)
  assert approx_equal(
    geometry_restraints.nonbonded_residuals(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies),
    [0.0824764182859]*2)
  assert approx_equal(
    geometry_restraints.nonbonded_residuals(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies,
      function=geometry_restraints.repulsion_function()),
    [0.0824764182859]*2)
  assert approx_equal(
    geometry_restraints.nonbonded_residual_sum(
      sites_cart=sites_cart,
      sorted_asu_proxies=sorted_asu_proxies,
      gradient_array=None),
    0.0824764182859*2)
  for disable_cache in [00000, 0001]:
    gradient_array = flex.vec3_double(2, [0,0,0])
    assert approx_equal(
      geometry_restraints.nonbonded_residual_sum(
        sites_cart=sites_cart,
        sorted_asu_proxies=sorted_asu_proxies,
        gradient_array=gradient_array,
        function=geometry_restraints.repulsion_function(),
        disable_cache=disable_cache),
      0.0824764182859*2)
    assert approx_equal(gradient_array,
      [(1.4216958630745458, 1.4216958630745458, 1.4216958630745458),
       (-1.4216958630745458, -1.4216958630745458, -1.4216958630745458)])
  #
  sorted_asu_proxies = geometry_restraints.nonbonded_sorted_asu_proxies(
    asu_mappings=asu_mappings)
  assert sorted_asu_proxies.asu_mappings().is_locked()
  assert not sorted_asu_proxies.process(proxy=proxies[0])
  assert not sorted_asu_proxies.process(proxy=sym_proxies[0])
  assert sorted_asu_proxies.simple.size() == 2
  assert sorted_asu_proxies.asu.size() == 0
  assert sorted_asu_proxies.n_total() == 2
  residual_0 = geometry_restraints.nonbonded(
    sites_cart=sites_cart,
    proxy=proxies[0]).residual()
  residual_1 = geometry_restraints.nonbonded(
    sites_cart=sites_cart,
    asu_mappings=asu_mappings,
    proxy=sym_proxies[0]).residual()
  gradient_array = flex.vec3_double(2, [0,0,0])
  assert approx_equal(geometry_restraints.nonbonded_residual_sum(
    sites_cart=sites_cart,
    sorted_asu_proxies=sorted_asu_proxies,
    gradient_array=gradient_array), residual_0+residual_1)
  assert approx_equal(gradient_array,
    [(1290.2817767146657, 1290.2817767146657, 1290.2817767146657),
     (-1290.2817767146657, -1290.2817767146657, -1290.2817767146657)])
  sorted_asu_proxies.process(sorted_asu_proxies.simple.deep_copy())
  assert sorted_asu_proxies.simple.size() == 4
  sorted_asu_proxies.process(sorted_asu_proxies.asu.deep_copy())
  assert sorted_asu_proxies.asu.size() == 0

def exercise_angle():
  p = geometry_restraints.angle_proxy(
    i_seqs=[2,1,0],
    angle_ideal=95,
    weight=1)
  assert p.i_seqs == (2,1,0)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2)
  assert approx_equal(p.angle_ideal, 95)
  assert approx_equal(p.weight, 1)
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

def exercise_dihedral():
  p = geometry_restraints.dihedral_proxy(
    i_seqs=[3,2,1,0],
    angle_ideal=-40,
    weight=1,
    periodicity=2)
  assert p.i_seqs == (3,2,1,0)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2,3)
  assert approx_equal(p.angle_ideal, -40)
  assert approx_equal(p.weight, 1)
  assert p.periodicity == 2
  p.angle_ideal = 50
  p.weight = 2
  p.periodicity = 3
  assert approx_equal(p.angle_ideal, 50)
  assert approx_equal(p.weight, 2)
  assert p.periodicity == 3
  p.angle_ideal = -40
  p.weight = 1
  p.periodicity = 2
  assert approx_equal(p.angle_ideal, -40)
  assert approx_equal(p.weight, 1)
  assert p.periodicity == 2
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
  assert d.periodicity == 2
  assert d.have_angle_model
  assert approx_equal(d.angle_model, 135)
  assert approx_equal(d.delta, 5)
  proxies = geometry_restraints.shared_dihedral_proxy([p,p])
  for proxy in proxies:
    assert proxy.periodicity == 2
  proxies[1].periodicity = 3
  assert proxies[1].periodicity == 3
  assert proxies[0].periodicity == 2
  proxies[1].periodicity = 2
  assert proxies[1].periodicity == 2
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
      if (approx_equal(p.angle_ideal, r_orig.angle_model)):
        n_equiv_direct += 1
      p_sorted = p.sort_i_seqs()
      assert p_sorted.i_seqs == (0,1,2,3)
      assert approx_equal(p_sorted.angle_ideal, r_orig.angle_model)
    n_perms += 1
    if (not perm.next_permutation()):
      break
  assert n_perms == 24
  assert n_equiv == 4
  assert n_equiv_direct == 2

def exercise_chirality():
  p = geometry_restraints.chirality_proxy(
    i_seqs=[0,2,3,1],
    volume_ideal=4,
    both_signs=00000,
    weight=1)
  assert p.i_seqs == (0,2,3,1)
  assert approx_equal(p.volume_ideal, 4)
  assert not p.both_signs
  assert approx_equal(p.weight, 1)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2,3)
  assert approx_equal(p.volume_ideal, 4)
  p = geometry_restraints.chirality_proxy(
    i_seqs=[0,2,1,3],
    volume_ideal=-4,
    both_signs=00000,
    weight=1)
  assert approx_equal(p.volume_ideal, -4)
  p = p.sort_i_seqs()
  assert p.i_seqs == (0,1,2,3)
  assert approx_equal(p.volume_ideal, 4)
  c = geometry_restraints.chirality(
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
  p = geometry_restraints.planarity_proxy(
    i_seqs=flex.size_t([0,1,2,3]),
    weights=weights)
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

def exercise():
  exercise_bond()
  exercise_nonbonded()
  exercise_angle()
  exercise_dihedral()
  exercise_chirality()
  exercise_planarity()
  print "OK"

if (__name__ == "__main__"):
  exercise()
