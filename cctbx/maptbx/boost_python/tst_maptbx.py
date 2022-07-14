from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex
from cctbx import crystal
from scitbx import matrix
from libtbx.test_utils import Exception_expected, approx_equal, \
  not_approx_equal
from libtbx.utils import n_dim_index_from_one_dim
import itertools
import time
import sys
import random
import math
from six.moves import range
from six.moves import zip

def flex_types():
  return (flex.float, flex.double)

def exercise_copy():
  for flex_type in flex_types():
    m = flex_type((1,2,3,4))
    m.resize(flex.grid(2,2))
    c = maptbx.copy(map=m, result_grid=m.accessor())
    assert tuple(m) == tuple(c)
    c = maptbx.copy(map=m, result_grid=flex.grid(2,3).set_focus(2,2))
    assert approx_equal(tuple(c), (1,2,0,3,4,0))
    n = maptbx.copy(c, result_grid=m.accessor())
    assert approx_equal(tuple(m), tuple(n))
    c = maptbx.copy(m, flex.grid(3,2).set_focus(2,2))
    assert approx_equal(tuple(c), (1,2,3,4,0,0))
    n = maptbx.copy(c, m.accessor())
    assert approx_equal(tuple(m), tuple(n))
    m = flex_type((1,2,3,4,5,6))
    m.resize(flex.grid((1,2),(3,5)))
    c = maptbx.copy(m, m.accessor())
    assert approx_equal(tuple(m), tuple(c))
    c = maptbx.copy(m, flex.grid((1,2),(3,6)).set_focus(3,5))
    assert approx_equal(tuple(c), (1,2,3,0,4,5,6,0))
    n = maptbx.copy(c, m.accessor())
    assert approx_equal(tuple(m), tuple(n))
    c = maptbx.copy(m, flex.grid((1,2),(4,5)).set_focus(3,5))
    assert approx_equal(tuple(c), (1,2,3,4,5,6,0,0,0))
    n = maptbx.copy(c, m.accessor())
    assert approx_equal(tuple(m), tuple(n))
    #
    m = flex_type()
    for i in range(2):
      for j in range(3):
        for k in range(5):
          m.append(i*100+j*10+k)
    m.resize(flex.grid(2,3,5).set_focus((2,3,4)))
    for i in range(-5,5):
      for j in range(-5,5):
        for k in range(-5,5):
          c = maptbx.copy(map_unit_cell=m, first=(i,j,k), last=(i,j,k))
          assert c.size() == 1
          assert c[(i,j,k)] == m[(i%2,j%3,k%4)]
    c = maptbx.copy(map_unit_cell=m, first=(-1,1,-2), last=(1,2,0))
    assert list(c) == [112, 113, 110, 122, 123, 120,  12,  13,  10,
                        22,  23,  20, 112, 113, 110, 122, 123, 120]
    #
    m2 = m.deep_copy()
    grid = flex.grid( (-1,-1,-1), (1,2,4) ).set_focus( (1,2,3) )
    m2.resize(grid)
    for i in range(-1,1):
      for j in range(-1,2):
        for k in range(-1,3):
          # aperiodic copy
          c = maptbx.copy_box(map=m2, first=(i,j,k), last=(i,j,k))
          assert c.size() == 1
          ind = ((i+1)%2-1,(j+1)%3-1,(k+1)%4-1)
          assert c[(i,j,k)] == m2[ind]
    c = maptbx.copy_box(map=m2, first=(-1,0,-1), last=(0,1,2))
    assert list(c) == [10, 11, 12, 13, 20, 21,  22,  23,  110,
                       111, 112, 113, 120, 121, 122, 123]
    #
    for n0 in range(4):
      for n1 in range(4):
        for n2 in range(4):
          for d2 in range(3):
            g = flex.grid((n0,n1,n2+d2)).set_focus((n0,n1,n2))
            map1 = flex_type(range(1,1+g.size_1d()))
            map1.resize(g)
            map2 = map1.deep_copy()
            maptbx.unpad_in_place(map=map2)
            assert map2.all() == (n0,n1,n2)
            assert not map2.is_padded()
            if (n0*n1*n2 != 0):
              for i in flex.nested_loop((n0,n1,n2)):
                assert map2[i] == map1[i]
    n0,n1,n2,d2 = 2,3,4,1
    g = flex.grid((n0,n1,n2+d2)).set_focus((n0,n1,n2))
    map1 = flex_type(range(1,1+g.size_1d()))
    map1.resize(g)
    map2 = map1.deep_copy()
    maptbx.unpad_in_place(map=map2)
    assert map2.all() == (n0,n1,n2)
    assert not map2.is_padded()
    assert list(map2) == [
       1, 2, 3, 4,
       6, 7, 8, 9,
      11,12,13,14,
      16,17,18,19,
      21,22,23,24,
      26,27,28,29]

def exercise_statistics():
  import scitbx.math
  for flex_type in flex_types():
    a = flex_type(flex.grid((3,5)))
    s = maptbx.statistics(a)
    assert s.min() == 0
    assert s.max() == 0
    assert s.mean() == 0
    assert s.mean_sq() == 0
    assert s.sigma() == 0
    a = flex_type([random.random() for i in range(3*5)])
    a.resize(flex.grid((3,5)))
    s = maptbx.statistics(a)
    assert approx_equal(flex.min(a), s.min())
    assert approx_equal(flex.max(a), s.max())
    assert approx_equal(flex.mean(a), s.mean())
    assert approx_equal(flex.mean_sq(a), s.mean_sq())
    assert approx_equal(flex.mean_sq(a)-flex.mean(a)**2, s.sigma()**2)
    b = flex_type(flex.grid((4,6)).set_focus((3,5)))
    for i in range(3):
      for j in range(5):
        b[(i,j)] = a[(i,j)]
    b[(3,5)] = -1
    b[(2,5)] = 2
    b.resize(flex.grid((-2,3), (2,9)).set_focus((1,8)))
    t = maptbx.statistics(b)
    assert not_approx_equal(flex.min(b), t.min())
    assert not_approx_equal(flex.max(b), t.max())
    assert not_approx_equal(flex.mean(b), t.mean())
    assert not_approx_equal(flex.mean_sq(b), t.mean_sq())
    assert not_approx_equal(flex.mean_sq(b)-flex.mean(b)**2, t.sigma()**2)
    assert approx_equal(s.min(), t.min())
    assert approx_equal(s.max(), t.max())
    assert approx_equal(s.mean(), t.mean())
    assert approx_equal(s.mean_sq(), t.mean_sq())
    assert approx_equal(s.sigma(), t.sigma())
  a = flex.double(flex.grid(5,3))
  s = maptbx.more_statistics(a)
  assert s.min() == 0
  assert s.max() == 0
  assert s.mean() == 0
  assert s.mean_sq() == 0
  assert s.sigma() == 0
  assert s.skewness() == 0
  assert s.kurtosis() == 0
  a = flex.random_double(5*3)
  reference = scitbx.math.basic_statistics(a)
  a.resize(flex.grid(5,3))
  s = maptbx.more_statistics(a)
  assert approx_equal(s.min(), reference.min)
  assert approx_equal(s.max(), reference.max)
  assert approx_equal(s.mean(), reference.mean)
  assert approx_equal(s.sigma(), reference.biased_standard_deviation)
  assert approx_equal(s.skewness(), reference.skew)
  assert approx_equal(s.kurtosis(), reference.kurtosis)
  b = flex.double(flex.grid((6,4)).set_focus((5,3)))
  for i in range(5):
    for j in range(3):
      b[(i,j)] = a[(i,j)]
  b[(5,3)] = -1
  b[(5,2)] = 2
  b.resize(flex.grid((3,-2), (9,2)).set_focus((8,1)))
  t = maptbx.statistics(b)
  assert approx_equal(s.min(), reference.min)
  assert approx_equal(s.max(), reference.max)
  assert approx_equal(s.mean(), reference.mean)
  assert approx_equal(s.sigma(), reference.biased_standard_deviation)
  assert approx_equal(s.skewness(), reference.skew)
  assert approx_equal(s.kurtosis(), reference.kurtosis)
  m = flex.double(flex.grid((6,4,8)).set_focus((5,3,7)))

def exercise_grid_tags():
  t = maptbx.grid_tags((8,10,12))
  assert not t.is_valid()
  assert t.tag_array().all() == (8,10,12)
  s = sgtbx.space_group_info("P 21")
  for i_flags in range(8):
    f = sgtbx.search_symmetry_flags(
      use_space_group_symmetry=i_flags % 2 != 0,
      use_space_group_ltr=0,
      use_seminvariants=(i_flags//4) % 2 != 0,
      use_normalizer_k2l=(i_flags//2) % 2 != 0,
      use_normalizer_l2n=False)
    t.build(s.type(), f)
    assert t.is_valid()
    assert t.space_group_type().group() == s.group()
    assert t.symmetry_flags() == f
    if (f.use_seminvariants()):
      assert [(vm.v, vm.m) for vm in t.grid_ss_continuous()] \
          == [((0, 1, 0), 10)]
    assert t.n_grid_misses() == 0
    assert t.n_independent() == (960, 480, 484, 242, 24, 14, 14, 14)[i_flags]
    assert t.n_independent() + t.n_dependent() == t.tag_array().size()
    for flex_type in flex_types():
      d = flex_type(t.tag_array().accessor())
      assert t.n_dependent() == 0 \
          or approx_equal(t.dependent_correlation(d, 1.e-10).coefficient(),0)
      assert t.verify(d, 0)
      t.sum_sym_equiv_points(d)
    if (i_flags == 0):
      assert t.n_independent() == t.tag_array().size()
    else:
      assert t.n_independent() < t.tag_array().size()
      for flex_type in flex_types():
        d = flex_type([random.random() for x in range(t.tag_array().size())])
        d.resize(t.tag_array().accessor())
        assert not t.verify(d)
        t.sum_sym_equiv_points(d)
        assert t.verify(d)

def exercise_peak_search():
  t = flex.long(flex.grid((3,4,5)))
  for flex_type in flex_types():
    d = flex_type(flex.grid((3,4,5)))
    l = maptbx.peak_list(d, t, peak_search_level=0, interpolate=False)
    assert l.gridding() == d.focus()
    assert l.grid_indices(0) == (0,0,0)
    assert list(l.grid_heights()) == [0]
    assert list(l.sites()) == [(0,0,0)]
    assert list(l.heights()) == [0]
    l = maptbx.peak_list(
      d, t, peak_search_level=0,peak_cutoff=-1, interpolate=False)
    assert l.gridding() == d.focus()
    assert l.grid_indices(0) == (0,0,0)
    assert list(l.grid_heights()) == [0]
    assert list(l.sites()) == [(0,0,0)]
    assert list(l.heights()) == [0]

def exercise_pymol_interface():
  for flex_type in flex_types():
    m = flex_type(flex.grid(3,4,6).set_focus(3,4,5))
    o = maptbx.as_CObjectZYX(m, first=(0,0,0), last=(4,5,6))

def exercise_structure_factors():
  uc = uctbx.unit_cell((11,13,17))
  sg = sgtbx.space_group_info("P 31")
  mi = flex.miller_index(((1,2,3),(2,3,4)))
  d = flex.complex_double((1+2j, 2+3j))
  for anomalous_flag in (False, True):
    for conjugate_flag in (False, True):
      t = maptbx.structure_factors.to_map(
        space_group=sg.group(),
        anomalous_flag=anomalous_flag,
        miller_indices=mi,
        structure_factors=d,
        n_real=(11,11,9),
        map_grid=flex.grid(11,11,9),
        conjugate_flag=conjugate_flag)
      assert t.complex_map().focus() == (11,11,9)
      t = maptbx.structure_factors.to_map(
        space_group=sg.group(),
        anomalous_flag=anomalous_flag,
        miller_indices=mi,
        structure_factors=d,
        n_real=(11,11,9),
        map_grid=flex.grid(11,11,9),
        conjugate_flag=conjugate_flag,
        treat_restricted=False)
      assert t.complex_map().focus() == (11,11,9)
      f = maptbx.structure_factors.from_map(
        unit_cell=uc,
        space_group_type=sg.type(),
        anomalous_flag=anomalous_flag,
        d_min=5.,
        complex_map=t.complex_map(),
        conjugate_flag=conjugate_flag)
      assert f.miller_indices().size() > 0
      assert f.miller_indices().size() == f.data().size()
      f = maptbx.structure_factors.from_map(
        anomalous_flag=anomalous_flag,
        miller_indices=mi,
        complex_map=t.complex_map(),
        conjugate_flag=conjugate_flag)
      assert f.miller_indices().size() == 0
      assert f.data().size() == mi.size()
      f = maptbx.structure_factors.from_map(
        anomalous_flag=anomalous_flag,
        miller_indices=mi,
        complex_map=t.complex_map(),
        conjugate_flag=conjugate_flag,
        allow_miller_indices_outside_map=True)
      assert f.miller_indices().size() == 0
      assert f.data().size() == mi.size()
      assert f.n_indices_affected_by_aliasing() == 0
      assert f.outside_map().size() == 0
      f = maptbx.structure_factors.from_map(
        space_group=sg.group(),
        anomalous_flag=anomalous_flag,
        miller_indices=mi,
        complex_map=t.complex_map(),
        conjugate_flag=conjugate_flag)
      assert f.miller_indices().size() == 0
      assert f.data().size() == mi.size()
      assert f.n_indices_affected_by_aliasing() == 0
      assert f.outside_map().size() == 0

def exercise_fft():
  sg = sgtbx.space_group_info("P 31").group()
  mi = flex.miller_index(((1,2,3),(2,3,4)))
  d = flex.complex_double((1+2j, 2+3j))
  map = maptbx.fft_to_real_map_unpadded(
    space_group=sg,
    n_real=(10,13,17),
    miller_indices=mi,
    data=d)
  assert map.is_0_based()
  assert not map.is_padded()
  assert map.focus() == (10,13,17)
  assert approx_equal(
    map[:5], [6,13.4163896,5.9603989,-8.1028328,-13.1838857])
  assert approx_equal(
    map[1000:1005], [-17.5217557,-1.4971115,20.1455825,-4.0350021,-2.7312275])
  assert approx_equal(
    map[-5:], [8.7337234,-0.7930404,-6.6343761,-1.9521735,2.6725642])

def exercise_gridding():
  u = uctbx.unit_cell((4,6,7))
  assert maptbx.ext.determine_gridding(u, 2, 1/3., (1,1,1), 5, True) \
      == (8,9,12)
  f = sgtbx.search_symmetry_flags(
    use_space_group_symmetry=True,
    use_space_group_ltr=0,
    use_seminvariants=False,
    use_normalizer_k2l=True,
    use_normalizer_l2n=False)
  t = sgtbx.space_group_info("F 2 2 2").primitive_setting().type()
  assert maptbx.ext.determine_gridding(u, 2, 1/3., f, t, (1,1,1), 5, True) \
      == (12, 12, 12)

def exercise_misc():
  for flex_type in flex_types():
    m = flex_type([1,2,-3,4,-5,6])
    maptbx.set_if_less_than(m, 0, 0)
    assert approx_equal(tuple(m), (1,2,0,4,0,6))
    maptbx.set_if_less_than(m, 2, 9)
    assert approx_equal(tuple(m), (9,2,9,4,9,6))
  from cctbx import xray
  structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(10.0,10.0,10.0,90,90,90),
        space_group_symbol="P1")),
    scatterers=flex.xray_scatterer([
      xray.scatterer(
        label="O",
        site=(0.5,0.5,0.5),
        u=0.2)]))
  fc = structure.structure_factors(d_min=2).f_calc()
  fc_map = fc.fft_map(resolution_factor=1/4.)
  fc_map.apply_sigma_scaling()
  real_map = fc_map.real_map_unpadded()
  stats = maptbx.spherical_variance_around_point(
    real_map=real_map,
    unit_cell=structure.unit_cell(),
    site_cart=(5.,5.,5.),
    radius=1.)
  # XXX exact numbers are *not* consistent across platforms!
  assert (approx_equal(stats.min, 8.3, eps=0.2))
  assert (approx_equal(stats.mean, 8.5, eps=0.1))
  assert (stats.standard_deviation < 0.15)
  stats = maptbx.spherical_variance_around_point(
    real_map=real_map,
    unit_cell=structure.unit_cell(),
    site_cart=(6.,6.,6.),
    radius=1.)
  assert (approx_equal(stats.min, -0.75, eps=0.15))
  assert (approx_equal(stats.mean, 1.35, eps=0.1))
  assert (approx_equal(stats.standard_deviation, 3.25, eps=0.1))
  # test principal_axes_of_inertia
  # XXX d_min=2.01 to get consistent behavior across platforms
  fc = structure.structure_factors(d_min=2.01).f_calc()
  assert (fc.indices().size() == 242)
  fc_map = fc.fft_map(resolution_factor=1/4.)
  fc_map.apply_sigma_scaling()
  real_map = fc_map.real_map_unpadded()
  pai = maptbx.principal_axes_of_inertia(
    real_map=real_map,
    unit_cell=structure.unit_cell(),
    site_cart=(5.,5.,5.),
    radius=2.0)
  assert approx_equal(pai.center_of_mass(), (5.0,5.0,5.0), 0.1)
  assert (approx_equal(pai.inertia_tensor(), (44.89,44.89,44.89,1.87,1.87,1.87),
    eps=0.5))
  assert (approx_equal(list(pai.eigensystem().values()), (48,43,43),
    eps=1))
  # and now with anisotropy
  structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal.symmetry(
        unit_cell=(10.0,10.0,10.0,90,90,90),
        space_group_symbol="P1")),
    scatterers=flex.xray_scatterer([
      xray.scatterer(
        label="O",
        site=(0.5,0.5,0.5),
        u=(0.2,0.2,0.2,0.1,0.0,0.0))]))
  fc = structure.structure_factors(d_min=2.01).f_calc()
  assert (fc.indices().size() == 242)
  fc_map = fc.fft_map(resolution_factor=1/4.)
  fc_map.apply_sigma_scaling()
  fc_map.as_ccp4_map("aniso.ccp4")
  real_map = fc_map.real_map_unpadded()
  pai = maptbx.principal_axes_of_inertia(
    real_map=real_map,
    unit_cell=structure.unit_cell(),
    site_cart=(5.,5.,5.),
    radius=2.0)
  assert (approx_equal(pai.center_of_mass(), (5.0,5.0,5.0), eps=0.2))

def exercise_eight_point_interpolation():
  map = flex.double(flex.grid(2,3,5), 10)
  for shift in [0,1,-1]:
    for index in flex.nested_loop(map.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, map.focus())]
      assert approx_equal(maptbx.eight_point_interpolation(map, x_frac), 10)
      assert approx_equal(
        maptbx.eight_point_interpolation_with_gradients(map, x_frac,[1,1,1])[0], 10)
      assert maptbx.closest_grid_point(map.accessor(), x_frac) == index
  for i in range(100):
    x_frac = [3*random.random()-1 for i in range(3)]
    assert approx_equal(map.eight_point_interpolation(x_frac), 10)
    assert approx_equal(
      map.eight_point_interpolation_with_gradients(x_frac,[1,1,1])[0], 10)
  map = flex.double(range(30))
  map.resize(flex.grid(2,3,5))
  for shift in [0,1,-1]:
    v = 0
    for index in flex.nested_loop(map.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, map.focus())]
      assert approx_equal(map.eight_point_interpolation(x_frac), v)
      assert approx_equal(
        map[maptbx.closest_grid_point(map.accessor(), x_frac)], v)
      assert approx_equal(map.value_at_closest_grid_point(x_frac), v)
      v += 1
  map = flex.double()
  for i in range(48): map.append(i%2)
  map.resize(flex.grid(2,4,6))
  for shift in [0,1,-1]:
    for offs in [.0,.5,.25,.75]:
      v = offs
      for index in flex.nested_loop(map.focus()):
        x_frac = [(i+offs)/n+shift for i,n in zip(index, map.focus())]
        assert approx_equal(map.eight_point_interpolation(x_frac), v)
        if (offs != .5):
          assert maptbx.closest_grid_point(map.accessor(), x_frac) == tuple(
            [int(i+offs+.5)%n for i,n in zip(index,map.focus())])
        v = 1-v

def exercise_real_space_gradients_simple(timing):
  uc = uctbx.unit_cell((11,13,17))
  def check():
    map = flex.double(flex.grid(22,26,36).set_focus(22,26,34))
    site_frac = [i/n for i,n in zip(grid_point, map.focus())]
    sites_cart = flex.vec3_double([uc.orthogonalize(site_frac)])
    target = maptbx.real_space_target_simple(
      unit_cell=uc, density_map=map, sites_cart=sites_cart,
      selection=flex.bool(sites_cart.size(), True))
    assert approx_equal(target, 0)
    terms = maptbx.real_space_target_simple_per_site(
      unit_cell=uc, density_map=map, sites_cart=sites_cart)
    assert approx_equal(terms, [0])
    grads = maptbx.real_space_gradients_simple(
      unit_cell=uc, density_map=map, sites_cart=sites_cart, delta=0.1,
      selection=flex.bool(sites_cart.size(), True))
    assert approx_equal(grads, [(0,0,0)])
    grid_point_mod = [i%n for i,n in zip(grid_point, map.focus())]
    map[grid_point_mod] = 1
    target = maptbx.real_space_target_simple(
      unit_cell=uc, density_map=map, sites_cart=sites_cart,
      selection=flex.bool(sites_cart.size(), True))
    assert approx_equal(target, 1)
    terms = maptbx.real_space_target_simple_per_site(
      unit_cell=uc, density_map=map, sites_cart=sites_cart)
    assert approx_equal(terms, [1])
    grads = maptbx.real_space_gradients_simple(
      unit_cell=uc, density_map=map, sites_cart=sites_cart, delta=0.1,
      selection=flex.bool(sites_cart.size(), True))
    assert approx_equal(grads, [(0,0,0)])
    i,j,k = grid_point_mod
    u,v,w = map.focus()
    map[((i+1)%u,j,k)] = 0.3
    map[(i,(j+1)%v,k)] = 0.5
    map[(i,j,(k+1)%w)] = 0.7
    target = maptbx.real_space_target_simple(
      unit_cell=uc, density_map=map, sites_cart=sites_cart,
      selection=flex.bool(sites_cart.size(), True))
    assert approx_equal(target, 1)
    for delta in [0.1, 0.2]:
      grads = maptbx.real_space_gradients_simple(
        unit_cell=uc, density_map=map, sites_cart=sites_cart, delta=delta,
        selection=flex.bool(sites_cart.size(), True))
      assert approx_equal(grads, [(0.3,0.5,0.7)])
  for grid_point in [(0,0,0), (3,4,5), (-3,15,20)]:
    check()
  for i_trial in range(10):
    grid_point = [random.randrange(-100,100) for i in [0,1,2]]
    check()
  if (timing): n = 1000000
  else:        n = 10
  sites_cart = flex.vec3_double(flex.random_double(size=n*3)*40-20)
  map = flex.double(flex.grid(22,26,36).set_focus(22,26,34), 1)
  target = maptbx.real_space_target_simple(
    unit_cell=uc, density_map=map, sites_cart=sites_cart,
    selection=flex.bool(sites_cart.size(), True))
  assert approx_equal(target, n)
  t0 = time.time()
  maptbx.real_space_gradients_simple(
    unit_cell=uc, density_map=map, sites_cart=sites_cart, delta=0.1,
    selection=flex.bool(sites_cart.size(), True))
  tm = time.time() - t0
  msg = "real_space_gradients_simple: %.2f s / %d sites" % (tm, n)
  if (tm >= 0.01): msg += ", %.0f sites / s" % (n / tm)
  if (timing): print(msg)

test_map  = flex.double([
  -0.069785, -0.109740, -0.172220, -0.209010, -0.255220, -0.285670,
  -0.303130, -0.221400, -0.136640, -0.121530, -0.215260, -0.292640,
  -0.498500, -0.371540, -0.180660, -0.093766, -0.200360, -0.334720,
  -0.356690, -0.330580, -0.249670, -0.204200, -0.264490, -0.320590,
  -0.190610, -0.302730, -0.375040, -0.377540, -0.327030, -0.219010,
  0.060113, -0.023043, -0.185520, -0.311580, -0.395500, -0.435890,
  -0.181560, -0.157460, -0.223360, -0.318560, -0.405230, -0.414470,
  -0.479920, -0.355250, -0.260400, -0.264970, -0.357940, -0.414640,
  -0.362330, -0.292680, -0.244540, -0.294960, -0.416090, -0.486070,
  -0.111790, -0.173040, -0.295180, -0.440520, -0.506410, -0.444720,
  0.077020, 0.038778, -0.054072, -0.178180, -0.323440, -0.434550,
  -0.015310, -0.030401, -0.141110, -0.275890, -0.399010, -0.461390,
  -0.241520, -0.209370, -0.248690, -0.322660, -0.384480, -0.383160,
  -0.209260, -0.185350, -0.221370, -0.317920, -0.396670, -0.400010,
  -0.001700, -0.066179, -0.224020, -0.427730, -0.515120, -0.456980,
  0.015884, 0.042340, 0.061897, -0.020660, -0.193310, -0.338580,
  0.034724, 0.063173, 0.054727, -0.019280, -0.169690, -0.341410,
  0.026087, -0.022744, -0.092896, -0.141160, -0.186520, -0.249360,
  0.052272, -0.015214, -0.111030, -0.175570, -0.180900, -0.163140,
  0.093106, 0.010392, -0.141450, -0.299110, -0.328730, -0.249500])

def exercise_asu_eight_point_interpolation():
  map = flex.double(flex.grid(2,3,5), 10)
  cs = crystal.symmetry(
    unit_cell=(1,1,1,90,90,90),
    space_group="P1")
  asu_mappings=cs.asu_mappings(buffer_thickness=0)
  for shift in [0,1,-1]:
    for index in flex.nested_loop(map.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, map.focus())]
      assert approx_equal(
        maptbx.asu_eight_point_interpolation(map, asu_mappings, x_frac), 10)
  assert approx_equal(
    maptbx.asu_eight_point_interpolation(map, asu_mappings, (10,11,12)), 10)

def exercise_transformers():
  unit_cell=uctbx.unit_cell((10,10,10,90,90,90))
  sites_cart = ((9.5,10.7,3.2),(-.5,1.7,13.2))
  sites_frac = ((1.0,1.9,-.21),(3.6,0.7,-100.7))
  sites_grid = ((5,8,51),(-12,6,105))
  extents = (100,100,100)
  c2f = maptbx.cart2frac(unit_cell.fractionalization_matrix())
  c2g = maptbx.cart2grid(unit_cell.fractionalization_matrix(),extents)
  c2c = maptbx.cart2cart()
  f2f = maptbx.frac2frac()
  f2c = maptbx.frac2cart(unit_cell.orthogonalization_matrix())
  f2g = maptbx.frac2grid(extents)
  g2f = maptbx.grid2frac(extents)
  g2g = maptbx.grid2grid()
  g2c = maptbx.grid2cart(extents,unit_cell.orthogonalization_matrix())
  for site_cart in sites_cart:
    assert approx_equal( c2f(site_cart), unit_cell.fractionalize(site_cart) )
    frac_pt = unit_cell.fractionalize(site_cart)
    grid_pt = (frac_pt[0]*extents[0],frac_pt[1]*extents[1],frac_pt[2]*extents[2])
    assert approx_equal( c2g(site_cart), grid_pt )
    assert approx_equal( c2c(site_cart), site_cart )
  for site_frac in sites_frac:
    assert approx_equal( f2f(site_frac), site_frac )
    assert approx_equal( f2c(site_frac), unit_cell.orthogonalize(site_frac) )
    grid_pt = (site_frac[0]*extents[0],site_frac[1]*extents[1],site_frac[2]*extents[2])
    assert approx_equal( f2g(site_frac), grid_pt )
  for site_grid in sites_grid:
    frac_pt = (site_grid[0]/float(extents[0]),site_grid[1]/float(extents[1]),site_grid[2]/float(extents[2]))
    assert approx_equal( g2f(site_grid), frac_pt )
    assert approx_equal( g2g(site_grid), site_grid )
    assert approx_equal( g2c(site_grid), unit_cell.orthogonalize(frac_pt) )

def exercise_non_crystallographic_eight_point_interpolation():
  unit_cell=130.45,130.245,288.405,90,90,120
  unit_cell_gridding_n=144,144,360
  grid_cell=uctbx.unit_cell((130.45/144,130.245/144,388.405/360,90,90,120))
  grid_mat = grid_cell.fractionalization_matrix()
  map = test_map.deep_copy()
  map.resize(flex.grid((-1,-2,-1),(3,3,5)))
  for site_cart,expected_result in ([(0.468661,-1.549268,3.352108),-0.333095],
                                    [(0.624992,1.553980,1.205578),-0.187556],
                                    [(0.278175,0.968454,2.578265),-0.375068],
                                    [(0.265198,-1.476055,0.704381),-0.147061],
                                    [(1.296042,0.002101,3.459270),-0.304401],
                                    [(0.296189,-1.346603,2.935777),-0.296395],
                                    [(0.551586,-1.284371,3.202145),-0.363263],
                                    [(0.856542,-0.782700,-0.985020),-0.106925],
                                    [(0.154407,1.078936,-0.917551),-0.151128]):
    assert approx_equal(maptbx.non_crystallographic_eight_point_interpolation(
      map=map,
      gridding_matrix=grid_mat,
      site_cart=site_cart,
      allow_out_of_bounds=False,
      out_of_bounds_substitute_value=0), expected_result)
  for x in range(0,2):
    for y in range(-1,2):
      for z in range(0,4):
        assert approx_equal(
          maptbx.non_crystallographic_eight_point_interpolation(
            map,
            grid_mat,
            grid_cell.orthogonalize((x,y,z))), map[x,y,z])
  try:
    val = maptbx.non_crystallographic_eight_point_interpolation(
      map, grid_mat, (5,5,5))
  except RuntimeError as e:
    assert str(e) == \
      "cctbx Error: non_crystallographic_eight_point_interpolation:" \
      " point required for interpolation is out of bounds."
  else: raise Exception_expected
  assert approx_equal(maptbx.non_crystallographic_eight_point_interpolation(
    map, grid_mat, (5,5,5), True, -123), -123)

def exercise_average_density():
  map = test_map.deep_copy()
  map.resize(flex.grid(4,5,6))
  sites_frac = flex.vec3_double([
    (-0.8492400683111605, 0.49159543530354166, 0.55624239788303198),
    (0.10631567870879444, -0.38726326483005269, -0.13581656178827783),
    (0.1895918946688977, -0.25027164520003642, -0.61981792226895172),
    (-0.88980846616667897, -0.79492758628794169, 0.015347308715653485)])
  unit_cell = uctbx.unit_cell((130,130,288,90,90,120))
  densities = maptbx.average_densities(
    unit_cell=unit_cell,
    data=map,
    sites_frac=sites_frac,
    radius=5)
  assert approx_equal(densities, [0,0,0,0])
  densities = maptbx.average_densities(
    unit_cell=unit_cell,
    data=map,
    sites_frac=sites_frac,
    radius=50)
  assert approx_equal(densities, [
    -0.27089094117647061,
    -0.34313799999999994,
    -0.2644232307692308,
    -0.20403226666666666])

def exercise_grid_indices_around_sites():
  unit_cell = uctbx.unit_cell((5,5,5))
  fft_n_real = (5,5,5)
  fft_m_real = (5,5,5)
  site_radii = flex.double([0.5*3**0.5+1e-6])
  def get():
    grid_indices = maptbx.grid_indices_around_sites(
      unit_cell=unit_cell, fft_n_real=fft_n_real, fft_m_real=fft_m_real,
      sites_cart=sites_cart, site_radii=site_radii)
    return list(grid_indices)
  sites_cart = flex.vec3_double([(0.5,0.5,0.5)])
  assert get() == [0, 1, 5, 6, 25, 26, 30, 31]
  sites_cart = flex.vec3_double([(1.5,1.5,1.5)])
  assert get() == [31, 32, 36, 37, 56, 57, 61, 62]
  def sample():
    for i in range(-2,7):
      for j in range(-2,7):
        for k in range(-2,7):
          sites_cart = flex.vec3_double([(i+.5,j+.5,k+.5)])
          assert len(get()) == 8
  sample()
  #
  unit_cell = uctbx.unit_cell((5,6,7))
  fft_n_real = (5,6,7)
  fft_m_real = (5,6,7)
  sites_cart = flex.vec3_double([(0.5,0.5,0.5)])
  assert get() == [0, 1, 7, 8, 42, 43, 49, 50]
  fft_m_real = (5,6,8)
  assert get() == [0, 1, 8, 9, 48, 49, 56, 57]
  fft_m_real = (5,7,8)
  assert get() == [0, 1, 8, 9, 56, 57, 64, 65]
  sample()
  #
  site_radii = flex.double([2])
  assert len(get()) == 8 + 6*4
  site_radii = flex.double([1000])
  assert len(get()) == 5*6*7
  #
  unit_cell = uctbx.unit_cell((18,26,27))
  fft_n_real = (18,26,27)
  fft_m_real = (18,27,28)
  for ish in range(5):
    x = 2*ish+.5
    sites_cart = flex.vec3_double([[x]*3])
    sh = 3**0.5*(ish+0.5)
    site_radii = flex.double([sh-1e-6])
    s1 = set(get())
    site_radii = flex.double([sh+1e-6])
    s2 = set(get())
    for gi in sorted(s2-s1):
      i,j,k = n_dim_index_from_one_dim(gi, fft_m_real)
      assert approx_equal(abs(matrix.col((i-x,j-x,k-x))), sh)
    assert len(s1) == [0, 56, 304, 912, 1904][ish]
    assert len(s2) == [8, 88, 360, 968, 2008][ish]
  #
  unit_cell = uctbx.unit_cell((8,9,7,80,100,110))
  fft_n_real = (11,13,15)
  fft_m_real = (18,26,19)
  sites_cart = flex.vec3_double([(3,11,5)])
  ls = []
  prev = 0
  for r in itertools.count(1):
    site_radii = flex.double([r])
    l = len(get())
    assert l > prev
    ls.append(l)
    if (l == 11*13*15):
      break
    assert r < 7
    prev = l
  assert ls == [18, 155, 524, 1225, 1940, 2139, 2145]
  #
  fft_m_real = (1073741824, 1073741824, 1073741824)
  try:
    maptbx.grid_indices_around_sites(
      unit_cell=unit_cell, fft_n_real=fft_n_real, fft_m_real=fft_m_real,
      sites_cart=sites_cart, site_radii=site_radii)
  except RuntimeError as e:
    assert str(e).startswith("product of fft_m_real")
  else: raise Exception_expected

def exercise_standard_devations_around_sites():
  unit_cell = uctbx.unit_cell((5,5,5))
  fft_n_real = (5,5,5)
  fft_m_real = (5,5,6)
  density_map = flex.double(flex.grid(fft_m_real).set_focus(fft_n_real), 0)
  sites_cart = flex.vec3_double([(2.5,2.5,2.5)])
  site_radii = flex.double([1.2])
  def get():
    return maptbx.standard_deviations_around_sites(
      unit_cell=unit_cell, density_map=density_map,
      sites_cart=sites_cart, site_radii=site_radii)
  assert approx_equal(get(), [0])
  density_map[(2,2,2)] = 1
  assert approx_equal(get(), [0.35355339059327379])

def exercise_region_density_correlation():
  sites_frac = flex.vec3_double([
    (0.02,0.10,0.02),
    (0.10,0.02,0.40),
    (0.98,0.10,0.60),
    (0.10,0.98,0.80),
    (0.20,0.50,0.98)])
  from cctbx import xray
  xray_structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(30,30,50,90,90,120),
      space_group_symbol="P1"),
    scatterers=flex.xray_scatterer([
      xray.scatterer(label=str(i), scattering_type="Si", site=site_frac)
        for i,site_frac in enumerate(sites_frac)]))
  d_min = 2
  f_calc = xray_structure.structure_factors(d_min=d_min).f_calc()
  density_map = f_calc.fft_map().real_map_unpadded()
  def get(region_sel):
    return maptbx.region_density_correlation(
      large_unit_cell=xray_structure.unit_cell(),
      large_d_min=d_min,
      large_density_map=density_map,
      sites_cart=xray_structure.sites_cart().select(region_sel),
      site_radii=flex.double(sites_frac.size(), 1).select(region_sel),
      work_scatterers=xray_structure.scatterers().select(region_sel))
  cc = get(region_sel=flex.bool(5, False))
  assert cc is None
  cc = get(region_sel=flex.bool(5, True))
  assert approx_equal(cc, 1)
  cc = get(region_sel=flex.size_t([4]))
  assert approx_equal(cc, 0.999923364584, eps=1.e-4)
  cc = get(region_sel=flex.size_t([0,2]))
  assert approx_equal(cc, 0.998640554144, eps=1.e-4)
  cc = get(region_sel=flex.size_t([1,3]))
  assert approx_equal(cc, 0.999324555256, eps=1.e-4)
  cc = get(region_sel=flex.size_t([0,4]))
  assert approx_equal(cc, 0.999570252441, eps=1.e-4)
  xray_structure.scatterers()[4].site = (0.205,0.503,0.974)
  cc = get(region_sel=flex.size_t([4]))
  assert approx_equal(cc, 0.6756590336, eps=1.e-3)

def exercise_boxing():
  n_real = (60, 100, 160)
  cs=crystal.symmetry(
    unit_cell=(21,37,58,80,111,117),
    space_group_symbol="P1")
  maptbx.boxes(n_real = n_real, fraction=0.1)

def exercise_hoppe_gassman_modification__and__convert_to_non_negative():
  values = [-2,-1,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,3,4]
  random.choice([0,1,2,3,4,5,6,7,8,9,10,11])
  av = [values[random.choice([0,1,2,3,4,5,6,7])] for i in range(10*20*30)]
  a = flex.double(av)
  a.resize(flex.grid((10,20,30)))
  # inefficient, but transparent way
  a1 = a.deep_copy()
  maptbx.convert_to_non_negative(data=a1, substitute_value=0)
  assert flex.min(a1)>=0
  ave = flex.mean(a1.as_1d().select(a1.as_1d()>0))
  a1.set_selected(a1>ave*2, ave*2)
  a1 = a1/flex.max(a1)
  a1 = 3*a1*a1-2*a1*a1*a1
  assert flex.min(a1)>=0
  assert flex.max(a1)<=1
  # unobvious but efficient way
  a2 = a.deep_copy()
  maptbx.hoppe_gassman_modification(data=a2, mean_scale=2, n_iterations=1)
  assert flex.min(a2)>=0
  assert flex.max(a2)<=1
  #
  assert approx_equal(flex.mean(a1), flex.mean(a2))

def exercise_set_box():
  n_real = (60, 100, 160)
  n = n_real[0]*n_real[1]*n_real[2]
  cs=crystal.symmetry(
    unit_cell=(21,37,58,80,111,117),
    space_group_symbol="P1")
  be = maptbx.boxes(n_real = n_real, fraction=0.1)
  #
  m1 = flex.double([-1 for i in range(n)])
  m1.resize(flex.grid(n_real))
  m2 = flex.double([1 for i in range(n)])
  m2.resize(flex.grid(n_real))
  #
  for s,e in zip(be.starts, be.ends):
    box = maptbx.copy(m2, s, e)
    box.reshape(flex.grid(box.all()))
    maptbx.set_box(
      map_data_from = box,
      map_data_to   = m1,
      start         = s,
      end           = e)
  assert m2.as_1d().min_max_mean().as_tuple() == (1.,1.,1.)


def exercise_set_box_0():
  # Create a grid of size 10x10x10 having value 0 everywhere
  box = flex.double(flex.grid(10,10,10), 0)
  # test 0: same start and end
  b1 = box.deep_copy()
  try:
    maptbx.set_box(
      value       = -1,
      map_data_to = b1,
      start       = b1.all(),
      end         = b1.all())
  except RuntimeError as e:
    assert str(e).endswith("CCTBX_ASSERT(end[i] > start[i]) failure.")
  else: raise Exception_expected
  # test 1: transform entire unit cell
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = -1,
    map_data_to = b1,
    start       = [0,0,0],
    end         = b1.all())
  assert approx_equal(b1.as_1d().min_max_mean().as_tuple(), [-1.0, -1.0, -1.0])
  # test 2: transform entire unit cell, this time translated
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = 1,
    map_data_to = b1,
    start       = [-30,-30,-30],
    end         = [-20,-20,-20])
  assert approx_equal(b1.as_1d().min_max_mean().as_tuple(), [1.0, 1.0, 1.0])
  # test 3: start in neighboring cell and end at 0,0,0
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = 1,
    map_data_to = b1,
    start       = [-5,-5,-5],
    end         = [0,0,0])
  assert approx_equal(b1.as_1d().min_max_mean().as_tuple(), [0.0, 1.0, 0.125])
  # test 4: slice instead of a box
  b1 = box.deep_copy()
  try:
    maptbx.set_box(
      value       = -1,
      map_data_to = b1,
      start       = [1,2,3],
      end         = [2,2,3])
  except RuntimeError as e:
    assert str(e).endswith("CCTBX_ASSERT(end[i] > start[i]) failure.")
  else: raise Exception_expected
  # test 5: another slice
  b1 = box.deep_copy()
  try:
    maptbx.set_box(
      value       = -1,
      map_data_to = b1,
      start       = [-1,0,2],
      end         = [0,0,3])
  except RuntimeError as e:
    assert str(e).endswith("CCTBX_ASSERT(end[i] > start[i]) failure.")
  else: raise Exception_expected
  # test 6: one point changed
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = 1,
    map_data_to = b1,
    start       = [0,0,0],
    end         = [1,1,1])
  assert (b1==0).count(True)==999
  assert (b1==1).count(True)==1
  # test 7: change 1/8 of the unit cell
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = 1,
    map_data_to = b1,
    start       = [0,0,0],
    end         = [5,5,5])
  assert approx_equal(b1.as_1d().min_max_mean().as_tuple(), [0.0, 1.0, 0.125])
  # test 8: change one point
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = 1,
    map_data_to = b1,
    start       = [-1,-1,-1],
    end         = [0,0,0])
  assert (b1==0).count(True)==999
  assert (b1==1).count(True)==1
  # test 9: entire cell, end point is 0,0,0
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = 1,
    map_data_to = b1,
    start       = [-10,-10,-10],
    end         = [0,0,0])
  assert approx_equal(b1.as_1d().min_max_mean().as_tuple(), [1.0, 1.0, 1.0])
  # test 10: slice of a box
  b1 = box.deep_copy()
  try:
    maptbx.set_box(
      value       = 1,
      map_data_to = b1,
      start       = [0,0,0],
      end         = [9,0,0])
  except RuntimeError as e:
    assert str(e).endswith("CCTBX_ASSERT(end[i] > start[i]) failure.")
  else: raise Exception_expected
  # test 11: box within unit cell one period apart
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = 1,
    map_data_to = b1,
    start       = [14,14,14],
    end         = [19,19,19])
  assert approx_equal(b1.as_1d().min_max_mean().as_tuple(), [0.0, 1.0, 0.125])
  # test 12 box between cells, translated by a period
  b1 = box.deep_copy()
  maptbx.set_box(
    value       = 1,
    map_data_to = b1,
    start       = [-14,-14,-14],
    end         = [-9,-9,-9])
  assert approx_equal(b1.as_1d().min_max_mean().as_tuple(), [0.0, 1.0, 0.125])
  # TEST 13: Reset map values in a box within a unit cell
  n_real = (100, 60, 80)
  n = n_real[0]*n_real[1]*n_real[2]
  m1 = flex.double([1 for i in range(n)])
  m1.resize(flex.grid(n_real))
  maptbx.set_box(
    value       = -1,
    map_data_to = m1,
    start       = [20,30,40],
    end         = [80,40,60])
  assert m1.as_1d().min_max_mean().as_tuple() == (-1.,1.,0.95)
  # TEST 14: reset map values in a box crossing the border of the unit cell
  n_real = (60, 100, 80)
  n = n_real[0]*n_real[1]*n_real[2]
  m2 = flex.double([1 for i in range(n)])
  m2.resize(flex.grid(n_real))
  maptbx.set_box(
    value       = -1,
    map_data_to = m2,
    start       = [-10,-20,20],
    end         = [30,40,50])
  assert m2.as_1d().min_max_mean().as_tuple() == (-1.,1.,0.7)

def exercise_median_filter():
  values = [-2,-1,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,3,4]
  av = [values[random.choice([0,1,2,3,4,5,6,7])] for i in range(10*20*30)]
  a = flex.double(av)
  a.resize(flex.grid((10,20,30)))
  maptbx.median_filter(map_data=a, index_span=1)

def exercise_kuwahara_filter():
  values = [-2,-1,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,3,4]
  av = [values[random.choice([0,1,2,3,4,5,6,7])] for i in range(10*20*30)]
  a = flex.double(av)
  a.resize(flex.grid((10,20,30)))
  maptbx.kuwahara_filter(map_data=a, index_span=2)

def exercise_intersection():
  thresholds = flex.double([0,0.1,0.2,0.3,0.4,0.5, 0.6,0.7,0.8,0.8, 1.0])
  def get_map():
    av = [random.random() for i in range(10*20*30)]
    m = flex.double(av)
    m.resize(flex.grid((10,20,30)))
    return m
  m1 = get_map()
  m2 = get_map()
  for average in [True, False]:
    maptbx.intersection(
      map_data_1 = m1,
      map_data_2 = m2,
      thresholds = thresholds,
      average    = average)

def exercise_binarize():
  thresholds = flex.double([0,0.1,0.2,0.3,0.4,0.5, 0.6,0.7,0.8,0.8, 1.0])
  def get_map():
    av = [random.random() for i in range(10*20*30)]
    m = flex.double(av)
    m.resize(flex.grid((10,20,30)))
    return m
  m = get_map()
  maptbx.binarize(map_data=m, threshold=0.5, substitute_value_below=0,
    substitute_value_above=1)
  assert m.as_1d().min_max_mean().as_tuple()[:2] == (0,1)

def exercise_intersection():
  sites_frac = flex.vec3_double([
    (0.02,0.10,0.02),
    (0.10,0.02,0.40),
    (0.98,0.10,0.60),
    (0.10,0.98,0.80),
    (0.20,0.50,0.98)])
  from cctbx import xray
  xray_structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(30,30,50,90,90,120),
      space_group_symbol="P1"),
    scatterers=flex.xray_scatterer([
      xray.scatterer(label=str(i), scattering_type="Si", site=site_frac)
        for i,site_frac in enumerate(sites_frac)]))
  d_min = 0.7
  f_calc = xray_structure.structure_factors(d_min=d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=1/6.)
  fft_map.apply_sigma_scaling()
  density_map = fft_map.real_map_unpadded()
  #
  cm1 = xray_structure.center_of_mass()
  cm2 = maptbx.center_of_mass(map_data=density_map, unit_cell=xray_structure.unit_cell(),
    cutoff=20) #large cutoff to make map look like point scattereres
  assert approx_equal(cm1, cm2, 0.1)

def exercise_map_accumulator(n1=2, n2=2, n3=2):
  ### Python prototype START
  def py_exercise(points, show=False, b=1):
    def smear(x, a, b):
      return math.exp(-(x-a)**2 / (2*b**2))
    def to_int(f):
      p0 = 0
      if(f<=p0): return 0
      return min(int(256*(f-p0)/(1.-p0))+1, 255)
    As = flex.int()
    for i, t_ in enumerate(points):
      a = to_int(t_)
      As.append(a)
      if(show): print("%2d: %8.4f %3d"%(i, t_, a))
    if(show): print(list(As))
    #
    R  = flex.double([0,]*256)
    Rx = flex.int(range(256))
    assert R.size()==Rx.size()
    #
    hit_l = False
    hit_r = False
    for a in As:
      for i in range(-10,11):
        x = a+i
        if(x>=0 and x<=255):
          R[x] += smear(x=x, a=a, b=b)
    #
    if(show):
      for i, rx in enumerate(R):
        print("%4d %10.7f"%(i, rx))
    #
    return R
  ### Python prototype END
  def get_ma(n1,n2,n3, points):
    ma = maptbx.map_accumulator(n_real = (n1,n2,n3), use_max_map=False)
    for value in points:
      m = [value for i in range(n1*n2*n3)]
      m = flex.double(m)
      m.resize(flex.grid((n1,n2,n3)))
      ma.add(map_data=m)
    return ma
  def mmm(m): return m.as_1d().min_max_mean().as_tuple()
  #
  # case 1
  points = flex.double([0,]*16)
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(0,0,0))
  # case 2
  points = flex.double([1,]*16)
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(255,255,255))
  # case 3
  points = flex.double([0.5,]*16)
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(129,129,129))
  # case 4
  points = flex.double([i/16. for i in range(16)])
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(0,0,0))
  # case 5
  points = flex.double([
    0,
    0.49375, 0.4875, 0.48125, 0.475, 0.46875, 0.4625, 0.45625,
    0.5,
    0.50625, 0.5125, 0.51875, 0.525, 0.53125, 0.5375, 0.54375,
    1])
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(128.7,128.7,128.7), 0.01)
  # case 6
  points = flex.double([0,]*8+[1,]*8)
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(0,0,0))
  # case 7
  points = flex.double([0,]*5+[0.5,]*6+[1,]*5)
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(0,0,0))
  # case 8
  points = flex.double([0,]*4+[0.5,]*8+[1,]*4)
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(0,0,0))
  # case 9
  points = flex.double([0,]*3+[0.5,]*10+[1,]*3)
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(129,129,129))
  # case 10
  points = flex.double([
    0,0,
    0.36,0.37,0.38,0.39,0.4,0.41,0.42,
    0.67,0.68,0.69,0.7,0.71,
    1,1])
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(100.17,100.17,100.17), 0.01)
  # case 11
  points = flex.double([0, 0.1,0.11, 0.44, 0.51,0.5101,0.515,0.534,0.54,0.55,
    0.577, 0.78,0.789, 0.77,0.79,0.8, 1])
  ma = get_ma(n1,n2,n3, points)
  # should be 27 or 28 if handles plateaus correctlu (see c++ code for a comment)
  assert approx_equal(mmm(ma.as_median_map()),(134.67,134.67,134.67), 0.01)
  # case 12
  points = flex.double([0.6423, 0.0000, 0.6346, 0.0000, 0.7042, 0.7037, 0.7092,
    0.0067, 0.0000, 0.6796, 0.7073, 0.7900, 0.8083, 0.7582, 0.6609, 0.6851])
  ma = get_ma(n1,n2,n3, points)
  assert approx_equal(mmm(ma.as_median_map()),(0,0,0))
  ###
  # useful for debugging
  #py_exercise(points=points, show=True, b=2)
  #print mmm(ma.as_median_map())
  ###
  # case xxx
  table = """
 0 0.8425 0.5849 0.7631 0.7929 0.7979 0.8319 0.6423 0.8151 0.6436 0.8064
 1 0.7954 0.0000 0.8245 0.8124 0.8245 0.8356 0.0000 0.7857 0.7891 0.8261
 2 0.7798 0.5470 0.8324 0.8453 0.8334 0.8256 0.6346 0.8167 0.5750 0.8150
 3 0.6782 0.7564 0.8275 0.7687 0.8300 0.8496 0.0000 0.6892 0.7873 0.8129
 4 0.7508 0.7976 0.8347 0.8284 0.8325 0.8334 0.7042 0.8344 0.6594 0.8312
 5 0.7192 0.5254 0.8341 0.8376 0.8336 0.8283 0.7037 0.8147 0.0693 0.7969
 6 0.7715 0.5572 0.8238 0.7971 0.8219 0.8312 0.7092 0.8138 0.7656 0.8116
 7 0.7813 0.5117 0.8071 0.8205 0.7839 0.7857 0.0067 0.8208 0.8069 0.7658
 8 0.8228 0.7991 0.8310 0.8057 0.8289 0.8343 0.0000 0.8295 0.7911 0.8072
 9 0.7975 0.6647 0.8411 0.8421 0.8437 0.8383 0.6796 0.8463 0.4164 0.8147
10 0.7570 0.4042 0.8310 0.8263 0.8357 0.8306 0.7073 0.8263 0.7749 0.8061
11 0.7565 0.7787 0.8140 0.8176 0.8302 0.8381 0.7900 0.8109 0.5895 0.8068
12 0.8105 0.3355 0.8273 0.8390 0.7878 0.8159 0.8083 0.8421 0.5090 0.8254
13 0.7748 0.5478 0.8075 0.8365 0.8278 0.8230 0.7582 0.8336 0.4904 0.8077
14 0.8250 0.7528 0.8053 0.8416 0.8251 0.8518 0.6609 0.8097 0.7752 0.7933
15 0.7847 0.4245 0.8193 0.8146 0.8283 0.8435 0.6851 0.8297 0.8027 0.8125
"""
  d = {}
  for l in table.splitlines():
    for i, v in enumerate(l.split()):
      if(i>0): d.setdefault(i, []).append(float(v))
  #
  Rs = []
  results = []
  for points in d.values():
    ma = get_ma(n1,n2,n3, points)
    r = mmm(ma.as_median_map())
    #print r
    results.append(r[0])
    Rs.append(py_exercise(points=points, show=False, b=5))
  assert approx_equal(results,
    (200.70, 140.54, 211.58, 212.27, 212.71, 213.86, 0.0, 211.27, 201.98, 208.10),
    0.1)
  # good to plot frequency distribution
  #for i in range(Rs[0].size()):
  #  print " ".join(["%10.7f"%r[i] for r in Rs])

def exercise_cc_peak():
  def get_map():
    av = [random.random() for i in range(10*20*30)]
    m = flex.double(av)
    m = m-flex.min(m)
    m = m/flex.max(m)
    m.resize(flex.grid((10,20,30)))
    return m
  m1 = get_map()
  m2 = get_map()
  for t in range(0,11):
    t=t/10.
    ccp=maptbx.cc_peak(map_1=m1, map_2=m2, cutoff=t)
  #
  sites_frac = flex.vec3_double([
    (0.50,0.50,0.50)])
  from cctbx import xray
  xray_structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(5,5,5,90,90,90),
      space_group_symbol="P1"),
    scatterers=flex.xray_scatterer([
      xray.scatterer(label=str(i), scattering_type="C", site=site_frac)
        for i,site_frac in enumerate(sites_frac)]))
  fc1 = xray_structure.structure_factors(d_min=1.6).f_calc()
  fc2 = xray_structure.structure_factors(d_min=1.7).f_calc()
  for t in range(0,11):
    t=t/10.
    ccp=maptbx.cc_peak(map_coeffs_1=fc1, map_coeffs_2=fc2, cutoff=t)
  #
  m1_he = maptbx.volume_scale(map = m1,  n_bins = 10000).map_data()
  m2_he = maptbx.volume_scale(map = m2,  n_bins = 10000).map_data()
  cutoffs = flex.double([i/20. for i in range(1,20)])
  df = maptbx.discrepancy_function(map_1=m1_he, map_2=m2_he, cutoffs=cutoffs)
  #
  fc1 = xray_structure.structure_factors(d_min=2.2).f_calc()
  fc2 = xray_structure.structure_factors(d_min=2.2).f_calc()
  for t in range(0,10):
    t=t/10.
    ccp=maptbx.cc_peak(map_coeffs_1=fc1, map_coeffs_2=fc2, cutoff=t)
    assert approx_equal(ccp, 1)
  # 1D case
  m1_he_1d = maptbx.volume_scale_1d(map = m1.as_1d(),  n_bins = 10000).map_data()
  m2_he_1d = maptbx.volume_scale_1d(map = m2.as_1d(),  n_bins = 10000).map_data()
  df_1d = maptbx.discrepancy_function(
    map_1=m1_he_1d, map_2=m2_he_1d, cutoffs=cutoffs)
  assert approx_equal(df, df_1d)

def exercise_gamma_compression():
  def get_map():
    av = [random.random() for i in range(10*20*30)]
    m = flex.double(av)
    m = (m-flex.min(m))*10
    m.resize(flex.grid((10,20,30)))
    return m
  m = get_map()
  maptbx.gamma_compression(map_data=m, gamma=0.5)


def exercise_sample_all_mask_regions():
  cmap = flex.double(flex.grid(30,30,30))
  cmap.fill(1)
  for i in range(0,10):
    for j in range(0,10):
      for k in range(0,10):
        cmap[i,j,k] = 10
  for i in range(15,25):
    for j in range(15,25):
      for k in range(15,25):
        cmap[i,j,k] = 20
  co = maptbx.connectivity(map_data=cmap, threshold=5, wrapping=False)
  uc = uctbx.unit_cell((10,10,10))
  mask_result = co.result()

  sample_regs_obj = maptbx.sample_all_mask_regions(
      mask=mask_result,
      volumes=flex.int([0, 1000,1000]),
      sampling_rates=flex.int([0, 10,10]),
      unit_cell=uc)
  a = sample_regs_obj.get_array(1)
  b = sample_regs_obj.get_array(2)

  assert a.size() == b.size() == 101
  assert approx_equal(a[0], (0,0,0))
  assert approx_equal(b[0], (5,5,5))

def exercise_map_values_along_line_connecting_two_points():
  pdb_str= """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1       1.000   2.000   3.000  1.00 10.00           O
HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O
END
"""
  import iotbx.pdb
  xrs = iotbx.pdb.input(lines = pdb_str, source_info=None).xray_structure_simple()
  fc = xrs.structure_factors(d_min=1).f_calc()
  fft_map = fc.fft_map(resolution_factor=1/4.)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  r = maptbx.map_values_along_line_connecting_two_points(
    map_data=map_data, points_cart=xrs.sites_cart(), step=0.001,
    unit_cell=xrs.unit_cell(), interpolation="tricubic")
  #
  sites_frac = xrs.sites_frac()
  m1 = map_data.tricubic_interpolation(sites_frac[0])
  m2 = map_data.tricubic_interpolation(sites_frac[1])
  #
  assert approx_equal(m1, r.vals[0])
  assert approx_equal(m2, r.vals[-1])

def run(args):
  assert args in [[], ["--timing"]]
  timing = len(args) != 0
  exercise_map_accumulator()
  exercise_gamma_compression()
  exercise_cc_peak()
  exercise_binarize()
  exercise_intersection()
  exercise_boxing()
  exercise_kuwahara_filter()
  exercise_median_filter()
  exercise_set_box()
  exercise_set_box_0()
  exercise_copy()
  exercise_statistics()
  exercise_grid_tags()
  exercise_gridding()
  exercise_misc()
  exercise_peak_search()
  exercise_pymol_interface()
  exercise_structure_factors()
  exercise_fft()
  exercise_transformers()
  exercise_eight_point_interpolation()
  exercise_real_space_gradients_simple(timing=timing)
  exercise_non_crystallographic_eight_point_interpolation()
  exercise_asu_eight_point_interpolation()
  exercise_average_density()
  exercise_grid_indices_around_sites()
  exercise_standard_devations_around_sites()
  exercise_region_density_correlation()
  exercise_hoppe_gassman_modification__and__convert_to_non_negative()
  exercise_sample_all_mask_regions()
  exercise_map_values_along_line_connecting_two_points()
  print("OK")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
