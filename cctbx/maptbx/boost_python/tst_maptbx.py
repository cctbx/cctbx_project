from cctbx import maptbx
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys
import random

def flex_types():
  return (flex.float, flex.double)

def exercise_copy():
  for flex_type in flex_types():
    m = flex_type((1,2,3,4))
    m.resize(flex.grid(2,2))
    c = maptbx.copy(m, m.accessor())
    assert tuple(m) == tuple(c)
    c = maptbx.copy(m, flex.grid(2,3).set_focus(2,2))
    assert approx_equal(tuple(c), (1,2,0,3,4,0))
    n = maptbx.copy(c, m.accessor())
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

def exercise_statistics():
  for flex_type in flex_types():
    a = flex_type(flex.grid((3,5)))
    s = maptbx.statistics(a)
    assert s.min() == 0
    assert s.max() == 0
    assert s.mean() == 0
    assert s.mean_sq() == 0
    assert s.sigma() == 0
    a = flex_type([random.random() for i in xrange(3*5)])
    a.resize(flex.grid((3,5)))
    s = maptbx.statistics(a)
    assert approx_equal(flex.min(a), s.min())
    assert approx_equal(flex.max(a), s.max())
    assert approx_equal(flex.mean(a), s.mean())
    assert approx_equal(flex.mean_sq(a), s.mean_sq())
    assert approx_equal(flex.mean_sq(a)-flex.mean(a)**2, s.sigma()**2)
    b = flex_type(flex.grid((4,6)).set_focus((3,5)))
    for i in xrange(3):
      for j in xrange(5):
        b[(i,j)] = a[(i,j)]
    b[(3,5)] = -1
    b[(2,5)] = 2
    b.resize(flex.grid((-2,3), (2,9)).set_focus((1,8)))
    t = maptbx.statistics(b)
    assert not approx_equal(flex.min(b), t.min())
    assert not approx_equal(flex.max(b), t.max())
    assert not approx_equal(flex.mean(b), t.mean())
    assert not approx_equal(flex.mean_sq(b), t.mean_sq())
    assert not approx_equal(flex.mean_sq(b)-flex.mean(b)**2, t.sigma()**2)
    assert approx_equal(s.min(), t.min())
    assert approx_equal(s.max(), t.max())
    assert approx_equal(s.mean(), t.mean())
    assert approx_equal(s.mean_sq(), t.mean_sq())
    assert approx_equal(s.sigma(), t.sigma())

def exercise_symmetry_flags():
  f = maptbx.symmetry_flags(0001)
  f = maptbx.symmetry_flags(0001, 0001, 0001)
  assert f.use_space_group_symmetry()
  assert f.use_normalizer_k2l()
  assert f.use_structure_seminvariants()
  sg_info = sgtbx.space_group_info("P 3 1 2")
  assert f.select_sub_space_group(sg_info.type()).type().lookup_symbol() \
      == "P -3 1 m"
  assert f == f
  assert not f != f
  assert f == maptbx.symmetry_flags(0001, 0001, 0001)
  assert not f != maptbx.symmetry_flags(0001, 0001, 0001)
  assert f != maptbx.symmetry_flags(0001)
  assert not f == maptbx.symmetry_flags(0001)

def exercise_grid_tags():
  t = maptbx.grid_tags((8,10,12))
  assert not t.is_valid()
  assert t.tag_array().all() == (8,10,12)
  s = sgtbx.space_group_info("P 21")
  for i_flags in xrange(8):
    f = maptbx.symmetry_flags(i_flags % 2 != 0,
                              (i_flags/2) % 2 != 0,
                              (i_flags/4) % 2 != 0)
    t.build(s.type(), f)
    assert t.is_valid()
    assert t.space_group_type().group() == s.group()
    assert t.symmetry_flags() == f
    if (f.use_structure_seminvariants()):
      assert [(vm.v, vm.m) for vm in t.grid_ss()] \
          == [((1, 0, 0), 2), ((0, 1, 0), 10), ((0, 0, 1), 2)]
    assert t.n_grid_misses() == 0
    assert t.n_independent() == (960, 480, 484, 242, 24, 14, 14, 14)[i_flags]
    assert t.n_independent() + t.n_dependent() == t.tag_array().size()
    for flex_type in flex_types():
      d = flex_type(t.tag_array().accessor())
      assert t.n_dependent() == 0 \
          or t.dependent_correlation(d, 1.e-10).coefficient() > 0.99
      assert t.verify(d, 0.999)
      t.sum_sym_equiv_points(d)
    if (i_flags == 0):
      assert t.n_independent() == t.tag_array().size()
    else:
      assert t.n_independent() < t.tag_array().size()
      for flex_type in flex_types():
        d = flex_type([random.random() for x in xrange(t.tag_array().size())])
        d.resize(t.tag_array().accessor())
        assert not t.verify(d)
        t.sum_sym_equiv_points(d)
        assert t.verify(d)

def exercise_peak_search():
  t = flex.long(flex.grid((3,4,5)))
  for flex_type in flex_types():
    d = flex_type(flex.grid((3,4,5)))
    l = maptbx.peak_list(d, t, peak_search_level=0, interpolate=00000)
    assert l.gridding() == d.focus()
    assert l.grid_indices(0) == (0,0,0)
    assert list(l.grid_heights()) == [0]
    assert list(l.sites()) == [(0,0,0)]
    assert list(l.heights()) == [0]
    l = maptbx.peak_list(
      d, t, peak_search_level=0,peak_cutoff=-1, interpolate=00000)
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
  for anomalous_flag in (00000, 0001):
    for conjugate_flag in (00000, 0001):
      t = maptbx.structure_factors.to_map(
        sg.group(), anomalous_flag, mi, d,
        (11,11,9), flex.grid(11,11,9), conjugate_flag)
      assert t.complex_map().focus() == (11,11,9)
      f = maptbx.structure_factors.from_map(
        uc, sg.type(), anomalous_flag, 5., t.complex_map(), conjugate_flag)
      assert f.miller_indices().size() > 0
      assert f.miller_indices().size() == f.data().size()
      f = maptbx.structure_factors.from_map(
        anomalous_flag, mi, t.complex_map(), conjugate_flag)
      assert f.miller_indices().size() == 0
      assert f.data().size() == mi.size()
      f = maptbx.structure_factors.from_map(
        anomalous_flag, mi, t.complex_map(), conjugate_flag, 0001)
      assert f.miller_indices().size() == 0
      assert f.data().size() == mi.size()
      assert f.n_indices_affected_by_aliasing() == 0
      assert f.outside_map().size() == 0
      f = maptbx.structure_factors.from_map(
        sg.group(), anomalous_flag, mi, t.complex_map(), conjugate_flag)
      assert f.miller_indices().size() == 0
      assert f.data().size() == mi.size()
      assert f.n_indices_affected_by_aliasing() == 0
      assert f.outside_map().size() == 0

def exercise_gridding():
  u = uctbx.unit_cell((4,6,7))
  assert maptbx.ext.determine_gridding(u, 2, 1/3., (1,1,1), 5, 0001) \
      == (8,9,12)
  f = maptbx.symmetry_flags(0001, 0001)
  t = sgtbx.space_group_info("F 2 2 2").primitive_setting().type()
  assert maptbx.ext.determine_gridding(u, 2, 1/3., f, t, (1,1,1), 5, 0001) \
      == (12, 12, 12)

def exercise_misc():
  for flex_type in flex_types():
    m = flex_type([1,2,-3,4,-5,6])
    maptbx.set_if_less_than(m, 0, 0)
    assert approx_equal(tuple(m), (1,2,0,4,0,6))
    maptbx.set_if_less_than(m, 2, 9)
    assert approx_equal(tuple(m), (9,2,9,4,9,6))

def exercise_eight_point_interpolation():
  map = flex.double(flex.grid(2,3,5), 10)
  for shift in [0,1,-1]:
    for index in flex.nested_loop(map.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, map.focus())]
      assert approx_equal(maptbx.eight_point_interpolation(map, x_frac), 10)
      assert maptbx.closest_grid_point(map, x_frac) == index
  for i in xrange(100):
    x_frac = [3*random.random()-1 for i in xrange(3)]
    assert approx_equal(map.eight_point_interpolation(x_frac), 10)
  map = flex.double(range(30))
  map.resize(flex.grid(2,3,5))
  for shift in [0,1,-1]:
    v = 0
    for index in flex.nested_loop(map.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, map.focus())]
      assert approx_equal(map.eight_point_interpolation(x_frac), v)
      assert approx_equal(map[map.closest_grid_point(x_frac)], v)
      assert approx_equal(map.value_at_closest_grid_point(x_frac), v)
      v += 1
  map = flex.double()
  for i in xrange(48): map.append(i%2)
  map.resize(flex.grid(2,4,6))
  for shift in [0,1,-1]:
    for offs in [.0,.5,.25,.75]:
      v = offs
      for index in flex.nested_loop(map.focus()):
        x_frac = [(i+offs)/n+shift for i,n in zip(index, map.focus())]
        assert approx_equal(map.eight_point_interpolation(x_frac), v)
        if (offs != .5):
          assert map.closest_grid_point(x_frac) == tuple(
            [int(i+offs+.5)%n for i,n in zip(index,map.focus())])
        v = 1-v

def run():
  exercise_copy()
  exercise_statistics()
  exercise_symmetry_flags()
  exercise_grid_tags()
  exercise_gridding()
  exercise_misc()
  exercise_peak_search()
  exercise_pymol_interface()
  exercise_structure_factors()
  exercise_eight_point_interpolation()
  print "OK"

if (__name__ == "__main__"):
  run()
