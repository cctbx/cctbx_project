from cctbx import maptbx
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal, not_approx_equal
import sys
import random

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
    for i in xrange(2):
     for j in xrange(3):
       for k in xrange(5):
         m.append(i*100+j*10+k)
    m.resize(flex.grid(2,3,5).set_focus((2,3,4)))
    for i in xrange(-5,5):
     for j in xrange(-5,5):
       for k in xrange(-5,5):
         c = maptbx.copy(map_unit_cell=m, first=(i,j,k), last=(i,j,k))
         assert c.size() == 1
         assert c[(i,j,k)] == m[(i%2,j%3,k%4)]
    c = maptbx.copy(map_unit_cell=m, first=(-1,1,-2), last=(1,2,0))
    assert list(c) == [112, 113, 110, 122, 123, 120,  12,  13,  10,
                        22,  23,  20, 112, 113, 110, 122, 123, 120]

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

def exercise_grid_tags():
  t = maptbx.grid_tags((8,10,12))
  assert not t.is_valid()
  assert t.tag_array().all() == (8,10,12)
  s = sgtbx.space_group_info("P 21")
  for i_flags in xrange(8):
    f = sgtbx.search_symmetry_flags(
      use_space_group_symmetry=i_flags % 2 != 0,
      use_space_group_ltr=0,
      use_seminvariants=(i_flags/4) % 2 != 0,
      use_normalizer_k2l=(i_flags/2) % 2 != 0,
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
        sg.group(), anomalous_flag, mi, d,
        (11,11,9), flex.grid(11,11,9), conjugate_flag)
      assert t.complex_map().focus() == (11,11,9)
      t = maptbx.structure_factors.to_map(
        sg.group(), anomalous_flag, mi, d,
        (11,11,9), flex.grid(11,11,9), conjugate_flag, False)
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
        anomalous_flag, mi, t.complex_map(), conjugate_flag, True)
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

def exercise_eight_point_interpolation():
  map = flex.double(flex.grid(2,3,5), 10)
  for shift in [0,1,-1]:
    for index in flex.nested_loop(map.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, map.focus())]
      assert approx_equal(maptbx.eight_point_interpolation(map, x_frac), 10)
      assert maptbx.closest_grid_point(map.accessor(), x_frac) == index
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
      assert approx_equal(
        map[maptbx.closest_grid_point(map.accessor(), x_frac)], v)
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
          assert maptbx.closest_grid_point(map.accessor(), x_frac) == tuple(
            [int(i+offs+.5)%n for i,n in zip(index,map.focus())])
        v = 1-v

def exercise_real_space_refinement():
  unit_cell=130.45,130.245,288.405,90,90,120
  unit_cell_gridding_n=144,144,360
  grid_cell=uctbx.unit_cell((130.45/144,130.245/144,388.405/360,90,90,120))
  grid_mat = grid_cell.fractionalization_matrix()
  map=flex.double([
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
  map.resize(flex.grid((-1,-2,-1),(3,3,5)))
  sites_cart = flex.vec3_double()
  sites_cart.append((0.468661,-1.549268,3.352108))
  sites_cart.append((0.624992,1.553980,1.205578))
  weights=flex.double(sites_cart.size(),1.0)
  assert approx_equal(maptbx.real_space_refinement_residual(map=map,
                                      gridding_matrix=grid_mat,
                                      sites_cart=sites_cart,
                                      weights=weights),
                      0.260325417539)
  sites_cart = flex.vec3_double()
  sites_cart.append((0.5,0.5,0.5))
  sites_cart.append((0.25,-0.25,0.25))
  sites_cart.append((0.25,0.0,0.5))
  sites_cart.append((0.5,-0.25,0.0))
  sites_cart.append((0,0.25,0.5))
  expected_grads=[(-0.11210999738405786,
                   -0.054035130900231883,
                   0.019751974451977447),
                  (-0.077865472603990016,
                   0.042424083871381954,
                   -0.041516193903453458),
                  (-0.062901734945740237,
                   0.012726058468002122,
                   -0.04725412710938201),
                  (-0.1200936645006829,
                   0.042460230048265907,
                   -0.030121635675126573),
                  (-0.05517232417432183,
                   -0.0043710248947356201,
                   -0.046833897401769242)]
  for grad, correct in zip(maptbx.real_space_refinement_gradients(map=map,
                               gridding_matrix=grid_mat,
                               sites_cart=sites_cart),
                           expected_grads):
    assert approx_equal(grad,correct)

def exercise_non_crystallographic_eight_point_interpolation():
  unit_cell=130.45,130.245,288.405,90,90,120
  unit_cell_gridding_n=144,144,360
  grid_cell=uctbx.unit_cell((130.45/144,130.245/144,388.405/360,90,90,120))
  grid_mat = grid_cell.fractionalization_matrix()
  map=flex.double([
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
  except RuntimeError, e:
    assert str(e) == \
      "cctbx Error: non_crystallographic_eight_point_interpolation:" \
      " point required for interpolation is out of bounds."
  else: raise RuntimeError("Exception expected.")
  assert approx_equal(maptbx.non_crystallographic_eight_point_interpolation(
    map, grid_mat, (5,5,5), True, -123), -123)

def run():
  exercise_copy()
  exercise_statistics()
  exercise_grid_tags()
  exercise_gridding()
  exercise_misc()
  exercise_peak_search()
  exercise_pymol_interface()
  exercise_structure_factors()
  exercise_eight_point_interpolation()
  exercise_non_crystallographic_eight_point_interpolation()
  exercise_real_space_refinement()
  print "OK"

if (__name__ == "__main__"):
  run()
