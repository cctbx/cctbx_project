from cctbx import maptbx
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex
from cctbx import crystal
from libtbx.test_utils import approx_equal, not_approx_equal
import sys
import random
from cPickle import dumps

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

def exercise_real_space_refinement():
  ## test non_symmetric
  unit_cell_parameters=130.45,130.245,388.405,90,90,120
  unit_cell_gridding_n=144,144,360
  unit_cell=uctbx.unit_cell(unit_cell_parameters)
  emap = test_map.deep_copy()
  emap.resize(flex.grid((-1,-2,-1),(3,3,5)))
  out_of_bounds_clamp = maptbx.out_of_bounds_clamp(0)
  basic_map = maptbx.basic_map(maptbx.basic_map_non_symmetric_flag(),
                               emap,
                               unit_cell_gridding_n,
                               unit_cell.orthogonalization_matrix(),
                               out_of_bounds_clamp.as_handle(),
                               unit_cell)

  sites_cart = flex.vec3_double()
  sites_cart.append((0.468661,-1.549268,3.352108))
  sites_cart.append((0.624992,1.553980,1.205578))
  weights=flex.double(sites_cart.size(),1.0)
  assert approx_equal(maptbx.real_space_refinement_residual(
                                      basic_map=basic_map,
                                      sites=sites_cart,
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
  for grad, correct in zip(maptbx.real_space_refinement_gradients(
                               basic_map=basic_map,
                               sites=sites_cart),
                           expected_grads):
    assert approx_equal(grad,correct)
  ## test unit_cell
  emap = flex.double(flex.grid(2,3,5), 10)
  unit_cell = uctbx.unit_cell((1,1,1,90,90,90))
  out_of_bounds_raise = maptbx.out_of_bounds_raise()
  basic_map = maptbx.basic_map(maptbx.basic_map_unit_cell_flag(),
                               emap,
                               emap.focus(),
                               unit_cell.orthogonalization_matrix(),
                               out_of_bounds_raise.as_handle(),
                               unit_cell)
  sites = flex.vec3_double()
  sites.append( (21,-3.4E8,2.6) )
  sites.append( (1.01,2.34,2.184) )
  sites.append( (3,14,15) )
  sites.append( (1,61,8) )
  weights=flex.double(sites.size(),1.0)
  assert approx_equal( maptbx.real_space_refinement_residual(
                                     basic_map=basic_map,
                                     sites=sites,
                                     weights=weights), -10, eps=1.e-4 )
  expected_gradients=flex.double(12,0.0)
  for grad, correct in zip(maptbx.real_space_refinement_gradients(
                                basic_map=basic_map,
                                sites=sites).as_double(),list(expected_gradients)):
    assert approx_equal( grad, correct, eps=1e-4 )
  ## test asu
  cs = crystal.symmetry(
    unit_cell=unit_cell,
    space_group="P1")
  basic_map.as_asu(emap,cs.space_group(),cs.direct_space_asu().as_float_asu(),emap.focus(),0.5,True)
  assert approx_equal( maptbx.real_space_refinement_residual(
                                     basic_map=basic_map,
                                     sites=sites,
                                     weights=weights), -10, eps=1.e-4 )
  for grad, correct in zip(maptbx.real_space_refinement_gradients(
                                basic_map=basic_map,
                                sites=sites).as_double(),list(expected_gradients)):
    assert approx_equal( grad, correct, eps=1e-4 )

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

def exercise_mappers():
  unit_cell=uctbx.unit_cell((10,10,10,90,90,90))
  cs = crystal.symmetry(unit_cell,"P2")
  sites_frac = ((1.0,1.9,-.21),(3.6,0.7,-100.7))
  ns_sites_frac = ((1.0,1.9,-.21),(3.6,0.7,-100.7))
  uc_sites_frac = ((0.0,0.9,0.79),(0.6,0.7,0.3))
  as_sites_frac = ((0,.9,.21),(.6,0.7,.3))

  nsf = maptbx.non_symmetric_factory()
  ucf = maptbx.unit_cell_factory()
  asf = maptbx.asu_factory(cs.space_group(),cs.direct_space_asu().as_float_asu(),0.5,True)

  for idx in xrange(len(sites_frac)):
    assert approx_equal( maptbx.get_non_symmetric_mapper(sites_frac[idx]).mapped_coordinate, ns_sites_frac[idx] )
    assert approx_equal( nsf.map(sites_frac[idx]).mapped_coordinate, ns_sites_frac[idx] )

    assert approx_equal( maptbx.get_unit_cell_mapper(sites_frac[idx]).mapped_coordinate, uc_sites_frac[idx] )
    assert approx_equal( ucf.map(sites_frac[idx]).mapped_coordinate, uc_sites_frac[idx] )

    assert approx_equal( maptbx.get_asu_mapper(sites_frac[idx],
        cs.space_group(),cs.direct_space_asu().as_float_asu(),0.5,True).mapped_coordinate, as_sites_frac[idx] )
    assert approx_equal( asf.map(sites_frac[idx]).mapped_coordinate, as_sites_frac[idx] )

def exercise_basic_map():
  #cs = crystal.symmetry(unit_cell,"P2")
  #### non-symmetric test ####
  unit_cell_parameters=130.45,130.245,388.405,90,90,120
  unit_cell_gridding_n=144,144,360
  unit_cell=uctbx.unit_cell(unit_cell_parameters)
  emap = test_map.deep_copy()
  emap.resize(flex.grid((-1,-2,-1),(3,3,5)))
  out_of_bounds_raise = maptbx.out_of_bounds_raise()
  out_of_bounds_clamp = maptbx.out_of_bounds_clamp(-123)
  basic_map = maptbx.basic_map(maptbx.basic_map_non_symmetric_flag(),
                               emap,
                               unit_cell_gridding_n,
                               unit_cell.orthogonalization_matrix(),
                               out_of_bounds_raise.as_handle(),
                               unit_cell)

  for site_cart,expected_result in ([(0.468661,-1.549268,3.352108),-0.333095],
                                    [(0.624992,1.553980,1.205578),-0.187556],
                                    [(0.278175,0.968454,2.578265),-0.375068],
                                    [(0.265198,-1.476055,0.704381),-0.147061],
                                    [(1.296042,0.002101,3.459270),-0.304401],
                                    [(0.296189,-1.346603,2.935777),-0.296395],
                                    [(0.551586,-1.284371,3.202145),-0.363263],
                                    [(0.856542,-0.782700,-0.985020),-0.106925],
                                    [(0.154407,1.078936,-0.917551),-0.151128]):
    assert approx_equal(basic_map.get_cart_value(site_cart), expected_result)
  for x in range(0,2):
    for y in range(-1,2):
      for z in range(0,4):
        assert approx_equal(
          basic_map.get_grid_value((x,y,z)), emap[x,y,z])
  try:
    val = basic_map.get_cart_value((5,5,5))
  except RuntimeError, e:
    assert str(e) == \
      "cctbx Error: basic_map<T>: the coordinate is out of bounds."
  else: raise RuntimeError("Exception expected.")
  basic_map.set_out_of_bounds_handle(out_of_bounds_clamp.as_handle())
  assert approx_equal(basic_map.get_cart_value((5,5,5)),-123)
  #### unit_cell test ####
  unit_cell = uctbx.unit_cell((1,1,1,90,90,90))
  emap = flex.double(flex.grid(2,3,5), 10)
  basic_map = maptbx.basic_map(maptbx.basic_map_unit_cell_flag(),
                               emap,
                               emap.focus(),
                               unit_cell.orthogonalization_matrix(),
                               out_of_bounds_raise.as_handle(),
                               unit_cell)

  assert basic_map.get_grid_value((0,0,0))==10
  emap[(0,0,0)] = -100
  assert basic_map.get_grid_value((0,0,0))==-100
  emap[(0,0,0)] = 10
  assert basic_map.get_grid_value((0,0,0))==10

  for shift in [0,1,-1]:
    for index in flex.nested_loop(emap.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, emap.focus())]
      assert approx_equal(basic_map.get_frac_value(x_frac), 10)
      assert basic_map.nearest_grid_point_fractional(x_frac) == index
  for i in xrange(100):
    x_frac = [3*random.random()-1 for i in xrange(3)]
    assert approx_equal(basic_map.get_frac_value(x_frac), 10)
  emap = flex.double(range(30))
  emap.resize(flex.grid(2,3,5))
  basic_map.as_unit_cell(emap)
  for shift in [0,1,-1]:
    v = 0
    for index in flex.nested_loop(emap.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, emap.focus())]
      assert approx_equal(basic_map.get_frac_value(x_frac), v)
      assert approx_equal(
        emap[basic_map.nearest_grid_point_fractional(x_frac)], v)
      assert approx_equal(basic_map.value_at_nearest_grid_point_fractional(x_frac), v)
      v += 1
  emap = flex.double()
  for i in xrange(48): emap.append(i%2)
  emap.resize(flex.grid(2,4,6))
  basic_map.as_unit_cell(emap)
  basic_map.rebuild_transformers(emap.focus(),unit_cell.orthogonalization_matrix())
  for shift in [0,1,-1]:
    for offs in [.0,.5,.25,.75]:
      v = offs
      for index in flex.nested_loop(emap.focus()):
        x_frac = [(i+offs)/n+shift for i,n in zip(index, emap.focus())]
        assert approx_equal(basic_map.get_frac_value(x_frac), v)
        if (offs != .5):
          assert basic_map.nearest_grid_point_fractional(x_frac) == tuple(
            [int(i+offs+.5)%n for i,n in zip(index,emap.focus())])
        v = 1-v
  #### asu test ####
  emap = flex.double(flex.grid(2,3,5), 10)
  cs = crystal.symmetry(
    unit_cell=(1,1,1,90,90,90),
    space_group="P1")
  basic_map.as_asu(emap,cs.space_group(),cs.direct_space_asu().as_float_asu(),emap.focus(),0.5,True)
  basic_map.rebuild_transformers(emap.focus(),cs.unit_cell().orthogonalization_matrix())
  for shift in [0,1,-1]:
    for index in flex.nested_loop(emap.focus()):
      x_frac = [float(i)/n+shift for i,n in zip(index, emap.focus())]
      assert approx_equal( basic_map.get_frac_value(x_frac), 10)
  assert approx_equal( basic_map.get_frac_value( (10,11,12)), 10)

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
  except RuntimeError, e:
    assert str(e) == \
      "cctbx Error: non_crystallographic_eight_point_interpolation:" \
      " point required for interpolation is out of bounds."
  else: raise RuntimeError("Exception expected.")
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

def run():
  exercise_copy()
  exercise_statistics()
  exercise_grid_tags()
  exercise_gridding()
  exercise_misc()
  exercise_peak_search()
  exercise_pymol_interface()
  exercise_structure_factors()
  exercise_transformers()
  exercise_mappers()
  exercise_basic_map()
  exercise_eight_point_interpolation()
  exercise_non_crystallographic_eight_point_interpolation()
  exercise_asu_eight_point_interpolation()
  exercise_real_space_refinement()
  exercise_average_density()
  print "OK"

if (__name__ == "__main__"):
  run()
