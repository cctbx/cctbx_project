from cctbx import uctbx
from cctbx import sgtbx
from cctbx import miller
from cctbx.array_family import flex
import scitbx.math
from libtbx.complex_math import polar
from libtbx.test_utils import Exception_expected, approx_equal
import pickle
import random
import math

def exercise_sym_equiv():
  s = sgtbx.space_group("P 31")
  e = miller.sym_equiv_indices(s, (0,0,0))
  assert len(e.indices()) == 1
  assert e.is_centric()
  h = (3,5,2)
  e = miller.sym_equiv_indices(s, h)
  i = e.indices()
  assert len(i) == 3
  assert i[0].h() == h
  assert i[0].hr() == h
  assert i[0].ht() == 0
  assert i[0].t_den() == s.t_den()
  assert i[0].ht_angle() == 0
  assert i[0].ht_angle(False) == 0
  assert i[0].ht_angle(True) == 0
  assert not i[0].friedel_flag()
  m = i[0].mate()
  assert m.h() == tuple([-x for x in h])
  assert m.friedel_flag()
  m = i[0].mate(1)
  assert m.h() == tuple([-x for x in h])
  assert m.friedel_flag()
  m = i[0].mate(0)
  assert m.h() == h
  assert not m.friedel_flag()
  for i_mate in (0,1):
    m = i[1].mate(i_mate)
    assert m.ht() == 8
    assert m.friedel_flag() == (i_mate != 0)
    assert approx_equal(m.phase_in(m.phase_eq(30)), 30)
    assert approx_equal(m.phase_in(m.phase_eq(30, False), False), 30)
    assert approx_equal(m.phase_in(m.phase_eq(30, True), True), 30)
    assert approx_equal(m.phase_eq(30*math.pi/180),
                        m.phase_eq(30, True)*math.pi/180)
    c = m.complex_in(m.complex_eq(1+2j))
    assert approx_equal(c.real, 1)
    assert approx_equal(c.imag, 2)
    h = m.hendrickson_lattman_in(m.hendrickson_lattman_eq((1,2,3,4)))
    for j in xrange(4):
      assert approx_equal(h[j], j+1)
  r = e.phase_restriction()
  assert not r.sys_abs_was_tested()
  assert r.ht() < 0
  assert not e.is_centric()
  assert e.multiplicity(False) == 6
  assert e.multiplicity(True) == 3
  assert e.f_mates(False) == 2
  assert e.f_mates(True) == 1
  assert e.epsilon() == 1
  for j in xrange(e.multiplicity(False)):
    if (j < e.multiplicity(True)):
      assert e(j).h() == e.indices()[j].h()
    else:
      assert e(j).friedel_flag()
  assert e.is_valid_phase(10)
  assert e.is_valid_phase(10, False)
  assert e.is_valid_phase(10, True)
  assert e.is_valid_phase(10, True, 1.e-5)
  for anomalous_flag in (False,True):
    j = e.p1_listing(anomalous_flag)
    assert len(j) == len(i)
  s = sgtbx.space_group("P 41")
  h = (3,5,0)
  e = miller.sym_equiv_indices(s, h)
  i = e.indices()
  assert e.is_centric()
  assert e.multiplicity(False) == 4
  assert e.multiplicity(True) == 4
  assert e.f_mates(False) == 1
  assert e.f_mates(True) == 1
  r = e.phase_restriction()
  assert r.ht() == 0
  j = e.p1_listing(False)
  assert len(j) == len(i)//2

def exercise_map_to_asu(sg_symbol):
  sg_type = sgtbx.space_group_type(sg_symbol)
  index_abs_range = (4,4,4)
  for anomalous_flag in (False,True):
    m = miller.index_generator(
      sg_type, anomalous_flag, index_abs_range).to_array()
    a = flex.double()
    p = flex.double()
    c = flex.hendrickson_lattman()
    for i in xrange(m.size()):
      a.append(random.random())
      p.append(random.random() * 2)
      c.append([random.random() for j in xrange(4)])
    f = flex.polar(a, p)
    p = [p, p*(180/math.pi)]
  m_random = flex.miller_index()
  p_random = [flex.double(), flex.double()]
  c_random = flex.hendrickson_lattman()
  f_random = flex.complex_double()
  for i,h_asym in enumerate(m):
    h_eq = miller.sym_equiv_indices(sg_type.group(), h_asym)
    i_eq = random.randrange(h_eq.multiplicity(anomalous_flag))
    h_i = h_eq(i_eq)
    m_random.append(h_i.h())
    for deg in (False,True):
      p_random[deg].append(h_i.phase_eq(p[deg][i], deg))
    f_random.append(h_i.complex_eq(f[i]))
    c_random.append(h_i.hendrickson_lattman_eq(c[i]))
  m_random_copy = m_random.deep_copy()
  miller.map_to_asu(sg_type, anomalous_flag, m_random_copy)
  for i,h_asym in enumerate(m):
    assert h_asym == m_random_copy[i]
  m_random_copy = m_random.deep_copy()
  miller.map_to_asu(sg_type, anomalous_flag, m_random_copy, f_random)
  for i,h_asym in enumerate(m):
    assert h_asym == m_random_copy[i]
  for i,f_asym in enumerate(f):
    assert abs(f_asym - f_random[i]) < 1.e-6
  m_random_copy = m_random.deep_copy()
  a_random = a.deep_copy()
  miller.map_to_asu(sg_type, anomalous_flag, m_random_copy, a_random)
  for i,h_asym in enumerate(m):
    assert h_asym == m_random_copy[i]
  for i,a_asym in enumerate(a):
    assert a_asym == a_random[i]
  for deg in (False,True):
    m_random_copy = m_random.deep_copy()
    miller.map_to_asu(
      sg_type, anomalous_flag, m_random_copy, p_random[deg], deg)
    for i,h_asym in enumerate(m):
      assert h_asym == m_random_copy[i]
    for i,p_asym in enumerate(p[deg]):
      assert scitbx.math.phase_error(p_asym, p_random[deg][i], deg) < 1.e-5
  m_random_copy = m_random.deep_copy()
  miller.map_to_asu(sg_type, anomalous_flag, m_random_copy, c_random)
  for i,h_asym in enumerate(m):
    assert h_asym == m_random_copy[i]
  for i,c_asym in enumerate(c):
    for j in xrange(4):
      assert abs(c_asym[j] - c_random[i][j]) < 1.e-5

def exercise_asu():
  sg_type = sgtbx.space_group_type("P 41")
  asu = sgtbx.reciprocal_space_asu(sg_type)
  miller_indices = flex.miller_index(((1,2,3), (3,5,0)))
  for h in miller_indices:
    h_eq = miller.sym_equiv_indices(sg_type.group(), h)
    for i_eq in xrange(h_eq.multiplicity(False)):
      h_i = h_eq(i_eq)
      for anomalous_flag in (False,True):
        a = miller.asym_index(sg_type.group(), asu, h_i.h())
        assert a.h() == h
        o = a.one_column(anomalous_flag)
        assert o.i_column() == 0
        t = a.two_column(anomalous_flag)
        assert t.h() == h
        assert (o.h() != h) == (t.i_column() == 1)
        assert not anomalous_flag or (t.i_column() != 0) == h_i.friedel_flag()
        assert anomalous_flag or t.i_column() == 0
  miller.map_to_asu(sg_type, False, miller_indices)
  data = flex.double((0,0))
  miller.map_to_asu(sg_type, False, miller_indices, data)
  miller.map_to_asu(sg_type, False, miller_indices, data, False)
  miller.map_to_asu(sg_type, False, miller_indices, data, True)
  data = flex.complex_double((0,0))
  miller.map_to_asu(sg_type, False, miller_indices, data)
  data = flex.hendrickson_lattman(((1,2,3,4),(2,3,4,5)))
  miller.map_to_asu(sg_type, False, miller_indices, data)
  for sg_symbol in ("P 41", "P 31 1 2"):
    exercise_map_to_asu(sg_symbol)
  #
  sg_type = sgtbx.space_group_type("P 2")
  miller_indices = flex.miller_index(((1,2,3), (-1,-2,-3)))
  assert not miller.is_unique_set_under_symmetry(
    space_group_type=sg_type,
    anomalous_flag=False,
    miller_indices=miller_indices)
  assert miller.is_unique_set_under_symmetry(
    space_group_type=sg_type,
    anomalous_flag=True,
    miller_indices=miller_indices)
  assert list(miller.unique_under_symmetry_selection(
    space_group_type=sg_type,
    anomalous_flag=False,
    miller_indices=miller_indices)) == [0]
  assert list(miller.unique_under_symmetry_selection(
    space_group_type=sg_type,
    anomalous_flag=True,
    miller_indices=miller_indices)) == [0,1]

def exercise_bins():
  uc = uctbx.unit_cell((11,11,13,90,90,120))
  sg_type = sgtbx.space_group_type("P 3 2 1")
  anomalous_flag = False
  d_min = 1
  m = miller.index_generator(uc, sg_type, anomalous_flag, d_min).to_array()
  f = flex.double()
  for i in xrange(m.size()): f.append(random.random())
  n_bins = 10
  b = miller.binning(uc, n_bins, 0, d_min)
  b = miller.binning(uc, n_bins, 0, d_min, 1.e-6)
  b = miller.binning(uc, n_bins, m)
  b = miller.binning(uc, n_bins, m, 0)
  b = miller.binning(uc, n_bins, m, 0, d_min)
  b = miller.binning(uc, n_bins, m, 0, d_min, 1.e-6)
  assert b.d_max() == -1
  assert approx_equal(b.d_min(), d_min)
  assert b.bin_d_range(0) == (-1,-1)
  assert approx_equal(b.bin_d_range(1), (-1,2.1544336))
  assert approx_equal(b.bin_d_range(b.n_bins_all()-1), (1,-1))
  d_star_sq = 0.5
  r = b.bin_d_range(b.get_i_bin(d_star_sq))
  d = 1/math.sqrt(d_star_sq)
  assert r[1] <= d <= r[0]
  h = (3,4,5)
  r = b.bin_d_range(b.get_i_bin(h))
  assert r[1] <= uc.d(h) <= r[0]
  # a quick test to excercise d-spacings on fractional Miller indices:
  assert approx_equal( uc.d((3,4,5)), uc.d_frac((3.001,4,5)), eps=0.001)
  binning1 = miller.binning(uc, n_bins, m)
  assert binning1.unit_cell().is_similar_to(uc)
  assert binning1.n_bins_used() == n_bins
  assert binning1.limits().size() == n_bins + 1
  assert binning1.n_bins_all() == n_bins + 2
  s = pickle.dumps(binning1)
  l = pickle.loads(s)
  assert str(l.unit_cell()) == "(11, 11, 13, 90, 90, 120)"
  assert approx_equal(l.limits(), binning1.limits())
  #
  binner1 = miller.ext.binner(binning1, m)
  assert binner1.miller_indices().id() == m.id()
  assert binner1.count(binner1.i_bin_d_too_large()) == 0
  assert binner1.count(binner1.i_bin_d_too_small()) == 0
  counts = binner1.counts()
  for i_bin in binner1.range_all():
    assert binner1.count(i_bin) == counts[i_bin]
    assert binner1.selection(i_bin).count(True) == counts[i_bin]
  assert list(binner1.range_all()) == range(binner1.n_bins_all())
  assert list(binner1.range_used()) == range(1, binner1.n_bins_used()+1)
  binning2 = miller.binning(uc, n_bins - 2,
    binning1.bin_d_min(2),
    binning1.bin_d_min(n_bins))
  binner2 = miller.ext.binner(binning2, m)
  assert tuple(binner1.counts())[1:-1] == tuple(binner2.counts())
  array_indices = flex.size_t(range(m.size()))
  perm_array_indices1 = flex.size_t()
  perm_array_indices2 = flex.size_t()
  for i_bin in binner1.range_all():
    perm_array_indices1.extend(array_indices.select(binner1.selection(i_bin)))
    perm_array_indices2.extend(binner1.array_indices(i_bin))
  assert perm_array_indices1.size() == m.size()
  assert perm_array_indices2.size() == m.size()
  assert tuple(perm_array_indices1) == tuple(perm_array_indices2)
  b = miller.ext.binner(miller.binning(uc, n_bins, m, 0, d_min), m)
  assert approx_equal(b.bin_centers(1),
    (0.23207956, 0.52448148, 0.62711856, 0.70311998, 0.7652538,
     0.818567, 0.86566877, 0.90811134, 0.94690405, 0.98274518))
  assert approx_equal(b.bin_centers(2),
    (0.10772184, 0.27871961, 0.39506823, 0.49551249, 0.58642261,
     0.67067026, 0.74987684, 0.82507452, 0.89697271, 0.96608584))
  assert approx_equal(b.bin_centers(3),
    (0.050000075, 0.15000023, 0.25000038, 0.35000053, 0.45000068,
     0.55000083, 0.65000098, 0.75000113, 0.85000128, 0.95000143))
  v = flex.double(xrange(b.n_bins_used()))
  i = b.interpolate(v, 0)
  for i_bin in b.range_used():
    assert i.select(b.selection(i_bin)).all_eq(v[i_bin-1])
  dss = uc.d_star_sq(m)
  for d_star_power in (1,2,3):
    j = b.interpolate(v, d_star_power)
    x = flex.pow(dss, (d_star_power/2.))
    r = flex.linear_correlation(x, j)
    assert r.is_well_defined()
    assert approx_equal(
      r.coefficient(), (0.946401,0.990764,1.0)[d_star_power-1],
      eps=1.e-4, multiplier=None)
  #
  s = pickle.dumps(binner2)
  l = pickle.loads(s)
  assert str(l.unit_cell()) == "(11, 11, 13, 90, 90, 120)"
  assert approx_equal(l.limits(), binner2.limits())
  assert l.miller_indices().all_eq(binner2.miller_indices())
  assert l.bin_indices().all_eq(binner2.bin_indices())
  #
  limits = flex.random_double(size=10)
  bng = miller.binning(uc, limits)
  assert bng.unit_cell().is_similar_to(uc)
  assert approx_equal(bng.limits(), limits)

def exercise_expand():
  sg = sgtbx.space_group("P 41 (1,-1,0)")
  h = flex.miller_index(((3,1,-2), (1,-2,0)))
  assert tuple(sg.is_centric(h)) == (0, 1)
  p1 = miller.expand_to_p1_iselection(
    space_group=sg, anomalous_flag=False, indices=h, build_iselection=False)
  p1_i0 = ((-3,-1,2), (-1, 3,2),(3,1,2),(1,-3,2),(1,-2, 0),(2,1,0))
  assert tuple(p1.indices) == p1_i0
  assert p1.iselection.size() == 0
  p1 = miller.expand_to_p1_iselection(
    space_group=sg, anomalous_flag=True, indices=h, build_iselection=False)
  assert tuple(p1.indices) \
      == ((3,1,-2), (1,-3,-2), (-3,-1,-2), (-1,3,-2),
          (1,-2,0), (-2,-1,0), (-1,2,0), (2,1,0))
  p1 = miller.expand_to_p1_iselection(
    space_group=sg, anomalous_flag=False, indices=h, build_iselection=True)
  assert tuple(p1.indices) == p1_i0
  assert tuple(p1.iselection) == (0,0,0,0,1,1)
  a = flex.double((1,2))
  p = flex.double((10,90))
  p1 = miller.expand_to_p1_phases(
    space_group=sg, anomalous_flag=False, indices=h, data=p, deg=True)
  assert approx_equal(tuple(p1.data), (-10,110,110,-10, 90,30))
  p1 = miller.expand_to_p1_phases(
    space_group=sg, anomalous_flag=True, indices=h, data=p, deg=True)
  assert approx_equal(tuple(p1.data), (10,-110,-110,10, 90,-30,-90,30))
  p = flex.double([x * math.pi/180 for x in p])
  v = [x * math.pi/180 for x in p1.data]
  p1 = miller.expand_to_p1_phases(
    space_group=sg, anomalous_flag=True, indices=h, data=p, deg=False)
  assert approx_equal(tuple(p1.data), v)
  f = flex.polar(a, p)
  p1 = miller.expand_to_p1_complex(
    space_group=sg, anomalous_flag=True, indices=h, data=f)
  assert approx_equal(tuple(flex.abs(p1.data)), (1,1,1,1,2,2,2,2))
  assert approx_equal(tuple(flex.arg(p1.data)), v)
  hl = flex.hendrickson_lattman([(1,2,3,4), (5,6,7,8)])
  p1 = miller.expand_to_p1_hendrickson_lattman(
    space_group=sg, anomalous_flag=True, indices=h, data=hl)
  assert approx_equal(p1.data, [
    [1,2,3,4],
    [1.232051,-1.866025,-4.964102,0.5980762],
    [1.232051,-1.866025,-4.964102,0.5980762],
    [1,2,3,4],
    [5,6,7,8],
    [2.696152,-7.330127,-10.4282,2.062178],
    [-5,-6,7,8],
    [7.696152,-1.330127,3.428203,-10.06218]])
  b = flex.bool([True,False])
  p1 = miller.expand_to_p1_iselection(
    space_group=sg, anomalous_flag=True, indices=h, build_iselection=True)
  assert b.select(p1.iselection).all_eq(
    flex.bool([True, True, True, True, False, False, False, False]))
  i = flex.int([13,17])
  p1 = miller.expand_to_p1_iselection(
    space_group=sg, anomalous_flag=True, indices=h, build_iselection=True)
  assert i.select(p1.iselection).all_eq(flex.int([13,13,13,13,17,17,17,17]))
  #
  assert approx_equal(miller.statistical_mean(sg, False, h, a), 4/3.)
  assert approx_equal(miller.statistical_mean(sg, True, h, a), 3/2.)

def exercise_index_generator():
  uc = uctbx.unit_cell((11,11,13,90,90,120))
  sg_type = sgtbx.space_group_type("P 3 1 2")
  for anomalous_flag in (False,True):
    mig = miller.index_generator(uc, sg_type, anomalous_flag, 8)
    assert mig.unit_cell().is_similar_to(uc)
    assert mig.space_group_type().group() == sg_type.group()
    assert mig.anomalous_flag() == anomalous_flag
    assert mig.asu().reference_as_string() == "h>=k and k>=0 and (k>0 or l>=0)"
    assert mig.next() == (0,0,1)
    if (not anomalous_flag):
      assert tuple(mig.to_array()) == ((1, 0, 0),)
    else:
      assert tuple(mig.to_array()) == ((1, 0, 0), (-1, 0, 0))
  assert tuple(miller.index_generator(uc, sg_type, False, 8)) \
         == ((0,0,1), (1, 0, 0))
  index_abs_range = (4,4,4)
  for sg_symbol in ("P 31 1 2", "P 31 2 1"):
    sg_type = sgtbx.space_group_type(sg_symbol)
    for anomalous_flag in (False,True):
      miller_indices = miller.index_generator(
        sg_type, anomalous_flag, index_abs_range)
      miller_dict = {}
      for h in miller_indices: miller_dict[h] = 0
      sg = sg_type.group()
      h = [0,0,0]
      for h[0] in range(-index_abs_range[0], index_abs_range[0]+1):
        for h[1] in range(-index_abs_range[1], index_abs_range[1]+1):
          for h[2] in range(-index_abs_range[2], index_abs_range[2]+1):
            if (sg.is_sys_absent(h) or h == [0,0,0]): continue
            h_eq = miller.sym_equiv_indices(sg, h)
            found_h_asu = 0
            for i_eq in xrange(h_eq.multiplicity(anomalous_flag)):
              h_i = h_eq(i_eq).h()
              if (h_i in miller_dict):
                assert found_h_asu == 0
                found_h_asu = 1
            assert found_h_asu != 0

def exercise_index_span():
  miller_indices = flex.miller_index(((1,-2,3), (-3,5,0)))
  s = miller.index_span(miller_indices)
  assert s.min() == (-3,-2,0)
  assert s.max() == (1,5,3)
  assert s.abs_range() == (4,6,4)
  assert s.map_grid() == (7,11,7)
  assert s.is_in_domain((-1,2,1))
  assert not s.is_in_domain((0,6,0))
  assert tuple(s.pack(miller_indices)) == (131, 28)

def exercise_match_bijvoet_mates():
  h0 = flex.miller_index(((1,2,3), (-1,-2,-3), (2,3,4), (-2,-3,-4), (3,4,5)))
  d0 = flex.double((1,2,3,4,5))
  bm = miller.match_bijvoet_mates(
    sgtbx.space_group_type(),
    h0)
  bm = miller.match_bijvoet_mates(
    sgtbx.reciprocal_space_asu(sgtbx.space_group_type()),
    h0)
  bm = miller.match_bijvoet_mates(
    h0)
  assert tuple(bm.pairs()) == ((0,1), (2,3))
  assert tuple(bm.singles("+")) == (4,)
  assert tuple(bm.singles("-")) == ()
  assert bm.n_singles() != 0
  assert tuple(bm.pairs_hemisphere_selection("+")) == (0, 2)
  assert tuple(bm.pairs_hemisphere_selection("-")) == (1, 3)
  assert tuple(bm.singles_hemisphere_selection("+")) == (4,)
  assert tuple(bm.singles_hemisphere_selection("-")) == ()
  assert tuple(bm.miller_indices_in_hemisphere("+")) == ((1,2,3), (2,3,4))
  assert tuple(bm.miller_indices_in_hemisphere("-")) == ((-1,-2,-3),(-2,-3,-4))
  assert approx_equal(tuple(bm.minus(d0)), (-1, -1))
  assert approx_equal(tuple(bm.additive_sigmas(d0)),
                      [math.sqrt(x*x+y*y) for x,y in ((1,2), (3,4))])
  assert approx_equal(tuple(bm.average(d0)), (3/2., 7/2.))
  h0.append((1,2,3))
  try: miller.match_bijvoet_mates(h0)
  except Exception: pass
  else: raise Exception_expected

def exercise_match_indices():
  h0 = flex.miller_index(((1,2,3), (-1,-2,-3), (2,3,4), (-2,-3,-4), (3,4,5)))
  d0 = flex.double((1,2,3,4,5))
  h1 = flex.miller_index(((-1,-2,-3), (-2,-3,-4), (1,2,3), (2,3,4)))
  d1 = flex.double((10,20,30,40))
  mi = miller.match_indices(h0, h0)
  assert mi.have_singles() == 0
  assert list(mi.pairs()) == zip(range(5), range(5))
  mi = miller.match_indices(h0, h1)
  assert tuple(mi.singles(0)) == (4,)
  assert tuple(mi.singles(1)) == ()
  assert tuple(mi.pairs()) == ((0,2), (1,0), (2,3), (3,1))
  assert tuple(mi.pair_selection(0)) == (1, 1, 1, 1, 0)
  assert tuple(mi.single_selection(0)) == (0, 0, 0, 0, 1)
  assert tuple(mi.pair_selection(1)) == (1, 1, 1, 1)
  assert tuple(mi.single_selection(1)) == (0, 0, 0, 0)
  assert tuple(mi.paired_miller_indices(0)) \
      == tuple(h0.select(mi.pair_selection(0)))
  l1 = list(mi.paired_miller_indices(1))
  l2 = list(h1.select(mi.pair_selection(1)))
  l1.sort()
  l2.sort()
  assert l1 == l2
  assert approx_equal(tuple(mi.plus(d0, d1)), (31, 12, 43, 24))
  assert approx_equal(tuple(mi.minus(d0, d1)), (-29,-8,-37,-16))
  assert approx_equal(tuple(mi.multiplies(d0, d1)), (30,20,120,80))
  assert approx_equal(tuple(mi.divides(d0, d1)), (1/30.,2/10.,3/40.,4/20.))
  assert approx_equal(tuple(mi.additive_sigmas(d0, d1)), [
    math.sqrt(x*x+y*y) for x,y in ((1,30), (2,10), (3,40), (4,20))])
  q = flex.size_t((3,2,0,4,1))
  h1 = h0.select(q)
  assert tuple(miller.match_indices(h1, h0).permutation()) == tuple(q)
  p = miller.match_indices(h0, h1).permutation()
  assert tuple(p) == (2,4,1,0,3)
  assert tuple(h1.select(p)) == tuple(h0)

def exercise_merge_equivalents():
  i = flex.miller_index(((1,2,3), (1,2,3), (3,0,3), (3,0,3), (3,0,3), (1,1,2)))
  d = flex.double((1,2,3,4,5,6))
  m = miller.ext.merge_equivalents_real(i, d)
  assert tuple(m.indices) == ((1,2,3), (3,0,3), (1,1,2))
  assert approx_equal(m.data, (3/2., 4, 6))
  assert tuple(m.redundancies) == (2,3,1)
  assert approx_equal(m.r_linear, (1/3., 1/6., 0))
  assert approx_equal(m.r_square, (0.1, 0.04, 0))
  assert approx_equal(m.r_int, (1.+2.)/(3.+12.))
  #
  s = flex.double((1/3.,1/2.,1/4.,1/6.,1/3.,1/5.))
  m = miller.ext.merge_equivalents_obs(i, d, s)
  assert tuple(m.indices) == ((1,2,3), (3,0,3), (1,1,2))
  assert approx_equal(m.data, (17/13., (16*3+36*4+9*5)/(16+36+9.), 6))
  assert approx_equal(m.sigmas, (math.sqrt(1/2./2),0.84077140277/3**0.5,1/5.))
  assert tuple(m.redundancies) == (2,3,1)
  assert approx_equal(m.r_linear, (1/3., 0.1762295, 0))
  assert approx_equal(m.r_square, (0.1147929, 0.0407901, 0))
  assert approx_equal(m.r_int, (abs(1-17/13.)+abs(2-17/13.)
                              + abs(3-237/61.)+abs(4-237/61.)+abs(5-237/61.)
                              ) / (1 + 2 + 3 + 4 + 5) )
  #
  d = flex.complex_double(
    [complex(-1.706478,  0.248638),
     complex( 1.097872, -0.983523),
     complex( 0.147183,  2.625064),
     complex(-0.933310,  2.496886),
     complex( 1.745500, -0.686761),
     complex(-0.620066,  2.097776)])
  m = miller.ext.merge_equivalents_complex(i, d)
  assert tuple(m.indices) == ((1,2,3), (3,0,3), (1,1,2))
  assert approx_equal(m.data, [
    complex(-0.304303,-0.367443),
    complex( 0.319791, 1.478396),
    complex(-0.620066, 2.097776)])
  assert tuple(m.redundancies) == (2,3,1)
  #
  d = flex.hendrickson_lattman(
    [(-1.706478,  0.248638,  1.653352, -2.411313),
     ( 1.097872, -0.983523, -2.756402,  0.294464),
     ( 0.147183,  2.625064,  1.003636,  2.563517),
     (-0.933310,  2.496886,  2.040418,  0.371885),
     ( 1.745500, -0.686761, -2.291345, -2.386650),
     (-0.620066,  2.097776,  0.099784,  0.268107)])
  m = miller.ext.merge_equivalents_hl(i, d)
  assert tuple(m.indices) == ((1,2,3), (3,0,3), (1,1,2))
  assert approx_equal(m.data, [
    (-0.3043030, -0.3674425, -0.5515250, -1.0584245),
    ( 0.3197910,  1.4783963,  0.2509030,  0.1829173),
    (-0.6200660,  2.0977760,  0.0997840,  0.2681070)])
  assert tuple(m.redundancies) == (2,3,1)
  #
  d = flex.bool((True,True,False,False,False,True))
  m = miller.ext.merge_equivalents_exact_bool(i, d)
  assert tuple(m.indices) == ((1,2,3), (3,0,3), (1,1,2))
  assert list(m.data) == [True, False, True]
  assert tuple(m.redundancies) == (2,3,1)
  d = flex.bool((True,True,False,True,False,True))
  try: m = miller.ext.merge_equivalents_exact_bool(i, d)
  except RuntimeError, e:
    assert str(e) == "cctbx Error: merge_equivalents_exact:"\
      " incompatible flags for hkl = (3, 0, 3)"
  else: raise Exception_expected
  d = flex.int((3,3,5,5,5,7))
  m = miller.ext.merge_equivalents_exact_int(i, d)
  assert list(m.data) == [3, 5, 7]
  #
  i = flex.miller_index(((1,2,3), (3,0,3), (1,1,2)))
  d = flex.double((1,2,3))
  m = miller.ext.merge_equivalents_real(i,d)
  assert m.r_int == 0

def exercise_phase_integral():
  sg = sgtbx.space_group_info("P 21 21 21").group()
  i = flex.miller_index([(1,2,3), (3,0,3)])
  hl = flex.hendrickson_lattman([(1,2,3,4),(-2,3,-4,-5)])
  integrator = miller.phase_integrator(n_steps=10)
  assert integrator.n_steps() == 10
  integrator = miller.phase_integrator()
  assert integrator.n_steps() == 360//5
  assert approx_equal(integrator(sg.phase_restriction(i[0]), hl[0]),
    0.78832161462+0.466993941292j)
  assert approx_equal(integrator(sg.phase_restriction(i[1]), hl[1]),
    6.09275186883e-17+0.995054753687j)
  assert approx_equal(integrator(
      space_group=sg, miller_indices=i, hendrickson_lattman_coefficients=hl),
    [(0.78832161462020822+0.46699394129187444j),
     (6.0927518688296534e-17+0.99505475368673046j)])

def exercise_phase_transfer():
  sg = sgtbx.space_group_info("P 21 21 21").group()
  i = flex.miller_index(((1,2,3), (3,0,3)))
  a = flex.double((-3.6,4.6))
  p = flex.complex_double((1+2j, 0))
  assert approx_equal(tuple(miller.phase_transfer(sg, i, a, p, 1.e-10)),
                      ((-1.6099689-3.2199379j), 0j))
  a = flex.complex_double((3.6,4.6))
  try:
    miller.phase_transfer(sg, i, a, p)
  except Exception, e:
    if (str(e.__class__).find("Boost.Python.ArgumentError") < 0):
      raise RuntimeError("Unexpected exception: %s" % str(e))
  else:
    raise Exception_expected

  a = flex.double((-3.6,4.6))
  p = flex.double((10,20))
  t = miller.phase_transfer(sg, i, a, p, True)
  assert approx_equal(tuple(flex.abs(t)), flex.abs(a))
  assert approx_equal(tuple(flex.arg(t, True)), (-170,90))
  p = p * (math.pi/180)
  t = miller.phase_transfer(sg, i, a, p, False)
  assert approx_equal(tuple(flex.abs(t)), flex.abs(a))
  assert approx_equal(tuple(flex.arg(t, True)), (-170,90))

def exercise_f_calc_map():
  i = flex.miller_index((   (1,1,0), (-1,-1,0), (1,2,3), (3,2,1) ))
  f = flex.complex_double((    1+1j,      2+2j,     3+3j,   4+4j ))
  f_map = miller.f_calc_map(i, f, anomalous_flag=True)
  assert f_map[(1,1,0)] == 1+1j
  assert f_map[(-1,-1,0)] == 2+2j
  assert f_map[(1,2,3)] == 3+3j
  assert f_map[(3,2,1)] == 4+4j
  assert f_map[(2,2,2)] == 0
  f_map = miller.f_calc_map(i, f, anomalous_flag=False)
  assert f_map[(1,1,0)] == 2-2j
  assert f_map[(-1,-1,0)] == 2+2j
  assert f_map[(1,2,3)] == 3+3j
  assert f_map[(3,2,1)] == 4+4j
  assert f_map[(-1,-3,-5)] == 0

def exercise_union_of_indices():
  u = miller.union_of_indices_registry()
  assert u.as_array().size() == 0
  u.update(indices=flex.miller_index([(1,2,0)]))
  assert list(u.as_array()) == [(1,2,0)]
  u.update(indices=flex.miller_index([(1,2,0)]))
  assert list(u.as_array()) == [(1,2,0)]
  u.update(indices=flex.miller_index([(3,2,0)]))
  assert sorted(u.as_array()) == [(1,2,0), (3,2,0)]

def exercise_slices () :
  i = flex.miller_index(((1,2,3), (3,0,3), (2,4,1),(0,1,2)))
  s = miller.simple_slice(
    indices=i,
    slice_axis=2,
    slice_index=3)
  assert (list(s) == [True, True, False, False])
  s = miller.multi_slice(
    indices=i,
    slice_axis=0,
    slice_start=1,
    slice_end=2)
  assert (list(s) == [True, False, True, False])

def run(args):
  assert len(args) == 0
  exercise_f_calc_map()
  exercise_sym_equiv()
  exercise_asu()
  exercise_bins()
  exercise_expand()
  exercise_index_generator()
  exercise_index_span()
  exercise_match_bijvoet_mates()
  exercise_merge_equivalents()
  exercise_match_indices()
  exercise_phase_integral()
  exercise_phase_transfer()
  exercise_union_of_indices()
  exercise_slices()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
