import math
import random
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import miller
from cctbx.array_family import flex
from cctbx import utils
from scitbx.python_utils.complex_math import polar
from scitbx.test_utils import approx_equal

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
  assert i[0].ht_angle(0) == 0
  assert i[0].ht_angle(1) == 0
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
    assert approx_equal(m.phase_in(m.phase_eq(30, 0), 0), 30)
    assert approx_equal(m.phase_in(m.phase_eq(30, 1), 1), 30)
    assert approx_equal(m.phase_eq(30*math.pi/180),
                        m.phase_eq(30, 1)*math.pi/180)
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
  assert e.multiplicity(0) == 6
  assert e.multiplicity(1) == 3
  assert e.f_mates(0) == 2
  assert e.f_mates(1) == 1
  assert e.epsilon() == 1
  for j in xrange(e.multiplicity(0)):
    if (j < e.multiplicity(1)):
      assert e(j).h() == e.indices()[j].h()
    else:
      assert e(j).friedel_flag()
  assert e.is_valid_phase(10)
  assert e.is_valid_phase(10, 0)
  assert e.is_valid_phase(10, 1)
  assert e.is_valid_phase(10, 1, 1.e-5)
  for anomalous_flag in (0,1):
    j = e.p1_listing(anomalous_flag)
    assert len(j) == len(i)
  s = sgtbx.space_group("P 41")
  h = (3,5,0)
  e = miller.sym_equiv_indices(s, h)
  i = e.indices()
  assert e.is_centric()
  assert e.multiplicity(0) == 4
  assert e.multiplicity(1) == 4
  assert e.f_mates(0) == 1
  assert e.f_mates(1) == 1
  r = e.phase_restriction()
  assert r.ht() == 0
  j = e.p1_listing(0)
  assert len(j) == len(i)/2

def exercise_map_to_asu(sg_symbol):
  sg_type = sgtbx.space_group_type(sg_symbol)
  index_abs_range = (4,4,4)
  for anomalous_flag in (0,1):
    m = miller.index_generator(
      sg_type, anomalous_flag, index_abs_range).to_array()
    a = flex.double()
    p = flex.double()
    c = flex.hendrickson_lattman()
    for i in m.indices():
      a.append(random.random())
      p.append(random.random() * 2)
      c.append([random.random() for i in xrange(4)])
    f = flex.polar(a, p)
    p = [p, p*(180/math.pi)]
  m_random = flex.miller_index()
  p_random = [flex.double(), flex.double()]
  c_random = flex.hendrickson_lattman()
  f_random = flex.complex_double()
  for i,h_asym in m.items():
    h_eq = miller.sym_equiv_indices(sg_type.group(), h_asym)
    i_eq = random.randrange(h_eq.multiplicity(anomalous_flag))
    h_i = h_eq(i_eq)
    m_random.append(h_i.h())
    for deg in (0,1):
      p_random[deg].append(h_i.phase_eq(p[deg][i], deg))
    f_random.append(h_i.complex_eq(f[i]))
    c_random.append(h_i.hendrickson_lattman_eq(c[i]))
  m_random_copy = m_random.deep_copy()
  miller.map_to_asu(sg_type, anomalous_flag, m_random_copy)
  for i,h_asym in m.items():
    assert h_asym == m_random_copy[i]
  m_random_copy = m_random.deep_copy()
  miller.map_to_asu(sg_type, anomalous_flag, m_random_copy, f_random)
  for i,h_asym in m.items():
    assert h_asym == m_random_copy[i]
  for i,f_asym in f.items():
    assert abs(f_asym - f_random[i]) < 1.e-6
  m_random_copy = m_random.deep_copy()
  a_random = a.deep_copy()
  miller.map_to_asu(sg_type, anomalous_flag, m_random_copy, a_random)
  for i,h_asym in m.items():
    assert h_asym == m_random_copy[i]
  for i,a_asym in a.items():
    assert a_asym == a_random[i]
  for deg in (0,1):
    m_random_copy = m_random.deep_copy()
    miller.map_to_asu(
      sg_type, anomalous_flag, m_random_copy, p_random[deg], deg)
    for i,h_asym in m.items():
      assert h_asym == m_random_copy[i]
    for i,p_asym in p[deg].items():
      assert utils.phase_error(p_asym, p_random[deg][i], deg) < 1.e-5
  m_random_copy = m_random.deep_copy()
  miller.map_to_asu(sg_type, anomalous_flag, m_random_copy, c_random)
  for i,h_asym in m.items():
    assert h_asym == m_random_copy[i]
  for i,c_asym in c.items():
    for j in xrange(4):
      assert abs(c_asym[j] - c_random[i][j]) < 1.e-5

def exercise_asu():
  sg_type = sgtbx.space_group_type("P 41")
  asu = sgtbx.reciprocal_space_asu(sg_type)
  miller_indices = flex.miller_index(((1,2,3), (3,5,0)))
  for h in miller_indices:
    h_eq = miller.sym_equiv_indices(sg_type.group(), h)
    for i_eq in xrange(h_eq.multiplicity(0)):
      h_i = h_eq(i_eq)
      for anomalous_flag in (0,1):
        a = miller.asym_index(sg_type.group(), asu, h_i.h())
        assert a.h() == h
        o = a.one_column(anomalous_flag)
        assert o.i_column() == 0
        t = a.two_column(anomalous_flag)
        assert t.h() == h
        assert (o.h() != h) == (t.i_column() == 1)
        assert not anomalous_flag or (t.i_column() != 0) == h_i.friedel_flag()
        assert anomalous_flag or t.i_column() == 0
  miller.map_to_asu(sg_type, 0, miller_indices)
  data = flex.double((0,0))
  miller.map_to_asu(sg_type, 0, miller_indices, data)
  miller.map_to_asu(sg_type, 0, miller_indices, data, 0)
  miller.map_to_asu(sg_type, 0, miller_indices, data, 1)
  data = flex.complex_double((0,0))
  miller.map_to_asu(sg_type, 0, miller_indices, data)
  data = flex.hendrickson_lattman(((1,2,3,4),(2,3,4,5)))
  miller.map_to_asu(sg_type, 0, miller_indices, data)
  for sg_symbol in ("P 41", "P 31 1 2"):
    exercise_map_to_asu(sg_symbol)

def exercise_bins():
  uc = uctbx.unit_cell((11,11,13,90,90,120))
  sg_type = sgtbx.space_group_type("P 3 2 1")
  anomalous_flag = 0
  d_min = 1
  m = miller.index_generator(uc, sg_type, anomalous_flag, d_min).to_array()
  f = flex.double()
  for i in m.indices(): f.append(random.random())
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
  binning1 = miller.binning(uc, n_bins, m)
  assert binning1.unit_cell().is_similar_to(uc)
  assert binning1.n_bins_used() == n_bins
  assert binning1.limits().size() == n_bins + 1
  assert binning1.n_bins_all() == n_bins + 2
  binner1 = miller.binner(binning1, m)
  assert binner1.count(binner1.i_bin_d_too_large()) == 0
  assert binner1.count(binner1.i_bin_d_too_small()) == 0
  counts = binner1.counts()
  for i_bin in binner1.range_all():
    assert binner1.count(i_bin) == counts[i_bin]
    assert binner1.selection(i_bin).count(1) == counts[i_bin]
  assert list(binner1.range_all()) == range(binner1.n_bins_all())
  assert list(binner1.range_used()) == range(1, binner1.n_bins_used()+1)
  binning2 = miller.binning(uc, n_bins - 2,
    binning1.bin_d_min(2),
    binning1.bin_d_min(n_bins))
  binner2 = miller.binner(binning2, m)
  assert tuple(binner1.counts())[1:-1] == tuple(binner2.counts())
  array_indices = flex.size_t(tuple(m.indices()))
  perm_array_indices1 = flex.size_t()
  perm_array_indices2 = flex.size_t()
  for i_bin in binner1.range_all():
    perm_array_indices1.append(array_indices.select(binner1.selection(i_bin)))
    perm_array_indices2.append(binner1.array_indices(i_bin))
  assert perm_array_indices1.size() == m.size()
  assert perm_array_indices2.size() == m.size()
  assert tuple(perm_array_indices1) == tuple(perm_array_indices2)

def exercise_expand():
  sg = sgtbx.space_group("P 41 (1,-1,0)")
  h = flex.miller_index(((3,1,-2), (1,-2,0)))
  assert tuple(sg.is_centric(h)) == (0, 1)
  p1 = miller.expand_to_p1(sg, 0, h)
  p1_i0 = ((-3,-1,2), (-1, 3,2),(3,1,2),(1,-3,2),(1,-2, 0),(2,1,0))
  assert tuple(p1.indices()) == p1_i0
  assert p1.amplitudes().size() == 0
  assert p1.phases().size() == 0
  assert p1.structure_factors().size() == 0
  p1 = miller.expand_to_p1(sg, 1, h)
  assert tuple(p1.indices()) \
      == ((3,1,-2), (1,-3,-2), (-3,-1,-2), (-1,3,-2),
          (1,-2,0), (-2,-1,0), (-1,2,0), (2,1,0))
  a = flex.double((1,2))
  p1 = miller.expand_to_p1(sg, 0, h, a)
  assert tuple(p1.indices()) == p1_i0
  assert tuple(p1.amplitudes()) == (1,1,1,1,2,2)
  p = flex.double((10,90))
  p1 = miller.expand_to_p1(sg, 0, h, p, 1)
  assert approx_equal(tuple(p1.phases()), (-10,110,110,-10, 90,30))
  p1 = miller.expand_to_p1(sg, 1, h, a, p, 1)
  assert tuple(p1.amplitudes()) == (1,1,1,1,2,2,2,2)
  assert approx_equal(tuple(p1.phases()), (10,-110,-110,10, 90,-30,-90,30))
  p = flex.double([x * math.pi/180 for x in p])
  v = [x * math.pi/180 for x in p1.phases()]
  p1 = miller.expand_to_p1(sg, 1, h, a, p)
  assert approx_equal(tuple(p1.phases()), v)
  p1 = miller.expand_to_p1(sg, 1, h, a, p, 0)
  assert approx_equal(tuple(p1.phases()), v)
  f = flex.polar(a, p)
  p1 = miller.expand_to_p1(sg, 1, h, f)
  assert approx_equal(tuple(flex.abs(p1.structure_factors())),
                      (1,1,1,1,2,2,2,2))
  assert approx_equal(tuple(flex.arg(p1.structure_factors())), v)
  assert approx_equal(miller.statistical_mean(sg, 0, h, a), 4/3.)
  assert approx_equal(miller.statistical_mean(sg, 1, h, a), 3/2.)

def exercise_index_generator():
  uc = uctbx.unit_cell((11,11,13,90,90,120))
  sg_type = sgtbx.space_group_type("P 3 1 2")
  for anomalous_flag in (0,1):
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
  assert tuple(miller.index_generator(uc, sg_type, 0, 8)) \
         == ((0,0,1), (1, 0, 0))
  index_abs_range = (4,4,4)
  for sg_symbol in ("P 31 1 2", "P 31 2 1"):
    sg_type = sgtbx.space_group_type(sg_symbol)
    for anomalous_flag in (0,1):
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
  assert tuple(bm.pairs_hemisphere_selection("+")) \
      == (0001,00000,0001,00000,00000)
  assert tuple(bm.pairs_hemisphere_selection("-")) \
      == (00000,0001,00000,0001,00000)
  assert tuple(bm.singles_hemisphere_selection("+")) \
      == (00000,00000,00000,00000,0001)
  assert tuple(bm.singles_hemisphere_selection("-")) \
      == (00000,00000,00000,00000,00000)
  assert tuple(bm.miller_indices_in_hemisphere("+")) == ((1,2,3), (2,3,4))
  assert tuple(bm.miller_indices_in_hemisphere("-")) == ((-1,-2,-3),(-2,-3,-4))
  assert approx_equal(tuple(bm.minus(d0)), (-1, -1))
  assert approx_equal(tuple(bm.additive_sigmas(d0)),
                      [math.sqrt(x*x+y*y) for x,y in ((1,2), (3,4))])
  assert approx_equal(tuple(bm.average(d0)), (3/2., 7/2.))

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
  s = flex.double((1/3.,1/2.,1/4.,1/6.,1/3.,1/5.))
  m = miller.ext.merge_equivalents(i, d)
  assert tuple(m.indices()) == ((1,2,3), (3,0,3), (1,1,2))
  assert approx_equal(m.data(), (3/2., 4, 6))
  assert m.sigmas().size() == 0
  assert tuple(m.redundancies()) == (2,3,1)
  m = miller.ext.merge_equivalents(i, d, 1./(s*s))
  assert tuple(m.indices()) == ((1,2,3), (3,0,3), (1,1,2))
  assert approx_equal(m.data(), (17/13., (16*3+36*4+9*5)/(16+36+9.), 6))
  assert approx_equal(m.sigmas(), (math.sqrt(1/2.), 0.84077140277, 1/5.))
  assert tuple(m.redundancies()) == (2,3,1)

def exercise_phase_transfer():
  sg = sgtbx.space_group_info("P 21 21 21").group()
  i = flex.miller_index(((1,2,3), (3,0,3)))
  a = flex.double((-3.6,4.6))
  p = flex.complex_double((1+2j, 0))
  assert approx_equal(tuple(miller.phase_transfer(sg, i, a, p, 1.e-10)),
                      ((1.6099689+3.2199379j), 0j))
  a = flex.complex_double((3.6,4.6))
  assert approx_equal(tuple(miller.phase_transfer(sg, i, a, p, 1.e-10)),
                      ((1.6099689+3.2199379j), 0j))
  for a in(flex.double((-3.6,4.6)), flex.complex_double((3.6,4.6))):
    p = flex.double((10,20))
    t = miller.phase_transfer(sg, i, a, p, 0001)
    assert approx_equal(tuple(flex.abs(t)), flex.abs(a))
    assert approx_equal(tuple(flex.arg(t, 1)), (10,90))
    p = p * (math.pi/180)
    t = miller.phase_transfer(sg, i, a, p, 00000)
    assert approx_equal(tuple(flex.abs(t)), flex.abs(a))
    assert approx_equal(tuple(flex.arg(t, 1)), (10,90))

def run():
  exercise_sym_equiv()
  exercise_asu()
  exercise_bins()
  exercise_expand()
  exercise_index_generator()
  exercise_index_span()
  exercise_match_bijvoet_mates()
  exercise_merge_equivalents()
  exercise_match_indices()
  exercise_phase_transfer()
  print "OK"

if (__name__ == "__main__"):
  run()
