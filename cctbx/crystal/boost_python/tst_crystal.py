from cctbx import crystal
from cctbx import sgtbx
import cctbx.crystal.direct_space_asu
from cctbx import uctbx
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.itertbx import count

def exercise_direct_space_asu():
  cp = crystal.direct_space_asu.float_cut_plane(n=[-1,0,0], c=1)
  assert approx_equal(cp.n, [-1,0,0])
  assert approx_equal(cp.c, 1)
  assert approx_equal(cp.evaluate(point=[0,2,3]), 1)
  assert approx_equal(cp.evaluate(point=[1,2,3]), 0)
  assert cp.is_inside(point=[0.99,0,0], epsilon=0)
  assert not cp.is_inside([1.01,0,0])
  assert approx_equal(cp.get_point_in_plane(), [1,0,0])
  cp.n = [0,-1,0]
  assert approx_equal(cp.n, [0,-1,0])
  cp.c = 2
  assert approx_equal(cp.c, 2)
  assert approx_equal(cp.get_point_in_plane(), [0,2,0])
  unit_cell = uctbx.unit_cell((1,1,1,90,90,90))
  cpb = cp.add_buffer(unit_cell=unit_cell, thickness=0.5)
  assert approx_equal(cpb.n, cp.n)
  assert approx_equal(cpb.c, 2.5)
  facets = []
  for i in xrange(3):
    n = [0,0,0]
    n[i] = -1
    facets.append(crystal.direct_space_asu.float_cut_plane(n=n, c=i+1))
  asu = crystal.direct_space_asu.float_asu(
    unit_cell=unit_cell,
    facets=facets,
    is_inside_epsilon=1.e-6)
  assert asu.unit_cell().is_similar_to(unit_cell)
  for i in xrange(3):
    n = [0,0,0]
    n[i] = -1
    assert approx_equal(asu.facets()[i].n, n)
    assert approx_equal(asu.facets()[i].c, i+1)
  assert asu.is_inside([0.99,0.49,0.32])
  eps = 0.02
  assert not asu.is_inside([0.99+eps,0.49+eps,0.32+eps])
  buf_asu = asu._add_buffer(0.2)
  assert buf_asu.is_inside([0.99+0.2,0.49+0.2,0.32+0.2])
  eps = 0.02
  assert not buf_asu.is_inside([0.99+0.2+eps,0.49+0.2+eps,0.32+0.2+eps])
  assert len(asu.volume_vertices()) == 1
  for cartesian in [00000,0001]:
    assert approx_equal(asu.volume_vertices(
      cartesian=cartesian, epsilon=1.e-6)[0], (1.0, 2.0, 3.0))
  asu = crystal.direct_space_asu.float_asu(
    unit_cell=unit_cell,
    facets=[crystal.direct_space_asu.float_cut_plane(n=n,c=c) for n,c in [
      [(0, 0, 1), -1/2.],
      [(-1, -1, 0), 1],
      [(0, 1, -1), 3/4.],
      [(1, 0, -1), 1/4.]]])
  assert approx_equal(asu.box_min(), [0.25, -0.25, 0.5])
  assert approx_equal(asu.box_max(), [1.25, 0.75, 1.0])
  asu_mappings = crystal.direct_space_asu.asu_mappings(
    space_group=sgtbx.space_group("P 2 3").change_basis(
      sgtbx.change_of_basis_op("x+1/4,y-1/4,z+1/2")),
    asu=asu,
    buffer_thickness=0.1,
    sym_equiv_epsilon=1.e-6)
  asu_mappings.reserve(n_sites_final=10)
  assert asu_mappings.space_group().order_z() == 12
  assert len(asu_mappings.asu().facets()) == 4
  assert asu_mappings.unit_cell().is_similar_to(unit_cell)
  assert approx_equal(asu_mappings.buffer_thickness(), 0.1)
  assert approx_equal(asu_mappings.asu_buffer().box_min(),
    [0.0085786, -0.4914214, 0.4])
  assert approx_equal(asu_mappings.sym_equiv_epsilon(), 1.e-6)
  assert approx_equal(asu_mappings.buffer_covering_sphere().radius(),0.8071081)
  sites_seq = [
    [3.1,-2.2,1.3],
    [-4.3,1.7,0.4]]
  assert asu_mappings.mappings().size() == 0
  asu_mappings.process(original_site=sites_seq[0])
  assert asu_mappings.mappings().size() == 1
  asu_mappings.process(original_site=sites_seq[1])
  assert asu_mappings.mappings().size() == 2
  assert not asu_mappings.is_locked()
  asu_mappings.lock()
  assert asu_mappings.is_locked()
  try: asu_mappings.process(original_site=[0,0,0])
  except RuntimeError, e: assert str(e).find("is_locked") > 0
  else: raise RuntimeError("Exception expected.")
  mappings = asu_mappings.mappings()[0]
  assert len(mappings) == 5
  am = mappings[0]
  assert am.i_sym_op() == 3
  assert am.unit_shifts() == (1,3,2)
  assert asu.is_inside(am.mapped_site())
  assert approx_equal(asu_mappings.mapped_sites_min(), [0.15,-0.4,0.4])
  assert approx_equal(asu_mappings.mapped_sites_max(), [1.05,0.6,0.65])
  assert approx_equal(asu_mappings.mapped_sites_span(), [0.9,1.0,0.25])
  for am in mappings:
    assert asu_mappings.asu_buffer().is_inside(am.mapped_site())
  o = matrix.sqr(asu_mappings.unit_cell().orthogonalization_matrix())
  f = matrix.sqr(asu_mappings.unit_cell().fractionalization_matrix())
  for i_seq,m_i_seq in asu_mappings.mappings().items():
    for i_sym in xrange(len(m_i_seq)):
      rt_mx = asu_mappings.get_rt_mx(i_seq, i_sym)
      site_frac = rt_mx * sites_seq[i_seq]
      site_cart = asu_mappings.unit_cell().orthogonalize(site_frac)
      assert approx_equal(m_i_seq[i_sym].mapped_site(), site_cart)
      r = matrix.sqr(float(rt_mx.r().inverse()))
      assert approx_equal(
        asu_mappings.r_inv_cart(i_seq=i_seq, i_sym=i_sym),
        (o*r*f).elems)
  pair_generator = crystal.neighbors_simple_pair_generator(asu_mappings)
  assert not pair_generator.at_end()
  assert len(asu_mappings.mappings()[1]) == 6
  index_pairs = []
  for index_pair in pair_generator:
    index_pairs.append((index_pair.i_seq, index_pair.j_seq, index_pair.j_sym))
    assert index_pair.dist_sq == -1
  assert pair_generator.at_end()
  assert index_pairs == [
    (0,0,1),(0,0,2),(0,0,3),(0,0,4),
    (0,1,0),(0,1,1),(0,1,2),(0,1,3),(0,1,4),(0,1,5),
    (1,1,1),(1,1,2),(1,1,3),(1,1,4),(1,1,5)]
  for two_flag,buffer_thickness,expected_index_pairs,expected_n_box in [
    (00000, 0.04, [], (1,1,1)),
    (00000, 0.1, [(0,0,1),(0,0,2),(0,0,3),(0,0,4)], (2,3,1)),
    (00001, 0, [(0, 1, 0)], (1,1,1)),
    (00001, 0.04, [(0, 1, 0), (0, 1, 1), (1, 1, 1)], (1,2,1))]:
    asu_mappings = crystal.direct_space_asu.asu_mappings(
      space_group=sgtbx.space_group("P 2 3").change_basis(
        sgtbx.change_of_basis_op("x+1/4,y-1/4,z+1/2")),
      asu=asu,
      buffer_thickness=buffer_thickness,
      sym_equiv_epsilon=1.e-6)
    asu_mappings.process(original_site=[3.1,-2.2,1.3])
    assert asu_mappings.mappings().size() == 1
    if (two_flag):
      asu_mappings.process(original_site=[-4.3,1.7,0.4])
    pair_generator = crystal.neighbors_simple_pair_generator(asu_mappings)
    index_pairs = []
    for index_pair in pair_generator:
      index_pairs.append((index_pair.i_seq,index_pair.j_seq,index_pair.j_sym))
      assert index_pair.dist_sq == -1
    assert pair_generator.at_end()
    assert index_pairs == expected_index_pairs
    pair_generator.restart()
    if (len(expected_index_pairs) == 0):
      assert pair_generator.at_end()
    else:
      assert not pair_generator.at_end()
    index_pairs = []
    for index_pair in pair_generator:
      index_pairs.append((index_pair.i_seq,index_pair.j_seq,index_pair.j_sym))
      assert index_pair.dist_sq == -1
    assert pair_generator.at_end()
    assert index_pairs == expected_index_pairs
    simple_pair_generator = crystal.neighbors_simple_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=100)
    assert approx_equal(simple_pair_generator.distance_cutoff_sq(), 100*100)
    fast_pair_generator = crystal.neighbors_fast_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=100,
      epsilon=1.e-6)
    assert approx_equal(fast_pair_generator.distance_cutoff_sq(), 100*100)
    assert approx_equal(fast_pair_generator.epsilon()/1.e-6, 1)
    assert fast_pair_generator.n_box() == (1,1,1)
    index_pairs = []
    dist_sq = flex.double()
    for index_pair in simple_pair_generator:
      index_pairs.append((index_pair.i_seq,index_pair.j_seq,index_pair.j_sym))
      assert index_pair.dist_sq > 0
      assert approx_equal(
        asu_mappings.diff_vec(pair=index_pair),
        index_pair.diff_vec)
      assert approx_equal(
        matrix.col(asu_mappings.diff_vec(pair=index_pair)).norm(),
        index_pair.dist_sq)
      dist_sq.append(index_pair.dist_sq)
    assert simple_pair_generator.at_end()
    assert index_pairs == expected_index_pairs
    if (len(index_pairs) == 0):
      assert fast_pair_generator.at_end()
    else:
      assert not fast_pair_generator.at_end()
    index_pairs = []
    for index_pair in fast_pair_generator:
      index_pairs.append((index_pair.i_seq,index_pair.j_seq,index_pair.j_sym))
    assert index_pairs == expected_index_pairs
    assert fast_pair_generator.at_end()
    distances = flex.sqrt(dist_sq)
    if (distances.size() > 0):
      cutoff = flex.mean(distances) + 1.e-5
    else:
      cutoff = 0
    short_dist_sq = dist_sq.select(distances <= cutoff)
    for pair_generator_type in [crystal.neighbors_simple_pair_generator,
                                crystal.neighbors_fast_pair_generator]:
      if (    pair_generator_type is crystal.neighbors_fast_pair_generator
          and cutoff == 0): continue
      pair_generator = pair_generator_type(
        asu_mappings=asu_mappings,
        distance_cutoff=cutoff)
      index_pairs = []
      dist_sq = flex.double()
      for index_pair in pair_generator:
        index_pairs.append((index_pair.i_seq,index_pair.j_seq,index_pair.j_sym))
        assert index_pair.dist_sq > 0
        assert approx_equal(
          asu_mappings.diff_vec(pair=index_pair),
          index_pair.diff_vec)
        assert approx_equal(
          matrix.col(asu_mappings.diff_vec(pair=index_pair)).norm(),
          index_pair.dist_sq)
        dist_sq.append(index_pair.dist_sq)
      assert pair_generator.at_end()
      if (pair_generator_type is crystal.neighbors_simple_pair_generator):
        assert approx_equal(dist_sq, short_dist_sq)
      else:
        assert pair_generator.n_box() == expected_n_box
        short_dist_sq_sorted = short_dist_sq.select(
          flex.sort_permutation(short_dist_sq))
        dist_sq_sorted = dist_sq.select(
          flex.sort_permutation(dist_sq))
        assert approx_equal(dist_sq_sorted, short_dist_sq_sorted)

def run():
  exercise_direct_space_asu()
  print "OK"

if (__name__ == "__main__"):
  run()
