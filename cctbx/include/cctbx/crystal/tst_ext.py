from cctbx import crystal
from cctbx import sgtbx
import cctbx.crystal.direct_space_asu
from cctbx import uctbx
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.itertbx import count
import sys

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
  asu = crystal.direct_space_asu_float_asu(
    unit_cell=unit_cell,
    facets=[])
  assert asu.volume_vertices().size() == 0
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
  assert approx_equal(asu.is_inside_epsilon(), 1.e-6)
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
  assert approx_equal(asu.box_min(cartesian=0001), [0.25, -0.25, 0.5])
  assert approx_equal(asu.box_max(cartesian=0001), [1.25, 0.75, 1.0])
  asu_mappings = crystal.direct_space_asu.asu_mappings(
    space_group=sgtbx.space_group("P 2 3").change_basis(
      sgtbx.change_of_basis_op("x+1/4,y-1/4,z+1/2")),
    asu=asu,
    buffer_thickness=0.1)
  asu_mappings.reserve(n_sites_final=10)
  assert asu_mappings.space_group().order_z() == 12
  assert len(asu_mappings.asu().facets()) == 4
  assert asu_mappings.unit_cell().is_similar_to(unit_cell)
  assert approx_equal(asu_mappings.buffer_thickness(), 0.1)
  assert approx_equal(asu_mappings.asu_buffer().box_min(),
    [0.0085786, -0.4914214, 0.4])
  assert approx_equal(asu_mappings.buffer_covering_sphere().radius(),0.8071081)
  sites_seq = [
    [3.1,-2.2,1.3],
    [-4.3,1.7,0.4]]
  assert asu_mappings.mappings().size() == 0
  asu_mappings.process(
    original_site=sites_seq[0],
    min_distance_sym_equiv=0.01)
  assert asu_mappings.mappings().size() == 1
  site_symmetry = sgtbx.site_symmetry(
    unit_cell=asu_mappings.asu().unit_cell(),
    space_group=asu_mappings.space_group(),
    original_site=sites_seq[1],
    min_distance_sym_equiv=0.01)
  asu_mappings.process(
    original_site=sites_seq[1],
    site_symmetry_ops=site_symmetry)
  assert asu_mappings.mappings().size() == 2
  assert asu_mappings.n_sites_in_asu_and_buffer() == 11
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
  assert list(asu_mappings.site_symmetry_table().indices()) == [0, 0]
  for am in mappings:
    assert asu_mappings.asu_buffer().is_inside(am.mapped_site())
  o = matrix.sqr(asu_mappings.unit_cell().orthogonalization_matrix())
  f = matrix.sqr(asu_mappings.unit_cell().fractionalization_matrix())
  for i_seq,m_i_seq in zip(count(), asu_mappings.mappings()):
    for i_sym in xrange(len(m_i_seq)):
      rt_mx = asu_mappings.get_rt_mx(i_seq, i_sym)
      site_frac = rt_mx * sites_seq[i_seq]
      site_cart = asu_mappings.unit_cell().orthogonalize(site_frac)
      assert approx_equal(m_i_seq[i_sym].mapped_site(), site_cart)
      assert approx_equal(
        asu_mappings.map_moved_site_to_asu(
          moved_original_site
            =asu_mappings.unit_cell().orthogonalize(sites_seq[i_seq]),
          i_seq=i_seq,
          i_sym=i_sym),
        site_cart)
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
    assert not pair_generator.is_simple_interaction(index_pair)
  assert pair_generator.at_end()
  assert index_pairs == [
    (0,0,1),(0,0,2),(0,0,3),(0,0,4),
    (0,1,0),(0,1,1),(0,1,2),(0,1,3),(0,1,4),(0,1,5),
    (1,0,1),(1,0,2),(1,0,3),(1,0,4),
    (1,1,1),(1,1,2),(1,1,3),(1,1,4),(1,1,5)]
  for two_flag,buffer_thickness,expected_index_pairs,expected_n_boxes in [
    (00000, 0.04, [], (1,1,1)),
    (00000, 0.1, [(0,0,1),(0,0,2),(0,0,3),(0,0,4)], (2,3,1)),
    (00001, 0, [(0, 1, 0)], (1,1,1)),
    (00001, 0.04, [(0, 1, 0), (0, 1, 1), (1, 1, 1)], (1,2,1))]:
    asu_mappings = crystal.direct_space_asu.asu_mappings(
      space_group=sgtbx.space_group("P 2 3").change_basis(
        sgtbx.change_of_basis_op("x+1/4,y-1/4,z+1/2")),
      asu=asu,
      buffer_thickness=buffer_thickness)
    asu_mappings.process_sites_frac(
      original_sites=flex.vec3_double([[3.1,-2.2,1.3]]),
      min_distance_sym_equiv=0.01)
    assert asu_mappings.mappings().size() == 1
    if (two_flag):
      asu_mappings.process_sites_cart(
        original_sites=flex.vec3_double([
          asu.unit_cell().orthogonalize([-4.3,1.7,0.4])]),
      min_distance_sym_equiv=0.01)
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
    assert pair_generator.count_pairs() == len(index_pairs)
    pair_generator.restart()
    index_pairs = []
    for index_pair in pair_generator:
      index_pairs.append((index_pair.i_seq,index_pair.j_seq,index_pair.j_sym))
      assert index_pair.dist_sq == -1
    assert pair_generator.at_end()
    assert index_pairs == expected_index_pairs
    simple_pair_generator = crystal.neighbors_simple_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=100,
      minimal=00000)
    assert simple_pair_generator.asu_mappings().is_locked()
    assert approx_equal(simple_pair_generator.distance_cutoff_sq(), 100*100)
    fast_pair_generator = crystal.neighbors_fast_pair_generator(
      asu_mappings=asu_mappings,
      distance_cutoff=100,
      minimal=00000,
      epsilon=1.e-6)
    assert fast_pair_generator.asu_mappings().is_locked()
    assert approx_equal(fast_pair_generator.distance_cutoff_sq(), 100*100)
    assert approx_equal(fast_pair_generator.epsilon()/1.e-6, 1)
    assert fast_pair_generator.n_boxes() == (1,1,1)
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
        assert not asu_mappings.is_simple_interaction(pair=index_pair)
        assert asu_mappings.interaction_type_id(pair=index_pair) == 0
        dist_sq.append(index_pair.dist_sq)
      assert pair_generator.at_end()
      if (pair_generator_type is crystal.neighbors_simple_pair_generator):
        assert approx_equal(dist_sq, short_dist_sq)
      else:
        assert pair_generator.n_boxes() == expected_n_boxes
        short_dist_sq_sorted = short_dist_sq.select(
          flex.sort_permutation(short_dist_sq))
        dist_sq_sorted = dist_sq.select(
          flex.sort_permutation(dist_sq))
        assert approx_equal(dist_sq_sorted, short_dist_sq_sorted)
      assert pair_generator.count_pairs() == 0
      pair_generator.restart()
      assert pair_generator.count_pairs() == len(index_pairs)
  pair = asu_mappings.make_trial_pair(i_seq=1, j_seq=0, j_sym=0)
  assert pair.i_seq == 1
  assert pair.j_seq == 0
  assert pair.j_sym == 0
  assert not pair.is_active()
  pair = asu_mappings.make_pair(i_seq=0, j_seq=1, j_sym=1)
  assert pair.i_seq == 0
  assert pair.j_seq == 1
  assert pair.j_sym == 1
  assert pair.is_active()
  from cctbx import xray
  structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell="12.548 12.548 20.789 90.000 90.000 120.000",
      space_group_symbol="P63/mmc"),
    scatterers=flex.xray_scatterer(
      [xray.scatterer(label="Si", site=site) for site in [
        (0.2466,0.9965,0.2500),
        (0.5817,0.6706,0.1254),
        (0.2478,0.0000,0.0000)]]))
  asu_mappings = structure.asu_mappings(buffer_thickness=3.5)
  assert list(asu_mappings.site_symmetry_table().indices()) == [1,0,2]
  assert asu_mappings.n_sites_in_asu_and_buffer() == 33
  pair = asu_mappings.make_pair(i_seq=1, j_seq=0, j_sym=1)
  assert pair.i_seq == 1
  assert pair.j_seq == 0
  assert pair.j_sym == 1
  assert pair.is_active(minimal=00000)
  assert not pair.is_active(minimal=0001)
  for i_seq,m in zip(count(), asu_mappings.mappings()):
    for i_sym in xrange(len(m)):
      rt_mx = asu_mappings.get_rt_mx(i_seq=i_seq, i_sym=i_sym)
      i_sym_found = asu_mappings.find_i_sym(i_seq=i_seq, rt_mx=rt_mx)
      assert i_sym_found == i_sym
      i_sym_found = asu_mappings.find_i_sym(
        i_seq=i_seq,
        rt_mx=sgtbx.rt_mx("0,0,0"))
      assert i_sym_found == -1
      if (i_sym != 0):
        pair_i_seq = max(0, i_seq-1)
        pair = asu_mappings.make_trial_pair(pair_i_seq, i_seq, i_sym)
        assert asu_mappings.get_rt_mx_i(pair=pair) \
            == asu_mappings.get_rt_mx(pair_i_seq, 0)
        assert asu_mappings.get_rt_mx_j(pair=pair) \
            == asu_mappings.get_rt_mx(i_seq, i_sym)
  assert str(asu_mappings.special_op(0)) == "x,y,1/4"
  assert str(asu_mappings.special_op(1)) == "x,y,z"
  assert str(asu_mappings.special_op(2)) == "x-1/2*y,0,0"
  asu_mappings = structure[:0].asu_mappings(buffer_thickness=3.5)
  asu_mappings.process_sites_frac(
    original_sites=structure.scatterers().extract_sites())
  assert asu_mappings.n_sites_in_asu_and_buffer() == 33
  asu_mappings = structure[:0].asu_mappings(buffer_thickness=3.5)
  asu_mappings.process_sites_frac(
    original_sites=structure.scatterers().extract_sites(),
    site_symmetry_table=structure.site_symmetry_table())
  assert asu_mappings.n_sites_in_asu_and_buffer() == 33
  asu_mappings = structure[:0].asu_mappings(buffer_thickness=3.5)
  asu_mappings.process_sites_cart(
    original_sites=structure.sites_cart())
  assert asu_mappings.n_sites_in_asu_and_buffer() == 33
  asu_mappings = structure[:0].asu_mappings(buffer_thickness=3.5)
  asu_mappings.process_sites_cart(
    original_sites=structure.sites_cart(),
    site_symmetry_table=structure.site_symmetry_table())
  assert asu_mappings.n_sites_in_asu_and_buffer() == 33

def check_pair_asu_table(asu_table, expected_asu_pairs):
  ip = count()
  for i_seq,asu_dict in enumerate(asu_table.table()):
    for j_seq,j_sym_groups in asu_dict.items():
      for j_sym_group in j_sym_groups:
        for j_sym in j_sym_group:
          if (0 or "--Verbose" in sys.argv[1:]):
            print str([i_seq, j_seq, j_sym]) + ","
          assert [i_seq, j_seq, j_sym] == expected_asu_pairs[ip.next()]

def exercise_pair_tables():
  d = crystal.pair_sym_dict()
  assert len(d) == 0
  sym_ops = sgtbx.space_group("P 41").all_ops()
  for i,j_sym in enumerate([10,18,13]):
    d[j_sym] = crystal.pair_sym_ops(sym_ops[:i])
    assert len(d) == i+1
    assert len(d[j_sym]) == i
    assert [str(s) for s in sym_ops[:i]] == [str(s) for s in d[j_sym]]
    d[j_sym] = sym_ops[:i]
    assert [str(s) for s in sym_ops[:i]] == [str(s) for s in d[j_sym]]
  assert [key for key in d] == [10,13,18]
  assert d[13].size() == 2
  d[13].append(sym_ops[-1])
  assert d[13].size() == 3
  del d[13][0]
  assert d[13].size() == 2
  d[13].clear()
  assert d[13].size() == 0
  t = crystal.pair_sym_table()
  t.append(d)
  assert t.size() == 1
  assert len(t[0][10]) == 0
  t.append(d)
  assert t.size() == 2
  assert len(t[1][18]) == 1
  t = crystal.pair_sym_table(3)
  for d in t:
    assert len(d) == 0
  t[1][10] = sym_ops[:2]
  assert len(t[1]) == 1
  assert len(t[1][10]) == 2
  #
  t = crystal.pair_asu_table_table(3)
  for d in t:
    assert len(d) == 0
  t[1][10] = crystal.pair_asu_j_sym_groups()
  assert t[1][10].size() == 0
  t[1][10].append(crystal.pair_asu_j_sym_group())
  assert t[1][10].size() == 1
  assert t[1][10][0].size() == 0
  t[1][10][0].insert(3)
  assert t[1][10][0].size() == 1
  t[1][10].append(crystal.pair_asu_j_sym_group())
  assert t[1][10][1].size() == 0
  t[1][10][1].insert([4,5,4])
  assert t[1][10][1].size() == 2
  #
  from cctbx import xray
  structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell="12.548 12.548 20.789 90.000 90.000 120.000",
      space_group_symbol="P63/mmc"),
    scatterers=flex.xray_scatterer(
      [xray.scatterer(label="Si", site=site) for site in [
        (0.2466,0.9965,0.2500),
        (0.5817,0.6706,0.1254),
        (0.2478,0.0000,0.0000)]]))
  asu_mappings = structure.asu_mappings(buffer_thickness=3.5)
  asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  assert asu_table.asu_mappings().is_locked()
  assert asu_table.table().size() == 3
  assert not asu_table.contains(i_seq=0, j_seq=1, j_sym=2)
  assert asu_table.add_all_pairs(distance_cutoff=3.5) is asu_table
  assert [d.size() for d in asu_table.table()] == [2,3,2]
  assert asu_table.add_all_pairs(3.5, epsilon=1.e-6) is asu_table
  assert [d.size() for d in asu_table.table()] == [2,3,2]
  expected_asu_pairs = [
    [0, 0, 2], [0, 0, 1], [0, 1, 0], [0, 1, 8],
    [1, 0, 0], [1, 1, 4], [1, 1, 1], [1, 2, 0],
    [2, 1, 0], [2, 1, 11], [2, 2, 1], [2, 2, 2]]
  check_pair_asu_table(asu_table, expected_asu_pairs)
  for p in expected_asu_pairs:
    assert asu_table.contains(i_seq=p[0], j_seq=p[1], j_sym=p[2])
    pair = asu_mappings.make_trial_pair(*p)
    assert pair in asu_table
  assert not asu_table.contains(i_seq=1, j_seq=1, j_sym=10)
  pair = asu_mappings.make_trial_pair(i_seq=1, j_seq=1, j_sym=10)
  assert not pair in asu_table
  sym_table = asu_table.extract_pair_sym_table()
  assert sym_table.size() == asu_table.table().size()
  expected_sym_pairs = [
    [0, 0, '-y+1,-x+1,-z+1/2'],
    [0, 0, 'x,x-y+2,-z+1/2'],
    [0, 1, '-y+1,x-y+1,-z+1/2'],
    [1, 1, 'x,x-y+1,z'],
    [1, 1, '-y+1,-x+1,z'],
    [1, 2, '-x+1,-x+y+1,-z'],
    [2, 2, 'y,-x+y,-z']]
  ip = count()
  for i_seq,sym_dict in enumerate(sym_table):
    for j_seq,rt_mx_list in sym_dict.items():
      for rt_mx in rt_mx_list:
        if (0 or "--Verbose" in sys.argv[1:]):
          print str([i_seq, j_seq, str(rt_mx)]) + ","
        assert [i_seq, j_seq, str(rt_mx)] == expected_sym_pairs[ip.next()]
  asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  assert asu_table.add_pair_sym_table(sym_table=sym_table) is asu_table
  check_pair_asu_table(asu_table, expected_asu_pairs)
  asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  for p in expected_sym_pairs:
    assert asu_table.add_pair(
      i_seq=p[0], j_seq=p[1], rt_mx_ji=sgtbx.rt_mx(p[2])) is asu_table
  check_pair_asu_table(asu_table, expected_asu_pairs)
  term_table = crystal.coordination_sequences_simple(
    pair_asu_table=asu_table,
    n_shells=5)
  assert [list(terms) for terms in term_table] \
      == [[1,4,10,20,34,54], [1,4,10,20,34,53], [1,4,10,20,34,54]]

def exercise_symmetry():
  symmetry = crystal.ext.symmetry(
    unit_cell=uctbx.unit_cell([1,2,3,80,90,100]),
    space_group=sgtbx.space_group("P 2"))
  assert symmetry.unit_cell().is_similar_to(
    uctbx.unit_cell([1,2,3,80,90,100]))
  assert symmetry.space_group().order_z() == 2
  symmetry_cb = symmetry.change_basis(
    change_of_basis_op=sgtbx.change_of_basis_op("z,x,y"))
  assert symmetry_cb.unit_cell().is_similar_to(
    uctbx.unit_cell([3,1,2,100,80,90]))
  assert str(sgtbx.space_group_info(
    group=symmetry_cb.space_group())) == "P 2 1 1"

def run():
  exercise_direct_space_asu()
  exercise_pair_tables()
  exercise_symmetry()
  print "OK"

if (__name__ == "__main__"):
  run()
