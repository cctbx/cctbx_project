from cctbx import crystal
import cctbx.crystal.direct_space_asu
from cctbx import xray
from cctbx import sgtbx
from cctbx.sgtbx.direct_space_asu import reference_table
from cctbx.sgtbx.direct_space_asu import facet_analysis
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.itertbx import count
from libtbx.test_utils import approx_equal
from boost import rational
import random
import copy
import sys

def exercise_cut_planes(cut_planes):
  for cut_plane in cut_planes:
    assert cut_plane.strip().is_inside(cut_plane.get_point_in_plane())

def exercise_volume_vertices(asu, unit_cell):
  volume_asu = asu.volume_only()
  asu_volume_vertices = asu.volume_vertices()
  for box_min,box_max in zip(asu.box_min(volume_vertices=asu_volume_vertices),
                             asu.box_max(volume_vertices=asu_volume_vertices)):
    assert box_min < box_max
  float_asu = asu.add_buffer(unit_cell=unit_cell, thickness=0)
  asu_shrunk = float_asu.add_buffer(relative_thickness=-1.e-5)
  mcs = asu.define_metric(unit_cell).minimum_covering_sphere()
  center = matrix.col(mcs.center())
  radius = mcs.radius()
  n_near_sphere_surface = 0
  for vertex in asu_volume_vertices:
    assert volume_asu.is_inside(vertex)
    float_vertex = [float(e) for e in vertex]
    assert float_asu.is_inside(float_vertex)
    assert not asu_shrunk.is_inside(float_vertex)
    r = abs(matrix.col(unit_cell.orthogonalize(float_vertex)) - center)
    assert r < radius + 1.e-5
    if (radius - r < radius * 1.e-2):
      n_near_sphere_surface += 1
  assert n_near_sphere_surface >= 2
  float_asu_volume_vertices = float_asu.volume_vertices()
  assert len(float_asu_volume_vertices) >= len(asu_volume_vertices)
  m_near_sphere_surface = 0
  for vertex in float_asu_volume_vertices:
    assert float_asu.is_inside(point=vertex)
    assert not asu_shrunk.is_inside(vertex)
    r = abs(matrix.col(unit_cell.orthogonalize(vertex)) - center)
    assert r < radius + 1.e-5
    if (radius - r < radius * 1.e-2):
      m_near_sphere_surface += 1
  assert m_near_sphere_surface >= n_near_sphere_surface
  line_asu = copy.copy(asu)
  line_asu.add_planes([(0,0,1),(1,1,1)], both_directions=0001)
  assert len(line_asu.facets) == len(asu.facets) + 4
  assert line_asu.facets[-2].n == (-line_asu.facets[-1]).n

def exercise_float_asu(space_group_info, n_grid=6):
  unit_cell = space_group_info.any_compatible_unit_cell(volume=1000)
  ref_asu = reference_table.get_asu(space_group_info.type().number())
  exercise_cut_planes(ref_asu.facets)
  inp_asu = space_group_info.direct_space_asu()
  exercise_cut_planes(inp_asu.facets)
  exercise_volume_vertices(inp_asu, unit_cell)
  assert sgtbx.space_group(inp_asu.hall_symbol) == space_group_info.group()
  float_asu = inp_asu.add_buffer(unit_cell=unit_cell, thickness=0.001)
  cb_mx_ref_inp = space_group_info.type().cb_op().c_inv().as_rational()
  n = n_grid
  for ref_n in flex.nested_loop((-n/2,-n/2,-n/2),(n,n,n),00000):
    # check correctness of space_group_info.direct_space_asu()
    ref_r = matrix.col([rational.int(g,n) for g in ref_n])
    inp_r = cb_mx_ref_inp * ref_r
    assert ref_asu.is_inside(ref_r.elems) == inp_asu.is_inside(inp_r.elems)
    # check correctness of cut_plane.add_buffer()
    inp_r = inp_r.elems
    inp_f = [float(r) for r in inp_r]
    for facet in inp_asu.facets:
      r_cut = facet.strip()
      r_inside = r_cut.is_inside(inp_r)
      for buffer_thickness in [0.001, 1, -1][:1]:
        f_cut = facet.as_float_cut_plane().add_buffer(
          unit_cell=unit_cell,
          thickness=buffer_thickness)
        f_inside = f_cut.is_inside(inp_f)
        if (buffer_thickness < 0):
          if (r_inside != f_inside):
            assert r_inside
        elif (buffer_thickness < 0.01):
          assert r_inside == f_inside
        elif (r_inside != f_inside):
          assert f_inside
    # check correctness of float_asu.add_buffer()
    assert float_asu.is_inside(inp_f) == inp_asu.volume_only().is_inside(inp_r)
  asu_with_metric = inp_asu.define_metric(unit_cell)
  assert asu_with_metric.hall_symbol is inp_asu.hall_symbol
  assert len(asu_with_metric.facets) == len(inp_asu.facets)
  assert asu_with_metric.unit_cell is unit_cell
  asu_tight = asu_with_metric.as_float_asu()
  asu_buffer = asu_with_metric.add_buffer(thickness=2)
  asu_shrunk = asu_with_metric.add_buffer(relative_thickness=-1.e-5)
  vertices = facet_analysis.volume_vertices(inp_asu)
  for vertex in vertices:
    assert inp_asu.volume_only().is_inside(vertex)
  for vertex in vertices:
    assert asu_tight.is_inside(float(matrix.col(vertex)).elems)
  for vertex in vertices:
    assert asu_buffer.is_inside(float(matrix.col(vertex)).elems)
  for vertex in vertices:
    assert not asu_shrunk.is_inside(float(matrix.col(vertex)).elems)

def exercise_asu_mappings(space_group_info, n_elements=10):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=1000,
    min_distance=3.,
    general_positions_only=00000)
  asu_mappings = crystal.direct_space_asu.asu_mappings(
    space_group=structure.space_group(),
    asu=structure.direct_space_asu().as_float_asu(),
    buffer_thickness=4)
  asu_mappings.reserve(structure.scatterers().size())
  for scatterer in structure.scatterers():
    asu_mappings.process(original_site=scatterer.site)
  assert asu_mappings.mappings().size() == structure.scatterers().size()
  frac = asu_mappings.unit_cell().fractionalize
  for mappings in asu_mappings.mappings():
    assert asu_mappings.asu().is_inside(frac(mappings[0].mapped_site()))
    for mapping in mappings:
      assert asu_mappings.asu_buffer().is_inside(frac(mapping.mapped_site()))

def exercise_neighbors_pair_generators(space_group_info, n_elements=10,
                                       verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=1000,
    min_distance=3.,
    general_positions_only=00000)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
    print
  for buffer_thickness in [1.e-5, 2, 4]:
    asu_mappings = crystal.direct_space_asu.asu_mappings(
      space_group=structure.space_group(),
      asu=structure.direct_space_asu().as_float_asu(),
      buffer_thickness=buffer_thickness)
    asu_mappings.reserve(structure.scatterers().size())
    for scatterer in structure.scatterers():
      asu_mappings.process(scatterer.site)
    array_of_array_of_mappings = asu_mappings.mappings()
    for full_matrix in [00000,0001]:
      pair_list = []
      for i_seq in xrange(array_of_array_of_mappings.size()):
        site_0 = matrix.col(array_of_array_of_mappings[i_seq][0].mapped_site())
        if (full_matrix):
          j_seq_start = 0
        else:
          j_seq_start = i_seq
        for j_seq in xrange(j_seq_start, array_of_array_of_mappings.size()):
          array_of_mappings = array_of_array_of_mappings[j_seq]
          if (i_seq == j_seq):
            j_start = 1
          else:
            j_start = 0
          for j_sym in xrange(j_start, len(array_of_mappings)):
            site_1 = matrix.col(array_of_mappings[j_sym].mapped_site())
            dist_sq = (site_1-site_0).norm()
            pair_list.append((i_seq,j_seq,j_sym,dist_sq))
      for pair_generator_type in [crystal.neighbors_simple_pair_generator,
                                  crystal.neighbors_fast_pair_generator]:
        pair_generator = pair_generator_type(
          asu_mappings=asu_mappings,
          distance_cutoff=1.e6,
          full_matrix=full_matrix)
        assert pair_generator.full_matrix() == full_matrix
        mps = asu_mappings.mappings()
        sc = structure.scatterers()
        uc = structure.unit_cell()
        sg = structure.space_group()
        col = matrix.col
        for pair_direct,pair in zip(pair_list, pair_generator):
          assert pair_direct[:3] == (pair.i_seq, pair.j_seq, pair.j_sym)
          assert approx_equal(pair_direct[3], pair.dist_sq)
          mp_i = mps[pair.i_seq][0]
          mp_j = mps[pair.j_seq][pair.j_sym]
          site_i = (  col(sg(mp_i.i_sym_op())*sc[pair.i_seq].site)
                    + col(mp_i.unit_shifts())).elems
          site_j = (  col(sg(mp_j.i_sym_op())*sc[pair.j_seq].site)
                    + col(mp_j.unit_shifts())).elems
          assert approx_equal(uc.orthogonalize(site_i), mp_i.mapped_site())
          assert approx_equal(uc.orthogonalize(site_j), mp_j.mapped_site())
          assert approx_equal(uc.distance(site_i, site_j)**2, pair.dist_sq)
          site_i = asu_mappings.get_rt_mx(pair.i_seq, 0) \
                 * sc[pair.i_seq].site
          site_j = asu_mappings.get_rt_mx(pair.j_seq, pair.j_sym) \
                 * sc[pair.j_seq].site
          assert approx_equal(uc.orthogonalize(site_i), mp_i.mapped_site())
          assert approx_equal(uc.orthogonalize(site_j), mp_j.mapped_site())
          j_frac = uc.fractionalize(mp_j.mapped_site())
          assert approx_equal(
            asu_mappings.get_rt_mx(pair.j_seq, pair.j_sym).inverse() * j_frac,
            sc[pair.j_seq].site)
          assert approx_equal(
            asu_mappings.map_moved_site_to_asu(
              moved_original_site
                =asu_mappings.unit_cell().orthogonalize(sc[pair.i_seq].site),
              i_seq=pair.i_seq,
              i_sym=0),
            mp_i.mapped_site())
          assert approx_equal(
            asu_mappings.map_moved_site_to_asu(
              moved_original_site
                =asu_mappings.unit_cell().orthogonalize(sc[pair.j_seq].site),
              i_seq=pair.j_seq,
              i_sym=pair.j_sym),
            mp_j.mapped_site())

def exercise_is_symmetry_interaction():
  for space_group_symbol in ["P1", "P4"]:
    for shifts in flex.nested_loop((-2,-2,-2),(2,2,2),00000):
      shifts = matrix.col(shifts)
      structure = xray.structure(
        crystal_symmetry=crystal.symmetry(
          unit_cell=(10,10,20,90,90,90),
          space_group_symbol=space_group_symbol),
        scatterers=flex.xray_scatterer([
          xray.scatterer(label="O", site=shifts+matrix.col((0,0,0))),
          xray.scatterer(label="N", site=shifts+matrix.col((0.5,0.5,0))),
          xray.scatterer(label="C", site=shifts+matrix.col((0.25,0.25,0)))]))
      asu_mappings = structure.asu_mappings(buffer_thickness=7)
      pair_generator = crystal.neighbors_simple_pair_generator(
        asu_mappings=asu_mappings,
        distance_cutoff=7)
      direct_interactions = {}
      for i_pair,pair in zip(count(),pair_generator):
        if (not asu_mappings.is_symmetry_interaction(pair)):
          key = (pair.i_seq,pair.j_seq)
          assert direct_interactions.get(key, None) is None
          direct_interactions[key] = 1
      assert len(direct_interactions) == 2
      assert direct_interactions[(0,2)] == 1
      assert direct_interactions[(1,2)] == 1

def exercise_non_crystallographic_asu_mappings():
  asu_mappings = crystal.direct_space_asu.non_crystallographic_asu_mappings(
    sites_cart=flex.vec3_double())
  assert approx_equal(asu_mappings.unit_cell().parameters(), (1,1,1,90,90,90))
  asu_mappings = crystal.direct_space_asu.non_crystallographic_asu_mappings(
    sites_cart=flex.vec3_double([(2,-3,4)]))
  assert approx_equal(asu_mappings.unit_cell().parameters(), (1,1,1,90,90,90))
  for i_trial in xrange(10):
    offs = random.uniform(-100,100)
    sites_cart = flex.vec3_double()
    for i_site in xrange(10):
      sites_cart.append([random.uniform(-10+offs,10+offs) for i in xrange(3)])
    asu_mappings = crystal.direct_space_asu.non_crystallographic_asu_mappings(
      sites_cart=sites_cart)
    for site,asu_mapping in zip(sites_cart,asu_mappings.mappings()):
      assert len(asu_mapping) == 1
      assert asu_mapping[0].i_sym_op() == 0
      assert asu_mapping[0].unit_shifts() == (0,0,0)
      assert approx_equal(asu_mapping[0].mapped_site(), site)
    assert approx_equal(asu_mappings.mapped_sites_min(), sites_cart.min())
    assert approx_equal(asu_mappings.mapped_sites_max(), sites_cart.max())

def exercise_all(flags, space_group_info):
  exercise_float_asu(space_group_info)
  exercise_asu_mappings(space_group_info)
  exercise_neighbors_pair_generators(
    space_group_info, verbose=flags.Verbose)

def run_call_back(flags, space_group_info):
  exercise_all(flags, space_group_info)
  if (space_group_info.group().n_ltr() != 1):
    exercise_all(flags, space_group_info.primitive_setting())

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  exercise_is_symmetry_interaction()
  exercise_non_crystallographic_asu_mappings()
  print "OK"

if (__name__ == "__main__"):
  run()
