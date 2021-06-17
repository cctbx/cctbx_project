from __future__ import absolute_import, division, print_function
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
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from boost_adaptbx.boost import rational
import random
import copy
import sys
from six.moves import range
from six.moves import zip

def exercise_reference_table():
  for space_group_number in range(1,230+1):
    space_group_info = sgtbx.space_group_info(number=space_group_number)
    asu = reference_table.get_asu(space_group_number)
    assert sgtbx.space_group(asu.hall_symbol) == space_group_info.group()
  #
  for space_group_number in range(1,230+1):
    asu = reference_table.get_asu(space_group_number)
    n_long_cuts = 0
    for cut in asu.cuts:
      s = str(cut)
      have_long_cut = (s.find("cut") >= 0)
      if (space_group_number != 213):
        assert not have_long_cut
      elif (have_long_cut):
        n_long_cuts += 1
    if (space_group_number == 213):
      assert n_long_cuts == 4 # done with change_basis

def exercise_cut_planes(cut_planes):
  for cut_plane in cut_planes:
    assert cut_plane.strip().is_inside(cut_plane.get_point_in_plane())

def exercise_shape_vertices(asu, unit_cell):
  volume_asu = asu.shape_only()
  asu_shape_vertices = asu.shape_vertices()
  assert len(asu_shape_vertices) >= 4
  facet_analysis_shape_vertices = facet_analysis.shape_vertices(asu)
  assert len(asu_shape_vertices) == len(facet_analysis_shape_vertices)
  asu_shape_vertices.sort()
  facet_analysis_shape_vertices.sort()
  assert asu_shape_vertices == facet_analysis_shape_vertices
  for box_min,box_max in zip(asu.box_min(shape_vertices=asu_shape_vertices),
                             asu.box_max(shape_vertices=asu_shape_vertices)):
    assert box_min < box_max
  if (unit_cell is None):
    return
  float_asu = asu.add_buffer(unit_cell=unit_cell, thickness=0)
  asu_shrunk = float_asu.add_buffer(relative_thickness=-1.e-5)
  mcs = asu.define_metric(unit_cell).minimum_covering_sphere()
  center = matrix.col(mcs.center())
  radius = mcs.radius()
  n_near_sphere_surface = 0
  for vertex in asu_shape_vertices:
    assert volume_asu.is_inside(vertex)
    float_vertex = [float(e) for e in vertex]
    assert float_asu.is_inside(float_vertex)
    assert not asu_shrunk.is_inside(float_vertex)
    r = abs(matrix.col(unit_cell.orthogonalize(float_vertex)) - center)
    assert r < radius + 1.e-5
    if (radius - r < radius * 1.e-2):
      n_near_sphere_surface += 1
  assert n_near_sphere_surface >= 2
  float_asu_shape_vertices = float_asu.shape_vertices()
  assert len(float_asu_shape_vertices) >= len(asu_shape_vertices)
  m_near_sphere_surface = 0
  for vertex in float_asu_shape_vertices:
    assert float_asu.is_inside(point=vertex)
    assert not asu_shrunk.is_inside(vertex)
    r = abs(matrix.col(unit_cell.orthogonalize(vertex)) - center)
    assert r < radius + 1.e-5
    if (radius - r < radius * 1.e-2):
      m_near_sphere_surface += 1
  assert m_near_sphere_surface >= n_near_sphere_surface
  line_asu = copy.copy(asu)
  line_asu.add_planes([(0,0,1),(1,1,1)], both_directions=True)
  assert len(line_asu.cuts) == len(asu.cuts) + 4
  assert line_asu.cuts[-2].n == (-line_asu.cuts[-1]).n

def exercise_float_asu(space_group_info, n_grid=6):
  unit_cell = space_group_info.any_compatible_unit_cell(volume=1000)
  ref_asu = reference_table.get_asu(space_group_info.type().number())
  exercise_cut_planes(ref_asu.cuts)
  inp_asu = space_group_info.direct_space_asu()
  assert sgtbx.space_group(inp_asu.hall_symbol) == space_group_info.group()
  exercise_cut_planes(inp_asu.cuts)
  exercise_shape_vertices(inp_asu, unit_cell)
  float_asu = inp_asu.add_buffer(unit_cell=unit_cell, thickness=0.001)
  cb_mx_ref_inp = space_group_info.type().cb_op().c_inv().as_rational()
  n = n_grid
  for ref_n in flex.nested_loop((-n//2,-n//2,-n//2),(n,n,n),False):
    # check correctness of space_group_info.direct_space_asu()
    ref_r = matrix.col([rational.int(g,n) for g in ref_n])
    inp_r = cb_mx_ref_inp * ref_r
    assert ref_asu.is_inside(ref_r.elems) == inp_asu.is_inside(inp_r.elems)
    # check correctness of cut_plane.add_buffer()
    inp_r = inp_r.elems
    inp_f = [float(r) for r in inp_r]
    for cut in inp_asu.cuts:
      r_cut = cut.strip()
      r_inside = r_cut.is_inside(inp_r)
      for buffer_thickness in [0.001, 1, -1][:1]:
        f_cut = cut.as_float_cut_plane().add_buffer(
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
    assert float_asu.is_inside(inp_f) == inp_asu.shape_only().is_inside(inp_r)
  asu_with_metric = inp_asu.define_metric(unit_cell)
  assert asu_with_metric.hall_symbol is inp_asu.hall_symbol
  assert len(asu_with_metric.cuts) == len(inp_asu.cuts)
  assert asu_with_metric.unit_cell is unit_cell
  asu_tight = asu_with_metric.as_float_asu()
  asu_buffer = asu_with_metric.add_buffer(thickness=2)
  asu_shrunk = asu_with_metric.add_buffer(relative_thickness=-1.e-5)
  vertices = facet_analysis.shape_vertices(inp_asu)
  for vertex in vertices:
    assert inp_asu.shape_only().is_inside(vertex)
  for vertex in vertices:
    assert asu_tight.is_inside(matrix.col(vertex).as_float().elems)
  for vertex in vertices:
    assert asu_buffer.is_inside(matrix.col(vertex).as_float().elems)
  for vertex in vertices:
    assert not asu_shrunk.is_inside(matrix.col(vertex).as_float().elems)

def exercise_asu_mappings(space_group_info, n_elements=10):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=1000,
    min_distance=3.,
    general_positions_only=False)
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

def exercise_neighbors_pair_generators(structure, verbose=0):
  if (0 or verbose):
    structure.show_summary().show_scatterers()
    print()
  for buffer_thickness in [1.e-5, 2, 4]:
    asu_mappings = crystal.direct_space_asu.asu_mappings(
      space_group=structure.space_group(),
      asu=structure.direct_space_asu().as_float_asu(),
      buffer_thickness=buffer_thickness)
    asu_mappings.reserve(structure.scatterers().size())
    for scatterer in structure.scatterers():
      asu_mappings.process(scatterer.site)
    array_of_array_of_mappings = asu_mappings.mappings()
    for minimal in [False, True]:
      pair_list = []
      for i_seq in range(array_of_array_of_mappings.size()):
        array_of_mappings_i = array_of_array_of_mappings[i_seq]
        site_0 = matrix.col(array_of_mappings_i[0].mapped_site())
        for j_seq in range(i_seq, array_of_array_of_mappings.size()):
          array_of_mappings_j = array_of_array_of_mappings[j_seq]
          if (i_seq == j_seq):
            j_sym_start = 1
          else:
            j_sym_start = 0
          for j_sym in range(j_sym_start, len(array_of_mappings_j)):
            site_1 = matrix.col(array_of_mappings_j[j_sym].mapped_site())
            dist_sq = (site_1-site_0).norm_sq()
            pair_list.append((i_seq,j_seq,j_sym,dist_sq))
          if (not minimal and i_seq != j_seq):
            site_1 = matrix.col(array_of_mappings_j[0].mapped_site())
            for i_sym in range(1, len(array_of_mappings_i)):
              site_2 = matrix.col(array_of_mappings_i[i_sym].mapped_site())
              dist_sq = (site_2-site_1).norm_sq()
              pair_list.append((j_seq,i_seq,i_sym,dist_sq))
      for pair_generator_type in [crystal.neighbors_simple_pair_generator,
                                  crystal.neighbors_fast_pair_generator]:
        pair_generator = pair_generator_type(
          asu_mappings=asu_mappings,
          distance_cutoff=1.e6,
          minimal=minimal)
        mps = asu_mappings.mappings()
        sc = structure.scatterers()
        uc = structure.unit_cell()
        sg = structure.space_group()
        col = matrix.col
        pair_dict = dict([(pair[:3], pair[3]) for pair in pair_list])
        for pair in pair_generator:
          assert pair.is_active(minimal=minimal)
          key = (pair.i_seq, pair.j_seq, pair.j_sym)
          assert approx_equal(pair_dict[key], pair.dist_sq)
          del pair_dict[key]
          mp_i = mps[pair.i_seq][0]
          mp_j = mps[pair.j_seq][pair.j_sym]
          site_i = (  col(sg(mp_i.i_sym_op())*sc[pair.i_seq].site)
                    + col(mp_i.unit_shifts())).elems
          site_j = (  col(sg(mp_j.i_sym_op())*sc[pair.j_seq].site)
                    + col(mp_j.unit_shifts())).elems
          assert approx_equal(uc.orthogonalize(site_i), mp_i.mapped_site())
          assert approx_equal(uc.orthogonalize(site_j), mp_j.mapped_site())
          assert approx_equal(uc.distance(site_i, site_j)**2, pair.dist_sq)
          site_i = asu_mappings.get_rt_mx_i(pair=pair) * sc[pair.i_seq].site
          site_j = asu_mappings.get_rt_mx_j(pair=pair) * sc[pair.j_seq].site
          assert approx_equal(uc.orthogonalize(site_i), mp_i.mapped_site())
          assert approx_equal(uc.orthogonalize(site_j), mp_j.mapped_site())
          j_frac = uc.fractionalize(mp_j.mapped_site())
          assert approx_equal(
            asu_mappings.get_rt_mx_j(pair=pair).inverse() * j_frac,
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
        assert len(pair_dict) == 0

def asu_mappings_is_simple_interaction_emulation(asu_mappings, pair):
  is_special_position = asu_mappings.site_symmetry_table().is_special_position
  if (is_special_position(i_seq=pair.i_seq)): return False
  if (is_special_position(i_seq=pair.j_seq)): return False
  return asu_mappings.get_rt_mx_i(pair) == asu_mappings.get_rt_mx_j(pair)

def exercise_is_simple_interaction():
  for space_group_symbol in ["P1", "P41"]:
    for shifts in flex.nested_loop((-2,-2,-2),(2,2,2),False):
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
      simple_interactions = {}
      for i_pair,pair in enumerate(pair_generator):
        if (asu_mappings.is_simple_interaction(pair)):
          assert asu_mappings_is_simple_interaction_emulation(
            asu_mappings, pair)
          key = (pair.i_seq,pair.j_seq)
          assert simple_interactions.get(key, None) is None
          simple_interactions[key] = 1
        else:
          assert not asu_mappings_is_simple_interaction_emulation(
            asu_mappings, pair)
      assert len(simple_interactions) == 2
      assert simple_interactions[(0,2)] == 1
      assert simple_interactions[(1,2)] == 1

def exercise_non_crystallographic_asu_mappings():
  asu_mappings = crystal.direct_space_asu.non_crystallographic_asu_mappings(
    sites_cart=flex.vec3_double())
  assert approx_equal(asu_mappings.unit_cell().parameters(), (1,1,1,90,90,90))
  asu_mappings = crystal.direct_space_asu.non_crystallographic_asu_mappings(
    sites_cart=flex.vec3_double([(2,-3,4)]))
  assert approx_equal(asu_mappings.unit_cell().parameters(), (1,1,1,90,90,90))
  for i_trial in range(10):
    offs = random.uniform(-100,100)
    sites_cart = flex.vec3_double()
    for i_site in range(10):
      sites_cart.append([random.uniform(-10+offs,10+offs) for i in range(3)])
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
    structure = random_structure.xray_structure(
      space_group_info,
      elements=["Si"]*5,
      volume_per_atom=100,
      min_distance=3.,
      general_positions_only=False),
    verbose=flags.Verbose)

def run_call_back(flags, space_group_info):
  exercise_all(flags, space_group_info)
  if (space_group_info.group().n_ltr() != 1):
    exercise_all(flags, space_group_info.primitive_setting())

def run():
  exercise_reference_table()
  for space_group_number in range(1,230+1):
    asu = reference_table.get_asu(space_group_number)
    exercise_shape_vertices(asu=asu, unit_cell=None)
  debug_utils.parse_options_loop_space_groups(
    sys.argv[1:], run_call_back, show_cpu_times=False)
  exercise_is_simple_interaction()
  exercise_non_crystallographic_asu_mappings()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
