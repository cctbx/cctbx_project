from cctbx import crystal
import cctbx.crystal.direct_space_asu
from cctbx import sgtbx
from cctbx.sgtbx.direct_space_asu import reference_table
from cctbx.sgtbx.direct_space_asu import facet_analysis
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx import matrix
from boost import rational
import sys

def exercise_cut_planes(cut_planes):
  for cut_plane in cut_planes:
    assert cut_plane.strip().is_inside(cut_plane.get_point_in_plane())

def exercise_volume_vertices(asu, unit_cell):
  volume_asu = asu.volume_only()
  asu_volume_vertices = asu.volume_vertices()
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
  for mappings in asu_mappings.mappings():
    assert asu_mappings.asu().is_inside(mappings[0].mapped_site())
    for mapping in mappings:
      assert asu_mappings.asu_buffer().is_inside(mapping.mapped_site())

def exercise_all(space_group_info):
  exercise_float_asu(space_group_info)
  exercise_asu_mappings(space_group_info)

def run_call_back(flags, space_group_info):
  exercise_all(space_group_info)
  if (space_group_info.group().n_ltr() != 1):
    exercise_all(space_group_info.primitive_setting())

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
