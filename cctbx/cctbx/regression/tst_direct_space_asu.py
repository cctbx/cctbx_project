from cctbx import sgtbx
from cctbx.sgtbx.direct_space_asu import reference_table
from cctbx.sgtbx.direct_space_asu import facet_analysis
from cctbx.development import debug_utils
from cctbx.array_family import flex
from cctbx import matrix
from boost import rational
import sys

def exercise_cut_planes(cut_planes):
  for cut_plane in cut_planes:
    assert cut_plane.strip().is_inside(cut_plane.get_point_in_plane())

def exercise_direct_space_asu(space_group_info, n_grid=6):
  unit_cell = space_group_info.any_compatible_unit_cell(volume=1000)
  ref_asu = reference_table.get_asu(space_group_info.type().number())
  exercise_cut_planes(ref_asu.facets)
  inp_asu = space_group_info.direct_space_asu()
  exercise_cut_planes(inp_asu.facets)
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
  asu_tight = asu_with_metric.add_buffer()
  asu_buffer = asu_with_metric.add_buffer(thickness=2)
  asu_shrunk = asu_with_metric.add_buffer(relative_thickness=-1.e-6)
  vertices = facet_analysis.volume_vertices(inp_asu)
  for vertex in vertices:
    assert inp_asu.volume_only().is_inside(vertex)
  for vertex in vertices:
    assert asu_tight.is_inside(float(matrix.col(vertex)).elems)
  for vertex in vertices:
    assert asu_buffer.is_inside(float(matrix.col(vertex)).elems)
  for vertex in vertices:
    assert not asu_shrunk.is_inside(float(matrix.col(vertex)).elems)

def run_call_back(flags, space_group_info):
  exercise_direct_space_asu(space_group_info)
  if (space_group_info.group().n_ltr() != 1):
    exercise_direct_space_asu(space_group_info.primitive_setting())

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
