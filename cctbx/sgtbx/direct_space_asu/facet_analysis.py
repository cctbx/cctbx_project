from __future__ import absolute_import, division, print_function
from operator import itemgetter
from cctbx.sgtbx.direct_space_asu import cut_plane
from cctbx.array_family import flex
from scitbx import matrix
import scitbx.math
from boost_adaptbx.boost import rational
from six.moves import range
from six.moves import zip

def intersection(cuts):
  assert len(cuts) == 3
  m = flex.int()
  t = flex.int()
  denominator = 1
  for cut in cuts:
    denominator = cut.lcm_of_denominators(start_lcm=denominator)
  for cut in cuts:
    for e in cut.n: m.append(int(e * denominator))
    t.append(-int(cut.c * denominator))
  m.reshape(flex.grid(3,3))
  t.reshape(flex.grid(3,1))
  r = scitbx.math.row_echelon_form_t(m, t)
  assert r in (2,3)
  if (r != 3): return None
  t.reshape(flex.grid(3))
  sol = flex.int(3)
  d = scitbx.math.row_echelon_back_substitution_int(m, t, sol)
  assert d > 0
  return tuple([rational.int(s,d) for s in sol])

def eliminate_outside_vertices(other_cuts, vertices):
  result = {}
  for vertex in vertices.keys():
    for cut in other_cuts:
      if (cut.evaluate(vertex) < 0):
        break
    else:
      result[vertex] = vertices[vertex]
  return result

def polygon_vertices(pivot, other_cuts):
  result = {}
  n = len(other_cuts)
  for i in range(n-1):
    for j in range(i+1,n):
      vertex = intersection(cuts=(pivot, other_cuts[i], other_cuts[j]))
      if (vertex is not None):
        result.setdefault(vertex, []).append((i,j))
  return eliminate_outside_vertices(other_cuts=other_cuts, vertices=result)

def depth1_cuts(expr, op_and="&", op_or="|"):
  if (isinstance(expr, cut_plane.cut)):
    return [[expr]]
  assert isinstance(expr, cut_plane.cut_expression)
  lhs = depth1_cuts(expr.lhs, op_and, op_or)
  rhs = depth1_cuts(expr.rhs, op_and, op_or)
  if (expr.op == op_or):
    return lhs + rhs
  assert expr.op == op_and
  result = []
  for l in lhs:
    for r in rhs:
      result.append(l + r)
  return result

def face_vertices(asu, i_pivot):
  pivot = asu.cuts[i_pivot]
  other_cuts = asu.cuts[:i_pivot] + asu.cuts[i_pivot+1:]
  if (not pivot.has_cuts()):
    return [(polygon_vertices(pivot, other_cuts), pivot.inclusive)]
  cuts_inside = depth1_cuts(pivot.cut_expr)
  cuts_outside = depth1_cuts(pivot.cut_expr, "|", "&")
  result = []
  for cut in cuts_inside:
    result.append((polygon_vertices(pivot, other_cuts + cut), True))
  for cut in cuts_outside:
    cut = [-f for f in cut]
    result.append((polygon_vertices(pivot, other_cuts + cut), False))
  return result

def sense_of_polygon(cut_n, polygon):
  assert len(polygon) >= 3
  row = matrix.row
  p00 = row(polygon[0][0])
  rows = [row(cut_n), row(polygon[-1][0])-p00, row(polygon[1][0])-p00]
  d = matrix.col(rows).resolve_partitions().determinant()
  assert d != 0
  return d

def trace_polygon(
      polygon,
      begin_cut_index,
      end_cut_index,
      unused_face_vertices):
  i_vertex = -1
  for vertex,list_of_other_cut_indices in unused_face_vertices:
    i_vertex += 1
    for cut_indices in list_of_other_cut_indices:
      if (cut_indices[0] != end_cut_index):
        cut_indices = list(cut_indices)
        cut_indices.reverse()
        cut_indices = tuple(cut_indices)
      if (cut_indices[0] == end_cut_index):
        full_polygon = trace_polygon(
          polygon=polygon+[(vertex, cut_indices)],
          begin_cut_index=begin_cut_index,
          end_cut_index=cut_indices[1],
          unused_face_vertices=
              unused_face_vertices[:i_vertex]
            + unused_face_vertices[i_vertex+1:])
        if (full_polygon is not None):
          return full_polygon
  if (len(unused_face_vertices) > 0): return None
  if (end_cut_index != begin_cut_index): return None
  return polygon

def face_polygons(asu, i_pivot):
  vertex_dicts = face_vertices(asu=asu, i_pivot=i_pivot)
  result = []
  for vertex_dict,flag in vertex_dicts:
    # FIXME ordering of lists in items changes in python2/3
    vertex_list = list(vertex_dict.items())
    vertex,list_of_other_cut_indices = vertex_list[0]
    for other_cut_indices in list_of_other_cut_indices:
      polygon = trace_polygon([(vertex, other_cut_indices)],
                              other_cut_indices[0],
                              other_cut_indices[1],
                              vertex_list[1:])
      if (polygon is not None): break
    if (polygon is not None):
      if (sense_of_polygon(asu.cuts[i_pivot].n, polygon) < 0):
        polygon.reverse()
      result.append((polygon, flag))
  return result

def asu_polygons(asu):
  result = []
  for i_pivot in range(len(asu.cuts)):
    result.append(face_polygons(asu=asu, i_pivot=i_pivot))
  verify_asu_polygons(asu=asu, list_of_polygons=result)
  return result

def extract_polygon_vertices(list_of_polygons):
  result = {}
  for polygons in list_of_polygons:
    for polygon,inclusive_flag in polygons:
      for vertex,cut_indices in polygon:
        result[vertex] = 1
  return list(result.keys())

def shape_vertices(asu):
  return extract_polygon_vertices(asu_polygons(asu.shape_only()))

def line_sample_point(a, b, f, gridding):
  return [a[i]+rational.int(f,gridding)*(b[i]-a[i]) for i in range(3)]

def verify_asu_polygons(asu, list_of_polygons, gridding=13):
  for polygons in list_of_polygons:
    for polygon,inclusive_flag in polygons:
      n = len(polygon)
      for i in range(n-1):
        for j in range(i+2,min(n,n+i-1)):
          a, b = polygon[i][0], polygon[j][0]
          for f in range(1, gridding):
            x = line_sample_point(a, b, f, gridding)
            assert asu.is_inside(x) == inclusive_flag

def collect_cuts(expr):
  if (expr is None): return []
  if (isinstance(expr, cut_plane.cut)):
    return [expr]
  assert isinstance(expr, cut_plane.cut_expression)
  lhs = collect_cuts(expr.lhs)
  rhs = collect_cuts(expr.rhs)
  for cut in rhs:
    if (not cut in lhs):
      lhs.append(cut)
  return lhs

def is_one_of(cut_list, addl_cut):
  for cut in cut_list:
    if (cut.n == addl_cut.n and cut.c == addl_cut.c):
      return True
  return False

def all_cut_points(asu):
  result = {}
  for pivot in asu.cuts:
    for first_cut in collect_cuts(pivot.cut_expr):
      for second_cut in collect_cuts(first_cut.cut_expr):
        all_cuts = (pivot, first_cut, second_cut)
        if (second_cut.has_cuts()):
          assert isinstance(second_cut.cut_expr, cut_plane.cut)
          assert second_cut.cut_expr.inclusive == False
          assert is_one_of(all_cuts, second_cut)
        point = intersection(cuts=all_cuts)
        assert point is not None
        if (asu.is_inside(point, shape_only=True)):
          result[point] = 1
  return list(result.keys())

def get_edge_vertices(list_of_polygons):
  result = {}
  for polygons in list_of_polygons:
    for polygon,inclusive_flag in polygons:
      n = len(polygon)
      for i in range(n):
        j = (i+1) % n
        v1 = polygon[i][0]
        v2 = polygon[j][0]
        key = (v1,v2)
        if (key not in result):
          key = (v2,v1)
        result[key] = 1
  return list(result.keys())

def edge_position(edge_end_points, other_point):
  a = edge_end_points[0]
  b = edge_end_points[1]
  x = other_point
  f = None
  for i in range(3):
    d_i = b[i] - a[i]
    if (f is None):
      if (d_i != 0):
        f = (x[i] - a[i]) / d_i
      else:
        if (x[i] != a[i]): return None
    else:
      y_i = a[i] + f * d_i
      if (y_i != x[i]): return None
  return f

class edge_with_cut_points(object):

  __slots__ = ["end_points", "cut_points"]

  def __init__(self, end_points, cut_points):
    self.end_points = tuple(end_points)
    self.cut_points = tuple(cut_points)

  def show_points(self):
    print("end points:", self.end_points)
    print("cut points:", self.cut_points)

  def cut_point_positions(self):
    return [edge_position(self.end_points, point) for point in self.cut_points]

  def sorted_cut_points(self):
    packed = list(zip(self.cut_points, self.cut_point_positions()))
    packed.sort(key=itemgetter(1))
    return [p[0] for p in packed]

  def sort_cut_points(self):
    return edge_with_cut_points(self.end_points, self.sorted_cut_points())

  def all_points(self):
    return (self.end_points[0],) + self.cut_points + (self.end_points[1],)

  def join_point(self, point):
    f = edge_position(self.end_points, point)
    if (f is None or f == 0 or f == 1 or point in self.cut_points):
      return f, self
    if (f < 0):
      return f, edge_with_cut_points(
        (point, self.end_points[1]),
        (self.end_points[0],) + self.cut_points)
    elif (f > 1):
      return f, edge_with_cut_points(
        (self.end_points[0], point),
        self.cut_points + (self.end_points[1],))
    return f, edge_with_cut_points(
      self.end_points,
      self.cut_points + (point,))

  def join_edge(self, edge):
    result = self
    linear_dependent = True
    for point in edge:
      f, result = result.join_point(point)
      if (f is None): linear_dependent = False
    return linear_dependent, result

class edge_segment(object):

  __slots__ = ["vertex", "vertex_inclusive_flag", "edge_inclusive_flag"]

  def __init__(self, vertex, vertex_inclusive_flag, edge_inclusive_flag):
    self.vertex = vertex
    self.vertex_inclusive_flag = vertex_inclusive_flag
    self.edge_inclusive_flag = edge_inclusive_flag

class consolidated_edges_with_cut_points(object):

  __slots__ = ["asu", "list"]

  def __init__(self, asu, list_of_polygons):
    self.asu = asu
    self.list = []
    addl_cut_points = all_cut_points(asu)
    for edge in get_edge_vertices(list_of_polygons):
      self.add(edge, addl_cut_points)
    self.list = [ec.sort_cut_points() for ec in self.list]

  def add(self, edge, addl_cut_points):
    for i in range(len(self.list)):
      linear_dependent, self.list[i] = self.list[i].join_edge(edge)
      if (linear_dependent): return
    edge_and_cuts = edge_with_cut_points(edge, ())
    for listed_edge_and_cuts in self.list:
      for point in listed_edge_and_cuts.all_points():
        f, edge_and_cuts = edge_and_cuts.join_point(point)
    for point in addl_cut_points:
      f, edge_and_cuts = edge_and_cuts.join_point(point)
    self.list.append(edge_and_cuts)

  def get_segments(self, edge_and_cuts):
    points = edge_and_cuts.all_points()
    result = []
    for i in range(len(points)-1):
      vertices = (points[i], points[i+1])
      mid_point = (
        (matrix.col(vertices[0]) + matrix.col(vertices[1])) / 2).elems
      result.append(edge_segment(
        vertex=vertices[0],
        vertex_inclusive_flag=self.asu.is_inside(vertices[0]),
        edge_inclusive_flag=self.asu.is_inside(mid_point)))
    result.append(edge_segment(
      vertex=points[-1],
      vertex_inclusive_flag=self.asu.is_inside(points[-1]),
      edge_inclusive_flag=None))
    return result

  def get_all_segments(self):
    result = []
    for edge_and_cuts in self.list:
      result.append(self.get_segments(edge_and_cuts))
    self.verify_edge_segments(result)
    return result

  def verify_edge_segments(self, all_edge_segments, gridding=13):
    for edge_segments in all_edge_segments:
      for i_segment in range(len(edge_segments)-1):
        a = edge_segments[i_segment].vertex
        b = edge_segments[i_segment+1].vertex
        is_inside = None
        for f in range(1, gridding):
          x = line_sample_point(a, b, f, gridding)
          if (is_inside is None):
            is_inside = self.asu.is_inside(x)
          else:
            assert is_inside == self.asu.is_inside(x)

def get_all_edge_segments(asu, list_of_polygons):
  return consolidated_edges_with_cut_points(
    asu, list_of_polygons).get_all_segments()

def get_all_vertices(all_edge_segments):
  result = {}
  for edge_segments in all_edge_segments:
    for segment in edge_segments:
      if (segment.vertex in result):
        assert result[segment.vertex] == segment.vertex_inclusive_flag
      else:
        result[segment.vertex] = segment.vertex_inclusive_flag
  return result
