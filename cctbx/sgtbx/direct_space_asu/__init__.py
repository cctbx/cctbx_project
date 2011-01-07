from cctbx import sgtbx
from cctbx.sgtbx.direct_space_asu import cut_plane
from cctbx.sgtbx.direct_space_asu.short_cuts import r1
from scitbx import matrix
import sys

class direct_space_asu(object):

  def __init__(self, hall_symbol, cuts=[]):
    self.hall_symbol = hall_symbol
    self.cuts = cuts[:]

  def __copy__(self):
    return direct_space_asu(
      hall_symbol=self.hall_symbol,
      cuts=self.cuts)

  def __and__(self, obj):
    self.cuts.append(obj)
    return self

  def show_summary(self, f=None):
    if (f == None): f = sys.stdout
    print >> f, "Hall symbol:", self.hall_symbol
    print >> f, "Number of cuts:", len(self.cuts)
    return self

  def show_comprehensive_summary(self, f=None):
    if (f == None): f = sys.stdout
    self.show_summary(f)
    for cut in self.cuts:
      print >> f, "    &", cut
    return self

  def is_inside(self, point, shape_only=False):
    if (shape_only):
      for cut in self.cuts:
        if (cut.evaluate(point) < 0): return False
    else:
      for cut in self.cuts:
        if (not cut.is_inside(point)): return False
    return True

  def in_which_cuts(self, point):
    result = []
    for cut in self.cuts:
      if (cut.evaluate(point) == 0):
        result.append(cut)
    return result

  def extract_all_cuts(self):
    result = []
    for cut in self.cuts:
      cut.extract_all_cuts(result)
    return result

  def shape_only(self):
    result = direct_space_asu(self.hall_symbol)
    for cut in self.cuts:
      result.cuts.append(cut.strip())
    return result

  def shape_vertices(self):
    result = set()
    cuts = self.cuts
    n_cuts = len(cuts)
    for i0 in xrange(0,n_cuts-2):
      for i1 in xrange(i0+1,n_cuts-1):
        for i2 in xrange(i1+1,n_cuts):
          m = matrix.rec(cuts[i0].n+cuts[i1].n+cuts[i2].n,(3,3))
          d = m.determinant()
          if (d != 0):
            m_inv = m.co_factor_matrix_transposed() * (r1/d)
            b = matrix.col([-cuts[i0].c,-cuts[i1].c,-cuts[i2].c])
            vertex = m_inv * b
            if (self.is_inside(vertex, shape_only=True)):
              result.add(vertex.elems)
    return sorted(result)

  def _box_corner(self, shape_vertices, min_or_max):
    if (shape_vertices is None):
      shape_vertices = self.shape_vertices()
    if (len(shape_vertices) == 0):
      return None
    result = list(shape_vertices[0])
    for vertex in shape_vertices[1:]:
      for i in xrange(3):
        result[i] = min_or_max(result[i], vertex[i])
    return result

  def box_min(self, shape_vertices=None):
    return self._box_corner(shape_vertices, min)

  def box_max(self, shape_vertices=None):
    return self._box_corner(shape_vertices, max)

  def add_plane(self, normal_direction, point=None):
    if (point is None):
      point = self.box_min()
    self.cuts.append(cut_plane.cut(
      n=normal_direction,
      c=-(matrix.col(normal_direction).dot(matrix.col(point)))))

  def add_planes(self, normal_directions, point=None, both_directions=False):
    if (point is None):
      point = self.box_min()
    for normal_direction in normal_directions:
      self.add_plane(normal_direction=normal_direction, point=point)
      if (both_directions):
        self.cuts.append(-self.cuts[-1])

  def change_basis(self, cb_op):
    if (not isinstance(cb_op, sgtbx.change_of_basis_op)):
      cb_op = sgtbx.change_of_basis_op(cb_op)
    cb_hall_symbol = None
    if (self.hall_symbol is not None):
      space_group_info = sgtbx.space_group_info("Hall: " + self.hall_symbol)
      cb_space_group_info = space_group_info.change_basis(cb_op)
      cb_hall_symbol = cb_space_group_info.type().hall_symbol()
    cb_asu = direct_space_asu(cb_hall_symbol)
    for cut in self.cuts:
      cb_asu.cuts.append(cut.change_basis(cb_op))
    return cb_asu

  def define_metric(self, unit_cell):
    import cctbx.crystal.direct_space_asu
    return cctbx.crystal.direct_space_asu.direct_space_asu(
      asu=self, unit_cell=unit_cell)

  def add_buffer(self, unit_cell, thickness=None, relative_thickness=None):
    return self.define_metric(unit_cell).add_buffer(
      thickness=thickness,
      relative_thickness=relative_thickness)
