from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.math import minimum_covering_sphere
import boost.python
import sys

float_cut_plane = sgtbx.direct_space_asu_float_cut_plane
float_asu = sgtbx.direct_space_asu_float_asu

class direct_space_asu:

  def __init__(self, hall_symbol, facets=[]):
    self.hall_symbol = hall_symbol
    self.facets = facets[:]

  def __and__(self, obj):
    self.facets.append(obj)
    return self

  def show_summary(self, f=None):
    if (f == None): f = sys.stdout
    print >> f, "Hall symbol:", self.hall_symbol
    print >> f, "Number of facets:", len(self.facets)
    return self

  def show_comprehensive_summary(self, f=None):
    if (f == None): f = sys.stdout
    self.show_summary(f)
    for facet in self.facets:
      print "    &", facet
    return self

  def is_inside(self, point, volume_only=00000):
    if (volume_only):
      for facet in self.facets:
        if (facet.evaluate(point) < 0): return 00000
    else:
      for facet in self.facets:
        if (not facet.is_inside(point)): return 00000
    return 0001

  def in_which_facets(self, point):
    result = []
    for facet in self.facets:
      if (facet.evaluate(point) == 0):
        result.append(facet)
    return result

  def volume_only(self):
    result = direct_space_asu(self.hall_symbol)
    for facet in self.facets:
      result.facets.append(facet.strip())
    return result

  def volume_vertices(self):
    from cctbx.sgtbx.direct_space_asu import facet_analysis
    return facet_analysis.volume_vertices(self)

  def change_basis(self, cb_op):
    if (not isinstance(cb_op, sgtbx.change_of_basis_op)):
      cb_op = sgtbx.change_of_basis_op(cb_op)
    cb_hall_symbol = None
    if (self.hall_symbol is not None):
      space_group_info = sgtbx.space_group_info("Hall: " + self.hall_symbol)
      cb_space_group_info = space_group_info.change_basis(cb_op)
      cb_hall_symbol = cb_space_group_info.type().hall_symbol()
    cb_asu = direct_space_asu(cb_hall_symbol)
    for facet in self.facets:
      cb_asu.facets.append(facet.change_basis(cb_op))
    return cb_asu

  def define_metric(self, unit_cell):
    return direct_space_asu_with_metric(asu=self, unit_cell=unit_cell)

  def add_buffer(self, unit_cell, thickness=None, relative_thickness=None):
    return self.define_metric(unit_cell).add_buffer(
      thickness=thickness,
      relative_thickness=relative_thickness)

class direct_space_asu_with_metric(direct_space_asu):

  def __init__(self, asu, unit_cell):
    direct_space_asu.__init__(self, asu.hall_symbol, asu.facets)
    self.unit_cell = unit_cell

  def minimum_covering_sphere(self, epsilon=None):
    if (epsilon is None):
      epsilon = self.unit_cell.volume()**(1/3.)*1.e-5
    points = flex.vec3_double()
    orth = self.unit_cell.orthogonalize
    for vertex in self.volume_vertices():
      points.append(orth([float(e) for e in vertex]))
    return minimum_covering_sphere(points=points, epsilon=epsilon)

  def as_float_asu(self):
    return float_asu(
      unit_cell=self.unit_cell,
      facets=[facet.as_float_cut_plane() for facet in self.facets])

  def add_buffer(self, thickness=None, relative_thickness=None):
    return self.as_float_asu().add_buffer(
      thickness=thickness,
      relative_thickness=relative_thickness)

class _float_asu(boost.python.injector, float_asu):

  def add_buffer(self, thickness=None, relative_thickness=None):
    assert [thickness, relative_thickness].count(None) > 0
    if (relative_thickness is None):
      relative_thickness = 1.e-6
    if (thickness is None):
      thickness = self.unit_cell().volume()**(1/3.)*relative_thickness
    return self._add_buffer(thickness)
