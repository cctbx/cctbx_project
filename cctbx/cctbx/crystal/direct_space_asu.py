from cctbx import crystal
from cctbx import sgtbx
import cctbx.sgtbx.direct_space_asu
from cctbx.array_family import flex
from scitbx.math import minimum_covering_sphere
import boost.python

float_cut_plane = crystal.direct_space_asu_float_cut_plane
float_asu = crystal.direct_space_asu_float_asu
asu_mapping = crystal.direct_space_asu_asu_mapping
asu_mappings = crystal.direct_space_asu_asu_mappings

class direct_space_asu(sgtbx.direct_space_asu.direct_space_asu):

  def __init__(self, asu, unit_cell):
    sgtbx.direct_space_asu.direct_space_asu.__init__(self,
      hall_symbol=asu.hall_symbol, facets=asu.facets)
    self.unit_cell = unit_cell

  def minimum_covering_sphere(self, epsilon=None):
    if (epsilon is None): epsilon = 1.e-3
    points = flex.vec3_double()
    orth = self.unit_cell.orthogonalize
    for vertex in self.volume_vertices():
      points.append(orth([float(e) for e in vertex]))
    return minimum_covering_sphere(points=points, epsilon=epsilon)

  def as_float_asu(self, is_inside_epsilon=None):
    if (is_inside_epsilon is None):
      is_inside_epsilon = 1.e-6
    return float_asu(
      unit_cell=self.unit_cell,
      facets=[facet.as_float_cut_plane() for facet in self.facets],
      is_inside_epsilon=is_inside_epsilon)

  def add_buffer(self, thickness=None, relative_thickness=None,
                       is_inside_epsilon=None):
    return self.as_float_asu(is_inside_epsilon=is_inside_epsilon).add_buffer(
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
