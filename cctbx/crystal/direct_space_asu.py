from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
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
      hall_symbol=asu.hall_symbol, cuts=asu.cuts)
    self.unit_cell = unit_cell

  def minimum_covering_sphere(self, epsilon=None):
    if (epsilon is None): epsilon = 1.e-3
    points = flex.vec3_double()
    orth = self.unit_cell.orthogonalize
    for vertex in self.shape_vertices():
      points.append(orth([float(e) for e in vertex]))
    return minimum_covering_sphere(points=points, epsilon=epsilon)

  def as_float_asu(self, is_inside_epsilon=None):
    if (is_inside_epsilon is None):
      is_inside_epsilon = 1.e-6
    return float_asu(
      unit_cell=self.unit_cell,
      cuts=[cut.as_float_cut_plane() for cut in self.cuts],
      is_inside_epsilon=is_inside_epsilon)

  def add_buffer(self, thickness=None, relative_thickness=None,
                       is_inside_epsilon=None):
    return self.as_float_asu(is_inside_epsilon=is_inside_epsilon).add_buffer(
      thickness=thickness,
      relative_thickness=relative_thickness)

class _(boost.python.injector, float_asu):

  def add_buffer(self, thickness=None, relative_thickness=None):
    assert [thickness, relative_thickness].count(None) > 0
    if (relative_thickness is None):
      relative_thickness = 1.e-6
    if (thickness is None):
      thickness = self.unit_cell().volume()**(1/3.)*relative_thickness
    return self._add_buffer(thickness)

class _(boost.python.injector, asu_mappings):

  def get_rt_mx_ji(self, pair):
    return self.get_rt_mx_i(pair).inverse().multiply(self.get_rt_mx_j(pair))

def non_crystallographic_asu_mappings(
      sites_cart,
      default_buffer_layer=0.5,
      min_unit_cell_length=0):
  sites_min = sites_cart.min()
  sites_max = sites_cart.max()
  crystal_symmetry = crystal.non_crystallographic_symmetry(
      sites_cart_min=sites_min,
      sites_cart_max=sites_max,
      default_buffer_layer=default_buffer_layer,
      min_unit_cell_length=min_unit_cell_length)
  buffer_layer = uctbx.non_crystallographic_buffer_layer(
      sites_cart_min=sites_min,
      sites_cart_max=sites_max,
      default_buffer_layer=default_buffer_layer)
  sites_min = crystal_symmetry.unit_cell().fractionalize(sites_min)
  sites_max = crystal_symmetry.unit_cell().fractionalize(sites_max)
  asu_cuts = [float_cut_plane(n=n,c=c) for n,c in [
    ([1,0,0],-sites_min[0]),
    ([-1,0,0],sites_max[0]),
    ([0,1,0],-sites_min[1]),
    ([0,-1,0],sites_max[1]),
    ([0,0,1],-sites_min[2]),
    ([0,0,-1],sites_max[2]),
  ]]
  result = asu_mappings(
    space_group=crystal_symmetry.space_group(),
    asu=float_asu(
      unit_cell=crystal_symmetry.unit_cell(),
      cuts=asu_cuts).add_buffer(thickness=buffer_layer),
    buffer_thickness=0)
  result.process_sites_cart(original_sites=sites_cart)
  return result
