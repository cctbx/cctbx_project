from cctbx import sgtbx
import sys

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
    import cctbx.crystal.direct_space_asu
    return cctbx.crystal.direct_space_asu.direct_space_asu(
      asu=self, unit_cell=unit_cell)

  def add_buffer(self, unit_cell, thickness=None, relative_thickness=None):
    return self.define_metric(unit_cell).add_buffer(
      thickness=thickness,
      relative_thickness=relative_thickness)
