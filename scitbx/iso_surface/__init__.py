from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("scitbx_iso_surface_ext")
import scitbx_iso_surface_ext
# the following is essential to load the binding
# of af::shared<tiny<size_t, 3> >
from scitbx.array_family import shared


class triangulation(scitbx_iso_surface_ext.triangulation):

  def __init__(self, map, iso_level, map_extent,
               from_here=None, to_there=None,
               periodic=False,
               lazy_normals=True, ascending_normal_direction=True):
    if from_here is None: from_here = (0,0,0)
    if to_there is None: to_there = map_extent
    super(triangulation, self).__init__(
      map, iso_level, map_extent,
      from_here, to_there, periodic,
      lazy_normals, ascending_normal_direction)

  def bounds(self):
    return self.from_here, self.to_there
  bounds = property(bounds)
