from __future__ import division
from __future__ import print_function
import cctbx.array_family.flex as flex# import dependency
import sys

import boost.python
ext = boost.python.import_ext("iotbx_ccp4_map_ext")
from iotbx_ccp4_map_ext import *
import iotbx_ccp4_map_ext as ext

class _(boost.python.injector, ext.map_reader) :

  def show_summary (self, out=None, prefix="") :
    if (out is None) : out = sys.stdout
    print(prefix + "header_min: ", self.header_min, file=out)
    print(prefix + "header_max: ", self.header_max, file=out)
    print(prefix + "header_mean:", self.header_mean, file=out)
    print(prefix + "header_rms: ", self.header_rms, file=out)
    print(prefix + "unit cell grid:", self.unit_cell_grid, file=out)
    print(prefix + "unit cell parameters:", self.unit_cell_parameters, file=out)
    print(prefix + "space group number:  ", self.space_group_number, file=out)
    print(prefix + "map origin:", self.data.origin(), file=out)
    print(prefix + "map grid:  ", self.data.all(), file=out)

  def crystal_symmetry(self):
    from cctbx import crystal
    return crystal.symmetry(self.unit_cell().parameters(),
      self.space_group_number)

  def unit_cell (self) :
    from cctbx import uctbx
    return uctbx.unit_cell(self.unit_cell_parameters)

  def statistics (self) :
    from cctbx import maptbx
    return maptbx.statistics(self.data)

  def map_data(self):
    return self.data.as_double()

  def grid_unit_cell (self) :
    """
    If we want to use maptbx.non_crystallographic_eight_point_interpolation,
    the "unit cell" is actually the original unit cell divided by the original
    grid size.
    """
    from cctbx import uctbx
    a = self.unit_cell_parameters[0] / self.unit_cell_grid[0]
    b = self.unit_cell_parameters[1] / self.unit_cell_grid[1]
    c = self.unit_cell_parameters[2] / self.unit_cell_grid[2]
    alpha,beta,gamma = self.unit_cell_parameters[3:6]
    return uctbx.unit_cell((a,b,c,alpha,beta,gamma))
