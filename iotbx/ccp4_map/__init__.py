from __future__ import division
import cctbx.array_family.flex as flex# import dependency
import sys

import boost.python
ext = boost.python.import_ext("iotbx_ccp4_map_ext")
from iotbx_ccp4_map_ext import *
import iotbx_ccp4_map_ext as ext

class _(boost.python.injector, ext.map_reader) :

  def show_summary (self, out=None, prefix="") :
    if (out is None) : out = sys.stdout
    print >> out, prefix + "header_min: ", self.header_min
    print >> out, prefix + "header_max: ", self.header_max
    print >> out, prefix + "header_mean:", self.header_mean
    print >> out, prefix + "header_rms: ", self.header_rms
    print >> out, prefix + "unit cell grid:", self.unit_cell_grid
    print >> out, prefix + "unit cell parameters:", self.unit_cell_parameters
    print >> out, prefix + "space group number:  ", self.space_group_number
    print >> out, prefix + "map origin:", self.data.origin()
    print >> out, prefix + "map grid:  ", self.data.all()

  def crystal_symmetry(self):
    from cctbx import crystal
    map_all = self.map_data().all()
    if(map_all != self.unit_cell_grid):
      a,b,c, al,be,ga = self.unit_cell().parameters()
      a = a * map_all[0]/self.unit_cell_grid[0]
      b = b * map_all[1]/self.unit_cell_grid[1]
      c = c * map_all[2]/self.unit_cell_grid[2]
      return crystal.symmetry((a,b,c, al,be,ga), self.space_group_number)
    else:
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

  def is_similar_map(self, other):
    f1 = self.crystal_symmetry().is_similar_symmetry(other.crystal_symmetry())
    s = self.map_data()
    o = other.map_data()
    f2 = s.focus()  == o.focus()
    f3 = s.origin() == o.origin()
    f4 = s.all()    == o.all()
    if([f1,f2,f3,f4].count(False)>0): return False
    else: return True

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
