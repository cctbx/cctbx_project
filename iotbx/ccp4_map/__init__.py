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

  def unit_cell (self) :
    from cctbx import uctbx
    return uctbx.unit_cell(self.unit_cell_parameters)

  def statistics (self) :
    from cctbx import maptbx
    return maptbx.statistics(self.data)

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

  def map_data(self):
    m_data = self.data.as_double()
    n_real = self.unit_cell_grid
    if(n_real == m_data.all()):
      return m_data
    else:
      # XXX hideously SLOW! MOVE TO C++
      map_new = flex.double(flex.grid(n_real), 0)
      o = m_data.origin()
      f = m_data.focus()
      for i in range(o[0],f[0]):
        for j in range(o[1],f[1]):
          for k in range(o[2],f[2]):
            map_new[i%n_real[0], j%n_real[1], k%n_real[2]] = m_data[i, j, k]
      return map_new
