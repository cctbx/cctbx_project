"""Outdated Low-level methods for reading ccp4 map files. (Use map_manager
instead).
"""

from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex as flex# import dependency
import sys

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("iotbx_ccp4_map_ext")
from iotbx_ccp4_map_ext import *
import iotbx_ccp4_map_ext as ext

class utils :  # These routines are used by ccp4_map

  def show_summary(self, out=None, prefix=""):
    if (out is None) : out = sys.stdout
    data=self.map_data()

    if hasattr(self,'header_min'):
      print(prefix + "header_min: ", self.header_min, file=out)
      print(prefix + "header_max: ", self.header_max, file=out)
      print(prefix + "header_mean:", self.header_mean, file=out)
      print(prefix + "header_rms: ", self.header_rms, file=out)
    print("\n"+prefix + "Information about FULL UNIT CELL:",file=out)
    print(prefix + "unit cell grid:", self.unit_cell_grid, file=out)
    print(prefix + "unit cell parameters:", self.unit_cell_parameters, file=out)
    print(prefix + "space group number:  ", self.space_group_number, file=out)

    if not data:
      print("No map data available")
    else:
      print("\n"+prefix + "Information about the PART OF MAP THAT IS PRESENT:",
       file=out)
      print(prefix + "map cell grid:", data.all(), file=out)
      print(prefix + "map cell parameters:",
        self.crystal_symmetry().unit_cell().parameters(), file=out)
      print(prefix + "map origin:", data.origin(), file=out)
      print(prefix + "pixel size: (%.4f, %.4f, %.4f) " %(
        self.pixel_sizes()), file=out)
    if hasattr(self,'origin_shift_grid_units'):
      print(prefix + "Shift (grid units) to place origin at original position:",
          self.origin_shift_grid_units, file=out)

    if hasattr(self,'_model') and self._model:
      print (prefix + "Associated model with",
          self._model.get_hierarchy().overall_counts().n_residues,"residues",
           file=out)

    if self.high_resolution():
      print (prefix + "High-resolution limit of map: ",self.high_resolution(),
            file=out)

  def pixel_sizes(self):
    # Return tuple with pixel size in each direction (normally all the same)
    data=self.map_data()
    if not data:
      return None
    cs=self.crystal_symmetry()
    cell_params=cs.unit_cell().parameters()[:3]
    map_all=data.all()
    pa=cell_params[0]/map_all[0]
    pb=cell_params[1]/map_all[1]
    pc=cell_params[2]/map_all[2]
    return (pa,pb,pc)

  def crystal_symmetry(self,sorry_message_if_incompatible=None):
    '''
      This is "crystal_symmetry" of a box the size of the map that is present
    '''

    data=self.map_data()
    if not data:
      return None
    from cctbx import crystal
    map_all=self.data.all()

    if(map_all != self.unit_cell_grid):
      # map that is present is not exactly one unit cell. Calculate cell params
      a,b,c, al,be,ga = self.unit_cell().parameters()
      a = a * map_all[0]/self.unit_cell_grid[0]
      b = b * map_all[1]/self.unit_cell_grid[1]
      c = c * map_all[2]/self.unit_cell_grid[2]
      try:
        return crystal.symmetry((a,b,c, al,be,ga),
           self.space_group_number)
      except Exception as e:
        from libtbx.utils import Sorry
        if str(e).find(
          "incompatible") and \
          sorry_message_if_incompatible:
          raise Sorry(sorry_message_if_incompatible)
        else:
          raise Sorry(str(e))
    else:
      # map that is present is exactly one unit cell. Use unit cell symmetry
      return self.unit_cell_crystal_symmetry()

  def unit_cell_crystal_symmetry(self):
    '''
     This is the cell dimensions and angles of the full unit_cell
    '''
    from cctbx import crystal
    return crystal.symmetry(self.unit_cell().parameters(),
      self.space_group_number)


  def unit_cell(self):
    '''
     This is the full unit_cell
    '''
    from cctbx import uctbx
    return uctbx.unit_cell(self.unit_cell_parameters)

  def statistics(self):
    from cctbx import maptbx
    return maptbx.statistics(self.map_data())

  def get_origin(self):
    data=self.map_data()
    if data:
      return data.origin()
    else:
      return None

  def map_data(self):

    '''
       Input data is converted to double and stored in self.data
       self.map_data() always returns self.data
    '''

    return self.data

  def high_resolution(self):
    if hasattr(self,'_high_resolution'):
      return self._high_resolution
    else:
      return None

  def is_similar_map(self, other):
    f1 = self.crystal_symmetry().is_similar_symmetry(other.crystal_symmetry())
    s = self.map_data()
    o = other.map_data()
    if not s or not o:
      return None

    f2 = s.focus()  == o.focus()
    f3 = s.origin() == o.origin()
    f4 = s.all()    == o.all()
    if([f1,f2,f3,f4].count(False)>0): return False
    else: return True

  def grid_unit_cell(self):
    """
    This is a unit cell describing one pixel of the map.
    It is used in maptbx.non_crystallographic_eight_point_interpolation.
    This grid_unit_cell is the original unit cell divided by the original
    grid size.
    """
    from cctbx import uctbx
    a = self.unit_cell_parameters[0] / self.unit_cell_grid[0]
    b = self.unit_cell_parameters[1] / self.unit_cell_grid[1]
    c = self.unit_cell_parameters[2] / self.unit_cell_grid[2]
    alpha,beta,gamma = self.unit_cell_parameters[3:6]
    return uctbx.unit_cell((a,b,c,alpha,beta,gamma))

bp.inject(ext.map_reader, utils) # A way to access these
@bp.inject_into(ext.map_reader) # A way to access these
class _():

  def dummy(self):
    pass # don't do anything
