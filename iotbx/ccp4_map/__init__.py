import cctbx.array_family.flex # import dependency
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
