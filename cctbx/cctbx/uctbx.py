import cctbx.array_family.flex

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.uctbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

import sys

def unit_cell(parameters, is_metrical_matrix=False):
  if (type(parameters) == type("")):
    parameters = [float(p) for p in parameters.replace(",", " ").split()]
  return ext.unit_cell(parameters, is_metrical_matrix)

def show_parameters(unit_cell, f=sys.stdout):
  print >> f, "Unit cell: (%.6g, %.6g, %.6g, %.6g, %.6g, %.6g)" \
              % unit_cell.parameters()
