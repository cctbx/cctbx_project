from scitbx.python_utils.misc import import_regular_symbols
from cctbx_boost import uctbx_ext as ext
import_regular_symbols(globals(), ext.__dict__)
del import_regular_symbols

def unit_cell(parameters, is_metrical_matrix=False):
  if (type(parameters) == type("")):
    parameters = [float(p) for p in parameters.replace(",", " ").split()]
  return ext.unit_cell(parameters, is_metrical_matrix)
