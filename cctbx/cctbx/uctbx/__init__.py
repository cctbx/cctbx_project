import cctbx.array_family.flex

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.uctbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from scitbx.boost_python_utils import injector
import sys

class unit_cell(ext.unit_cell):

  def __init__(self, parameters, is_metrical_matrix=00000):
    if (isinstance(parameters, str)):
      parameters = [float(p) for p in parameters.replace(",", " ").split()]
    ext.unit_cell.__init__(self, parameters, is_metrical_matrix)

class _unit_cell(injector, ext.unit_cell):

  def __str__(self):
    return "(%.6g, %.6g, %.6g, %.6g, %.6g, %.6g)" % self.parameters()

  def show_parameters(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "Unit cell:", str(self)
