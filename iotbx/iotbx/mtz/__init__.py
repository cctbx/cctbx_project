from scitbx.python_utils.misc import import_regular_symbols
from iotbx_boost import mtz as ext
import_regular_symbols(globals(), ext.__dict__)
del import_regular_symbols

from scitbx.array_family import flex as sciflex
from cctbx import miller
from cctbx.array_family import flex as ccflex

class Mtz (ext.Mtz):
  def __init__(self,s):
    ext.Mtz.__init__(self,s)

  def __getattr__(self,s):
    assert type(s) == type(str())
    return self.getColumn(s)
