from scitbx.python_utils.misc import import_regular_symbols
from scitbx.array_family import flex
from cctbx.array_family import flex #for exposing af_shared(miller_index)
from iotbx_boost import mtz
import_regular_symbols(globals(), mtz.__dict__)
del import_regular_symbols
del mtz

baseMtz = Mtz

class Mtz (baseMtz):
  def __init__(self,s):
    baseMtz.__init__(self,s)

  def __getattr__(self,s):
    return self.getShared(s)
