from scitbx.python_utils.misc import import_regular_symbols
from scitbx.array_family import flex as sciflex
from cctbx.array_family import flex as ccflex#for exposing af_shared(miller_index)
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
    
  def HL(self,hla_name,hlb_name,hlc_name,hld_name):
    hla = self.getColumn(hla_name)
    hlb = self.getColumn(hlb_name)
    hlc = self.getColumn(hlc_name)
    hld = self.getColumn(hld_name)
    assert hla.type()=='A'
    assert hlb.type()=='A'
    assert hlc.type()=='A'
    assert hld.type()=='A'
    return baseMtz.HL(self,hla,hlb,hlc,hld)

  def complex(self,amplitude_name,phi_name):
    amplitude = self.getColumn(amplitude_name)
    assert amplitude.type()=='F'
    phase = self.getColumn(phi_name)
    assert phase.type()=='P'
    amplitude = self.getShared(amplitude_name)
    phase = self.getShared(phi_name)
    return ccflex.polar(amplitude,phase,1) #degrees=1
