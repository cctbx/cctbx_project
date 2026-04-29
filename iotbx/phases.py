"""Tool to write miller array containing phases and amplitudes in XTALVIEW format"""
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import math
from six.moves import zip

def miller_array_as_phases_phs(self,
      out,
      scale_amplitudes=True,
      phases=None,
      phases_deg=None,
      figures_of_merit=None):
  """http://www.sdsc.edu/CCMS/Packages/XTALVIEW/xtalviewfaq.html"""
  if (phases is not None): assert phases_deg is False or phases_deg is True
  if (self.is_complex_array()):
    amplitudes = self.amplitudes().data()
  else:
    amplitudes = self.data()
  if (scale_amplitudes):
    amplitudes_max = flex.max(amplitudes)
    if (amplitudes_max > 0):
      amplitudes = (9999.99/amplitudes_max) * amplitudes
  assert len(" %7.2f" % flex.min(amplitudes)) == 8
  assert len(" %7.2f" % flex.max(amplitudes)) == 8
  if (phases is None):
    phases = self.phases(deg=True).data()
  else:
    if (hasattr(phases, "data")):
      phases = phases.data()
    if (not phases_deg):
      phases = phases * (180/math.pi)
  assert len(" %7.2f" % flex.min(phases)) == 8
  assert len(" %7.2f" % flex.max(phases)) == 8
  if (figures_of_merit is None):
    for h,a,p in zip(self.indices(),amplitudes,phases):
      print("%4d%4d%4d" % h + " %7.2f"%a + " %7.2f"%1 + " %7.2f"%p, file=out)
  else:
    if (hasattr(figures_of_merit, "data")):
      assert figures_of_merit.indices().all_eq(self.indices())
      figures_of_merit = figures_of_merit.data()
    assert len(" %7.2f" % flex.min(figures_of_merit)) == 8
    assert len(" %7.2f" % flex.max(figures_of_merit)) == 8
    for h,a,p,f in zip(self.indices(),amplitudes,phases,figures_of_merit):
      print("%4d%4d%4d" % h + " %7.2f"%a + " %7.2f"%f + " %7.2f"%p, file=out)
