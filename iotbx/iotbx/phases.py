from cctbx.array_family import flex
import math

def miller_array_as_phases_phs(self,
      out,
      scale_amplitudes=True,
      phases=None,
      phases_deg=None,
      figures_of_merit=None):
  if (phases is not None): assert phases_deg is False or phases_deg is True
  if (self.is_complex_array()):
    amplitudes = self.amplitudes().data()
  else:
    amplitudes = self.data()
  if (scale_amplitudes):
    amplitudes_max = flex.max(amplitudes)
    if (amplitudes_max > 0):
      amplitudes = (99999.99/amplitudes_max) * amplitudes
  assert len("%8.2f" % flex.min(amplitudes)) == 8
  assert len("%8.2f" % flex.max(amplitudes)) == 8
  if (phases is None):
    phases = self.phases(deg=True).data()
  else:
    if (hasattr(phases, "data")):
      phases = phases.data()
    if (not phases_deg):
      phases = phases * (180/math.pi)
  assert len("%8.2f" % flex.min(phases)) == 8
  assert len("%8.2f" % flex.max(phases)) == 8
  if (figures_of_merit is None):
    for h,a,p in zip(self.indices(),amplitudes,phases):
      print >> out, "%4d%4d%4d" % h + "%8.2f" % a + "%8.2f" % 1 + "%8.2f" % p
  else:
    if (hasattr(figures_of_merit, "data")):
      assert figures_of_merit.indices().all_eq(self.indices())
      figures_of_merit = figures_of_merit.data()
    assert len("%8.2f" % flex.min(figures_of_merit)) == 8
    assert len("%8.2f" % flex.max(figures_of_merit)) == 8
    for h,a,p,f in zip(self.indices(),amplitudes,phases,figures_of_merit):
      print >> out, "%4d%4d%4d" % h + "%8.2f" % a + "%8.2f" % f + "%8.2f" % p
