from __future__ import division
from scitbx.array_family import flex

from xfel.merging.algorithms.error_model.error_modeler_base import error_modeler_base

class errors_from_residuals(error_modeler_base):
  def __init__(self, scaler):
    self.scaler = scaler
    self.log = scaler.log

  def adjust_errors(self):
    """
    Use the distribution of intensities in a given miller index to compute the error for each merged reflection
    """
    print >> self.log, "Computing error estimates from sample residuals"
    self.scaler.summed_weight= flex.double(self.scaler.n_refl, 0.)
    self.scaler.summed_wt_I  = flex.double(self.scaler.n_refl, 0.)

    for hkl_id in xrange(self.scaler.n_refl):
      hkl = self.scaler.miller_set.indices()[hkl_id]
      if hkl not in self.scaler.ISIGI: continue

      n = len(self.scaler.ISIGI[hkl])
      if n > 1:
        variance = flex.mean_and_variance(flex.double([self.scaler.ISIGI[hkl][i][0] for i in xrange(n)])).unweighted_sample_variance()
      else:
        continue

      for i in xrange(n):
        Intensity = self.scaler.ISIGI[hkl][i][0] # scaled intensity
        self.scaler.summed_wt_I[hkl_id] += Intensity / variance
        self.scaler.summed_weight[hkl_id] += 1 / variance
    print >> self.log, "Done computing error estimates"
