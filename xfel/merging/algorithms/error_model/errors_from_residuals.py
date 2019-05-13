from __future__ import division
from __future__ import print_function
from six.moves import range
from scitbx.array_family import flex

from xfel.merging.algorithms.error_model.error_modeler_base import error_modeler_base

class errors_from_residuals(error_modeler_base):
  def adjust_errors(self):
    """
    Use the distribution of intensities in a given miller index to compute the error for each merged reflection
    """
    print("Computing error estimates from sample residuals", file=self.log)
    self.scaler.summed_weight= flex.double(self.scaler.n_refl, 0.)
    self.scaler.summed_wt_I  = flex.double(self.scaler.n_refl, 0.)

    for hkl_id in range(self.scaler.n_refl):
      hkl = self.scaler.miller_set.indices()[hkl_id]
      if hkl not in self.scaler.ISIGI: continue

      n = len(self.scaler.ISIGI[hkl])
      if n > 1:
        variance = flex.mean_and_variance(flex.double([self.scaler.ISIGI[hkl][i][0] for i in range(n)])).unweighted_sample_variance()
      else:
        continue

      for i in range(n):
        Intensity = self.scaler.ISIGI[hkl][i][0] # scaled intensity
        self.scaler.summed_wt_I[hkl_id] += Intensity / variance
        self.scaler.summed_weight[hkl_id] += 1 / variance
    print("Done computing error estimates", file=self.log)

class errors_from_residuals_refltable(error_modeler_base):
  def adjust_errors(self):
    """
    Use the distribution of intensities in a given miller index to compute the error for each merged reflection
    """
    print("Computing error estimates from sample residuals", file=self.log)
    self.scaler.summed_weight= flex.double(self.scaler.n_refl, 0.)
    self.scaler.summed_wt_I  = flex.double(self.scaler.n_refl, 0.)

    # 3-pass variance
    print("Variance step 1 of 3...", file=self.log)
    hkl_sumI = flex.double(self.scaler.n_refl, 0)
    n_obs = flex.double(self.scaler.n_refl, 0)

    for i in range(len(self.scaler.ISIGI)):
      idx = self.scaler.ISIGI['miller_id'][i]
      hkl_sumI[idx] += self.scaler.ISIGI['scaled_intensity'][i]

      n_obs[idx] += 1

    hkl_mean = flex.double(self.scaler.n_refl, 0)
    sel = n_obs > 0
    hkl_mean.set_selected(sel, hkl_sumI.select(sel)/n_obs.select(sel))

    print("Variance step 2 of 3...", file=self.log)
    hkl_sumsq = flex.double(self.scaler.n_refl, 0)
    for i in range(len(self.scaler.ISIGI)):
      idx = self.scaler.ISIGI['miller_id'][i]
      hkl_sumsq[idx] += (self.scaler.ISIGI['scaled_intensity'][i] - hkl_mean[idx])**2

    print("Variance step 3 of 3...", file=self.log)
    for i in range(len(self.scaler.ISIGI)):
      idx = self.scaler.ISIGI['miller_id'][i]
      n = n_obs[idx]
      if n <= 1:
        continue
      Intensity = self.scaler.ISIGI['scaled_intensity'][i]
      variance = hkl_sumsq[idx] / (n-1)

      self.scaler.summed_wt_I[idx] += Intensity / variance
      self.scaler.summed_weight[idx] += 1 / variance
    print("Done computing error estimates", file=self.log)
