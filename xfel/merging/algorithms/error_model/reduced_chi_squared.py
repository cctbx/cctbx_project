from __future__ import division
from six.moves import range
from scitbx.array_family import flex

from xfel.merging.algorithms.error_model.error_modeler_base import error_modeler_base

class reduced_chi_squared(error_modeler_base):
  def compute(self):
    """
    Computes the reduced chi squared parameter. Can be used to
    correct for under-estimation of experimental errors.
    See
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Correcting_for_over-_or_under-dispersion
    """
    print >> self.log, "Computing reduced chi squared"
    self.scaler.reduced_chi_squared = flex.double(self.scaler.n_refl, 1.)

    for hkl_id in range(self.scaler.n_refl):
      hkl = self.scaler.miller_set.indices()[hkl_id]
      if hkl not in self.scaler.ISIGI: continue

      n = len(self.scaler.ISIGI[hkl])
      if n <= 1:
        continue

      i = self.scaler.summed_wt_I[hkl_id] / self.scaler.summed_weight[hkl_id]
      x = flex.double([self.scaler.ISIGI[hkl][i][0] for i in range(n)])
      v = (x / flex.double([self.scaler.ISIGI[hkl][i][1] for i in range(n)]))**2

      self.scaler.reduced_chi_squared[hkl_id] = 1/(n-1) * flex.sum((x-i)**2/v)
    sel = self.scaler.reduced_chi_squared > 0

    print >> self.log, "Done computing reduced chi squared", flex.mean(self.scaler.reduced_chi_squared.select(sel))

    if False:
      from matplotlib import pyplot as plt
      plt.hist(self.scaler.reduced_chi_squared.select(sel).as_numpy_array(),
        bins = 100, range=(0,10))
      plt.show()
