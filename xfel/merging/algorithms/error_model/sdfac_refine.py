from __future__ import division
from scitbx.array_family import flex
import math

from xfel.merging.algorithms.error_model.error_modeler_base import error_modeler_base

from scitbx.simplex import simplex_opt
class simplex_minimizer(object):
  """Class for refining sdfac, sdb and sdadd"""
  def __init__(self, sdfac, sdb, sdadd, data, indices, bins, seed = None, log=None):
    """
    @param sdfac Initial value for sdfac
    @param sdfac Initial value for sdfac
    @param sdfac Initial value for sdfac
    @param data ISIGI dictionary of unmerged intensities
    @param indices array of miller indices to refine against
    @param bins array of flex.bool object specifying the bins to use to calculate the functional
    @param log Log to print to (none for stdout)
    """
    if log is None:
      log = sys.stdout
    self.log = log
    self.data = data
    self.intensity_bin_selections = bins
    self.indices = indices
    self.n = 3
    self.x = flex.double([sdfac, sdb, sdadd])
    self.starting_simplex = []
    if seed is None:
      random_func = flex.random_double
    else:
      print >> self.log, "Using random seed %d"%seed
      mt = flex.mersenne_twister(seed)
      random_func = mt.random_double

    for i in xrange(self.n+1):
      self.starting_simplex.append(random_func(self.n))

    self.optimizer = simplex_opt( dimension = self.n,
                                  matrix    = self.starting_simplex,
                                  evaluator = self,
                                  tolerance = 1e-1)
    self.x = self.optimizer.get_solution()

  def target(self, vector):
    """ Compute the functional by first applying the current values for the sd parameters
    to the input data, then computing the complete set of normalized deviations and finally
    using those normalized deviations to compute the functional."""
    from xfel import compute_normalized_deviations, apply_sd_error_params

    sdfac, sdb, sdadd = vector

    if sdfac < 0 or sdb < 0 or sdadd < 0:
      f = 1e6
    else:
      data = apply_sd_error_params(self.data, sdfac, sdb, sdadd)
      all_sigmas_normalized = compute_normalized_deviations(data, self.indices)

      f = 0
      for bin in self.intensity_bin_selections:
        binned_normalized_sigmas = all_sigmas_normalized.select(bin)
        n = len(binned_normalized_sigmas)
        if n == 0: continue
        # weighting scheme from Evans, 2011
        w = math.sqrt(n)
        # functional is weight * (1-rms(normalized_sigmas))^s summed over all intensitiy bins
        f += w * ((1-math.sqrt(flex.mean(binned_normalized_sigmas*binned_normalized_sigmas)))**2)

    print >> self.log, "f: % 12.1f, sdfac: %8.5f, sdb: %8.5f, sdadd: %8.5f"%(f, sdfac, sdb, sdadd)
    return f

class sdfac_refine(error_modeler_base):
  def get_overall_correlation_flex (self, data_a, data_b) :
    """
    Correlate any two sets of data.
    @param data_a data set a
    @param data_b data set b
    @return tuple containing correlation coefficent, slope and offset.
    """
    import math

    assert len(data_a) == len(data_b)
    corr = 0
    slope = 0
    offset = 0
    try:
      sum_xx = 0
      sum_xy = 0
      sum_yy = 0
      sum_x  = 0
      sum_y  = 0
      N      = 0
      for i in xrange(len(data_a)):
        I_r       = data_a[i]
        I_o       = data_b[i]
        N      += 1
        sum_xx += I_r**2
        sum_yy += I_o**2
        sum_xy += I_r * I_o
        sum_x  += I_r
        sum_y  += I_o
      slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
      offset = (sum_xx * sum_y - sum_x * sum_xy) / (N * sum_xx - sum_x**2)
      corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
                 math.sqrt(N * sum_yy - sum_y**2))
    except ZeroDivisionError:
      pass

    return corr, slope, offset

  def normal_probability_plot(self, data, rankits_sel=None, plot=False):
    """ Use normal probability analysis to determine if a set of data is normally distributed
    See https://en.wikipedia.org/wiki/Normal_probability_plot.
    Rankits are computed in the same way as qqnorm does in R.
    @param data flex array
    @param rankits_sel only use the rankits in a certain range. Useful for outlier rejection. Should be
    a tuple such as (-0.5,0.5).
    @param plot whether to show the normal probabilty plot
    """
    from scitbx.math import distributions
    import numpy as np
    norm = distributions.normal_distribution()

    n = len(data)
    if n <= 10:
      a = 3/8
    else:
      a = 0.5

    sorted_data = flex.sorted(data)
    rankits = flex.double([norm.quantile((i+1-a)/(n+1-(2*a))) for i in xrange(n)])

    if rankits_sel is None:
      corr, slope, offset = self.get_overall_correlation_flex(sorted_data, rankits)
    else:
      sel = (rankits >= rankits_sel[0]) & (rankits <= rankits_sel[1])
      corr, slope, offset = self.get_overall_correlation_flex(sorted_data.select(sel), rankits.select(sel))

    if plot:
      from matplotlib import pyplot as plt
      x = np.linspace(-500,500,100) # 100 linearly spaced numbers
      y = slope * x + offset
      plt.scatter(sorted_data, rankits)
      plt.plot(x,y)

      plt.figure()
      plt.hist(rankits, bins=100)
      plt.show()

    return corr, slope, offset

  def get_initial_sdparams_estimates(self):
    """
    Use normal probability analysis to compute intial sdfac and sdadd parameters.
    """
    from xfel import compute_normalized_deviations
    all_sigmas_normalized = compute_normalized_deviations(self.scaler.ISIGI, self.scaler.miller_set.indices())
    assert ((all_sigmas_normalized > 0) | (all_sigmas_normalized <= 0)).count(True) == len(all_sigmas_normalized) # no nans allowed

    corr, slope, offset = self.normal_probability_plot(all_sigmas_normalized, (-0.5, 0.5))
    sdfac = slope
    sdb = 0
    sdadd = offset

    return sdfac, sdb, sdadd


  def get_binned_intensities(self, n_bins=100):
    """
    Using self.ISIGI, bin the intensities using the following procedure:
    1) Find the minimum and maximum intensity values.
    2) Divide max-min by n_bins. This is the bin step size
    The effect is
    @param n_bins number of bins to use.
    @return a tuple with an array of selections for each bin and an array of median
    intensity values for each bin.
    """
    print >> self.log, "Computing intensity bins.",
    all_mean_Is = flex.double()
    only_means = flex.double()
    for hkl_id in xrange(self.scaler.n_refl):
      hkl = self.scaler.miller_set.indices()[hkl_id]
      if hkl not in self.scaler.ISIGI: continue
      n = len(self.scaler.ISIGI[hkl])
      # get scaled intensities
      intensities = flex.double([self.scaler.ISIGI[hkl][i][0] for i in xrange(n)])
      meanI = flex.mean(intensities)
      only_means.append(meanI)
      all_mean_Is.extend(flex.double([meanI]*n))
    step = (flex.max(only_means)-flex.min(only_means))/n_bins
    print >> self.log, "Bin size:", step

    sels = []
    binned_intensities = []
    min_all_mean_Is = flex.min(all_mean_Is)
    for i in xrange(n_bins):
      sel = (all_mean_Is > (min_all_mean_Is + step * i)) & (all_mean_Is < (min_all_mean_Is + step * (i+1)))
      if sel.all_eq(False): continue
      sels.append(sel)
      binned_intensities.append((step/2 + step*i)+min(only_means))

    for i, (sel, intensity) in enumerate(zip(sels, binned_intensities)):
      print >> self.log, "Bin %02d, number of observations: % 10d, midpoint intensity: %f"%(i, sel.count(True), intensity)

    return sels, binned_intensities

  def adjust_errors(self):
    """
    Adjust sigmas according to Evans, 2011 Acta D and Evans and Murshudov, 2013 Acta D
    """
    print >> self.log, "Starting adjust_errors"
    print >> self.log, "Computing initial estimates of sdfac, sdb and sdadd"
    sdfac, sdb, sdadd = self.get_initial_sdparams_estimates()

    print >> self.log, "Initial estimates:", sdfac, sdb, sdadd

    from xfel import compute_normalized_deviations, apply_sd_error_params

    print >> self.log, "Refining error correction parameters sdfac, sdb, and sdadd"
    sels, binned_intensities = self.get_binned_intensities()
    seed = self.scaler.params.raw_data.error_models.sdfac_refine.random_seed
    minimizer = simplex_minimizer(sdfac, sdb, sdadd, self.scaler.ISIGI, self.scaler.miller_set.indices(), sels, seed, self.log)
    sdfac, sdb, sdadd = minimizer.x
    print >> self.log, "Final sdadd: %8.5f, sdb: %8.5f, sdadd: %8.5f"%(sdfac, sdb, sdadd)

    print >> self.log, "Applying sdfac/sdb/sdadd 1"
    self.scaler.ISIGI = apply_sd_error_params(self.scaler.ISIGI, sdfac, sdb, sdadd)

    self.scaler.summed_weight= flex.double(self.scaler.n_refl, 0.)
    self.scaler.summed_wt_I  = flex.double(self.scaler.n_refl, 0.)

    print >> self.log, "Applying sdfac/sdb/sdadd 2"
    for hkl_id in xrange(self.scaler.n_refl):
      hkl = self.scaler.miller_set.indices()[hkl_id]
      if hkl not in self.scaler.ISIGI: continue

      n = len(self.scaler.ISIGI[hkl])

      for i in xrange(n):
        Intensity = self.scaler.ISIGI[hkl][i][0] # scaled intensity
        sigma = Intensity / self.scaler.ISIGI[hkl][i][1] # corrected sigma
        variance = sigma * sigma
        self.scaler.summed_wt_I[hkl_id] += Intensity / variance
        self.scaler.summed_weight[hkl_id] += 1 / variance

    if False:
      # validate using http://ccp4wiki.org/~ccp4wiki/wiki/index.php?title=Symmetry%2C_Scale%2C_Merge#Analysis_of_Standard_Deviations
      print >> self.log, "Validating"
      from matplotlib import pyplot as plt
      all_sigmas_normalized = compute_normalized_deviations(self.scaler.ISIGI, self.scaler.miller_set.indices())
      plt.hist(all_sigmas_normalized, bins=100)
      plt.figure()

      binned_rms_normalized_sigmas = []

      for i, sel in enumerate(sels):
        binned_rms_normalized_sigmas.append(math.sqrt(flex.mean(all_sigmas_normalized.select(sel)*all_sigmas_normalized.select(sel))))

      plt.plot(binned_intensities, binned_rms_normalized_sigmas, 'o')
      plt.show()

      self.normal_probability_plot(all_sigmas_normalized, (-0.5, 0.5), plot = True)

class simplex_minimizer_refltable(object):
  """Class for refining sdfac, sdb and sdadd"""
  def __init__(self, sdfac, sdb, sdadd, data, indices, bins, seed = None, log=None):
    """
    @param sdfac Initial value for sdfac
    @param sdfac Initial value for sdfac
    @param sdfac Initial value for sdfac
    @param data ISIGI dictionary of unmerged intensities
    @param indices array of miller indices to refine against
    @param bins array of flex.bool object specifying the bins to use to calculate the functional
    @param log Log to print to (none for stdout)
    """
    if log is None:
      log = sys.stdout
    self.log = log
    self.data = data
    self.intensity_bin_selections = bins
    self.indices = indices
    self.n = 3
    self.x = flex.double([sdfac, sdb, sdadd])
    self.starting_simplex = []
    if seed is None:
      random_func = flex.random_double
    else:
      print >> self.log, "Using random seed %d"%seed
      mt = flex.mersenne_twister(seed)
      random_func = mt.random_double

    for i in xrange(self.n+1):
      self.starting_simplex.append(random_func(self.n))

    self.optimizer = simplex_opt( dimension = self.n,
                                  matrix    = self.starting_simplex,
                                  evaluator = self,
                                  tolerance = 1e-1)
    self.x = self.optimizer.get_solution()

  def target(self, vector):
    """ Compute the functional by first applying the current values for the sd parameters
    to the input data, then computing the complete set of normalized deviations and finally
    using those normalized deviations to compute the functional."""
    from xfel import compute_normalized_deviations, apply_sd_error_params

    sdfac, sdb, sdadd = vector

    if sdfac < 0 or sdb < 0 or sdadd < 0:
      f = 1e6
    else:
      orig_isigi = self.data['isigi'] * 1
      apply_sd_error_params(self.data, sdfac, sdb, sdadd)
      all_sigmas_normalized = compute_normalized_deviations(self.data, self.indices)
      self.data['isigi'] = orig_isigi

      f = 0
      for bin in self.intensity_bin_selections:
        binned_normalized_sigmas = all_sigmas_normalized.select(bin)
        n = len(binned_normalized_sigmas)
        if n == 0: continue
        # weighting scheme from Evans, 2011
        w = math.sqrt(n)
        # functional is weight * (1-rms(normalized_sigmas))^s summed over all intensitiy bins
        f += w * ((1-math.sqrt(flex.mean(binned_normalized_sigmas*binned_normalized_sigmas)))**2)

    print >> self.log, "f: % 12.1f, sdfac: %8.5f, sdb: %8.5f, sdadd: %8.5f"%(f, sdfac, sdb, sdadd)
    return f

class sdfac_refine_refltable(sdfac_refine):
  def get_binned_intensities(self, n_bins=100):
    """
    Using self.ISIGI, bin the intensities using the following procedure:
    1) Find the minimum and maximum intensity values.
    2) Divide max-min by n_bins. This is the bin step size
    The effect is
    @param n_bins number of bins to use.
    @return a tuple with an array of selections for each bin and an array of median
    intensity values for each bin.
    """
    print >> self.log, "Computing intensity bins.",
    ISIGI = self.scaler.ISIGI
    sumI = flex.double(len(self.scaler.miller_set.indices()), 0)
    n_refl = flex.double(len(self.scaler.miller_set.indices()), 0)
    for i in xrange(len(ISIGI)):
      hkl_id = ISIGI['miller_id'][i]
      sumI[hkl_id] += ISIGI['scaled_intensity'][i]
      n_refl[hkl_id] += 1
    sel = n_refl > 0
    meanI = flex.double(len(sumI), 0)
    meanI.set_selected(sel, sumI.select(sel)/n_refl.select(sel))

    all_mean_Is = flex.double(len(ISIGI), 0)
    for i in xrange(len(ISIGI)):
      hkl_id = ISIGI['miller_id'][i]
      all_mean_Is[i] = sumI[hkl_id]/n_refl[hkl_id]

    min_meanI = flex.min(meanI)
    step = (flex.max(meanI)-min_meanI)/n_bins
    print >> self.log, "Bin size:", step

    sels = []
    binned_intensities = []
    for i in xrange(n_bins):
      sel = (all_mean_Is > (min_meanI + step * i)) & (all_mean_Is < (min_meanI + step * (i+1)))
      if sel.all_eq(False): continue
      sels.append(sel)
      binned_intensities.append((step/2 + step*i)+min(meanI))

    for i, (sel, intensity) in enumerate(zip(sels, binned_intensities)):
      print >> self.log, "Bin %02d, number of observations: % 10d, midpoint intensity: %f"%(i, sel.count(True), intensity)

    return sels, binned_intensities

  def adjust_errors(self):
    """
    Adjust sigmas according to Evans, 2011 Acta D and Evans and Murshudov, 2013 Acta D
    """
    print >> self.log, "Starting adjust_errors"
    print >> self.log, "Computing initial estimates of sdfac, sdb and sdadd"
    sdfac, sdb, sdadd = self.get_initial_sdparams_estimates()

    print >> self.log, "Initial estimates:", sdfac, sdb, sdadd

    from xfel import compute_normalized_deviations, apply_sd_error_params

    print >> self.log, "Refining error correction parameters sdfac, sdb, and sdadd"
    sels, binned_intensities = self.get_binned_intensities()
    seed = self.scaler.params.raw_data.error_models.sdfac_refine.random_seed
    minimizer = simplex_minimizer_refltable(sdfac, sdb, sdadd, self.scaler.ISIGI, self.scaler.miller_set.indices(), sels, seed, self.log)
    sdfac, sdb, sdadd = minimizer.x
    print >> self.log, "Final sdadd: %8.5f, sdb: %8.5f, sdadd: %8.5f"%(sdfac, sdb, sdadd)

    print >> self.log, "Applying sdfac/sdb/sdadd 1"
    apply_sd_error_params(self.scaler.ISIGI, sdfac, sdb, sdadd)

    self.scaler.summed_weight= flex.double(self.scaler.n_refl, 0.)
    self.scaler.summed_wt_I  = flex.double(self.scaler.n_refl, 0.)

    print >> self.log, "Applying sdfac/sdb/sdadd 2"
    for i in xrange(len(self.scaler.ISIGI)):
      hkl_id = self.scaler.ISIGI['miller_id'][i]
      Intensity = self.scaler.ISIGI['scaled_intensity'][i] # scaled intensity
      sigma = Intensity / self.scaler.ISIGI['isigi'][i] # corrected sigma
      variance = sigma * sigma
      self.scaler.summed_wt_I[hkl_id] += Intensity / variance
      self.scaler.summed_weight[hkl_id] += 1 / variance

    if False:
      # validate using http://ccp4wiki.org/~ccp4wiki/wiki/index.php?title=Symmetry%2C_Scale%2C_Merge#Analysis_of_Standard_Deviations
      print >> self.log, "Validating"
      from matplotlib import pyplot as plt
      all_sigmas_normalized = compute_normalized_deviations(self.scaler.ISIGI, self.scaler.miller_set.indices())
      plt.hist(all_sigmas_normalized, bins=100)
      plt.figure()

      binned_rms_normalized_sigmas = []

      for i, sel in enumerate(sels):
        binned_rms_normalized_sigmas.append(math.sqrt(flex.mean(all_sigmas_normalized.select(sel)*all_sigmas_normalized.select(sel))))

      plt.plot(binned_intensities, binned_rms_normalized_sigmas, 'o')
      plt.show()

      self.normal_probability_plot(all_sigmas_normalized, (-0.5, 0.5), plot = True)
