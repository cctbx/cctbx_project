from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex
import math

from xfel.merging.algorithms.error_model.sdfac_refine import sdfac_refine_refltable
from xfel.merging.algorithms.error_model import compute_normalized_deviations, apply_sd_error_params, df_dpsq as df_dpsq_cpp

def finite_difference(target, values, p, DELTA = 1.e-7):
  """ Compute finite difference given a target function """
  import copy
  functional = target(values)
  tempvals = copy.deepcopy(values)
  tempvals.reference[p] += DELTA

  dfunctional = target(tempvals)
  #calculate by finite_difference
  finite_g = (dfunctional-functional )/DELTA
  return finite_g

class sdfac_refine_refltable_lbfgs(sdfac_refine_refltable):
  def get_initial_sdparams_estimates(self):
    sdfac, sdb, sdadd = super(sdfac_refine_refltable_lbfgs, self).get_initial_sdparams_estimates()
    sdadd = math.sqrt(sdadd)
    sdb = math.sqrt(sdadd)
    return sdfac, sdb, sdadd

  def run_minimzer(self, values, sels, **kwargs):
    refinery = sdfac_refinery(self.scaler, self, self.scaler.miller_set.indices(), sels, self.log)
    return lbfgs_minimizer(values.reference, self.parameterization, refinery, self.log,
      show_finite_differences = self.scaler.params.raw_data.error_models.sdfac_refine.show_finite_differences)

from libtbx import adopt_init_args
class sdfac_refinery(object):
  def __init__(self, scaler, modeler, indices, bins, log):
    adopt_init_args(self, locals())
    self.ISIGI = self.scaler.ISIGI

    self.weights = flex.double()
    for bin in bins:
      # weighting scheme from Evans, 2011
      n = bin.count(True)
      self.weights.append(math.sqrt(n))

  def fvec_callable(self, values):
    """ Compute the functional by first applying the current values for the sd parameters
    to the input data, then computing the complete set of normalized deviations and finally
    using those normalized deviations to compute the functional."""

    all_sigmas_normalized, _ = self.get_normalized_sigmas(values)
    f = flex.double()
    for i, bin in enumerate(self.bins):
      binned_normalized_sigmas = all_sigmas_normalized.select(bin)
      n = len(binned_normalized_sigmas)
      if n == 0:
        f.append(0)
        continue
      # functional is weight * (1-rms(normalized_sigmas))^s summed over all intensitiy bins
      f.append(1-math.sqrt(flex.mean(binned_normalized_sigmas*binned_normalized_sigmas)))

    if self.scaler.params.raw_data.error_models.sdfac_refine.plot_refinement_steps:
      all_sigmas_normalized = all_sigmas_normalized.select(all_sigmas_normalized != 0)
      print(self.modeler.normal_probability_plot(all_sigmas_normalized, (-0.5, 0.5), plot = True))
    return f

  def functional(self, fvec):
    return flex.sum(self.weights * fvec**2)

  def jacobian_callable(self, values):
    all_sigmas_normalized, sigma_prime = self.get_normalized_sigmas(values)
    df_dsdfacsq = self.df_dsdfacsq(values, all_sigmas_normalized, sigma_prime)
    df_dsdbsq   = self.df_dsdbsq(values, all_sigmas_normalized, sigma_prime)
    df_dsaddbsq = self.df_dsaddbsq(values, all_sigmas_normalized, sigma_prime)
    return [df_dsdfacsq, df_dsdbsq, df_dsaddbsq]

  def gradients(self, values):
    g = flex.double()
    for d in self.jacobian_callable(values):
      g.append(flex.sum(self.weights * d))
    return g

  def apply_sd_error_params(self, data, values):
    apply_sd_error_params(data, values.SDFACSQ, values.SDBSQ, values.SDADDSQ, True)

  def get_normalized_sigmas(self, values):
    orig_isigi = self.ISIGI['isigi'] * 1
    # time split between these two functions
    self.apply_sd_error_params(self.ISIGI, values)
    all_sigmas_normalized = compute_normalized_deviations(self.ISIGI, self.indices)

    sigma_prime = self.ISIGI['scaled_intensity'] / self.ISIGI['isigi']

    self.ISIGI['isigi'] = orig_isigi
    return all_sigmas_normalized, sigma_prime

  def df_dpsq(self, all_sigmas_normalized, sigma_prime, dsigmasq_dpsq):
    g = df_dpsq_cpp(all_sigmas_normalized, sigma_prime, dsigmasq_dpsq, self.ISIGI, self.modeler.bin_indices, len(self.bins))
    return g

  def df_dsdfacsq(self, values, all_sigmas_normalized, sigma_prime):
    sigma = self.ISIGI['scaled_intensity'] / self.ISIGI['isigi']
    imean = self.ISIGI['mean_scaled_intensity']

    dsigmasq_dsdfac = 2 * values.SDFAC
    dsigmasq_dsdfacsq = (sigma**2 + values.SDBSQ * imean + values.SDADDSQ * imean**2) * dsigmasq_dsdfac

    return self.df_dpsq(all_sigmas_normalized, sigma_prime, dsigmasq_dsdfacsq)

  def df_dsdbsq(self, values, all_sigmas_normalized, sigma_prime):
    imean = self.ISIGI['mean_scaled_intensity']

    dsigmasq_dsdb = 2 * values.SDB
    dsigmasq_dsdbsq = values.SDFACSQ * imean * dsigmasq_dsdb

    return self.df_dpsq(all_sigmas_normalized, sigma_prime, dsigmasq_dsdbsq)

  def df_dsaddbsq(self, values, all_sigmas_normalized, sigma_prime):
    imean = self.ISIGI['mean_scaled_intensity']

    dsigmasq_dsdadd = 2 * values.SDADD
    dsigmasq_dsdaddsq = values.SDFACSQ * imean**2 * dsigmasq_dsdadd

    return self.df_dpsq(all_sigmas_normalized, sigma_prime, dsigmasq_dsdaddsq)

class lbfgs_minimizer(object):
  def __init__(self, current_x=None, parameterization=None, refinery=None,
               ISIGI = None, indices = None, bins = None, out=None,
               min_iterations=0, max_calls=1000, max_drop_eps=1.e-10,
               show_finite_differences = False):
    adopt_init_args(self, locals())
    self.n = current_x.size()
    self.x = current_x

    if False:
      self.diag_mode = "always"
      from scitbx import lbfgs
      self.minimizer = lbfgs.run(
        target_evaluator=self,
        termination_params=lbfgs.termination_parameters(
          traditional_convergence_test=False,
          drop_convergence_test_max_drop_eps=max_drop_eps,
          min_iterations=min_iterations,
          max_iterations = None,
          max_calls=max_calls),
        exception_handling_params=lbfgs.exception_handling_parameters(
           ignore_line_search_failed_rounding_errors=True,
           ignore_line_search_failed_step_at_lower_bound=True,#the only change from default
           ignore_line_search_failed_step_at_upper_bound=False,
           ignore_line_search_failed_maxfev=False,
           ignore_line_search_failed_xtol=False,
           ignore_search_direction_not_descent=False)
        )
    else:
      from scitbx import lbfgsb
      l = flex.double(self.n, 1e-8)

      if len(l) > 3:
        for p in range(7,len(l)):
          l[p] = 1e-15 # g*

      self.minimizer = lbfgsb.minimizer(
        n = self.n,
        l = l,
        u = flex.double(self.n, 0),
        nbd = flex.int(self.n, 1),
      )
      while True:
        self.compute_functional_and_gradients()
        if self.minimizer.process(self.x, self.f, self.g):
          pass
        elif self.minimizer.is_terminated():
          break

  def compute_functional_and_gradients(self):
    values = self.parameterization(self.x)
    fvec = self.refinery.fvec_callable(values)
    self.f = functional = self.refinery.functional(fvec)
    self.g = self.refinery.gradients(values)

    if self.show_finite_differences:
      finite_g = flex.double()
      for x in range(self.n):
        finite_g.append(finite_difference(
          lambda v: self.refinery.functional(self.refinery.fvec_callable(v)),
          values, x))

      for x in range(self.n):
        print("p%d finite % 20.10f analytical % 20.10f"%(x, finite_g[x], self.g[x]), file=self.out)

    print("functional value % 20.3f"%functional, end=' ', file=self.out)
    values.show(self.out)
    return self.f, self.g

  def compute_functional_gradients_diag(self):
    values = self.parameterization(self.x)
    fvec = self.refinery.fvec_callable(values)
    self.f = functional = self.refinery.functional(fvec)

    self.g = flex.double()
    self.c = flex.double()
    for d in self.refinery.jacobian_callable(values):
      self.g.append(flex.sum(self.refinery.weights * d))
      self.c.append(flex.sum(self.refinery.weights * d * d))

    diag = 1/self.c

    if self.show_finite_differences:
      finite_g = flex.double()
      for x in range(self.n):
        finite_g.append(finite_difference(
          lambda v: self.refinery.functional(self.refinery.fvec_callable(v)),
          values, x))

      for x in range(self.n):
        print("p%d finite % 20.10f analytical % 20.10f"%(x, finite_g[x], self.g[x]), file=self.out)

    print("functional value % 20.3f"%functional, end=' ', file=self.out)
    values.show(self.out)
    return self.f, self.g, diag


  def get_refined_params(self):
    return self.parameterization(self.x)

  def apply_sd_error_params(self, data, values):
    self.refinery.apply_sd_error_params(data, values)

  def __del__(self):
    values = self.parameterization(self.x)
    print("FINALMODEL", end=' ', file=self.out)
    print("functional value % 20.3f"%self.f, end=' ', file=self.out)
    values.show(self.out)
