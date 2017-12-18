from __future__ import division
from scitbx.array_family import flex
import math

from xfel.merging.algorithms.error_model.sdfac_refine import sdfac_refine_refltable
from xfel.merging.algorithms.error_model import compute_normalized_deviations, apply_sd_error_params

def finite_difference(target, values, p):
  """ Compute finite difference given a target function """
  import copy
  functional = target(values)
  DELTA = 1.E-7
  tempvals = copy.deepcopy(values)
  tempvals.reference[p] += DELTA

  dfunctional = target(tempvals)
  #calculate by finite_difference
  finite_g = (dfunctional-functional )/DELTA
  return finite_g

class sdfac_refine_refltable_lbfgs(sdfac_refine_refltable):
  def __init__(self, scaler):
    sdfac_refine_refltable.__init__(self, scaler)
    self.parameterization = sdfac_sq_parameterization

  def get_initial_sdparams_estimates(self):
    sdfac, sdb, sdadd = super(sdfac_refine_refltable_lbfgs, self).get_initial_sdparams_estimates()
    sdadd = math.sqrt(sdadd)
    sdb = math.sqrt(sdadd)
    return sdfac, sdb, sdadd

  def run_minimzer(self, values, sels, **kwargs):
    # base class uses non-squared values, but lbfgs version refines the squares.
    values = self.parameterization(values.reference**2)
    refinery = sdfac_refinery(self.scaler.ISIGI, self.scaler.miller_set.indices(), sels, self.log)
    return lbfgs_minimizer(values.reference, self.parameterization, refinery, self.log)

from libtbx import adopt_init_args
from xfel.cxi.postrefinement_legacy_rs import unpack_base
class sdfac_sq_parameterization(unpack_base):
  def __getattr__(YY,item):
    try:
      if item=="SDFAC" : return math.sqrt(YY.reference[0])
      if item=="SDB"   : return math.sqrt(YY.reference[1])
      if item=="SDADD" : return math.sqrt(YY.reference[2])
    except ValueError:
      return float("nan")
    if item=="SDFACSQ" : return YY.reference[0]
    if item=="SDBSQ"   : return YY.reference[1]
    if item=="SDADDSQ" : return YY.reference[2]
    raise AttributeError(item)

  def show(YY, out):
    print >> out, "Sdfac^2: % 15.12f"%YY.SDFACSQ,
    print >> out, "SdB^2: % 15.12f"%YY.SDBSQ,
    print >> out, "Sdadd^2: % 15.12f"%YY.SDADDSQ

class sdfac_refinery(object):
  def __init__(self, ISIGI, indices, bins, log):
    adopt_init_args(self, locals())

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
    self.apply_sd_error_params(self.ISIGI, values)
    all_sigmas_normalized = compute_normalized_deviations(self.ISIGI, self.indices)

    sigma_prime = self.ISIGI['scaled_intensity'] / self.ISIGI['isigi']

    self.ISIGI['isigi'] = orig_isigi
    return all_sigmas_normalized, sigma_prime

  def df_dpsq(self, values, all_sigmas_normalized, sigma_prime, dsigmasq_dpsq, p = None):
    c = self.ISIGI['nn']*((self.ISIGI['scaled_intensity']-self.ISIGI['meanprime_scaled_intensity'])**2)

    dsigmanormsq_dpsq = ( -c / ((sigma_prime**2)**2)) * dsigmasq_dpsq

    g = flex.double(len(self.bins), 0)
    for b, bin in enumerate(self.bins):
      binned_normalized_sigmas = all_sigmas_normalized.select(bin)
      n = len(binned_normalized_sigmas)
      if n == 0: continue
      if binned_normalized_sigmas.count(0) == n: continue

      bnssq = binned_normalized_sigmas * binned_normalized_sigmas

      t1 = 2*(1 - math.sqrt(flex.sum(bnssq)/n))
      t2 = -0.5 / math.sqrt(flex.sum(bnssq)/n)
      t3 = flex.sum(dsigmanormsq_dpsq.select(bin))/n
      g[b] = t1 * t2 * t3
    return g

  def df_dsdfacsq(self, values, all_sigmas_normalized, sigma_prime):
    sigma = self.ISIGI['scaled_intensity'] / self.ISIGI['isigi']
    imean = self.ISIGI['mean_scaled_intensity']

    dsigmasq_dsdfacsq = sigma**2 + values.SDBSQ * imean + values.SDADDSQ * imean**2

    return self.df_dpsq(values, all_sigmas_normalized, sigma_prime, dsigmasq_dsdfacsq, 0)

  def df_dsdbsq(self, values, all_sigmas_normalized, sigma_prime):
    imean = self.ISIGI['mean_scaled_intensity']

    dsigmasq_dsddbsq = values.SDFACSQ * imean

    return self.df_dpsq(values, all_sigmas_normalized, sigma_prime, dsigmasq_dsddbsq, 1)

  def df_dsaddbsq(self, values, all_sigmas_normalized, sigma_prime):
    imean = self.ISIGI['mean_scaled_intensity']

    dsigmasq_dsdaddsq = values.SDFACSQ * imean**2

    return self.df_dpsq(values, all_sigmas_normalized, sigma_prime, dsigmasq_dsdaddsq, 2)

class lbfgs_minimizer(object):
  def __init__(self, current_x=None, parameterization=None, refinery=None,
               ISIGI = None, indices = None, bins = None, out=None,
               min_iterations=0, max_calls=1000, max_drop_eps=1.e-5):
    adopt_init_args(self, locals())
    self.n = current_x.size()
    self.x = current_x
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

  def compute_functional_and_gradients(self):
    values = self.parameterization(self.x)
    fvec = self.refinery.fvec_callable(values)
    self.func = functional = self.refinery.functional(fvec)
    self.f = functional
    finite_g = flex.double()
    for x in xrange(self.n):
      finite_g.append(finite_difference(
        lambda v: self.refinery.functional(self.refinery.fvec_callable(v)),
        values, x))

    self.g = self.refinery.gradients(values)

    for x in xrange(self.n):
      print >> self.out, "p%d finite % 20.7f analytical % 20.7f"%(x, finite_g[x], self.g[x])

    print >> self.out, "functional value % 20.3f"%functional,
    values.show(self.out)
    return self.f, self.g

  def get_refined_params(self):
    return self.parameterization(self.x)

  def apply_sd_error_params(self, data, values):
    self.refinery.apply_sd_error_params(data, values)

  def __del__(self):
    values = self.parameterization(self.x)
    print >> self.out, "FINALMODEL",
    print >> self.out, "functional value % 20.3f"%self.func,
    values.show(self.out)
