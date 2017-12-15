from __future__ import division
from scitbx.array_family import flex
import math

from xfel.merging.algorithms.error_model.sdfac_refine import sdfac_refine_refltable

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
  def run_minimzer(self, sdfac, sdb, sdadd, sels, **kwargs):
    assert kwargs['squared']
    refinery = sdfac_refinery(self.scaler.ISIGI, self.scaler.miller_set.indices(), sels, self.log)
    # note the power of 2! We refine the square of the sd terms, not the sd terms themselves
    return lbfgs_minimizer(flex.double([sdfac, sdb, sdadd])**2, sdfac_parameterization, refinery, self.log)

from libtbx import adopt_init_args
from xfel.cxi.postrefinement_legacy_rs import unpack_base
class sdfac_parameterization(unpack_base):
  def __getattr__(YY,item):
    if item=="SDFAC" : return math.sqrt(YY.reference[0])
    if item=="SDB"   : return math.sqrt(YY.reference[1])
    if item=="SDADD" : return math.sqrt(YY.reference[2])
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

  def fvec_callable(self, values):
    """ Compute the functional by first applying the current values for the sd parameters
    to the input data, then computing the complete set of normalized deviations and finally
    using those normalized deviations to compute the functional."""
    from xfel import compute_normalized_deviations, apply_sd_error_params

    if values.SDFACSQ < 0 or values.SDBSQ < 0 or values.SDADDSQ < 0:
      f = 1e6
    else:
      all_sigmas_normalized, _ = self.get_normalized_sigmas(values)
      f = 0
      for bin in self.bins:
        binned_normalized_sigmas = all_sigmas_normalized.select(bin)
        n = len(binned_normalized_sigmas)
        if n == 0: continue
        # weighting scheme from Evans, 2011
        w = math.sqrt(n)
        # functional is weight * (1-rms(normalized_sigmas))^s summed over all intensitiy bins
        f += w * ((1-math.sqrt(flex.mean(binned_normalized_sigmas*binned_normalized_sigmas)))**2)

    return f

  def gradients(self, values):
    if values.SDFACSQ < 0 or values.SDBSQ < 0 or values.SDADDSQ < 0:
      return flex.double(3, 0)
    all_sigmas_normalized, sigma_prime = self.get_normalized_sigmas(values)
    df_dsdfacsq = self.df_dsdfacsq(values, all_sigmas_normalized, sigma_prime)
    df_dsdbsq   = self.df_dsdbsq(values, all_sigmas_normalized, sigma_prime)
    df_dsaddbsq = self.df_dsaddbsq(values, all_sigmas_normalized, sigma_prime)
    return flex.double([df_dsdfacsq, df_dsdbsq, df_dsaddbsq])

  def get_normalized_sigmas(self, values):
    from xfel import compute_normalized_deviations, apply_sd_error_params

    orig_isigi = self.ISIGI['isigi'] * 1
    apply_sd_error_params(self.ISIGI, values.SDFAC, values.SDB, values.SDADD, True)
    all_sigmas_normalized = compute_normalized_deviations(self.ISIGI, self.indices)

    sigma_prime = self.ISIGI['scaled_intensity'] / self.ISIGI['isigi']

    self.ISIGI['isigi'] = orig_isigi
    return all_sigmas_normalized, sigma_prime

  def df_dpsq(self, values, all_sigmas_normalized, sigma_prime, dsigmasq_dpsq, p = None):
    c = self.ISIGI['nn']*((self.ISIGI['scaled_intensity']-self.ISIGI['meanprime_scaled_intensity'])**2)

    dsigmanormsq_dpsq = ( -c / ((sigma_prime**2)**2)) * dsigmasq_dpsq

    g = 0
    for bin in self.bins:
      binned_normalized_sigmas = all_sigmas_normalized.select(bin)
      n = len(binned_normalized_sigmas)
      if n == 0: continue
      if binned_normalized_sigmas.count(0) == n: continue

      bnssq = binned_normalized_sigmas * binned_normalized_sigmas
      # weighting scheme from Evans, 2011
      w = math.sqrt(n)

      t1 = 2*(1 - math.sqrt(flex.sum(bnssq)/n))
      t2 = -0.5 / math.sqrt(flex.sum(bnssq)/n)
      t3 = flex.sum(dsigmanormsq_dpsq.select(bin))/n
      g += w * t1 * t2 * t3
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
    self.func = self.refinery.fvec_callable(values)
    functional = self.func
    self.f = functional
    finite_g = flex.double()
    for x in xrange(self.n):
      finite_g.append(finite_difference(self.refinery.fvec_callable, values, x))

    self.g = self.refinery.gradients(values)

    for x in xrange(self.n):
      print >> self.out, "p%d finite % 20.7f analytical % 20.7f"%(x, finite_g[x], self.g[x])

    print >> self.out, "functional value % 20.3f"%self.func,
    values.show(self.out)
    return self.f, self.g

  def get_refined_params(self):
    values = self.parameterization(self.x)
    return values.SDFAC, values.SDB, values.SDADD

  def __del__(self):
    values = self.parameterization(self.x)
    print >> self.out, "FINALMODEL",
    print >> self.out, "functional value % 20.3f"%self.func,
    values.show(self.out)
