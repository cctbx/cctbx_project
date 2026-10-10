""" ScatFrac(resolution), the model's own normalisation for the F-scale LLGI
target: the fraction of the expected intensity the model accounts for, so
that Ecalc = k Fcalc/(sqrt(TEPS ScatFrac) RESN) has unit mean square. It
enters the target's mean only (the variance is TEPS - (Dobs sigmaA)^2, as
in phasertng). Estimated from the model's moments (default) or fitted
against the F-scale target with sigmaA fixed (see mmtbx.refinement.
llgi_e_sigmaa.estimate_sigmaa_e_then_scatfrac_f). """
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import group_args
import iotbx.phil
import math
import mmtbx.refinement.llgi_tncs as llgi_tncs

llgi_scatfrac_params = iotbx.phil.parse("""\
  method = *moments likelihood
    .type = choice
    .short_caption = ScatFrac estimate
    .help = "moments: a smoothed <k^2 |F_model|^2/(TEPS RESN^2)> by "\
            "resolution. likelihood: ScatFrac_inf exp(-B_scatfrac ss) "\
            "fitted against the F-scale target on the working set."
  moments_n_bins = 20
    .type = int
    .expert_level = 3
    .help = "Resolution bins (equal counts) for the moments ScatFrac, "\
            "interpolated linearly in d*^2 between bin means."
  max_iterations = 100
    .type = int
    .expert_level = 3
    .help = "Maximum L-BFGS iterations for the likelihood ScatFrac fit."
  scatfrac_b_factor_restraint_sigma = 0.0
    .type = float
    .short_caption = ScatFrac B-factor restraint sigma (Angstrom^2)
    .help = "Likelihood ScatFrac fit: width of the restraint "\
            "0.5*(B_scatfrac/sigma)^2 on the B-factor of ScatFrac = "\
            "ScatFrac_inf*exp(-B_scatfrac*ss), added to the mean target "\
            "per reflection. 0 (the default) disables it."
""")

def scatfrac_ratio_of_sums(f_calc, teps, resn, scale_factor=1.0):
  """ Overall ScatFrac as k^2*sum(|Fcalc|^2)/sum(TEPS*RESN^2), the starting
  value for the likelihood fit. Summing before dividing keeps a single
  large |Fcalc| at small RESN from dominating, as it would in a mean of
  per-reflection ratios. Returns 1 if the estimate is not positive.
  """
  denominator = flex.sum(teps * resn * resn)
  numerator = scale_factor ** 2 * flex.sum(flex.norm(f_calc))
  if(not (denominator > 0 and numerator > 0)):
    return 1.0
  return numerator / denominator

def _b_factor_restraint_penalty_and_gradient(b_scatfrac, sigma):
  """ Restraint 0.5*(b_scatfrac/sigma)^2 toward 0; sigma <= 0 disables it.
  Returns (penalty, d(penalty)/d(b_scatfrac)).
  """
  if(sigma <= 0):
    return 0.0, 0.0
  return 0.5 * (b_scatfrac / sigma) ** 2, b_scatfrac / (sigma * sigma)

class llgi_scatfrac_b_factor_target_evaluator(object):
  """ scitbx.lbfgs target evaluator for ScatFrac = ScatFrac_inf*exp(-B*ss)
  (ss = d*^2/4, so B_scatfrac is on the scale of a crystallographic B, and
  may have either sign), against the mean F-scale LLGI over selection,
  with sigmaA fixed. ScatFrac_inf and B_scatfrac are strongly correlated
  over a finite resolution range (as scale and B are in a Wilson plot), so
  the parameters are x = [ln ScatFrac(ss_ref), B_scatfrac*ss_sd], with
  ss_ref and ss_sd the mean and standard deviation of ss over the
  selected reflections: this removes most of the correlation and puts
  both parameters on the same scale. Only the curve over the observed
  range is well determined.
  There is no upper bound: Feff need not be on an absolute scale, so
  ScatFrac can exceed 1. B_scatfrac is optionally restrained toward 0
  (restraint_sigma).
  """

  def __init__(self,
        f_eff,
        selection,
        f_calc,
        dobs,
        sigmaa,
        teps,
        resn,
        ss,
        centric_flags,
        scale_factor,
        scatfrac_inf_start,
        b_scatfrac_start,
        restraint_sigma=0.0,
        max_iterations=100,
        hybrid=None,
        tncs_u=None):
    self.hybrid = hybrid
    self.tncs_u = tncs_u
    self.f_eff = f_eff
    self.selection = selection
    self.f_calc = f_calc
    self.dobs = dobs
    self.sigmaa = sigmaa
    self.teps = teps
    self.resn = resn
    self.centric_flags = centric_flags
    self.scale_factor = scale_factor
    self.restraint_sigma = restraint_sigma
    assert ss.size() == f_eff.size()
    self.ss = ss
    ss_sel = ss.select(selection)
    self.ss_ref = flex.mean(ss_sel) if ss_sel.size() else 0.
    self.ss_sd = max(ss_sel.standard_deviation_of_the_sample()
      if ss_sel.size() > 1 else 0., 1.e-6)
    self.x = self.x_from(scatfrac_inf_start, b_scatfrac_start)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations))

  def compute_functional_and_gradients(self):
    scatfrac = self.scatfrac()
    result = llgi_tncs.sigmaa_scatfrac_functor(
      tncs_u=self.tncs_u,
      f_eff=self.f_eff,
      selection=self.selection,
      f_calc=self.f_calc,
      dobs=self.dobs,
      sigmaa=self.sigmaa,
      scatfrac=scatfrac,
      scale_factor=self.scale_factor,
      teps=self.teps,
      resn=self.resn,
      centric_flags=self.centric_flags, hybrid=self.hybrid)
    # d(target)/dz for each reflection, z = ln ScatFrac; dz/dx[0] = 1,
    # dz/dx[1] = -(ss - ss_ref)/ss_sd
    g_z = result.d_target_by_dscatfrac() * scatfrac
    b = self.x[1] / self.ss_sd
    penalty, penalty_grad = _b_factor_restraint_penalty_and_gradient(
      b, self.restraint_sigma)
    return result.target() + penalty, flex.double([
      flex.sum(g_z),
      (-flex.sum(g_z * (self.ss - self.ss_ref)) + penalty_grad) / self.ss_sd])

  def x_from(self, scatfrac_inf, b_scatfrac):
    """ the parameters for ScatFrac = scatfrac_inf*exp(-b_scatfrac*ss) """
    return flex.double([math.log(scatfrac_inf) - b_scatfrac * self.ss_ref,
      b_scatfrac * self.ss_sd])

  def scatfrac(self):
    """ The fitted ScatFrac for every reflection. """
    b = self.x[1] / self.ss_sd
    return flex.exp(self.x[0] - b * (self.ss - self.ss_ref))

  def scatfrac_inf_and_b(self):
    b = self.x[1] / self.ss_sd
    return math.exp(self.x[0] + b * self.ss_ref), b

def scatfrac_moments(f_calc, teps, resn, d_star_sq, n_bins=20,
      scale_factor=1.0):
  """ ScatFrac as the model's own normalisation: <k^2 |F_calc|^2/(TEPS
  RESN^2)> in n_bins equal-count resolution bins, interpolated linearly
  in d*^2 between bin means (constant beyond the outermost means). """
  import numpy as np
  k = scale_factor if scale_factor > 0 else 1.0
  x = np.array(flex.norm(f_calc) * (k * k) / (teps * resn * resn))
  dss = np.array(d_star_sq)
  order = np.argsort(dss)
  bins = [c for c in np.array_split(order, max(1, min(n_bins, x.size)))
          if c.size > 0]
  x_mean = np.array([x[c].mean() for c in bins])
  d_mean = np.array([dss[c].mean() for c in bins])
  return flex.double(np.interp(dss, d_mean, x_mean))

def estimate_llgi_scatfrac_moments(
      f_eff,
      working_selection,
      f_calc,
      dobs,
      sigmaa,
      teps,
      resn,
      centric_flags,
      d_star_sq,
      scale_factor=1.0,
      params=None,
      hybrid=None,
      tncs_u=None):
  """ ScatFrac from the model's moments (scatfrac_moments), with the
  F-scale target on working_selection. Returns a group_args like
  estimate_llgi_scatfrac_likelihood's, with scatfrac_inf and b_scatfrac
  None. """
  if(params is None):
    params = llgi_scatfrac_params.extract()
  scatfrac = scatfrac_moments(f_calc, teps, resn, d_star_sq,
    n_bins=params.moments_n_bins, scale_factor=scale_factor)
  final_result = llgi_tncs.sigmaa_scatfrac_functor(tncs_u=tncs_u,
    f_eff=f_eff, selection=working_selection, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, scatfrac=scatfrac, scale_factor=scale_factor, teps=teps,
    resn=resn, centric_flags=centric_flags, hybrid=hybrid)
  return group_args(
    scatfrac=scatfrac, target=final_result.target(),
    scatfrac_inf=None, b_scatfrac=None, lbfgs_error=None)

def estimate_llgi_scatfrac(params=None, **args):
  """ estimate_llgi_scatfrac_moments or estimate_llgi_scatfrac_likelihood,
  as params.method says """
  if(params is None):
    params = llgi_scatfrac_params.extract()
  if(params.method == "moments"):
    return estimate_llgi_scatfrac_moments(params=params, **args)
  return estimate_llgi_scatfrac_likelihood(params=params, **args)

def estimate_llgi_scatfrac_likelihood(
      f_eff,
      working_selection,
      f_calc,
      dobs,
      sigmaa,
      teps,
      resn,
      centric_flags,
      d_star_sq,
      scale_factor=1.0,
      params=None,
      hybrid=None,
      tncs_u=None):
  """ Fit ScatFrac = ScatFrac_inf*exp(-B_scatfrac*ss) against the F-scale
  LLGI on working_selection (pass ~r_free_flags), with sigmaA fixed (see
  llgi_scatfrac_b_factor_target_evaluator). Starts from B_scatfrac = 0
  and scatfrac_ratio_of_sums.

  Returns a group_args with .scatfrac (flex.double, every reflection),
  .target (final mean target on the working set), .scatfrac_inf and
  .b_scatfrac, and .lbfgs_error (None, or the message L-BFGS stopped
  with).

  tncs_u: the tNCS error-variance modulation (the fit is in the u-frame:
  f_eff = Eeff RESN sqrt(u) and the hybrid's observations normalised by
  u; teps must then be 1), or None.
  """
  if(hybrid is not None):
    # Exact likelihood wherever possible, as for the sigmaA fit.
    hybrid = hybrid.with_rice_kappa(0.0)
  if(params is None):
    params = llgi_scatfrac_params.extract()
  n_refl = f_eff.size()
  assert working_selection.size() == n_refl
  assert f_calc.size() == n_refl
  assert dobs.size() == n_refl
  assert sigmaa.size() == n_refl
  assert teps.size() == n_refl
  assert resn.size() == n_refl
  assert centric_flags.size() == n_refl
  assert d_star_sq.size() == n_refl
  if(working_selection.count(True) == 0):
    raise RuntimeError(
      "No working-set reflections available for the LLGI ScatFrac fit.")
  evaluator = llgi_scatfrac_b_factor_target_evaluator(
    f_eff=f_eff,
    selection=working_selection,
    f_calc=f_calc,
    dobs=dobs,
    sigmaa=sigmaa,
    teps=teps,
    resn=resn,
    ss=d_star_sq / 4.0,
    centric_flags=centric_flags,
    scale_factor=scale_factor,
    scatfrac_inf_start=scatfrac_ratio_of_sums(
      f_calc=f_calc, teps=teps, resn=resn, scale_factor=scale_factor),
    b_scatfrac_start=0.0,
    restraint_sigma=params.scatfrac_b_factor_restraint_sigma,
    max_iterations=params.max_iterations, hybrid=hybrid, tncs_u=tncs_u)
  scatfrac_inf, b_scatfrac = evaluator.scatfrac_inf_and_b()
  scatfrac = evaluator.scatfrac()
  final_result = llgi_tncs.sigmaa_scatfrac_functor(tncs_u=tncs_u,
    f_eff=f_eff, selection=working_selection, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, scatfrac=scatfrac, scale_factor=scale_factor,
    teps=teps, resn=resn, centric_flags=centric_flags, hybrid=hybrid)
  return group_args(
    scatfrac=scatfrac, target=final_result.target(),
    scatfrac_inf=scatfrac_inf, b_scatfrac=b_scatfrac,
    lbfgs_error=evaluator.minimizer.error)
