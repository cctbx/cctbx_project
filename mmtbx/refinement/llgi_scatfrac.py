""" ScatFrac(resolution), the fraction of the scattering accounted for by
the model, fitted against the F-scale LLGI target with sigmaA fixed (see
mmtbx.refinement.llgi_e_sigmaa.estimate_sigmaa_e_then_scatfrac_f). """
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.xray import ext
import scitbx.lbfgs
from libtbx import group_args
import iotbx.phil
import math

llgi_scatfrac_params = iotbx.phil.parse("""\
  max_iterations = 100
    .type = int
    .expert_level = 3
    .help = "Maximum L-BFGS iterations for the ScatFrac fit."
  scatfrac_b_factor_restraint_sigma = 10.0
    .type = float
    .short_caption = ScatFrac B-factor restraint sigma (Angstrom^2)
    .help = "Width of the restraint 0.5*(B_scatfrac/sigma)^2 on the "\
            "B-factor of ScatFrac = ScatFrac_inf*exp(-B_scatfrac*ss). "\
            "sigmaA and ScatFrac trade off against each other in the "\
            "F-scale target, so a weak restraint keeps B_scatfrac near 0 "\
            "unless the data call for a trend. 0 disables it."
""")

# The ScatFrac fit keeps the effective F-scale sigmaA, sigmaA*k/
# sqrt(ScatFrac), below A_EFF_MAX by a smooth floor on ln(ScatFrac) of
# width FLOOR_WIDTH (see llgi_scatfrac_b_factor_target_evaluator).
A_EFF_MAX = 0.995
FLOOR_WIDTH = 0.01

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

def ln_scatfrac_floor(sigmaa, scale_factor):
  """ ln of the ScatFrac at which sigmaA*k/sqrt(ScatFrac) = A_EFF_MAX, for
  each reflection; -1000 (no floor) where sigmaA <= 0. """
  k = scale_factor if scale_factor > 0 else 1.0  # as llgi::effective_sigmaa
  no_floor = ~(sigmaa > 0)
  sa = sigmaa.deep_copy()
  sa.set_selected(no_floor, 1.0)
  result = 2 * flex.log(sa * (k / A_EFF_MAX))
  result.set_selected(no_floor, -1.e3)
  return result

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
  with sigmaA fixed. Parameters x = [ln(ScatFrac_inf), B_scatfrac]. There
  is no upper bound: Feff need not be on an absolute scale, so ScatFrac can
  exceed 1. ScatFrac_inf and B_scatfrac are strongly correlated over a
  finite resolution range (as scale and B are in a Wilson plot), so only
  the curve over the observed range is well determined. B_scatfrac is
  restrained toward 0 (restraint_sigma).

  ScatFrac is kept above the floor at which the effective F-scale sigmaA,
  sigmaA*k/sqrt(ScatFrac), would reach A_EFF_MAX: an effective sigmaA
  above 1 is unphysical, and at 0.999 the hybrid LLGI switches from the
  exact to the Rice likelihood, which would make the objective
  discontinuous. The floor is smooth in z = ln(ScatFrac):
    z' = z_floor + w*ln(1 + exp((z - z_floor)/w)),  w = FLOOR_WIDTH,
  and changes ScatFrac by less than 1e-6 (relative) where the curve is
  more than 10% above the floor. Below the floor the target is flat in
  the curve, so the fit should start above it (see
  estimate_llgi_scatfrac_likelihood).
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
        hybrid=None):
    self.hybrid = hybrid
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
    self.z_floor = ln_scatfrac_floor(sigmaa, scale_factor)
    self.x = flex.double([math.log(scatfrac_inf_start), b_scatfrac_start])
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations))

  def _scatfrac_and_derivative(self):
    """ ScatFrac with the floor applied, and dz'/dz. """
    z = self.x[0] - self.x[1] * self.ss
    u = (z - self.z_floor) * (1. / FLOOR_WIDTH)
    e = flex.exp(-flex.abs(u))
    u_pos = u.deep_copy()
    u_pos.set_selected(u < 0, 0.0)
    z_eff = self.z_floor + FLOOR_WIDTH * (u_pos + flex.log(1 + e))
    negative = (u < 0).as_double()
    dz_eff_dz = (negative * e + (1 - negative)) / (1 + e)
    return flex.exp(z_eff), dz_eff_dz

  def compute_functional_and_gradients(self):
    scatfrac, dz_eff_dz = self._scatfrac_and_derivative()
    result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
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
    # d(target)/dz for each reflection; dz/d(ln ScatFrac_inf) = 1,
    # dz/dB = -ss
    g_z = result.d_target_by_dscatfrac() * scatfrac * dz_eff_dz
    penalty, penalty_grad = _b_factor_restraint_penalty_and_gradient(
      self.x[1], self.restraint_sigma)
    return result.target() + penalty, flex.double([
      flex.sum(g_z), -flex.sum(g_z * self.ss) + penalty_grad])

  def scatfrac(self):
    """ The fitted ScatFrac, floor applied, for every reflection. """
    return self._scatfrac_and_derivative()[0]

  def n_at_floor(self):
    """ Number of selected reflections whose curve is below the floor. """
    z = self.x[0] - self.x[1] * self.ss
    return ((z < self.z_floor) & self.selection).count(True)

  def scatfrac_inf_and_b(self):
    return math.exp(self.x[0]), self.x[1]

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
      hybrid=None):
  """ Fit ScatFrac = ScatFrac_inf*exp(-B_scatfrac*ss) against the F-scale
  LLGI on working_selection (pass ~r_free_flags), with sigmaA fixed (see
  llgi_scatfrac_b_factor_target_evaluator). Starts from B_scatfrac = 0
  and scatfrac_ratio_of_sums, raised if necessary to 1.1 times the highest
  floor among the working reflections, so that the fit does not start
  where the target is flat.

  Returns a group_args with .scatfrac (flex.double, every reflection),
  .target (final mean target on the working set), .scatfrac_inf and
  .b_scatfrac, .n_at_floor (working reflections held at the floor) and
  .lbfgs_error (None, or the message L-BFGS stopped with).
  """
  if(hybrid is not None):
    # Exact likelihood wherever possible. On the F scale the exact/Rice
    # switch depends on ScatFrac through the effective sigmaA, which would
    # make the objective discontinuous; with rice_kappa = 0 only the
    # effective sigmaA < 0.999 limit remains, which the floor in the
    # evaluator keeps clear of.
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
    scatfrac_inf_start=max(
      scatfrac_ratio_of_sums(
        f_calc=f_calc, teps=teps, resn=resn, scale_factor=scale_factor),
      1.1 * math.exp(flex.max(
        ln_scatfrac_floor(sigmaa, scale_factor).select(working_selection)))),
    b_scatfrac_start=0.0,
    restraint_sigma=params.scatfrac_b_factor_restraint_sigma,
    max_iterations=params.max_iterations, hybrid=hybrid)
  scatfrac_inf, b_scatfrac = evaluator.scatfrac_inf_and_b()
  scatfrac = evaluator.scatfrac()
  final_result = ext.llgi_sigmaa_scatfrac_target_and_gradients(
    f_eff=f_eff, selection=working_selection, f_calc=f_calc, dobs=dobs,
    sigmaa=sigmaa, scatfrac=scatfrac, scale_factor=scale_factor,
    teps=teps, resn=resn, centric_flags=centric_flags, hybrid=hybrid)
  return group_args(
    scatfrac=scatfrac, target=final_result.target(),
    scatfrac_inf=scatfrac_inf, b_scatfrac=b_scatfrac,
    n_at_floor=evaluator.n_at_floor(),
    lbfgs_error=evaluator.minimizer.error)
