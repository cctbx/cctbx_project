""" Fit of D_model(s; theta) (mmtbx.refinement.llgi_e_dmodel) as sigmaA
against the E-scale LLGI target on the R-free set, with Emodel fixed: the
default sigmaA model of mmtbx.refinement.llgi_e_sigmaa. """
from __future__ import absolute_import, division, print_function
import numpy as np
import scitbx.minimizers
from cctbx.array_family import flex
from cctbx.xray import ext as xray_ext
from libtbx import group_args
import iotbx.phil
import mmtbx.refinement.llgi_e_dmodel as dmodel

llgi_e_dmodel_params = iotbx.phil.parse("""\
  n_gaussian_terms = 2
    .type = int
    .short_caption = Number of coordinate-error Gaussian terms (K)
    .help = "Number of decaying Gaussian terms exp(-B_k*s^2) in D_model, "\
            "with B_k on a fixed log-spaced ladder over the resolution "\
            "range of the data; they stand for the spread of coordinate "\
            "errors in the model."
  max_iterations = 200
    .type = int
    .expert_level = 3
    .help = "Maximum L-BFGS-B iterations for the D_model fit."
  include_constant_term = True
    .type = bool
    .short_caption = Add a constant (B=0) term to the D_model ladder
    .help = "Prepend a B=0 rung (a resolution-independent amplitude) to "\
            "the n_gaussian_terms log-spaced B_k ladder. The ladder's "\
            "lowest decaying rung (4*d_min^2) has already fallen to 1/e "\
            "at d_min, so without this term D_model cannot stay flat "\
            "at high resolution, as sigmaA does for a well-refined "\
            "model. With it, D_model need not fall to 0 "\
            "at infinite resolution."
""")

def default_b_k_grid(k, s2):
  """ K decay constants B_k, log-spaced from 1/max(s2) to 1/min(s2): B_k
  smaller than this behave as constants over the data, larger ones are
  zero everywhere in it. s2 = d*^2/4 of the fitted reflections. Falls
  back to 10..300 if s2 has no spread. Returns an ascending numpy array.
  """
  if(k <= 0):
    return np.zeros(0, dtype=float)
  s2 = np.asarray(s2, dtype=float)
  s2_pos = s2[s2 > 0.0]
  if(s2_pos.size == 0 or float(np.min(s2_pos)) == float(np.max(s2_pos))):
    b_lo, b_hi = 10.0, 300.0
  else:
    s2_min, s2_max = float(np.min(s2_pos)), float(np.max(s2_pos))
    b_lo, b_hi = 1.0 / s2_max, 1.0 / s2_min
  if(k == 1):
    return np.array([np.sqrt(b_lo * b_hi)], dtype=float)
  return np.exp(np.linspace(np.log(b_lo), np.log(b_hi), k))

def target_and_gradient(theta, s2, e_eff, e_model, dobs, centric_flags,
      b_k_grid, hybrid=None):
  """ Mean E-scale LLGI target (minimize-me) over the given reflections at
  sigmaA = D_model(s2; theta), and its gradient with respect to theta.
  s2 is a numpy array; e_eff, e_model, dobs (flex.double), centric_flags
  (flex.bool) and hybrid (cctbx.xray.llgi_hybrid or None) are as for
  cctbx.xray.ext.llgi_e_sigmaa_target_and_gradients, which computes the
  target and its derivative with respect to each reflection's sigmaA.
  """
  result = xray_ext.llgi_e_sigmaa_target_and_gradients(
    e_eff=e_eff,
    selection=flex.bool(e_eff.size(), True),
    e_model=e_model,
    dobs=dobs,
    sigmaa=flex.double(dmodel.d_model(s2, theta, b_k_grid)),
    centric_flags=centric_flags,
    hybrid=hybrid)
  gradient = dmodel.d_model_gradient(s2, theta, b_k_grid).dot(
    result.d_target_by_dsigmaa().as_numpy_array())
  return result.target(), gradient

class d_model_target_evaluator(object):
  """ L-BFGS-B fit of theta for D_model against the mean E-scale LLGI on
  the R-free set, with Emodel fixed. Bounds: 0 <= a_k, b <= amplitude_max,
  b_defect_min <= B_defect <= b_defect_max. Fitting theta directly with
  bounds (rather than e.g. ln(theta) without) lets amplitudes leave 0.

  B_defect is fixed at b_sol_anchor (bss's B_sol estimate) when given;
  otherwise it is fitted. Free, it trades off against coordinate-error
  terms with a similar decay (e.g. a_1 = 233, b = 232), and the result
  depends on the start.
  """

  b_defect_min = 0.1
  b_defect_max = 1.e4
  # tanh(D_raw) is saturated well before D_raw reaches this, so larger
  # amplitudes would only drift along flat directions.
  amplitude_max = 100.

  def __init__(self,
        e_eff, r_free_flags, e_model, dobs, centric_flags, d_star_sq,
        n_gaussian_terms=2, max_iterations=200,
        b_sol_anchor=None, b_k_grid=None,
        include_constant_term=True,
        hybrid=None):
    n_refl = e_eff.size()
    assert r_free_flags.size() == n_refl
    assert e_model.size() == n_refl
    assert dobs.size() == n_refl
    assert centric_flags.size() == n_refl
    assert d_star_sq.size() == n_refl
    if(r_free_flags.count(True) == 0):
      raise RuntimeError(
        "d_model_target_evaluator: no R-free/test-set reflections "
        "available for the D_model(s) LLGI sigmaA fit.")
    self.s2 = (d_star_sq.select(r_free_flags) / 4).as_numpy_array()
    self.e_eff = e_eff.select(r_free_flags)
    self.e_model = e_model.select(r_free_flags)
    self.dobs = dobs.select(r_free_flags)
    self.centric_flags = centric_flags.select(r_free_flags)
    self.hybrid = None
    if(hybrid is not None):
      self.hybrid = hybrid.select(r_free_flags)
    if(b_k_grid is None):
      b_k_grid = default_b_k_grid(n_gaussian_terms, self.s2)
      if(include_constant_term):
        b_k_grid = np.concatenate([[0.0], b_k_grid])
    self.b_k_grid = np.asarray(b_k_grid, dtype=float)
    n_terms = n_gaussian_terms + int(include_constant_term)
    assert self.b_k_grid.shape == (n_terms,)
    self.b_defect_fixed = None
    if(b_sol_anchor is not None and b_sol_anchor > 0):
      self.b_defect_fixed = float(b_sol_anchor)
    theta_start = self._default_theta_start(n_terms)
    # Fitted parameters: all of theta, or all but B_defect when it is fixed
    self.n_fit = n_terms + 2 - int(self.b_defect_fixed is not None)
    lower = np.zeros(self.n_fit)
    upper = np.full(self.n_fit, self.amplitude_max)
    if(self.b_defect_fixed is None):
      lower[-1] = self.b_defect_min
      upper[-1] = self.b_defect_max
    self.x = flex.double(np.clip(theta_start[:self.n_fit], lower, upper))
    self.bound_flags = flex.int(self.n_fit, 2)  # lower and upper bounds
    self.lower_bound = flex.double(lower)
    self.upper_bound = flex.double(upper)
    self.minimizer = scitbx.minimizers.lbfgs(
      mode="lbfgsb", calculator=self, max_iterations=max_iterations)
    self.update(self.minimizer.x)

  @staticmethod
  def _default_theta_start(k):
    """ Each a_k = 0.5/K, b = 0.05, B_defect = 40 (used only when B_defect
    is fitted).
    """
    theta = np.empty(k + 2, dtype=float)
    if(k > 0):
      theta[:k] = 0.5 / k
    theta[-2] = 0.05
    theta[-1] = 40.0
    return theta

  def _full_theta(self, theta_fit):
    theta_fit = np.asarray(theta_fit, dtype=float)
    if(self.b_defect_fixed is None): return theta_fit
    return np.concatenate([theta_fit, [self.b_defect_fixed]])

  # calculator interface for scitbx.minimizers.lbfgs
  def update(self, x):
    self.x = x
    self._f, g = target_and_gradient(
      self.theta(), self.s2, self.e_eff, self.e_model, self.dobs,
      self.centric_flags, self.b_k_grid, hybrid=self.hybrid)
    self._g = g[:self.n_fit]

  def target(self):
    return self._f

  def gradients(self):
    return flex.double(self._g)

  def theta(self):
    return self._full_theta(np.array(self.x))

  def lbfgs_error(self):
    """ None, or why L-BFGS-B stopped if it was not convergence or the
    iteration limit. """
    m = self.minimizer.minimizer
    if(m.error is not None): return m.error
    task = m.task()
    if(task.startswith("ABNORMAL") or task.startswith("ERROR")): return task
    return None

def estimate_d_model_sigmaa(e_eff, r_free_flags, e_model, dobs,
      centric_flags, d_star_sq, n_gaussian_terms=2, max_iterations=200,
      b_sol_anchor=None, b_k_grid=None, include_constant_term=True,
      hybrid=None):
  """ Fit D_model against the E-scale LLGI on the R-free set (see
  d_model_target_evaluator) and evaluate it at every reflection.

  b_k_grid: the fixed B_k ladder, or None to derive it from the R-free
  reflections (default_b_k_grid, with a leading B=0 if
  include_constant_term).

  Returns a group_args with .sigmaa (flex.double, every reflection),
  .theta ([a_1..a_K, b, B_defect]) and .b_k_grid (flex.double), .target
  (final mean target on the R-free set), .evaluate_at (the fitted curve
  at other d*^2 values) and .lbfgs_error (see
  d_model_target_evaluator.lbfgs_error).
  """
  evaluator = d_model_target_evaluator(
    e_eff=e_eff, r_free_flags=r_free_flags, e_model=e_model, dobs=dobs,
    centric_flags=centric_flags, d_star_sq=d_star_sq,
    n_gaussian_terms=n_gaussian_terms,
    max_iterations=max_iterations,
    b_sol_anchor=b_sol_anchor,
    b_k_grid=b_k_grid,
    include_constant_term=include_constant_term,
    hybrid=hybrid)
  theta = evaluator.theta()
  b_k_grid_used = evaluator.b_k_grid
  def evaluate_at(d_star_sq_new):
    """ The fitted D_model at other d*^2 values (e.g. missing
    reflections). D_model is defined at any resolution, so no clamping.
    """
    s2_new = np.asarray(d_star_sq_new, dtype=float) / 4.0
    return flex.double(dmodel.d_model(s2_new, theta, b_k_grid_used))
  return group_args(
    sigmaa=evaluate_at(d_star_sq),
    theta=flex.double(theta),
    b_k_grid=flex.double(b_k_grid_used),
    target=evaluator.target(),
    evaluate_at=evaluate_at,
    lbfgs_error=evaluator.lbfgs_error())
