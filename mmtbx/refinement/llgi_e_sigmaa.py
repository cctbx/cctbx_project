""" E-scale quantities for the LLGI target (Eeff, Emodel) and the fit of
sigmaA(resolution) against the E-scale LLGI on the R-free set, followed
by the F-scale ScatFrac fit (estimate_sigmaa_e_then_scatfrac_f). """
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.refinement.llgi_tncs as llgi_tncs
from mmtbx import scaling
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
import scitbx.lbfgs
from libtbx import group_args
import iotbx.phil
import mmtbx.refinement.llgi_e_dmodel_fit as llgi_e_dmodel_fit
import mmtbx.refinement.llgi_model_anisotropy as llgi_model_anisotropy

llgi_e_sigmaa_params = iotbx.phil.parse("""\
  sigmaa_model = spline *d_model
    .type = choice
    .short_caption = E-scale sigmaA(resolution) functional form
    .help = "Functional form of sigmaA(resolution) fitted on the E scale. "\
            "d_model: a sum of Gaussian coordinate-error terms minus a "\
            "bulk-solvent defect term (see d_model_params). spline: a "\
            "B-spline over a bounded sigmoid with a curvature restraint "\
            "(n_sigmaa_coeffs, spline_degree, sigmaa_max_iterations, "\
            "sigmaa_curvature_weight)."
  d_model_params {
    include scope mmtbx.refinement.llgi_e_dmodel_fit.llgi_e_dmodel_params
  }
  model_anisotropy = off ellipsoid resolution_dependent *auto
    .type = choice
    .short_caption = Anisotropic normalisation of the model (LLGI)
    .help = "Normalisation of the model's E (for the sigmaA fit and the "\
            "LLGI maps) and of ScatFrac (for the target): off, isotropic "\
            "(SigmaP(d*^2), and ScatFrac in resolution bins); ellipsoid, "\
            "also an anisotropic factor exp(-q^T B q/4) with B trace-free "\
            "and invariant under the point group; resolution_dependent, "\
            "with B changing linearly with d*^2 as nacelle's anisotropy "\
            "correction can (mmtbx.refinement.llgi_model_anisotropy). Each "\
            "is fitted to the model's own intensities. auto: the same "\
            "family as nacelle's correction of the data (resolution_"\
            "dependent if its ANISOBETA changes with resolution, otherwise "\
            "ellipsoid). ScatFrac is anisotropic only in its moments form "\
            "(llgi_scatfrac method=moments, the default)."
  n_sigmap_nodes = 15
    .type = int
    .expert_level = 3
    .help = "Number of Chebyshev nodes (and terms) for the smoothed "\
            "SigmaP(d*^2) trend of |F_model|^2 used to normalise Emodel."
  auto_kernel_number = 50
    .type = int
    .expert_level = 3
    .help = "Number of sorted reflections spanned by the auto-tuned "\
            "SigmaP kernel width (as kernel_normalisation's "\
            "number_of_sorted_reflections_for_auto_kernel)."
  n_sigmaa_coeffs = 8
    .type = int
    .short_caption = Number of sigmaA spline coefficients (E-scale)
    .help = "Number of B-spline coefficients for sigmaA(resolution) "\
            "(spline only)."
  spline_degree = 3
    .type = int
    .expert_level = 3
    .help = "Degree of the clamped B-spline basis (spline only)."
  sigmaa_max_iterations = 100
    .type = int
    .expert_level = 3
    .help = "Maximum L-BFGS iterations for the spline fit."
  sigmaa_curvature_weight = 0.02
    .type = float
    .short_caption = SigmaA spline curvature restraint weight (E-scale)
    .help = "Weight of the restraint on the second difference of the "\
            "spline's pre-sigmoid coefficients (spline only). 0 disables "\
            "it."
""", process_includes=True)

def b_spline_design_matrix(x, n_coeffs, degree, x_range=None):
  """ Clamped B-spline design matrix B[i,k] = B_k(x[i]), with interior
  knots evenly spaced over x_range (x mapped to [0,1]). Returns a numpy
  array of shape (len(x), n_coeffs). x outside x_range raises ValueError.

  x_range: (x_min, x_max), or None to use the range of x. When a curve is
  fitted on one set of x and evaluated on another, both calls must pass
  the same x_range.
  """
  import numpy as np
  from scipy.interpolate import BSpline
  x = np.asarray(x, dtype=float)
  if(x_range is None):
    x_min = float(x.min())
    x_max = float(x.max())
  else:
    x_min, x_max = x_range
  if(x_max <= x_min):
    # Degenerate resolution range (e.g. a single reflection, or all
    # reflections at identical resolution): fall back to a constant
    # basis (every reflection maps to the same single coefficient) so
    # this does not crash; the caller ends up fitting one overall value.
    x_max = x_min + 1.0
  x_norm = (x - x_min) / (x_max - x_min)
  n_knots = n_coeffs + degree + 1
  n_interior = n_knots - 2 * degree
  if(n_interior < 2):
    raise RuntimeError(
      "n_coeffs=%d is too small for spline_degree=%d (need at least "
      "%d coefficients)." % (n_coeffs, degree, 2 * degree - degree + 1))
  interior = np.linspace(0.0, 1.0, n_interior)
  knots = np.concatenate([[0.0] * degree, interior, [1.0] * degree])
  design = BSpline.design_matrix(
    x_norm, knots, degree, extrapolate=False).toarray()
  return design

def spline_curvature_penalty_and_gradient(coeffs, weight):
  """ Roughness restraint on a B-spline's pre-sigmoid coefficients c:
  R(c) = weight * sum_i (c[i-1] - 2*c[i] + c[i+1])^2, i.e. the squared
  second difference, a proxy for the curvature of z(x) = design(x).c
  with evenly spaced knots. It keeps the sigmaA curve from collapsing
  where the R-free set is too sparse to constrain it (the highest
  resolution shells). It acts on z rather than sigmaA: near the lower
  bound of the sigmoid, log(sigmaA - lower) ~ const + z, and the penalty
  stays quadratic in c. It is zero for c linear in the index.

  coeffs: numpy array of the coefficients (not sigmaA).
  weight: 0 (or fewer than 3 coefficients) gives (0.0, zeros).

  Returns (penalty, d(penalty)/d(coeffs)), as a float and a numpy array.
  """
  import numpy as np
  n = coeffs.shape[0]
  if(weight == 0 or n < 3):
    return 0.0, np.zeros_like(coeffs)
  d2 = coeffs[:-2] - 2.0 * coeffs[1:-1] + coeffs[2:]
  penalty = weight * float(np.sum(d2 * d2))
  grad = np.zeros_like(coeffs)
  grad[:-2]  += 2.0 * weight * d2
  grad[1:-1] += -4.0 * weight * d2
  grad[2:]   += 2.0 * weight * d2
  return penalty, grad

def _auto_kernel_width(d_star_sq, number=50):
  """ Kernel width as chosen by kernel_normalisation(auto_kernel=True):
  the d*^2 span of the lowest `number` reflections, widened while that
  is zero. Reimplemented here because the kernel_normalisation class
  divides by epsilon (see build_sigma_p). Falls back to the full d*^2
  range, or to 1 if d*^2 is constant, rather than failing.
  """
  n = d_star_sq.size()
  assert n > 1
  d_star_sq_low = flex.min(d_star_sq)
  d_star_sq_high = flex.max(d_star_sq)
  if(d_star_sq_high - d_star_sq_low < 1.e-12):
    return 1.0
  d_sorted = d_star_sq.select(flex.sort_permutation(d_star_sq))
  number = min(number, n - 1)
  kernel_width = d_sorted[number] - d_star_sq_low
  if(kernel_width == 0):
    original_number = number
    while number < n / 2:
      number = min(number + original_number, n - 1)
      kernel_width = d_sorted[number] - d_star_sq_low
      if(kernel_width > 0):
        break
  if(kernel_width <= 0):
    kernel_width = d_star_sq_high - d_star_sq_low
  return kernel_width

def build_sigma_p(f_model_no_aniso_scale, d_star_sq, n_nodes=15,
      auto_kernel_number=50, d_star_sq_eval=None):
  """ SigmaP(d*^2): smooth resolution trend of |f_model_no_aniso_scale|^2,
  without epsilon (build_e_model applies epsilon per reflection). Follows
  kernel_normalisation: Gaussian-kernel average at Chebyshev nodes in
  d*^2, Chebyshev fit to the log of those averages. Calls the
  kernel_normalisation extension directly with epsilon = 1, because the
  Python class divides by each reflection's epsilon before averaging.

  Returns a flex.double of SigmaP at each reflection's d*^2, or at
  d_star_sq_eval if given (e.g. reflections without data). d_star_sq_eval
  is clamped to the fitted range, since the polynomial is not constrained
  outside it.
  """
  d_star_sq_low = flex.min(d_star_sq)
  d_star_sq_high = flex.max(d_star_sq)
  intensity = flex.norm(f_model_no_aniso_scale)
  kernel_width = _auto_kernel_width(d_star_sq, number=auto_kernel_number)
  nodes = chebyshev_lsq_fit.chebyshev_nodes(
    n=n_nodes, low=d_star_sq_low, high=d_star_sq_high, include_limits=True)
  mean_intensity_at_nodes = scaling.kernel_normalisation(
    d_star_sq_hkl=d_star_sq,
    I_hkl=intensity,
    epsilon=flex.double(intensity.size(), 1.0),
    d_star_sq_array=nodes,
    kernel_width=kernel_width)
  floor = 1.e-16
  mean_intensity_at_nodes.set_selected(mean_intensity_at_nodes < floor, floor)
  fit = chebyshev_lsq_fit.chebyshev_lsq_fit(
    n_nodes, nodes, flex.log(mean_intensity_at_nodes))
  poly = chebyshev_polynome(
    n_nodes, d_star_sq_low, d_star_sq_high, fit.coefs)
  if(d_star_sq_eval is None):
    eval_at = d_star_sq
  else:
    eval_at = d_star_sq_eval.deep_copy()
    eval_at.set_selected(eval_at < d_star_sq_low, d_star_sq_low)
    eval_at.set_selected(eval_at > d_star_sq_high, d_star_sq_high)
  return flex.exp(poly.f(eval_at))

def build_e_model(f_model_no_aniso_scale, epsilons, d_star_sq,
      n_sigmap_nodes=15, auto_kernel_number=50, miller_set=None,
      model_anisotropy="off"):
  """ Emodel = f_model_no_aniso_scale / sqrt(EPS * SigmaP * A), with SigmaP
  from build_sigma_p and A = 1, or with model_anisotropy ("ellipsoid" or
  "resolution_dependent", which need miller_set, the reflections of
  f_model_no_aniso_scale) the anisotropic factor fitted with SigmaP to
  |F|^2/EPS (mmtbx.refinement.llgi_model_anisotropy.fit); SigmaP is then
  the trend of |F|^2/A.

  epsilons: per-reflection epsilon factors, as flex.double.

  Returns a group_args with .e_model (flex.complex_double), .sigma_p
  (flex.double), .anisotropy (A, flex.double), .factor
  (llgi_model_anisotropy.anisotropic_factor or None) and
  .sigma_p_at(d_star_sq_eval), SigmaP at other resolutions (as
  build_sigma_p's d_star_sq_eval).
  """
  llgi_model_anisotropy.check_mode(model_anisotropy)
  n = f_model_no_aniso_scale.size()
  assert epsilons.size() == n
  assert d_star_sq.size() == n
  def sigma_p_of(f, d_star_sq_eval=None):
    return build_sigma_p(f, d_star_sq, n_nodes=n_sigmap_nodes,
      auto_kernel_number=auto_kernel_number, d_star_sq_eval=d_star_sq_eval)
  if(model_anisotropy == "off" or miller_set is None):
    sigma_p = sigma_p_of(f_model_no_aniso_scale)
    a = flex.double(n, 1)
    factor = None
  else:
    assert miller_set.size() == n
    # A is fitted to |F|^2/EPS against the trend of |F|^2/A (SigmaP excludes
    # epsilon, as above); A is trace-free, so the 1/<EPS> this leaves in
    # shells with EPS > 1 reflections hardly affects it
    def iso_fn(i_over_a):  # SigmaP of |F|^2/A (without epsilon, as above)
      return sigma_p_of(flex.complex_double(flex.sqrt(i_over_a * epsilons),
        flex.double(n, 0)))
    r = llgi_model_anisotropy.fit(
      intensity=flex.norm(f_model_no_aniso_scale) / epsilons,
      indices=miller_set.indices(),
      centric_flags=miller_set.centric_flags().data(),
      unit_cell=miller_set.unit_cell(), space_group=miller_set.space_group(),
      iso_fn=iso_fn, mode=model_anisotropy)
    sigma_p, a, factor = r.sigma_iso, r.a, r.factor
  e_model = f_model_no_aniso_scale * (1.0 / flex.sqrt(epsilons * sigma_p * a))
  f_over_root_a = f_model_no_aniso_scale * (1.0 / flex.sqrt(a))
  return group_args(e_model=e_model, sigma_p=sigma_p, anisotropy=a,
    factor=factor,
    sigma_p_at=lambda d_star_sq_eval: sigma_p_of(f_over_root_a, d_star_sq_eval))

def build_e_model_as_fitted(llgi_data, f_model_no_aniso_scale, miller_set):
  """ build_e_model for the reflections of miller_set with the settings of
  the sigmaA fit that made llgi_data's sigmaA (llgi_data.fit_params.
  e_scale, or the defaults; the model anisotropy as llgi_model_anisotropy.
  mode_for), for the maps """
  fp = getattr(llgi_data, "fit_params", None)
  e_params = getattr(fp, "e_scale", None) if fp is not None else None
  if(e_params is None):
    e_params = llgi_e_sigmaa_params.extract()
  return build_e_model(f_model_no_aniso_scale,
    miller_set.epsilons().data().as_double(), miller_set.d_star_sq().data(),
    n_sigmap_nodes=e_params.n_sigmap_nodes,
    auto_kernel_number=e_params.auto_kernel_number, miller_set=miller_set,
    model_anisotropy=llgi_model_anisotropy.mode_for(llgi_data))

def isotropic_scatfrac_warnings(scatfrac_params, model_anisotropy):
  """ A warning (list of messages) if ScatFrac stays isotropic although the
  model's E is normalised anisotropically (the likelihood form) """
  if(model_anisotropy == "off" or scatfrac_params is None or
     scatfrac_params.method == "moments"):
    return []
  return ["ScatFrac (method=%s) is isotropic; the model anisotropy (%s) "
    "applies to the sigmaA fit and the maps only" % (scatfrac_params.method,
      model_anisotropy)]

def build_e_eff(feff, resn):
  """ Eeff = Feff / RESN. nacelle's RESN (root-epsilon-SigmaN) already
  includes epsilon and the expected-intensity trend, so no further
  normalisation is applied. Returns a flex.double.
  """
  assert resn.size() == feff.size()
  return feff / resn

def _sigmoid(z, lower=0.01, upper=0.99):
  """ Bounded sigmoid with sigmaa_point_estimator's bounds. Returns
  (value, d(value)/dz) as numpy arrays.
  """
  import numpy as np
  z = np.asarray(z, dtype=float)
  exp_neg_z = np.exp(-z)
  value = lower + (upper - lower) / (1.0 + exp_neg_z)
  dvalue_dz = (upper - lower) * exp_neg_z / (1.0 + exp_neg_z) ** 2
  return value, dvalue_dz

def bss_k_sol_b_sol(fmodel, k_sol_default=0.35, b_sol_default=46.0):
  """ Scalar (k_sol, b_sol) summary of bss's per-reflection k_mask, via
  fmodel.k_sol_b_sol_from_k_mask() (clipped to [0, 0.6]/[0, 150]); the
  defaults if bss has no estimate. Used for logging and to anchor
  D_model's B_defect; Emodel uses k_mask itself.
  """
  k_masks = fmodel.k_masks()
  assert len(k_masks) == 1, (
    "bss_k_sol_b_sol: only a single bulk-solvent mask is supported; "
    "got %d." % len(k_masks))
  k_sol, b_sol = fmodel.k_sol_b_sol_from_k_mask()
  if(k_sol is None or b_sol is None):
    return k_sol_default, b_sol_default
  return float(k_sol), float(b_sol)

class spline_sigmaa_curve(object):
  """ A fitted spline sigmaA as a function of d*^2: curve(d_star_sq)
  returns a flex.double, with d*^2 clamped to the fitted range, outside
  which the B-spline basis is not defined. Holds only the coefficients
  and basis, so it can be copied and pickled with llgi_data.
  """

  def __init__(self, coeffs, x_range, spline_degree):
    import numpy as np
    self.coeffs = np.array(coeffs, dtype=float)
    self.x_range = tuple(x_range)
    self.spline_degree = spline_degree

  def __call__(self, d_star_sq):
    import numpy as np
    d_star_sq_clamped = np.clip(
      np.asarray(d_star_sq, dtype=float), self.x_range[0], self.x_range[1])
    design = b_spline_design_matrix(
      d_star_sq_clamped, self.coeffs.size, self.spline_degree,
      x_range=self.x_range)
    sigmaa, _ = _sigmoid(design.dot(self.coeffs))
    return flex.double(sigmaa)

class e_sigmaa_target_evaluator(object):
  """ scitbx.lbfgs target evaluator for the spline sigmaA(d*^2): the
  B-spline coefficients of z(d*^2), sigmaA = _sigmoid(z), fitted against
  the mean E-scale LLGI over test_selection with Emodel fixed, plus the
  curvature restraint spline_curvature_penalty_and_gradient.
  """

  def __init__(self,
        e_eff, test_selection, e_model, dobs, centric_flags, d_star_sq,
        n_sigmaa_coeffs, spline_degree=3, max_iterations=100,
        curvature_weight=0.0, hybrid=None, tncs_u=None):
    self.hybrid = hybrid
    self.tncs_u = tncs_u
    self.e_eff = e_eff
    self.test_selection = test_selection
    self.e_model = e_model
    self.dobs = dobs
    self.centric_flags = centric_flags
    self.n_sigmaa_coeffs = n_sigmaa_coeffs
    self.spline_degree = spline_degree
    self.curvature_weight = curvature_weight
    self.x_range = (flex.min(d_star_sq), flex.max(d_star_sq))
    self.sigmaa_design = b_spline_design_matrix(
      d_star_sq.as_numpy_array(), n_sigmaa_coeffs, spline_degree,
      x_range=self.x_range)
    # z = 0 everywhere: sigmaA = 0.5
    self.x = flex.double(n_sigmaa_coeffs, 0.0)
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations=max_iterations)
    exception_handling_parameters = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound=True,
      ignore_line_search_failed_step_at_upper_bound=True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=term_parameters,
      exception_handling_params=exception_handling_parameters)

  def _current_sigmaa(self):
    z = self.sigmaa_design.dot(self.x.as_numpy_array())
    return _sigmoid(z)

  def compute_functional_and_gradients(self):
    sigmaa, dsigmaa_dz = self._current_sigmaa()
    result = llgi_tncs.e_sigmaa_functor(
      e_eff=self.e_eff,
      selection=self.test_selection,
      e_model=self.e_model,
      dobs=self.dobs,
      sigmaa=flex.double(sigmaa),
      centric_flags=self.centric_flags, hybrid=self.hybrid,
      tncs_u=self.tncs_u)
    f = result.target()
    d_target_by_dsigmaa = result.d_target_by_dsigmaa().as_numpy_array()
    g = self.sigmaa_design.T.dot(d_target_by_dsigmaa * dsigmaa_dz)
    penalty, penalty_grad = spline_curvature_penalty_and_gradient(
      self.x.as_numpy_array(), self.curvature_weight)
    return f + penalty, flex.double(g + penalty_grad)

  def sigmaa(self):
    sigmaa, _ = self._current_sigmaa()
    return flex.double(sigmaa)

  def curve(self):
    """ The current sigmaA curve, as a spline_sigmaa_curve. """
    return spline_sigmaa_curve(
      self.x.as_numpy_array(), self.x_range, self.spline_degree)

def estimate_e_sigmaa(e_eff, r_free_flags, e_model, dobs, centric_flags,
      d_star_sq, n_coeffs=8, spline_degree=3, max_iterations=100,
      curvature_weight=0.0, hybrid=None, tncs_u=None):
  """ Spline fit of sigmaA(d*^2) against the E-scale LLGI on the R-free
  set, with Emodel fixed.

  Returns a group_args with .sigmaa (flex.double, every reflection),
  .target (final mean target on the R-free set), .evaluate_at (the
  fitted curve, a spline_sigmaa_curve) and .lbfgs_error (None, or the
  message L-BFGS stopped with).
  """
  n_refl = e_eff.size()
  assert r_free_flags.size() == n_refl
  assert e_model.size() == n_refl
  assert dobs.size() == n_refl
  assert centric_flags.size() == n_refl
  assert d_star_sq.size() == n_refl
  if(r_free_flags.count(True) == 0):
    raise RuntimeError(
      "No R-free/test-set reflections available for the E-scale LLGI "
      "sigmaA fit.")
  evaluator = e_sigmaa_target_evaluator(
    e_eff=e_eff,
    test_selection=r_free_flags,
    e_model=e_model,
    dobs=dobs,
    centric_flags=centric_flags,
    d_star_sq=d_star_sq,
    n_sigmaa_coeffs=n_coeffs,
    spline_degree=spline_degree,
    max_iterations=max_iterations,
    curvature_weight=curvature_weight,
    hybrid=hybrid, tncs_u=tncs_u)
  sigmaa = evaluator.sigmaa()
  final_result = llgi_tncs.e_sigmaa_functor(
    e_eff=e_eff, selection=r_free_flags, e_model=e_model, dobs=dobs,
    sigmaa=sigmaa, centric_flags=centric_flags, hybrid=hybrid, tncs_u=tncs_u)
  return group_args(
    sigmaa=sigmaa, target=final_result.target(),
    evaluate_at=evaluator.curve(),
    lbfgs_error=evaluator.minimizer.error)

def estimate_e_sigmaa_for_fmodel(fmodel, dobs, feff, resn, params=None,
      hybrid=None, tncs_u=None):
  """ Fit sigmaA(resolution) on the E scale against fmodel's current model,
  with bulk solvent as bss left it (fmodel is not modified).

  fmodel: an mmtbx.f_model.manager, already scaled (e.g. after
  update_all_scales()).

  dobs, feff, resn: nacelle DOBS/FEFF/RESN on fmodel.f_obs()'s index set.

  params: extracted llgi_e_sigmaa_params, or None for defaults.

  hybrid: the llgi_hybrid to use instead of fmodel's (tNCS: the u-frame
  data, with dobs and feff in that frame), or None.

  tncs_u: the tNCS error-variance modulation (the fit is in the u-frame:
  feff/resn and the hybrid's observations normalised by u), or None.

  Returns a group_args with .sigmaa (flex.double, every reflection),
  .target (final mean target on the R-free set), .evaluate_at (the fitted
  curve at other d*^2 values), .k_sol/.b_sol (bss_k_sol_b_sol; b_sol
  anchors D_model's B_defect), .e_model (|E_model| as used by the fit) and
  .lbfgs_error (None, or why the fit's L-BFGS stopped early).
  """
  if(params is None):
    params = llgi_e_sigmaa_params.extract()
  f_obs = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags().data()
  centric_flags = f_obs.centric_flags().data()
  epsilons = f_obs.epsilons().data().as_double()
  d_star_sq = f_obs.d_star_sq().data()
  e_eff = build_e_eff(feff, resn)
  e_model = flex.abs(build_e_model(
    fmodel.f_model_no_aniso_scale().data(), epsilons, d_star_sq,
    n_sigmap_nodes=params.n_sigmap_nodes,
    auto_kernel_number=params.auto_kernel_number, miller_set=f_obs,
    model_anisotropy=llgi_model_anisotropy.resolve(params.model_anisotropy,
      fmodel.llgi_data())).e_model)
  k_sol, b_sol = bss_k_sol_b_sol(fmodel)
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  if(hybrid is None):
    hybrid = llgi_hybrid.get_hybrid(fmodel.llgi_data())
  if(hybrid is not None):
    # Exact likelihood wherever possible: a sigmaA-dependent switch would
    # make the fitted objective discontinuous (see llgi_exact.h, hybrid).
    hybrid = hybrid.with_rice_kappa(0.0)
  if(params.sigmaa_model == "d_model"):
    dp = params.d_model_params
    result = llgi_e_dmodel_fit.estimate_d_model_sigmaa(
      e_eff=e_eff, r_free_flags=r_free_flags, e_model=e_model, dobs=dobs,
      centric_flags=centric_flags, d_star_sq=d_star_sq,
      n_gaussian_terms=dp.n_gaussian_terms,
      max_iterations=dp.max_iterations,
      b_sol_anchor=b_sol,
      include_constant_term=dp.include_constant_term, hybrid=hybrid,
      tncs_u=tncs_u)
  else:
    result = estimate_e_sigmaa(
      e_eff=e_eff, r_free_flags=r_free_flags, e_model=e_model, dobs=dobs,
      centric_flags=centric_flags, d_star_sq=d_star_sq,
      n_coeffs=params.n_sigmaa_coeffs, spline_degree=params.spline_degree,
      max_iterations=params.sigmaa_max_iterations,
      curvature_weight=params.sigmaa_curvature_weight, hybrid=hybrid,
      tncs_u=tncs_u)
  return group_args(
    sigmaa=result.sigmaa,
    target=result.target,
    evaluate_at=result.evaluate_at,
    k_sol=k_sol, b_sol=b_sol,
    e_model=e_model,
    lbfgs_error=result.lbfgs_error)

def estimate_sigmaa_e_then_scatfrac_f(
      fmodel, dobs, feff, resn, e_params=None, scatfrac_params=None):
  """ Fit sigmaA on the E scale, then ScatFrac (the F-scale model's
  normalisation, which the E-scale target does not need) with sigmaA
  fixed.

  Step 1: estimate_e_sigmaa_for_fmodel (R-free set).
  Step 2: llgi_scatfrac.estimate_llgi_scatfrac (working set), against
  fmodel.f_model(), which includes k_isotropic as Feff does: the
  moments ScatFrac or the likelihood fit (scatfrac_params.method).

  fmodel: an mmtbx.f_model.manager, already scaled.

  dobs, feff, resn: nacelle DOBS/FEFF/RESN on fmodel.f_obs()'s index set.

  e_params, scatfrac_params: extracted llgi_e_sigmaa_params and
  llgi_scatfrac.llgi_scatfrac_params, or None for defaults.

  Returns a group_args with .sigmaa and .scatfrac (flex.double, every
  reflection), .sigmaa_curve (the fitted sigmaA as a function of d*^2,
  see estimate_e_sigmaa_for_fmodel's .evaluate_at), .target (final
  ScatFrac target on the working set),
  .scatfrac_inf/.b_scatfrac (ScatFrac = scatfrac_inf*exp(-b_scatfrac*ss),
  likelihood method; None for moments) and .warnings (list of
  messages about fits that stopped early, and isotropic_scatfrac_warnings).
  """
  import mmtbx.refinement.llgi_scatfrac as llgi_scatfrac
  if(e_params is None):
    e_params = llgi_e_sigmaa_params.extract()
  if(scatfrac_params is None):
    scatfrac_params = llgi_scatfrac.llgi_scatfrac_params.extract()
  sigmaa_result = estimate_e_sigmaa_for_fmodel(
    fmodel, dobs=dobs, feff=feff, resn=resn, params=e_params)
  sigmaa = sigmaa_result.sigmaa
  f_obs = fmodel.f_obs()
  llgi_data = fmodel.llgi_data()
  model_anisotropy = llgi_model_anisotropy.resolve(e_params.model_anisotropy,
    llgi_data)
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  scatfrac_result = llgi_scatfrac.estimate_llgi_scatfrac(
    f_eff=feff,
    working_selection=~fmodel.r_free_flags().data(),
    f_calc=fmodel.f_model().data(),
    dobs=dobs,
    sigmaa=sigmaa,
    teps=llgi_data.teps.data(),
    resn=resn,
    centric_flags=f_obs.centric_flags().data(),
    d_star_sq=f_obs.d_star_sq().data(),
    scale_factor=fmodel.scale_ml_wrapper(),
    params=scatfrac_params, hybrid=llgi_hybrid.get_hybrid(llgi_data),
    miller_set=f_obs, model_anisotropy=model_anisotropy)
  return group_args(
    sigmaa=sigmaa, sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=scatfrac_result.scatfrac,
    target=scatfrac_result.target,
    scatfrac_inf=scatfrac_result.scatfrac_inf,
    b_scatfrac=scatfrac_result.b_scatfrac,
    warnings=[
      "%s fit: L-BFGS stopped early: %s" % (name, error)
      for name, error in (("sigmaA", sigmaa_result.lbfgs_error),
                          ("ScatFrac", scatfrac_result.lbfgs_error))
      if error is not None] + isotropic_scatfrac_warnings(scatfrac_params,
        model_anisotropy))
