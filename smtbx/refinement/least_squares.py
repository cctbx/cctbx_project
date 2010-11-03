import boost.python
ext = boost.python.import_ext("smtbx_refinement_least_squares_ext")
from smtbx_refinement_least_squares_ext import *

import libtbx
from libtbx import adopt_optional_init_args
from scitbx import linalg, lstbx
from scitbx.array_family import flex
from cctbx import xray
from smtbx.structure_factors import direct

from stdlib import math

class normal_equations(object):

  default_weighting_scheme = mainstream_shelx_weighting
  weighting_scheme = "default"
  floating_origin_restraint_relative_weight = 1e3
  scale_factor = None
  f_mask = None
  restraints_manager=None
  n_restraints = None

  def __init__(self, fo_sq, reparametrisation, **kwds):
    self.fo_sq = fo_sq
    self.reparametrisation = reparametrisation
    adopt_optional_init_args(self, kwds)
    self.one_h_linearisation = direct.linearisation_of_f_calc_modulus_squared(
      self.xray_structure)
    if self.weighting_scheme == "default":
      self.weighting_scheme = self.default_weighting_scheme()
    self.floating_origin_restraints = floating_origin_restraints(
      self.xray_structure.space_group(),
      reparametrisation.asu_scatterer_parameters,
      reparametrisation.jacobian_transpose_matching_grad_fc(),
      self.floating_origin_restraint_relative_weight)
    self._core_normal_eqns = lstbx.normal_equations_separating_scale_factor(
      self.reparametrisation.n_independent_params)
    self.reduced = None
    self.shifts = None

  class xray_structure(libtbx.property):
    def fget(self):
      return self.reparametrisation.structure

  def compute_quick_scale_factor_approximation(self):
    self.fo_sq.set_observation_type_xray_intensity()
    f_calc = xray.structure_factors.from_scatterers_direct(
      self.xray_structure, self.fo_sq).f_calc()
    self.scale_factor = self.fo_sq.scale_factor(f_calc, cutoff_factor=0.99)

  def build_up(self):
    if self.scale_factor is None:
      self.compute_quick_scale_factor_approximation()
    if self.reduced is not None:
      self._core_normal_eqns.reset()
    if self.f_mask is not None:
      f_mask = self.f_mask.data()
    else:
      f_mask = flex.complex_double()
    if self.reduced is None:
      self.reparametrisation.linearise()
      self.reparametrisation.store()
    result = ext.build_normal_equations(
      self._core_normal_eqns,
      self.fo_sq.indices(),
      self.fo_sq.data(),
      self.fo_sq.sigmas(),
      f_mask,
      self.weighting_scheme,
      self.scale_factor,
      self.one_h_linearisation,
      self.reparametrisation.jacobian_transpose_matching_grad_fc())
    self.f_calc = self.fo_sq.array(data=result.f_calc(), sigmas=None)
    self.weights = result.weights()
    self.reduced = self._core_normal_eqns.reduced_equations()
    self.scale_factor = self._core_normal_eqns.optimal_scale_factor()
    self.objective_data_only = self._core_normal_eqns.objective()
    self.sum_w_yo_sq = self._core_normal_eqns.sum_w_yo_sq()
    if self.restraints_manager is not None:
      # Here we determine a normalisation factor to place the restraints on the
      # same scale as the observations. This is the normalisation factor
      # suggested in Giacovazzo. In contrast, shelxl simply uses the mean
      # value of the deltas (shelx manual, page 5-1).
      dof = self.fo_sq.size() - self.reparametrisation.n_independent_params
      normalisation_factor = self.objective_data_only/dof
      linearised_eqns = self.restraints_manager.build_linearised_eqns(
        self.xray_structure)
      jacobian = \
        self.reparametrisation.jacobian_transpose_matching_grad_fc().transpose()
      self.reduced.add_equations(linearised_eqns.deltas,
                                 linearised_eqns.design_matrix * jacobian,
                                 linearised_eqns.weights * normalisation_factor,
                                 negate_right_hand_side=True)
      self.n_restraints = linearised_eqns.n_restraints()
    self.parameter_vector_norm = \
        self.reparametrisation.norm_of_independent_parameter_vector
    self.objective = self._core_normal_eqns.objective()
    self.gradient = self._core_normal_eqns.gradient()
    self.floating_origin_restraints.add_to(self.reduced)

  def solve(self):
    self.reduced.solve()
    self.shifts = self.reduced.solution

  def apply_shifts(self):
    assert self.shifts is not None
    self.reparametrisation.apply_shifts(self.shifts)
    self.reparametrisation.linearise()
    self.reparametrisation.store()

  def solve_and_apply_shifts(self):
    self.solve()
    self.apply_shifts()

  def goof(self):
    dof = self.fo_sq.size() - self.reparametrisation.n_independent_params
    return math.sqrt(self.objective_data_only*self.sum_w_yo_sq/dof)

  def restrained_goof(self):
    if self.n_restraints is not None: n_restraints = self.n_restraints
    else: n_restraints = 0
    dof = (self.fo_sq.size() + n_restraints
           - self.reparametrisation.n_independent_params)
    return math.sqrt(self.objective*self.sum_w_yo_sq/dof)

  def wR2(self):
    return math.sqrt(self.objective_data_only)

  def r1_factor(self, cutoff_factor=None):
    f_obs = self.fo_sq.f_sq_as_f()
    if cutoff_factor is not None:
      strong = f_obs.data() > cutoff_factor*f_obs.sigmas()
      f_obs = f_obs.select(strong)
      f_calc = self.f_calc.select(strong)
    else:
      f_calc = self.f_calc
    R1 = f_obs.r1_factor(f_calc, scale_factor=math.sqrt(self.scale_factor))
    return R1, f_obs.size()

  def covariance_matrix(self,
                        independent_params=False,
                        normalised_by_goof=True):
    cov = linalg.inverse_of_u_transpose_u(
      self.reduced.cholesky_factor_packed_u)
    cov /= self.sum_w_yo_sq
    if not independent_params:
      jac_tr = self.reparametrisation.jacobian_transpose_matching_grad_fc()
      cov = jac_tr.self_transpose_times_symmetric_times_self(cov)
    if normalised_by_goof: cov *= self.restrained_goof()**2
    return cov

  def covariance_matrix_and_annotations(self):
    return (self.covariance_matrix(),
            self.reparametrisation.component_annotations)


class _mainstream_shelx_weighting(boost.python.injector,
                                  mainstream_shelx_weighting):

  def __str__(self):
    """ A string representation of the weighting scheme in a format that is
        appropriate for the CIF item _refine_ls_weighting_details.
    """
    if round(self.a, 4) in (0.1, 0.2):
      a = "%.1f" %self.a
    else:
      a = "%.4f" %self.a
    if round(self.b, 4) == 0: b_part=""
    else: b_part = "+%.4fP" %self.b
    return ("w=1/[\s^2^(Fo^2^)+(%sP)^2^%s]"
            " where P=(Fo^2^+2Fc^2^)/3" %(a, b_part))

  def type(self):
    return "calc"

  def optimise_parameters(self, normal_eqns):
    """ Find optimal values of a and b that give a flat analysis of the variance
        when binned by Fc/max(Fc), and a goodness of fit close to 1.

        This is done in a grid search fashion similar to Shelxl.

        self is not modified in place; instead a new instance of the weighting
        scheme is returned.
    """
    weighting = mainstream_shelx_weighting(a=self.a, b=self.b)

    def compute_chi_sq(fo_sq, fc_sq, a,b):
      weighting.a = a
      weighting.b = b
      weights = weighting(
        fo_sq.data(), fo_sq.sigmas(), fc_sq.data(), scale_factor)
      return (flex.sum(
        weights * flex.pow2(fo_sq.data() - scale_factor * fc_sq.data())))

    scale_factor = normal_eqns.scale_factor
    fo_sq = normal_eqns.fo_sq.deep_copy()
    fo_sq.data().set_selected(fo_sq.data() < 0, 0)
    fc_sq = normal_eqns.f_calc.as_intensity_array()

    fo2 = fo_sq.data().deep_copy()
    fo2 /= scale_factor
    sigmas = fo_sq.sigmas() / scale_factor
    sigmas_sq = flex.pow2(sigmas)
    fc2 = fc_sq.data()

    # determine starting values for a and b, formulae taken from shelxl code
    p = (fo2 + 2 * fc2)/3
    p_sq = flex.pow2(p)
    x = flex.sum((flex.pow2(fo2-fc2)-sigmas) * (p_sq/sigmas_sq))
    y = flex.sum( flex.pow2(p_sq)/sigmas_sq)
    z = flex.sum(p)
    start_a = math.sqrt(max(0.0001, 0.64*x/max(1e-8, y)))
    start_b = 0.5 * z * start_a**2 /fo_sq.size()
    a_step = 0.2 * start_a
    b_step = 0.4 * start_b

    # sort data and setup binning by fc/fc_max
    fc_sq_over_fc_sq_max = fc_sq.data()/flex.max(fc_sq.data())
    permutation = flex.sort_permutation(fc_sq_over_fc_sq_max)
    fc_sq_over_fc_sq_max = fc_sq.customized_copy(
      data=fc_sq_over_fc_sq_max).select(permutation)
    fc_sq = fc_sq.select(permutation)
    fo_sq = fo_sq.select(permutation)
    n_bins = 10
    bins = []
    bin_max = 0
    bin_limits = flex.size_t(1, 0)
    bin_count = flex.size_t()
    for i in range(n_bins):
      bin_limits.append(int(math.ceil((i+1) * fc_sq.size()/n_bins)))
      bin_count.append(bin_limits[i+1] - bin_limits[i])

    n = fo_sq.size() \
      /(fo_sq.size()-normal_eqns.reparametrisation.n_independent_params)

    # search on a 9x9 grid to determine best values of a and b
    gridding = flex.grid(9,9)
    while (a_step > 1e-4 and b_step > 5e-3):
      tmp = flex.double(gridding, 0)
      binned_chi_sq = [tmp.deep_copy() for i in range(n_bins)]
      start_a = max(start_a, 4*a_step) - 4*a_step
      start_b = max(start_b, 4*b_step) - 4*b_step
      for i_bin in range(n_bins):
        sel = flex.size_t_range(bin_limits[i_bin], bin_limits[i_bin+1])
        fc2 = fc_sq.select(sel)
        fo2 = fo_sq.select(sel)
        b = start_b
        for j in range(9):
          a = start_a
          b += b_step
          for k in range(9):
            a += a_step
            binned_chi_sq[i_bin][j,k] += compute_chi_sq(fo2, fc2, a, b)
      min_variance = 9e9
      j_min, k_min = (0, 0)
      for j in range(9):
        for k in range(9):
          variance = 0
          for i_bin in range(n_bins):
            if bin_count[i_bin] == 0: continue
            goof = math.sqrt(binned_chi_sq[i_bin][j,k]*n/bin_count[i_bin])
            variance += (goof-1)**2
          min_variance = min(variance, min_variance)
          if variance == min_variance:
            j_min = j
            k_min = k
      start_a += k_min*a_step
      start_b += j_min*b_step
      if k_min == 8:
        a_step *= 2
        continue
      elif k_min != 0:
        a_step /= 4
      if j_min == 8:
        b_step *= 2
        continue
      elif j_min != 0:
        b_step /=4
      if start_a <= 1e-4: a_step /= 4
      if start_b <= 1e-3: b_step /= 4
    if start_a > 0.2:
      start_a = 0.2
      start_b = 0
    weighting.a = start_a
    weighting.b = start_b
    return weighting

class _unit_weighting(boost.python.injector,
                      unit_weighting):

  def __str__(self):
    return "w=1"

  def type(self):
    return "unit"

  def optimise_parameters(self, normal_eqns):
    # no parameters to optimise!
    return self
