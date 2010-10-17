import boost.python
ext = boost.python.import_ext("smtbx_refinement_least_squares_ext")
from smtbx_refinement_least_squares_ext import *

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

  def __init__(self, xray_structure, fo_sq, reparametrisation, **kwds):
    self.xray_structure = xray_structure
    self.fo_sq = fo_sq
    self.reparametrisation = reparametrisation
    adopt_optional_init_args(self, kwds)
    self.one_h_linearisation = direct.linearisation_of_f_calc_modulus_squared(
      self.xray_structure)
    if self.weighting_scheme == "default":
      self.weighting_scheme = self.default_weighting_scheme()
    self.floating_origin_restraints = floating_origin_restraints(
      xray_structure.space_group(),
      reparametrisation.asu_scatterer_parameters,
      reparametrisation.jacobian_transpose_matching_grad_fc(),
      self.floating_origin_restraint_relative_weight)
    self._core_normal_eqns = lstbx.normal_equations_separating_scale_factor(
      self.reparametrisation.n_independent_params)
    self.reduced = None
    self.shifts = None

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
    if self.restraints_manager is not None:
      # Here we determine a normalisation factor to place the restraints on the
      # same scale as the observations. This is the normalisation factor
      # suggested in Giacovazzo. In contrast, shelxl simply uses the mean
      # value of the deltas (shelx manual, page 5-1).
      normalisation_factor = flex.sum(self.weights * flex.pow2(
        self.fo_sq.data() - self.scale_factor * flex.norm(self.f_calc.data())))\
          / (self.fo_sq.size() - self.reparametrisation.n_independent_params)
      linearised_eqns = self.restraints_manager.build_linearised_eqns(
        self.xray_structure)
      jacobian = \
        self.reparametrisation.jacobian_transpose_matching_grad_fc().transpose()
      self.reduced.add_equations(linearised_eqns.deltas,
                                 linearised_eqns.design_matrix * jacobian,
                                 linearised_eqns.weights * normalisation_factor,
                                 negate_right_hand_side=True)
      self.n_restraints = linearised_eqns.n_restraints()
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
    return math.sqrt(
      self.objective_data_only
      /(self.fo_sq.size() - self.reparametrisation.n_independent_params))

  def restrained_goof(self):
    if self.n_restraints is not None: n_restraints = self.n_restraints
    else: n_restraints = 0
    return math.sqrt(
      self.objective
      /(self.fo_sq.size() + n_restraints
        - self.reparametrisation.n_independent_params))

  def wR2(self):
    return math.sqrt(flex.sum(self.weights * flex.pow2(
      self.fo_sq.data() - self.scale_factor * flex.norm(self.f_calc.data())))
                     /flex.sum(self.weights * flex.pow2(self.fo_sq.data())))

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

  def covariance_matrix(self, normalised_by_goof=True):
    cov_ind_params = linalg.inverse_of_u_transpose_u(
      self.reduced.cholesky_factor_packed_u)
    jac_tr = self.reparametrisation.jacobian_transpose_matching_grad_fc()
    cov = jac_tr.self_transpose_times_symmetric_times_self(cov_ind_params)
    if normalised_by_goof: cov *= self.restrained_goof()**2
    return cov

  def covariance_matrix_and_annotations(self):
    return (self.covariance_matrix(),
            self.reparametrisation.component_annotations)
