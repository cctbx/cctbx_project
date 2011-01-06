import boost.python
ext = boost.python.import_ext("smtbx_refinement_least_squares_ext")
from smtbx_refinement_least_squares_ext import *

import smtbx.refinement.weighting_schemes # import dependency

import libtbx
from libtbx import adopt_optional_init_args
from scitbx import linalg
from scitbx.lstbx import normal_eqns
from scitbx.array_family import flex
from smtbx.structure_factors import direct

from stdlib import math

class crystallographic_ls(
  normal_eqns.non_linear_ls_with_separable_scale_factor):

  default_weighting_scheme = mainstream_shelx_weighting
  weighting_scheme = "default"
  floating_origin_restraint_relative_weight = 1e3
  f_mask = None
  restraints_manager=None
  n_restraints = None
  initial_scale_factor = None

  def __init__(self, fo_sq, reparametrisation, **kwds):
    super(crystallographic_ls, self).__init__(
      reparametrisation.n_independent_params)
    self.fo_sq = fo_sq
    self.reparametrisation = reparametrisation
    adopt_optional_init_args(self, kwds)
    if self.f_mask is not None:
      assert self.f_mask.size() == self.fo_sq.size()
    self.one_h_linearisation = direct.f_calc_modulus_squared(
      self.xray_structure)
    if self.weighting_scheme == "default":
      self.weighting_scheme = self.default_weighting_scheme()
    self.floating_origin_restraints = floating_origin_restraints(
      self.xray_structure.space_group(),
      reparametrisation.asu_scatterer_parameters,
      reparametrisation.jacobian_transpose_matching_grad_fc(),
      self.floating_origin_restraint_relative_weight)
    self.taken_step = None

  class xray_structure(libtbx.property):
    def fget(self):
      return self.reparametrisation.structure

  class twin_components(libtbx.property):
    def fget(self):
      return self.reparametrisation.twin_components

  def build_up(self, objective_only=False):
    if self.f_mask is not None:
      f_mask = self.f_mask.data()
    else:
      f_mask = flex.complex_double()
    twin_components = self.reparametrisation.twin_components
    if twin_components is None:
      twin_components = ()

    def args(scale_factor, weighting_scheme, objective_only):
      args = (self,
              self.fo_sq.indices(),
              self.fo_sq.data(),
              self.fo_sq.sigmas(),
              f_mask,
              weighting_scheme,
              scale_factor,
              self.one_h_linearisation,
              self.reparametrisation.jacobian_transpose_matching_grad_fc(),
              twin_components)
      if objective_only:
        args += (True,)
      return args

    if not self.finalised: #i.e. never been called
      self.reparametrisation.linearise()
      self.reparametrisation.store()
      scale_factor = self.initial_scale_factor
      if scale_factor is None: # we haven't got one from previous refinement
        result = ext.build_normal_equations(
          *args(scale_factor=None, weighting_scheme=sigma_weighting(),
                objective_only=True))
        scale_factor = self.scale_factor()
    else: # use scale factor from the previous step
      scale_factor = self.scale_factor()

    self.reset()
    result = ext.build_normal_equations(*args(scale_factor,
                                              self.weighting_scheme,
                                              objective_only))
    self.f_calc = self.fo_sq.array(data=result.f_calc(), sigmas=None)
    self.fc_sq = self.fo_sq.array(data=result.observables(), sigmas=None)\
        .set_observation_type_xray_intensity()
    self.weights = result.weights()
    self.objective_data_only = self.objective()
    self.chi_sq_data_only = self.chi_sq()
    if self.restraints_manager is not None:
      # Here we determine a normalisation factor to place the restraints on the
      # same scale as the observations. This is the normalisation factor
      # suggested in Giacovazzo. In contrast, shelxl simply uses the mean
      # value of the deltas (shelx manual, page 5-1).
      # The factor 2 comes from the fact that we minimize 1/2 sum w delta^2
      normalisation_factor = self.chi_sq_data_only/2
      linearised_eqns = self.restraints_manager.build_linearised_eqns(
        self.xray_structure)
      jacobian = \
        self.reparametrisation.jacobian_transpose_matching_grad_fc().transpose()
      self.reduced_problem().add_equations(
        linearised_eqns.deltas,
        linearised_eqns.design_matrix * jacobian,
        linearised_eqns.weights * normalisation_factor)
      self.n_restraints = linearised_eqns.n_restraints()
      self.chi_sq_data_and_restraints = self.chi_sq()
    if not objective_only:
      self.floating_origin_restraints.add_to(self.step_equations())

  def parameter_vector_norm(self):
    return self.reparametrisation.norm_of_independent_parameter_vector

  def scale_factor(self): return self.optimal_scale_factor()

  def step_forward(self):
    self.reparametrisation.apply_shifts(self.step())
    self.reparametrisation.linearise()
    self.reparametrisation.store()
    self.taken_step = self.step().deep_copy()

  def step_backward(self):
    self.reparametrisation.apply_shifts(-self.taken_step)
    self.reparametrisation.linearise()
    self.reparametrisation.store()
    self.taken_step = None

  def goof(self):
    return math.sqrt(self.chi_sq_data_only)

  def restrained_goof(self):
    return math.sqrt(self.chi_sq_data_and_restraints)

  def wR2(self):
    return math.sqrt(self.objective_data_only)

  def r1_factor(self, cutoff_factor=None):
    f_obs = self.fo_sq.f_sq_as_f()
    f_calc = self.fc_sq.f_sq_as_f()
    if cutoff_factor is not None:
      strong = f_obs.data() > cutoff_factor*f_obs.sigmas()
      f_obs = f_obs.select(strong)
      f_calc = f_calc.select(strong)
    else:
      f_calc = f_calc
    R1 = f_obs.r1_factor(f_calc, scale_factor=math.sqrt(self.scale_factor()))
    return R1, f_obs.size()

  def covariance_matrix(self,
                        jacobian_transpose=None,
                        normalised_by_goof=True):
    """ The columns of the jacobian_transpose determine which crystallographic
        parameters appear in the covariance matrix.
        If jacobian_transpose is None, then the covariance matrix returned will
        be that for the independent L.S. parameters.
    """
    if not self.step_equations().solved:
      self.solve()
    cov = linalg.inverse_of_u_transpose_u(
      self.step_equations().cholesky_factor_packed_u())
    cov /= self.sum_w_yo_sq()
    if jacobian_transpose is not None:
      cov = jacobian_transpose.self_transpose_times_symmetric_times_self(cov)
    if normalised_by_goof: cov *= self.restrained_goof()**2
    return cov

  def covariance_matrix_and_annotations(self):
    jac_tr = self.reparametrisation.jacobian_transpose_matching_grad_fc()
    return covariance_matrix_and_annotations(
      self.covariance_matrix(jacobian_transpose=jac_tr),
      self.reparametrisation.component_annotations)


class covariance_matrix_and_annotations(object):

  def __init__(self, covariance_matrix, annotations):
    """ The covariance matrix is assumed to be a symmetric matrix stored as a
        packed upper diagonal matrix.
    """
    self.matrix = covariance_matrix
    self.annotations = annotations
    self._2_n_minus_1 = 2*len(self.annotations)-1 # precompute for efficiency

  def __call__(self, i, j):
    return self.matrix[i*0.5*(self._2_n_minus_1-i)+j]

  def variance_of(self, annotation):
    i = self.annotations.index(annotation)
    return self(i, i)

  def covariance_of(self, annotation_1, annotation_2):
    i = self.annotations.index(annotation_1)
    j = self.annotations.index(annotation_2)
    if j > i:
      i, j = j, i
    return self(i, j)

  def diagonal(self):
    return self.matrix.matrix_packed_u_diagonal()
