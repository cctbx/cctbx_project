import boost.python
ext = boost.python.import_ext("smtbx_refinement_least_squares_ext")
from smtbx_refinement_least_squares_ext import *

from libtbx import adopt_optional_init_args
from scitbx import linalg, lstbx
from scitbx.array_family import flex
from cctbx import xray
from smtbx.structure_factors import direct

class normal_equations(object):

  default_weighting_scheme = mainstream_shelx_weighting
  weighting_scheme = "default"
  floating_origin_restraint_relative_weight = 1e3
  scale_factor = None
  f_mask = None

  def __init__(self, xray_structure, fo_sq, **kwds):
    self.xray_structure = xray_structure
    self.fo_sq = fo_sq
    adopt_optional_init_args(self, kwds)
    self.one_h_linearisation = direct.linearisation_of_f_calc_modulus_squared(
      self.xray_structure)
    if self.weighting_scheme == "default":
      self.weighting_scheme = self.default_weighting_scheme()
    self.floating_origin_restraints = floating_origin_restraints(
      xray_structure.space_group(),
      xray_structure.site_symmetry_table(),
      xray_structure.scatterers(),
      self.floating_origin_restraint_relative_weight)
    self.special_position_constraints = special_position_constraints(
      xray_structure.unit_cell(),
      xray_structure.site_symmetry_table(),
      xray_structure.scatterers())
    self._core_normal_eqns = lstbx.normal_equations_separating_scale_factor(
      self.special_position_constraints.n_independent_params)
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
    ext.build_normal_equations(
      self._core_normal_eqns,
      self.fo_sq.indices(),
      self.fo_sq.data(),
      self.fo_sq.sigmas(),
      f_mask,
      self.weighting_scheme,
      self.scale_factor,
      self.one_h_linearisation,
      self.special_position_constraints)
    self.reduced = self._core_normal_eqns.reduced_equations()
    self.scale_factor = self._core_normal_eqns.optimal_scale_factor()
    self.objective = self._core_normal_eqns.objective()
    self.gradient = self._core_normal_eqns.gradient()
    self.floating_origin_restraints.add_to(self.reduced)

  def solve(self):
    self.reduced.solve()
    self.shifts = self.reduced.solution

  def apply_shifts(self):
    assert self.shifts is not None
    self.special_position_constraints.apply_shifts(self.shifts)

  def solve_and_apply_shifts(self):
    self.solve()
    self.apply_shifts()

  def covariance_matrix(self):
    """ Covariance matrix for crystallographic parameters
    They are ordered scatterer by scatterer, in the order
    they are stored in self.xray_structure, as follow:
       x, y, z, u_iso, u_11, u_22, u_33, u_12, u_13, u_23, occupancy
    If a parameter is not refined, it is taken out from the list.
    The upper diagonal of the covariance matrix is returned packed by rows
    (packed-u format).
    """
    from scitbx import sparse

    # compute jacobian matrix (sparse)
    jac = sparse.matrix(
      self.special_position_constraints.n_crystallographic_params,
      self.special_position_constraints.n_independent_params)
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    i = 0 # crystallographic param index
    j = 0 # independent param index
    for i_sc, sc in enumerate(self.xray_structure.scatterers()):
      site_symm = site_symmetry_table.get(i_sc)
      if sc.flags.grad_site():
        if site_symm.is_point_group_1():
          jac[i, j] = jac[i+1, j+1] = jac[i+2, j+2] = 1
          i += 3; j += 3
        else:
          site_constraints = site_symm.site_constraints()
          site_jac_tr = site_constraints.gradient_sum_matrix()
          for k in xrange(3):
            for l in xrange(site_constraints.n_independent_params()):
              jac[i+k, j+l] = site_jac_tr[l, k]
          i += 3; j += site_constraints.n_independent_params()
      if sc.flags.use_u_iso() and sc.flags.grad_u_iso():
        jac[i,j] = 1
        i += 1; j += 1
      if sc.flags.use_u_aniso() and sc.flags.grad_u_aniso():
        if site_symm.is_point_group_1():
          for l in xrange(6): jac[i+l, j+l] = 1
          i += 6; j += 6
        else:
          adp_constraints = site_symm.cartesian_adp_constraints(
            self.xray_structure.unit_cell())
          adp_jac = adp_constraints.jacobian()
          for k in xrange(6):
            for l in xrange(adp_constraints.n_independent_params()):
              jac[i+k, j+l] = adp_jac[k, l]
          i += 6; j += adp_constraints.n_independent_params()
      if sc.flags.grad_occupancy():
        jac[i,j] = 1
        i += 1; j += 1

    # compute covariance matrix for crystallographic parameters
    cov_ind_params = linalg.inverse_of_u_transpose_u(
      self.reduced.cholesky_factor_packed_u)
    cov = jac.self_times_symmetric_times_self_transpose(cov_ind_params)
    return cov
