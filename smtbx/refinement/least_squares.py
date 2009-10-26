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

  def compute_quick_scale_factor_approximation(self):
    sel = self.fo_sq.data() > 0.99*flex.max(self.fo_sq.data())
    fo_sq = self.fo_sq.select(sel)
    fc_sq = xray.structure_factors.from_scatterers_direct(
      self.xray_structure, fo_sq).f_calc().norm()
    self.scale_factor = (
      flex.mean_weighted(fc_sq.data()*fo_sq.data(), fo_sq.sigmas())
     /flex.mean_weighted(flex.pow2(fc_sq.data()), fo_sq.sigmas()))

  def build_up(self):
    if self.scale_factor is None:
      self.compute_quick_scale_factor_approximation()
    if self.reduced is not None:
      self._core_normal_eqns.reset()
    ext.build_normal_equations(
      self._core_normal_eqns,
      self.fo_sq.indices(),
      self.fo_sq.data(),
      self.fo_sq.sigmas(),
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
    return self.reduced.solution

  def apply_shifts(self, shifts):
    self.special_position_constraints.apply_shifts(shifts)

  def solve_and_apply_shifts(self):
    shifts = self.solve()
    self.apply_shifts(shifts)
