from __future__ import absolute_import, division, print_function
import math
import os

from cctbx.array_family import flex
from cctbx import xray
from libtbx import adopt_init_args
import libtbx.phil
import libtbx.utils
import scitbx
import scitbx.lbfgs
import scitbx.lstbx
import scitbx.lstbx.normal_eqns_solving
from smtbx.refinement import least_squares
from smtbx.refinement import constraints
import smtbx.utils
from smtbx.refinement.constraints.tests import tst_constrained_structure
from six.moves import range


class lbfgs(object):

  def __init__(self, target_functor,
                     xray_structure,
                     reparametrisation,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None,
                     cos_sin_table=True,
                     structure_factor_algorithm=None,
                     verbose=0,
                     reference_structure=None):
    adopt_init_args(self, locals())
    self.scatterer_grad_flags_counts = xray.minimization.ext.scatterer_grad_flags_counts(
                                              self.xray_structure.scatterers())
    self.grad_flags_counts = \
      xray.minimization.ext.scatterer_grad_flags_counts(self.xray_structure.scatterers())
    self.structure_factors_from_scatterers = \
      xray.structure_factors.from_scatterers(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.structure_factor_gradients = \
      xray.structure_factors.gradients(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.x = flex.double(reparametrisation.n_independents, 0)
    xray_structure.tidy_us(u_min=1.e-6)
    self.last_shift = None
    import sys
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params,
      #log=sys.stdout
    )
    self.apply_shifts()
    self.compute_target(compute_gradients=False)
    self.final_target_value = self.target_result.target()

  def apply_shifts(self):
    shifts = self.x.deep_copy()
    if self.last_shift is not None:
      # beware: self.x are the shifts from the starting parameters, not from the
      # current parameters!
      shifts -= self.last_shift
      self.reparametrisation.apply_shifts(shifts)
    self.reparametrisation.linearise()
    self.reparametrisation.store()
    self.last_shift = self.x.deep_copy()
    return flex.double() # XXX

  def compute_target(self, compute_gradients):
    self.f_calc = self.structure_factors_from_scatterers(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      algorithm=self.structure_factor_algorithm).f_calc()
    self.target_result = self.target_functor(
      self.f_calc,
      compute_gradients)
    assert self.target_result.target() is not None

  def compute_functional_and_gradients(self):
    u_iso_refinable_params = self.apply_shifts()
    self.compute_target(compute_gradients=True)
    self.f = self.target_result.target()
    gradients = self.structure_factor_gradients(
      xray_structure=self.xray_structure,
      u_iso_refinable_params=u_iso_refinable_params,
      miller_set=self.target_functor.f_obs(),
      d_target_d_f_calc=self.target_result.derivatives(),
      n_parameters=0, # so the gradients aren't packed
      algorithm=self.structure_factor_algorithm)
    if self.scatterer_grad_flags_counts.site:
      d_target_d_site_frac = gradients.d_target_d_site_frac()
    if self.scatterer_grad_flags_counts.u_iso:
      d_target_d_u_iso = gradients.d_target_d_u_iso()
    if self.scatterer_grad_flags_counts.u_aniso:
      d_target_d_u_star = gradients.d_target_d_u_star()
    if self.scatterer_grad_flags_counts.fp:
      d_target_d_fp = gradients.d_target_d_fp()
    if self.scatterer_grad_flags_counts.fdp:
      d_target_d_fdp = gradients.d_target_d_fdp()
    if self.scatterer_grad_flags_counts.occupancy:
      d_target_d_occupancy = gradients.d_target_d_occupancy()
    # pack the gradients ourselves - currently the constraints system assumes
    # we are refining fractional coordinates and u_star
    self.g = flex.double()
    for i, sc in enumerate(self.xray_structure.scatterers()):
      if sc.flags.grad_site():
        for j in range(3):
          self.g.append(d_target_d_site_frac[i][j])
      if sc.flags.use_u_iso() and sc.flags.grad_u_iso():
        self.g.append(d_target_d_u_iso[i])
      if sc.flags.use_u_aniso() and sc.flags.grad_u_aniso():
        for j in range(6):
          self.g.append(d_target_d_u_star[i][j])
      if sc.flags.grad_occupancy():
        self.g.append(d_target_d_occupancy[i])
      if sc.flags.grad_fp():
        self.g.append(d_target_d_fp[i])
      if sc.flags.grad_fdp():
        self.g.append(d_target_d_fdp[i])
    if self.verbose > 0:
      print("xray.minimization line search: f,rms(g):", end=' ')
      print(self.f, math.sqrt(flex.mean_sq(self.g)))
    jacobian = self.reparametrisation.jacobian_transpose_matching_grad_fc()
    self.g = jacobian * self.g
    return self.f, self.g

  def callback_after_step(self, minimizer):
    if self.verbose > 0:
      print("xray.minimization step: f,iter,nfun:", end=' ')
      print(self.f,minimizer.iter(),minimizer.nfun())
    if self.verbose > 1 and self.reference_structure is not None:
      xray.meaningful_site_cart_differences(self.xray_structure, self.reference_structure).show()


def run(args):
  master_phil = libtbx.phil.parse("""
    d_min = 0.5
      .type = float
    constrained_refinement = True
      .type = bool
    random_seed = 1
      .type = int
    shake_sites_rmsd = 0.5
      .type = float
    shake_adp_spread = 20
      .type = float
    grad_site=True
      .type = bool
    grad_u_iso=False
      .type = bool
    grad_u_aniso=False
      .type = bool
    grad_occupancy=False
      .type = bool
    grad_fp_fdp=False
      .type = bool
    lbfgs_m = 5
      .type = int
    lbfgs_max_iterations = 1000
      .type = int
    verbose = 0
      .type = int
""")

  argument_interpreter = master_phil.command_line_argument_interpreter()
  phil_objects = []
  remaining_args = []
  for arg in args:
    if (arg.find("=") >= 0):
      phil_objects.append(argument_interpreter.process(arg=arg))
    else:
      remaining_args.append(arg)
  work_phil = master_phil.fetch(sources=phil_objects)
  work_phil.show()
  params = work_phil.extract()

  if params.random_seed is not None:
    import scitbx.random
    import random
    scitbx.random.set_random_seed(params.random_seed)
    flex.set_random_seed(params.random_seed)
    random.seed(params.random_seed)

  if len(remaining_args):
    assert len(remaining_args) == 1
    file_path = remaining_args[0]
    root, ext = os.path.splitext(file_path)

    if ext == ".cif":
      xs_dict = xray.structure.from_cif(file_path=file_path)
      assert len(xs_dict) == 1, "CIF should contain only one xray structure"
      xs = list(xs_dict.values())[0]
      xs.show_summary().show_scatterers()
      print()
      constraints_list = None
      t_celsius = 20
    else:
      raise RuntimeError("Only CIF format currently supported!")

  else:
    test_case = tst_constrained_structure.sucrose_test_case(None)
    t_celsius = test_case.t_celsius
    xs = test_case.xray_structure
    constraints_list = test_case.constraints

  #from cctbx import adptbx
  #for sc in xs.scatterers():
    #if sc.flags.use_u_aniso():
      #sc.u_iso = adptbx.u_star_as_u_iso(xs.unit_cell(), sc.u_star)
      #sc.set_use_u_iso_only()

  if not params.constrained_refinement:
    constraints_list = []

  exercise_constrained_lbfgs(xray_structure=xs,
                             constraints_list=constraints_list,
                             t_celsius=t_celsius,
                             d_min=params.d_min,
                             shake_sites_rmsd=params.shake_sites_rmsd,
                             shake_u_iso_spread=params.shake_adp_spread,
                             grad_site=params.grad_site,
                             grad_u_iso=params.grad_u_iso,
                             grad_u_aniso=params.grad_u_aniso,
                             grad_occupancy=params.grad_occupancy,
                             grad_fp_fdp=params.grad_fp_fdp,
                             lbfgs_m=params.lbfgs_m,
                             lbfgs_max_iterations=params.lbfgs_max_iterations,
                             verbose=params.verbose)

def exercise_constrained_lbfgs(xray_structure,
                               constraints_list,
                               t_celsius,
                               d_min=0.5,
                               shake_sites_rmsd=0.5,
                               shake_u_iso_spread=0,
                               shake_u_aniso_spread=0,
                               grad_site=True,
                               grad_u_iso=True,
                               grad_u_aniso=False,
                               grad_occupancy=False,
                               grad_fp_fdp=False,
                               lbfgs_m=5,
                               lbfgs_max_iterations=1000,
                               verbose=0):

  xs = xray_structure
  xray.set_scatterer_grad_flags(scatterers=xs.scatterers(),
                                site=grad_site,
                                u_iso=grad_u_iso,
                                u_aniso=grad_u_aniso,
                                occupancy=grad_occupancy,
                                fp=grad_fp_fdp,
                                fdp=grad_fp_fdp,
                                tan_u_iso=False,
                                param=0)

  xs0 = xs.deep_copy_scatterers()
  mi = xs0.build_miller_set(anomalous_flag=False, d_min=d_min)
  fo_sq = mi.structure_factors_from_scatterers(
    xs0, algorithm="direct").f_calc().norm()
  fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.size(), 1))
  fo_sq.set_observation_type_xray_intensity()
  y_obs = fo_sq
  #y_obs = fo_sq.f_sq_as_f()
  if grad_site:
    xs.shake_sites_in_place(rms_difference=shake_sites_rmsd)
  if not grad_u_aniso: shake_u_aniso_spread = 0
  if not grad_u_iso: shake_u_iso_spread = 0
  if grad_u_aniso or grad_u_iso:
    xs.shake_adp(spread=shake_u_iso_spread, aniso_spread=shake_u_aniso_spread)
  xs1 = xs.deep_copy_scatterers()

  core_params = scitbx.lbfgs.core_parameters(m=lbfgs_m, maxfev=100, xtol=1e-5)

  connectivity_table = smtbx.utils.connectivity_table(xs0)

  if constraints_list is None:
    from smtbx.development import generate_hydrogen_constraints
    constraints_list = generate_hydrogen_constraints(
      structure=xs0, connectivity_table=connectivity_table)

  reparametrisation = constraints.reparametrisation(
    xs,
    constraints_list,
    connectivity_table,
    temperature=t_celsius)

  lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
    traditional_convergence_test=False,
    drop_convergence_test_max_drop_eps=1.e-20,
    drop_convergence_test_iteration_coefficient=1,
    min_iterations=500,
    max_iterations=lbfgs_max_iterations)

  minimizer = lbfgs(
    target_functor=xray.target_functors.unified_least_squares_residual(y_obs),
    xray_structure=xs,
    reparametrisation=reparametrisation,
    structure_factor_algorithm="direct",
    lbfgs_termination_params=lbfgs_termination_params,
    lbfgs_core_params=core_params,
    reference_structure=xs0,
    verbose=verbose)

  if verbose > 0:
    print("Total parameters: ", xs.n_parameters())
    print("Independent parameters: ", reparametrisation.n_independents)

    print("Reference model: ")
    xs0.show_angles(distance_cutoff=1.5)
    print()
    print("Starting model: ")
    xs1.show_angles(distance_cutoff=1.5)
    print()
    print("Refined model: ")
    xs.show_angles(distance_cutoff=1.5)
    print()

    print("n_iter, n_fun: ", minimizer.minimizer.iter(), minimizer.minimizer.nfun())

  h_selection = xs.element_selection('H')

  diff = xray.meaningful_site_cart_differences(
    xs0.select(h_selection, negate=True),
    xs.select(h_selection, negate=True))
  #diff = xray.meaningful_site_cart_differences(xs0, xs)
  #assert diff.max_absolute() < 1e-3
  if verbose > 0:
    diff.show()
    print()
  assert diff.max_absolute() < 2e-2, diff.max_absolute()

  diff = xray.meaningful_site_cart_differences(xs0, xs)
  if verbose > 0:
    diff.show()
    print()
  # XXX why does this tolerance have to be so high?
  assert diff.max_absolute() < 0.5, diff.max_absolute()

  #ls = least_squares.crystallographic_ls(
    #fo_sq.as_xray_observations(),
    #reparametrisation,
    #restraints_manager=None,
    #weighting_scheme=least_squares.sigma_weighting())

  #steps = scitbx.lstbx.normal_eqns_solving.naive_iterations(
    #non_linear_ls=ls,
    #n_max_iterations=100,
  #)
  #steps.do()

  #diff = xray.meaningful_site_cart_differences(xs0, xs)
  #assert diff.max_absolute() < 1e-4


if __name__ == '__main__':
  import sys
  libtbx.utils.show_times_at_exit()
  run(sys.argv[1:])
