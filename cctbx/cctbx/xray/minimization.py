import cctbx.xray.structure_factors
from cctbx.xray import ext
from cctbx.xray.structure import structure as cctbx_xray_structure
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.lbfgs
import scitbx.math
from libtbx import adopt_init_args
from stdlib import math
from cctbx import adptbx

def add_gradients(
      scatterers,
      xray_gradients,
      site_gradients      = None,
      u_iso_gradients     = None,
      u_aniso_gradients   = None,
      occupancy_gradients = None):
  ext.minimization_add_gradients(scatterers          = scatterers,
                                 xray_gradients      = xray_gradients,
                                 site_gradients      = site_gradients,
                                 u_iso_gradients     = u_iso_gradients,
                                 u_aniso_gradients   = u_aniso_gradients,
                                 occupancy_gradients = occupancy_gradients)

class occupancy_penalty_exp(object):

  def __init__(self,
        penalty_factor=1,
        penalty_scale=100,
        min_functional=1.e-10):
    adopt_init_args(self, locals())
    self.occupancy_max = -math.log(min_functional) \
                       / (penalty_factor * penalty_scale)

  def functional(self, occupancy):
    if (occupancy > self.occupancy_max): return 0
    if (occupancy < -1): occupancy = -1 - math.log(-occupancy)
    s = self.penalty_scale
    return math.exp(-s*occupancy*self.penalty_factor)

  def gradient(self, occupancy):
    if (occupancy > self.occupancy_max): return 0
    if (occupancy < -1): occupancy = -1 - math.log(-occupancy)
    s = self.penalty_scale
    return -s*self.penalty_factor*math.exp(-s*occupancy*self.penalty_factor)

class lbfgs(object):

  def __init__(self, target_functor, xray_structure,
                     occupancy_penalty=None,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None,
                     correct_special_position_tolerance=1.e-2,
                     cos_sin_table=True,
                     structure_factor_algorithm=None,
                     verbose=0):
    adopt_init_args(self, locals())
    self.scatterer_grad_flags_counts = ext.scatterer_grad_flags_counts(
                                              self.xray_structure.scatterers())
    self.grad_flags_counts = \
              ext.scatterer_grad_flags_counts(self.xray_structure.scatterers())
    self.structure_factors_from_scatterers = \
      cctbx.xray.structure_factors.from_scatterers(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.structure_factor_gradients = \
      cctbx.xray.structure_factors.gradients(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.x = flex.double(xray_structure.n_parameters_XXX(), 0)
    xray_structure.tidy_us(u_min=1.e-6)
    self._scatterers_start = xray_structure.scatterers()
    self._d_min = self.target_functor.f_obs().d_min()
    self.first_target_value = None
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    del self._scatterers_start
    del self._d_min
    self.compute_target(compute_gradients=False)
    self.final_target_value = self.target_result.target()

  def apply_shifts(self):
    apply_shifts_result = ext.minimization_apply_shifts(
      unit_cell=self.xray_structure.unit_cell(),
      scatterers=self._scatterers_start,
      shifts=self.x)
    shifted_scatterers = apply_shifts_result.shifted_scatterers
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      shifted_scatterers[i_seq].site = crystal.correct_special_position(
        crystal_symmetry=self.xray_structure,
        special_op=site_symmetry_table.get(i_seq).special_op(),
        site_frac=shifted_scatterers[i_seq].site,
        tolerance=self.correct_special_position_tolerance)
    self.xray_structure.replace_scatterers(scatterers=shifted_scatterers)
    return apply_shifts_result.u_iso_reinable_params

  def compute_target(self, compute_gradients):
    self.f_calc = self.structure_factors_from_scatterers(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      algorithm=self.structure_factor_algorithm).f_calc()
    self.target_result = self.target_functor(
      self.f_calc,
      compute_gradients)

  def compute_functional_and_gradients(self):
    u_iso_reinable_params = self.apply_shifts()
    self.compute_target(compute_gradients=True)
    self.f = self.target_result.target()
    if (self.first_target_value is None):
      self.first_target_value = self.f
    if (self.occupancy_penalty is not None
        and self.grad_flags_counts != 0):
      occupancies = self.xray_structure.scatterers().extract_occupancies()
      for occupancy in occupancies:
        self.f += self.occupancy_penalty.functional(occupancy=occupancy)
    self.g = self.structure_factor_gradients(
      xray_structure=self.xray_structure,
      u_iso_reinable_params=u_iso_reinable_params,
      miller_set=self.target_functor.f_obs(),
      d_target_d_f_calc=self.target_result.derivatives(),
      n_parameters=self.x.size(),
      algorithm=self.structure_factor_algorithm).packed()
    if (self.occupancy_penalty is not None
        and self.grad_flags_counts != 0):
      g = flex.double()
      for occupancy in occupancies:
        g.append(self.occupancy_penalty.gradient(occupancy=occupancy))
      del occupancies
      add_gradients(
        scatterers=self.xray_structure.scatterers(),
        xray_gradients=self.g,
        occupancy_gradients=g)
      del g
    if (self.verbose > 1):
      print "xray.minimization line search: f,rms(g):",
      print self.f, math.sqrt(flex.mean_sq(self.g))
    return self.f, self.g

  def callback_after_step(self, minimizer):
    if (self.verbose > 0):
      print "xray.minimization step: f,iter,nfun:",
      print self.f,minimizer.iter(),minimizer.nfun()
