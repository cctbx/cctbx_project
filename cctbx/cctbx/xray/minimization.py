import cctbx.xray.structure_factors
from cctbx.xray import ext
from cctbx.xray.structure import structure as cctbx_xray_structure
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.lbfgs
import scitbx.math
from scitbx.python_utils.misc import adopt_init_args
from libtbx import introspection
from stdlib import math

class u_penalty_singular_at_zero:

  def __init__(self,
        penalty_factor=1,
        penalty_scale=8*math.pi**2,
        u_min=1.e-6,
        min_functional=1.e-10):
    introspection.adopt_init_args()
    self.u_max = scitbx.math.lambertw(penalty_factor/min_functional) \
               / (penalty_factor * penalty_scale)

  def functional(self, u):
    if (u > self.u_max): return 0
    u = max(u, self.u_min)
    s = self.penalty_scale
    return 1/(s*u)*math.exp(-(s*u)*self.penalty_factor)

  def gradient(self, u):
    if (u > self.u_max): return 0
    u = max(u, self.u_min)
    s = self.penalty_scale
    return -1/(math.exp(s*self.penalty_factor*u)*s*u**2) \
           - self.penalty_factor/(math.exp(s*self.penalty_factor*u)*u)

class u_penalty_exp:

  def __init__(self,
        penalty_factor=1,
        penalty_scale=10*8*math.pi**2,
        min_functional=1.e-10):
    introspection.adopt_init_args()
    self.u_max = -math.log(min_functional) / (penalty_factor * penalty_scale)

  def functional(self, u):
    if (u > self.u_max): return 0
    s = self.penalty_scale
    return math.exp(-s*u*self.penalty_factor)

  def gradient(self, u):
    if (u > self.u_max): return 0
    s = self.penalty_scale
    return -s*self.penalty_factor*math.exp(-s*u*self.penalty_factor)

class lbfgs:

  def __init__(self, target_functor, gradient_flags, xray_structure,
                     u_penalty=None,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None,
                     cos_sin_table=True,
                     structure_factor_algorithm=None,
                     verbose=0):
    adopt_init_args(self, locals())
    self.structure_factors_from_scatterers = \
      cctbx.xray.structure_factors.from_scatterers(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.structure_factor_gradients = \
      cctbx.xray.structure_factors.gradients(
        miller_set=self.target_functor.f_obs(),
        cos_sin_table=cos_sin_table)
    self.n = xray_structure.n_parameters(gradient_flags)
    self.x = flex.double(self.n, 0)
    xray_structure.tidy_us(u_min=1.e-6)
    self._scatterers_start = xray_structure.scatterers()
    self._scattering_dict = xray_structure.scattering_dict()
    self._d_min = self.target_functor.f_obs().d_min()
    self.first_target_value = None
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    del self._scatterers_start
    del self._scattering_dict
    del self._d_min
    self.compute_target(compute_gradients=False)
    self.final_target_value = self.target_result.target()

  def apply_shifts(self):
    unit_cell = self.xray_structure.unit_cell()
    apply_shifts_result = ext.minimization_apply_shifts(
      unit_cell=unit_cell,
      scatterers=self._scatterers_start,
      gradient_flags=self.gradient_flags,
      shifts=self.x)
    shifted_scatterers = apply_shifts_result.shifted_scatterers
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      shifted_scatterers[i_seq].site = crystal.correct_special_position(
        unit_cell=unit_cell,
        special_op=site_symmetry_table.get(i_seq).special_op(),
        site_frac=shifted_scatterers[i_seq].site)
    self.xray_structure.replace_scatterers(scatterers=shifted_scatterers)
    return apply_shifts_result.mean_displacements

  def compute_target(self, compute_gradients):
    self.f_calc = self.structure_factors_from_scatterers(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      algorithm=self.structure_factor_algorithm).f_calc()
    self.target_result = self.target_functor(
      self.f_calc,
      compute_gradients)

  def __call__(self):
    mean_displacements = self.apply_shifts()
    self.compute_target(compute_gradients=True)
    self.f = self.target_result.target()
    if (self.first_target_value is None):
      self.first_target_value = self.f
    if (self.u_penalty is not None and self.gradient_flags.u_iso):
      u_isos = self.xray_structure.scatterers().extract_u_iso(
        unit_cell=self.xray_structure.unit_cell())
      for u_iso in u_isos:
        self.f += self.u_penalty.functional(u=u_iso)
    self.g = self.structure_factor_gradients(
      xray_structure=self.xray_structure,
      mean_displacements=mean_displacements,
      miller_set=self.target_functor.f_obs(),
      d_target_d_f_calc=self.target_result.derivatives(),
      gradient_flags=self.gradient_flags,
      n_parameters=self.x.size(),
      algorithm=self.structure_factor_algorithm).packed()
    if (self.u_penalty is not None and self.gradient_flags.u_iso):
      g = flex.double()
      if (self.gradient_flags.sqrt_u_iso):
        for mean_displacement in mean_displacements:
          g.append(2*mean_displacement
                  *self.u_penalty.gradient(u=mean_displacement**2))
      else:
        for u_iso in u_isos:
          g.append(self.u_penalty.gradient(u=u_iso))
      del u_isos
      ext.minimization_add_u_iso_gradients(
        scatterers=self.xray_structure.scatterers(),
        gradient_flags=self.gradient_flags,
        xray_gradients=self.g,
        u_iso_gradients=g)
    if (self.verbose > 1):
      print "xray.minimization line search: f,rms(g):",
      print self.f, math.sqrt(flex.mean_sq(self.g))
    return self.x, self.f, self.g

  def callback_after_step(self, minimizer):
    if (self.verbose > 0):
      print "xray.minimization step: f,iter,nfun:",
      print self.f,minimizer.iter(),minimizer.nfun()
