import cctbx.xray.structure_factors
from cctbx.xray import ext
from cctbx.xray.structure import structure as cctbx_xray_structure
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.lbfgs
from scitbx.python_utils.misc import adopt_init_args
from stdlib import math

class lbfgs:

  def __init__(self, target_functor, gradient_flags, xray_structure,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None,
                     cos_sin_table=0001,
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
    self.compute_target(compute_gradients=00000)
    self.final_target_value = self.target_result.target()

  def apply_shifts(self):
    unit_cell = self.xray_structure.unit_cell()
    scatterers_shifted = ext.minimization_apply_shifts(
      unit_cell=unit_cell,
      space_group_type=self.xray_structure.space_group_info().type(),
      scatterers=self._scatterers_start,
      scattering_dict=self._scattering_dict,
      gradient_flags=self.gradient_flags,
      shifts=self.x,
      d_min=self._d_min)
    if (0): # XXX ASSERTION FAILURE
      site_symmetry_table = self.xray_structure.site_symmetry_table()
      for i_seq in site_symmetry_table.special_position_indices():
        scatterers_shifted[i_seq].site = crystal.correct_special_position(
          unit_cell=unit_cell,
          special_op=site_symmetry_table.get(i_seq).special_op(),
          site_frac=scatterers_shifted[i_seq].site)
    self.xray_structure.replace_scatterers(scatterers_shifted)

  def compute_target(self, compute_gradients):
    self.f_calc = self.structure_factors_from_scatterers(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      algorithm=self.structure_factor_algorithm).f_calc()
    self.target_result = self.target_functor(
      self.f_calc,
      compute_gradients)

  def __call__(self):
    if (self.first_target_value is None):
      assert self.x.all_eq(0)
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    self.f = self.target_result.target()
    if (self.first_target_value is None):
      self.first_target_value = self.f
    sf = self.structure_factor_gradients(
      xray_structure=self.xray_structure,
      miller_set=self.target_functor.f_obs(),
      d_target_d_f_calc=self.target_result.derivatives(),
      gradient_flags=self.gradient_flags,
      n_parameters=self.x.size(),
      algorithm=self.structure_factor_algorithm)
    self.g = sf.packed()
    if (self.verbose > 1):
      print "xray.minimization line search: f,rms(g):",
      print self.f, math.sqrt(flex.mean_sq(self.g))
    return self.x, self.f, self.g

  def callback_after_step(self, minimizer):
    if (self.verbose > 0):
      print "xray.minimization step: f,iter,nfun:",
      print self.f,minimizer.iter(),minimizer.nfun()
