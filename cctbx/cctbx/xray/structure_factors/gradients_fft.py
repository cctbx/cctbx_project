from cctbx.xray.structure_factors.gradients_base import gradients_base
from cctbx.xray.structure_factors.misc import expensive_function_call_message
from cctbx.xray import ext
from cctbx import miller
from cctbx import adptbx
from cctbx import matrix
from cctbx.array_family import flex
from scitbx.python_utils.misc import user_plus_sys_time

class gradients_fft(gradients_base):

  def __init__(self, manager,
                     xray_structure,
                     miller_set,
                     d_target_d_f_calc,
                     gradient_flags):
    gradients_base.__init__(self, manager, xray_structure, miller_set)
    self._d_target_d_f_calc = d_target_d_f_calc
    manager.setup_fft() # before timing
    time_process_coefficients = user_plus_sys_time()
    coeff = self.gradient_map_coeff()
    time_process_coefficients = time_process_coefficients.elapsed()
    time_fft = user_plus_sys_time()
    gradient_map = self.fft_gradient_map_coeff(coeff)
    if (not gradient_map.anomalous_flag()):
      gradient_map_real = gradient_map.real_map()
      gradient_map_complex = flex.complex_double(flex.grid(0,0,0))
    else:
      gradient_map_real = flex.double(flex.grid(0,0,0))
      gradient_map_complex = gradient_map.complex_map()
    time_fft = time_fft.elapsed()
    time_sampling = user_plus_sys_time()
    self._results = ext.fast_gradients(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      gradient_map_real,
      gradient_map_complex,
      gradient_flags,
      manager.u_extra(),
      manager.wing_cutoff(),
      manager.exp_table_one_over_step_size(),
      manager.electron_density_must_be_positive())
    time_sampling = time_sampling.elapsed()
    manager.estimate_time_fft.register(
      n_scatterers=xray_structure.scatterers().size(),
      n_miller_indices=miller_set.indices().size(),
      time_sampling=time_sampling,
      time_symmetry_mapping=0,
      time_fft=time_fft,
      time_process_coefficients=time_process_coefficients)
    self.d_target_d_site_frac_was_used = 00000
    self.d_target_d_u_star_was_used = 00000

  def gradient_map_coeff(self):
    multiplier = (  self.manager().unit_cell().volume()
                  / matrix.row(self.manager().rfft().n_real()).product()
                  * self.manager().space_group().order_z()
                  / self.miller_set().multiplicities().data().as_double())
    coeff = self.d_target_d_f_calc().deep_copy()
    ext.apply_u_extra(
      self.manager().unit_cell(),
      self.manager().u_extra(),
      self.miller_set().indices(),
      coeff,
      multiplier)
    return coeff

  def fft_gradient_map_coeff(self, coeff):
    return miller.fft_map(
      crystal_gridding=self.manager().crystal_gridding(),
      fourier_coefficients=miller.array(
        miller_set=self.miller_set(),
        data=coeff))

  def d_target_d_site_cart(self):
    return self.check_size(self._results.d_target_d_site_cart())

  def d_target_d_site_frac(self):
    if (self.d_target_d_site_frac_was_used):
      raise RuntimeError(expensive_function_call_message)
    self.d_target_d_site_frac_was_used = 0001
    return self.d_target_d_site_cart() \
         * self.xray_structure().unit_cell().orthogonalization_matrix()

  def d_target_d_u_cart(self):
    return self.check_size(self._results.d_target_d_u_cart())

  def d_target_d_u_star(self):
    if (self.d_target_d_u_star_was_used):
      raise RuntimeError(expensive_function_call_message)
    self.d_target_d_u_star_was_used = 0001
    return adptbx.grad_u_cart_as_u_star(
      self.xray_structure().unit_cell(), self.d_target_d_u_cart())
