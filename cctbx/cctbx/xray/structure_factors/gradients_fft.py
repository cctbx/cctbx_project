from cctbx.xray.structure_factors.gradients_base import gradients_base
from cctbx.xray.structure_factors.misc import expensive_function_call_message
from cctbx.xray import ext
from cctbx import miller
from cctbx import adptbx
from cctbx import matrix
from scitbx.python_utils.misc import user_plus_sys_time

class gradients_fft(gradients_base):

  def __init__(self, manager,
                     xray_structure,
                     miller_set,
                     d_target_d_f_calc,
                     gradient_flags,
                     electron_density_must_be_positive=0001):
    gradients_base.__init__(self, manager, xray_structure, miller_set)
    self._d_target_d_f_calc = d_target_d_f_calc
    manager.setup_fft()
    # XXX timing
    gradient_map = self.fft_d_target_d_f_calc(
      d_target_d_f_calc=d_target_d_f_calc)
    self._results = ext.fast_gradients(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      gradient_map.complex_map(),
      gradient_flags,
      manager.u_extra(),
      manager.wing_cutoff(),
      manager.exp_table_one_over_step_size(),
      electron_density_must_be_positive)
    self.d_target_d_site_frac_was_used = 00000
    self.d_target_d_u_star_was_used = 00000

  def fft_d_target_d_f_calc(self, d_target_d_f_calc):
    multiplier = self.miller_set().epsilons().data().as_double() * (
                   self.manager().unit_cell().volume()
                 / matrix.row(self.manager().rfft().n_real()).product()
                 * self.manager().space_group().n_ltr())
    coeff = d_target_d_f_calc.deep_copy()
    ext.apply_u_extra(
      self.manager().unit_cell(),
      self.manager().u_extra(),
      self.miller_set().indices(),
      coeff,
      multiplier)
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
