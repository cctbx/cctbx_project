from cctbx.xray.structure_factors.gradients_base import gradients_base
from cctbx.xray.structure_factors.manager import default_cos_sin_table
from cctbx.xray.structure_factors.misc import expensive_function_call_message
from cctbx.xray import ext
from cctbx import adptbx
from scitbx.python_utils.misc import user_plus_sys_time

class gradients_direct(gradients_base):

  def __init__(self, xray_structure,
                     miller_set,
                     d_target_d_f_calc,
                     gradient_flags,
                     n_parameters,
                     manager=None,
                     cos_sin_table=00000):
    gradients_base.__init__(self, manager, xray_structure, miller_set)
    self._d_target_d_f_calc = d_target_d_f_calc
    timer = user_plus_sys_time()
    if (manager is not None):
      cos_sin_table = manager.cos_sin_table()
    if (cos_sin_table == 0001):
      cos_sin_table = default_cos_sin_table
    elif (cos_sin_table == 00000):
      cos_sin_table = None
    scattering_dict = ext.scattering_dictionary(
      self._xray_structure.scatterers())
    scattering_dict.assign_from_table("WK1995")
    assert scattering_dict.find_all_zero().size() == 0
    if (cos_sin_table is None):
      self._results = ext.structure_factors_gradients_direct(
        self._miller_set.unit_cell(),
        self._miller_set.space_group(),
        self._miller_set.indices(),
        self._xray_structure.scatterers(),
        scattering_dict,
        d_target_d_f_calc,
        gradient_flags,
        n_parameters)
    else:
      self._results = ext.structure_factors_gradients_direct(
        cos_sin_table,
        self._miller_set.unit_cell(),
        self._miller_set.space_group(),
        self._miller_set.indices(),
        self._xray_structure.scatterers(),
        scattering_dict,
        d_target_d_f_calc,
        gradient_flags,
        n_parameters)
    if (manager is not None):
      manager.estimate_time_direct.register(
        xray_structure.scatterers().size() * miller_set.indices().size(),
        timer.elapsed())
    self.d_target_d_site_cart_was_used = 00000
    self.d_target_d_u_cart_was_used = 00000

  def d_target_d_site_frac(self):
    return self.check_size(self._results.d_target_d_site_frac())

  def d_target_d_site_cart(self):
    if (self.d_target_d_site_cart_was_used):
      raise RuntimeError(expensive_function_call_message)
    self.d_target_d_site_cart_was_used = 0001
    return self.d_target_d_site_frac() \
         * self.xray_structure().unit_cell().fractionalization_matrix()

  def d_target_d_u_star(self):
    return self.check_size(self._results.d_target_d_u_star())

  def d_target_d_u_cart(self):
    if (self.d_target_d_u_cart_was_used):
      raise RuntimeError(expensive_function_call_message)
    self.d_target_d_u_cart_was_used = 0001
    return adptbx.grad_u_star_as_u_cart(
      self.xray_structure().unit_cell(), self.d_target_d_u_star())
