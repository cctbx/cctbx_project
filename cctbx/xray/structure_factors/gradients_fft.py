from cctbx.xray.structure_factors.gradients_base import gradients_base
from cctbx.xray.structure_factors.misc import expensive_function_call_message
from cctbx.xray.structure_factors import global_counters
from cctbx.xray import ext
from cctbx import maptbx
from cctbx import adptbx
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.utils import user_plus_sys_time
from libtbx import introspection

class gradients_fft(gradients_base):

  def __init__(self, manager,
                     xray_structure,
                     u_iso_refinable_params,
                     miller_set,
                     d_target_d_f_calc,
                     n_parameters):
    time_all = time_apply_u_extra = user_plus_sys_time()
    gradients_base.__init__(self,
      manager, xray_structure, miller_set, algorithm="fft")
    self._d_target_d_f_calc = d_target_d_f_calc
    manager.setup_fft() # before timing
    time_apply_u_extra = user_plus_sys_time()
    self._results = ext.fast_gradients(
      unit_cell=xray_structure.unit_cell(),
      scatterers=xray_structure.scatterers(),
      scattering_type_registry=xray_structure.scattering_type_registry(),
      u_base=manager.u_base(),
      wing_cutoff=manager.wing_cutoff(),
      exp_table_one_over_step_size=manager.exp_table_one_over_step_size(),
      tolerance_positive_definite=manager.tolerance_positive_definite())
    coeff = self._gradient_map_coeff()
    time_apply_u_extra = time_apply_u_extra.elapsed()
    time_from_or_to_map = user_plus_sys_time()
    coeff_map = self._gradient_map_coeff_to_map(coeff)
    time_from_or_to_map = time_from_or_to_map.elapsed()
    time_fft = user_plus_sys_time()
    if (not coeff.anomalous_flag()):
      gradient_map = manager.rfft().backward(coeff_map.complex_map())
    else:
      gradient_map = manager.cfft().backward(coeff_map.complex_map())
    time_fft = time_fft.elapsed()
    time_sampling = user_plus_sys_time()
    self._results.sampling(
      scatterers=xray_structure.scatterers(),
      u_iso_refinable_params=u_iso_refinable_params,
      scattering_type_registry=xray_structure.scattering_type_registry(),
      site_symmetry_table=xray_structure.site_symmetry_table(),
      ft_d_target_d_f_calc=gradient_map,
      n_parameters=n_parameters,
      sampled_density_must_be_positive=
        manager.sampled_density_must_be_positive())
    time_sampling = time_sampling.elapsed()
    introspection.virtual_memory_info().update_max()
    manager.estimate_time_fft.register(
      n_scatterers=xray_structure.scatterers().size(),
      n_miller_indices=miller_set.indices().size(),
      time_sampling=time_sampling,
      time_fft=time_fft,
      time_from_or_to_map=time_from_or_to_map,
      time_apply_u_extra=time_apply_u_extra)
    self.d_target_d_site_frac_was_used = False
    self.d_target_d_u_star_was_used = False
    global_counters.gradients_fft.process(time_all.elapsed())

  def _gradient_map_coeff(self):
    coeff = self.miller_set().array(data=self.d_target_d_f_calc().deep_copy())
    multiplier = (  self.manager().unit_cell().volume()
                  / matrix.row(self.manager().rfft().n_real()).product()
                  * self.manager().space_group().n_ltr())
    if (    not coeff.anomalous_flag()
        and not coeff.space_group().is_centric()):
      multiplier /= 2
    ext.apply_u_extra(
      self.manager().unit_cell(),
      self._results.u_extra(),
      coeff.indices(),
      coeff.data(),
      multiplier)
    return coeff

  def _gradient_map_coeff_to_map(self, coeff):
    if (not coeff.anomalous_flag()):
      n_complex = self.manager().rfft().n_complex()
    else:
      n_complex = self.manager().rfft().n_real()
    return maptbx.structure_factors.to_map(
      space_group=coeff.space_group(),
      anomalous_flag=coeff.anomalous_flag(),
      miller_indices=coeff.indices(),
      structure_factors=coeff.data(),
      n_real=self.manager().rfft().n_real(),
      map_grid=flex.grid(n_complex),
      conjugate_flag=True,
      treat_restricted=False)

  def d_target_d_site_cart(self):
    return self.check_size(self._results.d_target_d_site_cart())

  def d_target_d_site_frac(self):
    if (self.d_target_d_site_frac_was_used):
      raise RuntimeError(expensive_function_call_message)
    self.d_target_d_site_frac_was_used = True
    return self.d_target_d_site_cart() \
         * self.xray_structure().unit_cell().orthogonalization_matrix()

  def d_target_d_u_cart(self):
    return self.check_size(self._results.d_target_d_u_cart())

  def d_target_d_u_star(self):
    if (self.d_target_d_u_star_was_used):
      raise RuntimeError(expensive_function_call_message)
    self.d_target_d_u_star_was_used = True
    return adptbx.grad_u_cart_as_u_star(
      self.xray_structure().unit_cell(), self.d_target_d_u_cart())
