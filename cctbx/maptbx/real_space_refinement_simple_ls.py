from cctbx.array_family import flex
import scitbx.lbfgs
from cctbx import miller
from cctbx import maptbx
from libtbx.utils import Sorry, user_plus_sys_time

FFT_TIME = 0

class minimization(object):
  def __init__(self,
        xray_structure,
        miller_array,
        crystal_gridding,
        map_target,
        step,
        max_iterations=30,
        min_iterations=25,
        geometry_restraints_manager=None,
        real_space_target_weight=1):
    self.step = step
    self.geometry_restraints_manager = geometry_restraints_manager
    self.real_space_target_weight = real_space_target_weight
    self.sites_cart = xray_structure.sites_cart()
    self.x = self.sites_cart.as_double()
    self.xray_structure = xray_structure
    self.map_target = map_target
    self.miller_array = miller_array
    self.crystal_gridding = crystal_gridding
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations,
        min_iterations=min_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))
    self.sites_cart.clear()
    self.sites_cart.extend(flex.vec3_double(self.x))
    self.xray_structure = self.xray_structure.replace_sites_cart(
      self.sites_cart)

  def compute_functional_and_gradients(self):
    self.sites_cart = flex.vec3_double(self.x)
    self.xray_structure = self.xray_structure.replace_sites_cart(self.sites_cart)
    self.map_current = self.compute_map()
    self.mngr = maptbx.target_and_gradients(
      unit_cell   = self.xray_structure.unit_cell(),
      map_target  = self.map_target,
      map_current = self.map_current,
      step        = self.step,
      sites_frac  = self.xray_structure.sites_frac()) 
    rs_f = self.mngr.target()
    rs_g = self.mngr.gradients()
    if(self.geometry_restraints_manager):
      rs_f *= self.real_space_target_weight
      rs_g *= self.real_space_target_weight
      gr_e = self.geometry_restraints_manager.energies_sites(
        sites_cart = self.sites_cart, compute_gradients = True)
      f = rs_f + gr_e.target
      g = rs_g + gr_e.gradients
    else:
      f = rs_f
      g = rs_g
    return f, g.as_double()

  def compute_map(self):
    global FFT_TIME
    timer = user_plus_sys_time()
    map_coefficients = self.miller_array.structure_factors_from_scatterers(
      xray_structure = self.xray_structure).f_calc()
    fft_map = miller.fft_map(crystal_gridding     = self.crystal_gridding,
                             fourier_coefficients = map_coefficients)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    FFT_TIME += timer.elapsed()
    return result
