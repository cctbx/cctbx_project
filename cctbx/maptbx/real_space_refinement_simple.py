from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.lbfgs
from cctbx import miller

class lbfgs(object):

  def __init__(O,
        sites_cart,
        density_map,
        gradients_delta,
        real_space_weight_scale,
        unit_cell,
        geometry_restraints_manager = None,
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None):
    O.gr = geometry_restraints_manager
    O.unit_cell = unit_cell
    O.density_map = density_map
    O.gradients_delta = gradients_delta
    O.rs_weight_scale = real_space_weight_scale
    O.rs_weight = None
    O.x = sites_cart.as_double()
    O.minimizer = scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=lbfgs_termination_params,
      exception_handling_params=lbfgs_exception_handling_params)
    O.f, O.g = O.compute_functional_and_gradients()
    O.sites_cart = flex.vec3_double(O.x)
    del O.x

  def compute_functional_and_gradients(O):
    sites_cart = flex.vec3_double(O.x)
    gr_e, gr_f, gr_g = [None]*3
    if(O.gr is not None):
      gr_e = O.gr.energies_sites(sites_cart=sites_cart, compute_gradients=True)
      gr_f = gr_e.target
      gr_g = gr_e.gradients
    rs_f = -maptbx.real_space_target_simple(
      unit_cell=O.unit_cell,
      density_map=O.density_map,
      sites_cart=sites_cart)
    rs_g = maptbx.real_space_gradients_simple(
      unit_cell=O.unit_cell,
      density_map=O.density_map,
      sites_cart=sites_cart,
      delta=O.gradients_delta)
    rs_g *= -1
    if (O.rs_weight is None and O.gr is not None):
      rms_gr_g = flex.mean_sq(gr_g.as_double())**0.5
      rms_rs_g = flex.mean_sq(rs_g.as_double())**0.5
      if (rms_rs_g < 1):
        O.rs_weight = 1
      else:
        O.rs_weight = rms_gr_g / rms_rs_g * O.rs_weight_scale
    if(O.gr is not None):
      f = gr_f + rs_f * O.rs_weight
      g = gr_g + rs_g * O.rs_weight
    else:
      f = rs_f
      g = rs_g
    return f, g.as_double()

class minimization(object):
  def __init__(self,
        xray_structure,
        miller_array,
        crystal_gridding,
        map_target,
        max_iterations,
        min_iterations,
        step):
    self.step = step
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
    self.xray_structure = self.xray_structure.replace_sites_cart(self.sites_cart)

  def compute_functional_and_gradients(self):
    self.sites_cart = flex.vec3_double(self.x)
    self.xray_structure = self.xray_structure.replace_sites_cart(self.sites_cart)
    self.map_current = self.compute_map()
    f = 0
    g = flex.vec3_double(self.sites_cart.size(), (0,0,0))
    f += self.target()
    g = self.gradients()
    return f, g.as_double()

  def target(self):
    r = self.map_target-self.map_current
    r = r*r
    r = flex.sum(r)
    return r

  def gradients(self):
    g = flex.vec3_double()
    uc = self.xray_structure.unit_cell()
    r = 2 * (self.map_target - self.map_current)
    step_x = flex.double([self.step,0,0])
    step_y = flex.double([0,self.step,0])
    step_z = flex.double([0,0,self.step])
    for site_cart in self.xray_structure.sites_cart():
      site_cart = flex.double(site_cart)
      sxp = uc.fractionalize(tuple(site_cart+step_x))
      sxm = uc.fractionalize(tuple(site_cart-step_x))
      syp = uc.fractionalize(tuple(site_cart+step_y))
      sym = uc.fractionalize(tuple(site_cart-step_y))
      szp = uc.fractionalize(tuple(site_cart+step_z))
      szm = uc.fractionalize(tuple(site_cart-step_z))
      gx = (self.map_target.eight_point_interpolation(sxp) -
            self.map_target.eight_point_interpolation(sxm)) / (2*self.step)
      gy = (self.map_target.eight_point_interpolation(syp) -
            self.map_target.eight_point_interpolation(sym)) / (2*self.step)
      gz = (self.map_target.eight_point_interpolation(szp) -
            self.map_target.eight_point_interpolation(szm)) / (2*self.step)
      if 0: # Equivalent to lines below
        s = 2*(map_target.eight_point_interpolation(site)-
          map_current.eight_point_interpolation(site))
        g.append([gx*s,gy*s,gz*s])
      site_frac = uc.fractionalize(tuple(site_cart))
      x = (gx*r).eight_point_interpolation(site_frac)
      y = (gy*r).eight_point_interpolation(site_frac)
      z = (gz*r).eight_point_interpolation(site_frac)
      g.append([x,y,z])
    return g

  def compute_map(self):
    map_coefficients = self.miller_array.structure_factors_from_scatterers(
      xray_structure = self.xray_structure).f_calc()
    fft_map = miller.fft_map(crystal_gridding     = self.crystal_gridding,
                             fourier_coefficients = map_coefficients)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result
