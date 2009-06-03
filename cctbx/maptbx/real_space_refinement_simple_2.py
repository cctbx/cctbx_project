from cctbx.array_family import flex
import scitbx.lbfgs
from cctbx import miller
from cctbx import maptbx

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
    self.xray_structure = self.xray_structure.replace_sites_cart(
      self.sites_cart)

  def compute_functional_and_gradients(self):
    self.sites_cart = flex.vec3_double(self.x)
    self.xray_structure = self.xray_structure.replace_sites_cart(
      self.sites_cart)
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
    #zz = maptbx.real_space_target_simple_2(
    #  unit_cell   = self.xray_structure.unit_cell(),
    #  map_target  = self.map_target,
    #  map_current = self.map_current,
    #  box_size    = 100,
    #  sites_frac  = self.xray_structure.sites_frac())
    #print r, zz
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
      tmp = r.eight_point_interpolation(site_frac)
      #x = (gx*r).eight_point_interpolation(site_frac)
      #y = (gy*r).eight_point_interpolation(site_frac)
      #z = (gz*r).eight_point_interpolation(site_frac)
      x = gx*tmp
      y = gy*tmp
      z = gz*tmp
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
