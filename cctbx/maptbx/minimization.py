from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import scitbx.lbfgs
from cctbx import miller
from cctbx import maptbx
from libtbx import adopt_init_args
from cctbx.maptbx import real_space_refinement_simple # special import

def show(histogram):
  h_1 = histogram
  lc_1 = histogram.data_min()
  s_1 = enumerate(histogram.slots())
  for (i_1,n_1) in s_1:
    hc_1 = h_1.data_min() + h_1.slot_width() * (i_1+1)
    print("%8.3f - %8.3f: %5d" % (lc_1,hc_1,n_1))
    lc_1 = hc_1

class target_and_gradients(object):

  def __init__(self,
               unit_cell,
               map_target,
               real_space_gradients_delta,
               sites_frac,
               map_current = None,
               geometry_restraints_manager = None,
               real_space_target_weight = None,
               restraints_target_weight = None,
               sites_cart = None,
               target_type = "diffmap"):
    adopt_init_args(self, locals())
    assert target_type in ["simple", "diffmap"]
    if(geometry_restraints_manager is not None):
      assert [real_space_target_weight, restraints_target_weight, sites_cart
             ].count(None)==0
      self.geom_tg_obj = geometry_restraints_manager.energies_sites(
        sites_cart = sites_cart, compute_gradients = True)
    if(target_type == "diffmap"):
      assert map_current is not None
      self.rsr_tg_obj = maptbx.target_and_gradients_diffmap(
        unit_cell   = unit_cell,
        map_target  = map_target,
        map_current = map_current,
        step        = real_space_gradients_delta,
        sites_frac  = sites_frac)
    if(target_type == "simple"):
      assert sites_cart is not None
      self.rsr_tg_obj = maptbx.target_and_gradients_simple(
        unit_cell   = unit_cell,
        map_target  = map_target,
        sites_cart  = sites_cart,
        delta       = real_space_gradients_delta,
        selection   = flex.bool(sites_cart.size(),True))

  def target(self):
    rs_f = self.rsr_tg_obj.target()
    if(self.target_type=="simple"): rs_f *= -1.
    if(self.geometry_restraints_manager):
      result = self.real_space_target_weight * rs_f + \
               self.restraints_target_weight * self.geom_tg_obj.target
    else:
      result = rs_f
    return result

  def gradients(self):
    rs_g = self.rsr_tg_obj.gradients()
    if(self.target_type=="simple"): rs_g *= -1.
    if(self.geometry_restraints_manager):
      g = self.real_space_target_weight * rs_g + \
          self.restraints_target_weight * self.geom_tg_obj.gradients
    else:
      g = rs_g
    return g.as_double()

  def weight(self):
    gx = self.rsr_tg_obj.gradients()
    gc = self.geom_tg_obj.gradients
    tmp = flex.sqrt(gc.dot())
    sel = tmp < 3*flex.mean(tmp)
    gc = gc.select(sel)
    result = gx.norm()/gc.norm()
    return result


class run(object):
  def __init__(self,
        xray_structure,
        miller_array,
        crystal_gridding,
        map_target,
        step,
        target_type,
        max_iterations=30,
        min_iterations=25,
        geometry_restraints_manager=None,
        real_space_target_weight=1,
        restraints_target_weight=1):
    assert target_type in ["diffmap", "simple"]
    self.step = step
    self.geometry_restraints_manager = geometry_restraints_manager
    self.real_space_target_weight = real_space_target_weight
    self.sites_cart = xray_structure.sites_cart()
    self.x = self.sites_cart.as_double()
    self.xray_structure = xray_structure
    self.map_target = map_target
    self.miller_array = miller_array
    self.crystal_gridding = crystal_gridding
    self.restraints_target_weight = restraints_target_weight
    self.target_type = target_type
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
    self.map_current = None
    if(self.target_type == "diffmap"):
      self.map_current = self.compute_map()
    tg_obj = target_and_gradients(
      unit_cell                   = self.xray_structure.unit_cell(),
      map_target                  = self.map_target,
      map_current                 = self.map_current,
      real_space_gradients_delta  = self.step,
      sites_frac                  = self.xray_structure.sites_frac(),
      geometry_restraints_manager = self.geometry_restraints_manager,
      real_space_target_weight    = self.real_space_target_weight,
      restraints_target_weight    = self.restraints_target_weight,
      sites_cart                  = self.sites_cart,
      target_type                 = self.target_type)
    return tg_obj.target(), tg_obj.gradients()

  def compute_map(self):
    map_coefficients = self.miller_array.structure_factors_from_scatterers(
      xray_structure = self.xray_structure).f_calc()
    fft_map = miller.fft_map(crystal_gridding     = self.crystal_gridding,
                             fourier_coefficients = map_coefficients)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result
