from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args

class target_and_gradients(object):

  def __init__(self,
               unit_cell,
               density_map,
               sites_cart,
               real_space_gradients_delta):
    adopt_init_args(self, locals())

  def target(self):
    return -1.*maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.density_map,
      sites_cart  = self.sites_cart)

  def gradients(self):
    return -1.*maptbx.real_space_gradients_simple(
      unit_cell   =self.unit_cell,
      density_map =self.density_map,
      sites_cart  =self.sites_cart,
      delta       =self.real_space_gradients_delta)

def local_standard_deviations_target_per_site(
      unit_cell, density_map, weight_map, weight_map_scale_factor,
      sites_cart, site_radii):
  if (weight_map is None):
    return maptbx.standard_deviations_around_sites(
      unit_cell=unit_cell,
      density_map=density_map,
      sites_cart=sites_cart,
      site_radii=site_radii)
  d = maptbx.real_space_target_simple_per_site(
    unit_cell=unit_cell,
    density_map=density_map,
    sites_cart=sites_cart)
  w = flex.pow2(maptbx.standard_deviations_around_sites(
    unit_cell=unit_cell,
    density_map=weight_map,
    sites_cart=sites_cart,
    site_radii=site_radii))
  w_min = 0.01
  w.set_selected((w < w_min), w_min)
  w = 1. / w
  if (weight_map_scale_factor is not None):
    assert weight_map_scale_factor > 0
    w *= weight_map_scale_factor
  return d / w

def local_standard_deviations_target(
      unit_cell, density_map, weight_map, weight_map_scale_factor,
      sites_cart, site_radii):
  return flex.sum(local_standard_deviations_target_per_site(
    unit_cell, density_map, weight_map, weight_map_scale_factor,
    sites_cart, site_radii))

def local_standard_deviations_gradients(
      unit_cell, density_map, weight_map, weight_map_scale_factor,
      sites_cart, site_radii, delta):
  # XXX inefficient implementation:
  # inner-most loop (in maptbx.standard_deviations_around_sites)
  # should be outer-most
  grad_cols = []
  for i_dim in [0,1,2]:
    shift = [0,0,0]
    targets = []
    for signed_delta in [delta, -delta]:
      shift[i_dim] = signed_delta
      sites_cart_shifted = sites_cart + shift
      targets.append(local_standard_deviations_target_per_site(
        unit_cell=unit_cell,
        density_map=density_map,
        weight_map=weight_map,
        weight_map_scale_factor=weight_map_scale_factor,
        sites_cart=sites_cart_shifted,
        site_radii=site_radii))
    grad_cols.append((targets[0]-targets[1])/(2*delta))
  return flex.vec3_double(*grad_cols)

class lbfgs(object):

  def __init__(O,
        sites_cart,
        density_map,
        weight_map=None,
        unit_cell=None,
        selection_variable=None,
        geometry_restraints_manager=None,
        energies_sites_flags=None,
        real_space_target_weight=1,
        real_space_gradients_delta=None,
        local_standard_deviations_radius=None,
        weight_map_scale_factor=None,
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None):
    assert [unit_cell, geometry_restraints_manager].count(None) == 1
    assert real_space_gradients_delta is not None
    if (unit_cell is None):
      unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
    O.density_map = density_map
    O.weight_map = weight_map
    O.unit_cell = unit_cell
    O.sites_cart = sites_cart
    O.geometry_restraints_manager = geometry_restraints_manager
    O.energies_sites_flags = energies_sites_flags
    O.real_space_target_weight = real_space_target_weight
    O.real_space_gradients_delta = real_space_gradients_delta
    O.local_standard_deviations_radius = local_standard_deviations_radius
    if (O.local_standard_deviations_radius is None):
      O.site_radii = None
    else:
      O.site_radii = flex.double(
        O.sites_cart.size(), O.local_standard_deviations_radius)
    O.weight_map_scale_factor = weight_map_scale_factor
    O.selection_variable = selection_variable
    if (O.selection_variable is None):
      O.sites_cart = sites_cart
      O.x = sites_cart.as_double()
    else:
      O.sites_cart = sites_cart.deep_copy()
      O.x = sites_cart.select(O.selection_variable).as_double()
    O.number_of_function_evaluations = -1
    O.f_start, O.g_start = O.compute_functional_and_gradients()
    O.minimizer = scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=lbfgs_termination_params,
      exception_handling_params=lbfgs_exception_handling_params)
    O.f_final, O.g_final = O.compute_functional_and_gradients()
    del O.x
    del O.site_radii

  def compute_functional_and_gradients(O):
    if (O.number_of_function_evaluations == 0):
      O.number_of_function_evaluations += 1
      return O.f_start, O.g_start
    O.number_of_function_evaluations += 1
    O.sites_cart_variable = flex.vec3_double(O.x)
    if (O.real_space_target_weight == 0):
      rs_f = 0.
      rs_g = flex.vec3_double(O.sites_cart_variable.size(), (0,0,0))
    else:
      if (O.local_standard_deviations_radius is None):
        rs_f = maptbx.real_space_target_simple(
          unit_cell   = O.unit_cell,
          density_map = O.density_map,
          sites_cart  = O.sites_cart_variable)
        rs_g = maptbx.real_space_gradients_simple(
          unit_cell   = O.unit_cell,
          density_map = O.density_map,
          sites_cart  = O.sites_cart_variable,
          delta       = O.real_space_gradients_delta)
      else:
        rs_f = local_standard_deviations_target(
          unit_cell=O.unit_cell,
          density_map=O.density_map,
          weight_map=O.weight_map,
          weight_map_scale_factor=O.weight_map_scale_factor,
          sites_cart=O.sites_cart_variable,
          site_radii=O.site_radii)
        rs_g = local_standard_deviations_gradients(
          unit_cell=O.unit_cell,
          density_map=O.density_map,
          weight_map=O.weight_map,
          weight_map_scale_factor=O.weight_map_scale_factor,
          sites_cart=O.sites_cart_variable,
          site_radii=O.site_radii,
          delta=O.real_space_gradients_delta)
      rs_f *= -O.real_space_target_weight
      rs_g *= -O.real_space_target_weight
    if (O.geometry_restraints_manager is None):
      f = rs_f
      g = rs_g
    else:
      if (O.selection_variable is None):
        O.sites_cart = O.sites_cart_variable
      else:
        O.sites_cart.set_selected(O.selection_variable, O.sites_cart_variable)
      gr_e = O.geometry_restraints_manager.energies_sites(
        sites_cart=O.sites_cart,
        flags=O.energies_sites_flags,
        compute_gradients=True)
      gr_e_gradients = gr_e.gradients
      if (O.selection_variable is not None):
        gr_e_gradients = gr_e.gradients.select(O.selection_variable)
      f = rs_f + gr_e.target
      g = rs_g + gr_e_gradients
    return f, g.as_double()
