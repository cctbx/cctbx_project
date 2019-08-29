from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.lbfgs
from cctbx import xray
from scitbx import matrix

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

class magnification_anisotropic_minimization(object):

  def __init__(O,
        sites_cart,
        density_map,
        unit_cell,
        K, # magnification 3*3  matrix, or a triplet (diagonal only)
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None):
    assert [isinstance(K, matrix.sqr), isinstance(K, matrix.col)].count(True)==1
    O.density_map = density_map
    O.unit_cell = unit_cell
    O.sites_cart = sites_cart
    O.K = K
    O.x = flex.double(O.K.elems)
    O.number_of_function_evaluations = -1
    O.f_start, O.g_start = O.compute_functional_and_gradients()
    O.minimizer = scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=lbfgs_termination_params,
      exception_handling_params=lbfgs_exception_handling_params)
    O.f_final, O.g_final = O.compute_functional_and_gradients()
    del O.x

  def compute_functional_and_gradients(O):
    if (O.number_of_function_evaluations == 0):
      O.number_of_function_evaluations += 1
      return O.f_start, O.g_start
    O.number_of_function_evaluations += 1
    if(isinstance(O.K, matrix.sqr)):
      O.K = matrix.sqr(
        [O.x[0], O.x[1], O.x[2],
         O.x[3], O.x[4], O.x[5],
         O.x[6], O.x[7], O.x[8]])
    else:
      O.K = matrix.col([O.x[0], O.x[1], O.x[2]])
    o = maptbx.target_and_gradients_simple_magnification(
      unit_cell  = O.unit_cell,
      map_target = O.density_map,
      sites_cart = O.sites_cart,
      K          = O.K)
    return -1.* o.target(), o.gradients()*(-1.)

class lbfgs(object):

  def __init__(O,
        sites_cart,
        density_map,
        gradients_method,
        weight_map=None,
        unit_cell=None,
        selection_variable=None,
        selection_variable_real_space=None,
        geometry_restraints_manager=None,
        energies_sites_flags=None,
        gradient_only=False,
        line_search=True,
        real_space_target_weight=1,
        real_space_gradients_delta=None,
        local_standard_deviations_radius=None,
        weight_map_scale_factor=None,
        lbfgs_core_params=None,
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None,
        states_collector=None):
    #assert [unit_cell, geometry_restraints_manager].count(None) == 1
    assert real_space_gradients_delta is not None
    if (unit_cell is None):
      unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
    if(selection_variable_real_space is not None):
      assert selection_variable_real_space.size() == sites_cart.size()
    else:
      selection_variable_real_space = flex.bool(sites_cart.size(), True)
    O.lbfgs_core_params = lbfgs_core_params
    O.gradients_method = gradients_method
    O.x_previous = None
    O.states_collector = states_collector
    O.gradient_only=gradient_only
    O.line_search=line_search
    O.density_map = density_map
    O.weight_map = weight_map
    O.unit_cell = unit_cell
    O.sites_cart = sites_cart
    O.geometry_restraints_manager = geometry_restraints_manager
    O.energies_sites_flags = energies_sites_flags
    O.real_space_target_weight = real_space_target_weight
    O.real_space_gradients_delta = real_space_gradients_delta
    O.local_standard_deviations_radius = local_standard_deviations_radius
    O.selection_variable_real_space = selection_variable_real_space
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
      gradient_only=O.gradient_only,
      line_search=O.line_search,
      core_params=O.lbfgs_core_params,
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
    #
    x_current = O.x
    if(O.x_previous is None):
      O.x_previous = x_current.deep_copy()
    else:
      xray.ext.damp_shifts(previous=O.x_previous, current=x_current,
        max_value=10.)
      O.x_previous = x_current.deep_copy()
    #
    O.sites_cart_variable = flex.vec3_double(x_current)
    if (O.real_space_target_weight == 0):
      rs_f = 0.
      rs_g = flex.vec3_double(O.sites_cart_variable.size(), (0,0,0))
    else:
      if (O.local_standard_deviations_radius is None):
        if(O.gradients_method=="fd"):
          o = maptbx.target_and_gradients_simple(
            unit_cell   = O.unit_cell,
            map_target  = O.density_map,
            sites_cart  = O.sites_cart_variable,
            delta       = O.real_space_gradients_delta,
            selection   = O.selection_variable_real_space)
        else:
          o = maptbx.target_and_gradients_simple(
            unit_cell     = O.unit_cell,
            map_target    = O.density_map,
            sites_cart    = O.sites_cart_variable,
            selection     = O.selection_variable_real_space,
            interpolation = O.gradients_method)
        rs_f = o.target()
        rs_g = o.gradients()
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
        if(O.states_collector is not None):
          O.states_collector.add(sites_cart = O.sites_cart)
      gr_e = O.geometry_restraints_manager.energies_sites(
        sites_cart=O.sites_cart,
        compute_gradients=True)
      gr_e_gradients = gr_e.gradients
      if (O.selection_variable is not None):
        gr_e_gradients = gr_e.gradients.select(O.selection_variable)
      f = rs_f + gr_e.target
      g = rs_g + gr_e_gradients
    return f, g.as_double()
