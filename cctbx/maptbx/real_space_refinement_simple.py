from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.lbfgs

class lbfgs(object):

  def __init__(O,
        sites_cart,
        density_map,
        unit_cell=None,
        geometry_restraints_manager=None,
        real_space_target_weight=1,
        real_space_gradients_delta=None,
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None):
    assert [unit_cell, geometry_restraints_manager].count(None) == 1
    assert real_space_gradients_delta is not None
    if(unit_cell is None and geometry_restraints_manager is not None):
      unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
    O.density_map = density_map
    O.real_space_gradients_delta = real_space_gradients_delta
    O.unit_cell = unit_cell
    O.geometry_restraints_manager = geometry_restraints_manager
    O.real_space_target_weight = real_space_target_weight
    O.x = sites_cart.as_double()
    O.number_of_function_evaluations = -1
    O.f_start, O.g_start = O.compute_functional_and_gradients()
    O.minimizer = scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=lbfgs_termination_params,
      exception_handling_params=lbfgs_exception_handling_params)
    O.f_final, O.g_final = O.compute_functional_and_gradients()
    O.sites_cart = flex.vec3_double(O.x)
    del O.x

  def compute_functional_and_gradients(O):
    if (O.number_of_function_evaluations == 0):
      O.number_of_function_evaluations += 1
      return O.f_start, O.g_start
    O.number_of_function_evaluations += 1
    sites_cart = flex.vec3_double(O.x)
    rs_f = maptbx.real_space_target_simple(
      unit_cell=O.unit_cell,
      density_map=O.density_map,
      sites_cart=sites_cart)
    rs_g = maptbx.real_space_gradients_simple(
      unit_cell=O.unit_cell,
      density_map=O.density_map,
      sites_cart=sites_cart,
      delta=O.real_space_gradients_delta)
    rs_f *= -O.real_space_target_weight
    rs_g *= -O.real_space_target_weight
    if (O.geometry_restraints_manager is None):
      f = rs_f
      g = rs_g
    else:
      gr_e = O.geometry_restraints_manager.energies_sites(
        sites_cart=sites_cart, compute_gradients=True)
      f = rs_f + gr_e.target
      g = rs_g + gr_e.gradients
    return f, g.as_double()
