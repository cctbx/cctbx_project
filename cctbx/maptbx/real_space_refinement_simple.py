from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.lbfgs

class lbfgs(object):

  def __init__(O,
        sites_cart,
        density_map,
        gradients_delta,
        unit_cell=None,
        geometry_restraints_manager=None,
        real_space_weight=None,
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None):
    assert [unit_cell, geometry_restraints_manager].count(None) == 1
    if (geometry_restraints_manager is None):
      assert real_space_weight is None
    else:
      assert real_space_weight is not None
      unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
    O.density_map = density_map
    O.gradients_delta = gradients_delta
    O.unit_cell = unit_cell
    O.gr = geometry_restraints_manager
    O.rs_weight = real_space_weight
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
      delta=O.gradients_delta)
    if (O.gr is None):
      rs_f *= -1
      rs_g *= -1
      f = rs_f
      g = rs_g
    else:
      rs_f *= -O.rs_weight
      rs_g *= -O.rs_weight
      gr_e = O.gr.energies_sites(sites_cart=sites_cart, compute_gradients=True)
      f = rs_f + gr_e.target
      g = rs_g + gr_e.gradients
    return f, g.as_double()
