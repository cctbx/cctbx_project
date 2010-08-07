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

class lbfgs(object):

  def __init__(O,
        sites_cart,
        density_map,
        unit_cell=None,
        selection_variable=None,
        geometry_restraints_manager=None,
        energies_sites_flags=None,
        real_space_target_weight=1,
        real_space_gradients_delta=None,
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None):
    assert [unit_cell, geometry_restraints_manager].count(None) == 1
    assert real_space_gradients_delta is not None
    if (unit_cell is None):
      unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
    O.density_map = density_map
    O.unit_cell = unit_cell
    O.sites_cart = sites_cart
    O.geometry_restraints_manager = geometry_restraints_manager
    O.energies_sites_flags = energies_sites_flags
    O.real_space_gradients_delta = real_space_gradients_delta
    O.real_space_target_weight = real_space_target_weight
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
      rs_f = maptbx.real_space_target_simple(
        unit_cell   = O.unit_cell,
        density_map = O.density_map,
        sites_cart  = O.sites_cart_variable)
      rs_g = maptbx.real_space_gradients_simple(
        unit_cell   = O.unit_cell,
        density_map = O.density_map,
        sites_cart  = O.sites_cart_variable,
        delta       = O.real_space_gradients_delta)
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
