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
        geometry_restraints_manager=None,
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
    O.geometry_restraints_manager = geometry_restraints_manager
    O.real_space_gradients_delta = real_space_gradients_delta
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

class fragment_refine_restrained(object):

  def __init__(self,
        sites_cart_all,
        fragment_i_seqs,
        density_map,
        geometry_restraints_manager,
        real_space_target_weight,
        real_space_gradients_delta,
        lbfgs_termination_params,
        restraints_target_weight = 1.0):
    self.density_map = density_map
    self.geometry_restraints_manager = geometry_restraints_manager
    self.real_space_gradients_delta = real_space_gradients_delta
    self.real_space_target_weight = real_space_target_weight
    self.restraints_target_weight = restraints_target_weight
    #
    self.unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
    self.sites_cart_all = sites_cart_all
    self.fragment_i_seqs = fragment_i_seqs
    self.x = self.sites_cart_all.select(indices=self.fragment_i_seqs).as_double()
    #
    self.real_space_target = None
    self.number_of_function_evaluations = -1
    self.f_start, self.g_start = self.compute_functional_and_gradients()
    self.rs_f_start = self.rs_f
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params)
    self.f_final, self.g_final = self.compute_functional_and_gradients()
    self.rs_f_final = self.rs_f
    del self.rs_f
    del self.x
    del self.fragment_i_seqs
    #del self.sites_cart_all
    del self.unit_cell

  def compute_functional_and_gradients(self):
    if(self.number_of_function_evaluations == 0):
      self.number_of_function_evaluations += 1
      return self.f_start, self.g_start
    self.number_of_function_evaluations += 1
    self.sites_cart_fragment = flex.vec3_double(self.x)
    rs_f = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.density_map,
      sites_cart  = self.sites_cart_fragment)
    self.real_space_target = rs_f
    rs_g = maptbx.real_space_gradients_simple(
      unit_cell   = self.unit_cell,
      density_map = self.density_map,
      sites_cart  = self.sites_cart_fragment,
      delta       = self.real_space_gradients_delta)
    self.rs_f = rs_f
    rs_f *= -self.real_space_target_weight
    rs_g *= -self.real_space_target_weight
    if(self.geometry_restraints_manager is None or
       self.restraints_target_weight == 0):
      f = rs_f
      g = rs_g
    else:
      self.sites_cart_all.set_selected(self.fragment_i_seqs,
        self.sites_cart_fragment)
      gr_e = self.geometry_restraints_manager.energies_sites(
        sites_cart = self.sites_cart_all, compute_gradients=True)
      f = rs_f + gr_e.target * self.restraints_target_weight
      g = rs_g + gr_e.gradients.select(indices = self.fragment_i_seqs) * \
          self.restraints_target_weight
    return f, g.as_double()
