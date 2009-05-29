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
    O.minimizer = scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=lbfgs_termination_params,
      exception_handling_params=lbfgs_exception_handling_params)
    O.f, O.g = O.compute_functional_and_gradients()
    O.sites_cart = flex.vec3_double(O.x)
    del O.x

  def compute_functional_and_gradients(O):
    sites_cart = flex.vec3_double(O.x)
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
    if (O.gr is None):
      f = rs_f
      g = rs_g
      print "f:", rs_f
    else:
      gr_e = O.gr.energies_sites(sites_cart=sites_cart, compute_gradients=True)
      f = O.rs_weight * rs_f + gr_e.target
      g = O.rs_weight * rs_g + gr_e.gradients
      print "f:", O.rs_weight * rs_f, gr_e.target
    return f, g.as_double()
