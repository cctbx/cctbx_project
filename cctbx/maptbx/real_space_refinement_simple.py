from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.lbfgs

class lbfgs(object):

  def __init__(O,
        geometry_restraints_manager,
        sites_cart,
        density_map,
        gradients_delta,
        real_space_weight_scale,
        lbfgs_termination_params=None,
        lbfgs_exception_handling_params=None):
    O.gr = geometry_restraints_manager
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
    unit_cell = O.gr.crystal_symmetry.unit_cell()
    gr_e = O.gr.energies_sites(sites_cart=sites_cart, compute_gradients=True)
    gr_f = gr_e.target
    gr_g = gr_e.gradients
    rs_f = -maptbx.real_space_target_simple(
      unit_cell=unit_cell,
      density_map=O.density_map,
      sites_cart=sites_cart)
    rs_g = maptbx.real_space_gradients_simple(
      unit_cell=unit_cell,
      density_map=O.density_map,
      sites_cart=sites_cart,
      delta=O.gradients_delta)
    rs_g *= -1
    if (O.rs_weight is None):
      rms_gr_g = flex.mean_sq(gr_g.as_double())**0.5
      rms_rs_g = flex.mean_sq(rs_g.as_double())**0.5
      if (rms_rs_g < 1):
        O.rs_weight = 1
      else:
        O.rs_weight = rms_gr_g / rms_rs_g * O.rs_weight_scale
    f = gr_f + rs_f * O.rs_weight
    g = gr_g + rs_g * O.rs_weight
    return f, g.as_double()
