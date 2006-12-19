from cctbx.array_family import flex
import mmtbx.f_model
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import sgtbx
from cctbx import adptbx
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import random
import sys, math
from cctbx import xray

random.seed(0)
flex.set_random_seed(0)

class mask_params:
  def __init__(self):
    self.solvent_radius = 1
    self.shrink_truncation_radius = 1
    self.grid_step_factor = 4
    self.verbose = 1
    self.mean_shift_for_mask_update = 0.1

def finite_differences_site(cartesian_flag, fmodel, delta=0.00001):
  structure = fmodel.xray_structure
  unit_cell = structure.unit_cell()
  abc = unit_cell.parameters()[:3]
  derivatives = flex.vec3_double()
  for i_scatterer in xrange(structure.scatterers().size()):
    d_target_d_site = [0,0,0]
    for ix in xrange(3):
      target_values = []
      for d_sign in (-1, 1):
        modified_structure = structure.deep_copy_scatterers()
        ms = modified_structure.scatterers()[i_scatterer]
        site = list(ms.site)
        if (not cartesian_flag):
          site[ix] += d_sign * delta / abc[ix]
        else:
          site_cart = list(unit_cell.orthogonalize(site))
          site_cart[ix] += d_sign * delta
          site = unit_cell.fractionalize(site_cart)
        ms.site = site
        fmodel.update_xray_structure(xray_structure = modified_structure,
                                     update_f_calc = True)
        target_values.append(fmodel.target_w())
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      if (not cartesian_flag): derivative *= abc[ix]
      d_target_d_site[ix] = derivative
    derivatives.append(d_target_d_site)
  return derivatives


def exercise(space_group_info,
             n_elements       = 10,
             sf_cos_sin_table = False,
             sf_algorithm     = "direct",
             table            = "wk1995",
             d_min            = 2.0,
             k_sol            = 0.35,
             b_sol            = 45.0,
             b_cart           = None):
  xray_structure = random_structure.xray_structure(
         space_group_info       = space_group_info,
         elements               =(("O","N","C")*(n_elements/3+1))[:n_elements],
         volume_per_atom        = 100,
         min_distance           = 1.5,
         general_positions_only = True,
         random_u_iso           = False,
         random_occupancy       = False)
  xray_structure.scattering_type_registry(table = table)
  sg = xray_structure.space_group()
  uc = xray_structure.unit_cell()
  u_cart_1 = adptbx.random_u_cart(u_scale=5, u_min=5)
  u_star_1 = adptbx.u_cart_as_u_star(uc, u_cart_1)
  b_cart   = adptbx.u_star_as_u_cart(uc, sg.average_u_star(u_star = u_star_1))
  for anomalous_flag in [False, True]:
      f_obs = abs(xray_structure.structure_factors(
                                       d_min          = d_min,
                                       anomalous_flag = anomalous_flag,
                                       cos_sin_table  = sf_cos_sin_table,
                                       algorithm      = sf_algorithm).f_calc())
      f_obs_comp = f_obs.structure_factors_from_scatterers(
                                    xray_structure = xray_structure,
                                    algorithm      = sf_algorithm,
                                    cos_sin_table  = sf_cos_sin_table).f_calc()
      f_obs = abs(f_obs_comp)
      flags = f_obs.generate_r_free_flags(fraction = 0.1,
                                          max_free = 99999999)
      #flags = flags.array(data = flex.bool(f_obs.data().size(), False))
      xrs = xray_structure.deep_copy_scatterers()
      xrs.shake_sites_in_place(rms_difference=0.3)
      for target in mmtbx.f_model.target_names:
          #XXX Must find out why ml-tolerance is so BIG.
          if(target == "ml"): tolerance = 1.5
          else: tolerance = 1.e-6
          if(target != "mlhl"):
             print "  ",target
             xray.set_scatterer_grad_flags(
                                      scatterers = xrs.scatterers(),
                                      site       = True)
             fmodel = mmtbx.f_model.manager(
                                          xray_structure    = xrs,
                                          f_obs             = f_obs,
                                          r_free_flags      = flags,
                                          target_name       = target,
                                          sf_cos_sin_table  = sf_cos_sin_table,
                                          sf_algorithm      = sf_algorithm,
                                          k_sol             = k_sol,
                                          b_sol             = b_sol,
                                          b_cart            = b_cart,
                                          mask_params       = mask_params())
             fmodel.update_xray_structure(xray_structure = xrs,
                                          update_f_calc = True,
                                          update_f_mask = True)
             if(0):
                fmodel.show_essential()
                print  f_obs.data().size()
             if(0):
                fmodel.show_comprehensive()
                print  f_obs.data().size()
             gs = flex.vec3_double(
                   fmodel.gradient_wrt_atomic_parameters(site = True).packed())
             xray.set_scatterer_grad_flags(
                               scatterers = fmodel.xray_structure.scatterers(),
                               site       = True)
             gfd = finite_differences_site(cartesian_flag = True,
                                           fmodel = fmodel)
             diff = (gs - gfd).as_double()
             assert approx_equal(abs(flex.min(diff) ), 0.0, tolerance)
             assert approx_equal(abs(flex.mean(diff)), 0.0, tolerance)
             assert approx_equal(abs(flex.max(diff) ), 0.0, tolerance)
             fmodel.model_error_ml()


def run_call_back(flags, space_group_info):
  exercise(space_group_info = space_group_info)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
