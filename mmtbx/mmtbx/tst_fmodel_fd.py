from cctbx.array_family import flex
import mmtbx.f_model
from cctbx.regression.tst_miller import generate_random_hl
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

def finite_differences_site(fmodel, eps=1.e-5):
  structure = fmodel.xray_structure
  unit_cell = structure.unit_cell()
  if (fmodel.target_name != "ml_sad"):
    gs_old = flex.double()
  else:
    gs_old = None
  gs_new = flex.double()
  alpha, beta = fmodel.alpha_beta_w(only_if_required_by_target=True)
  target_functor = fmodel.target_functor()
  for i_scatterer in xrange(structure.scatterers().size()):
    sc = structure.scatterers()[i_scatterer]
    site_orig = sc.site
    d_target_d_site = []
    for ix in xrange(3):
      ts_old = []
      ts_new = []
      for signed_eps in [eps, -eps]:
        site_cart = list(unit_cell.orthogonalize(site_orig))
        site_cart[ix] += signed_eps
        sc.site = unit_cell.fractionalize(site_cart)
        fmodel.update_xray_structure(update_f_calc=True)
        if (gs_old is not None):
          ts_old.append(fmodel.target_w(alpha=alpha, beta=beta))
        ts_new.append(target_functor().target_work())
      if (gs_old is not None):
        gs_old.append((ts_old[0]-ts_old[1])/(2*eps))
      gs_new.append((ts_new[0]-ts_new[1])/(2*eps))
    sc.site = site_orig
  return gs_old, gs_new

def exercise(space_group_info,
             n_elements       = 10,
             sf_cos_sin_table = False,
             sf_algorithm     = "direct",
             table            = "wk1995",
             d_min            = 2.0,
             k_sol            = 0.35,
             b_sol            = 45.0,
             b_cart           = None,
             quick=False,
             verbose=0):
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
      for target in mmtbx.f_model.manager.target_names:
          if (quick):
            if (target not in ["ls_wunit_k1", "ml", "mlhl", "ml_sad"]):
              continue
          if (target == "mlhl"):
            experimental_phases = generate_random_hl(miller_set=f_obs)
          else:
            experimental_phases = None
          if (target == "ml_sad"
                and (not anomalous_flag or mmtbx.f_model.phaser is None)):
            continue
          print "  ",target
          xray.set_scatterer_grad_flags(
                                   scatterers = xrs.scatterers(),
                                   site       = True)
          fmodel = mmtbx.f_model.manager(
                                       xray_structure    = xrs,
                                       f_obs             = f_obs,
                                       r_free_flags      = flags,
                                       target_name       = target,
                                       abcd=experimental_phases,
                                       sf_cos_sin_table  = sf_cos_sin_table,
                                       sf_algorithm      = sf_algorithm,
                                       k_sol             = k_sol,
                                       b_sol             = b_sol,
                                       b_cart            = b_cart,
                                       mask_params       = mask_params())
          fmodel.update_xray_structure(
            xray_structure=xrs,
            update_f_calc=True,
            update_f_mask=True)
          if ((0 or verbose) and fmodel.target_name != "ml_sad"):
            fmodel.show_essential()
            print  f_obs.data().size()
          if ((0 or verbose) and fmodel.target_name != "ml_sad"):
            fmodel.show_comprehensive()
            print  f_obs.data().size()
          xray.set_scatterer_grad_flags(
            scatterers=fmodel.xray_structure.scatterers(),
            site=True)
          if (fmodel.target_name != "ml_sad"):
            gs_old = fmodel.gradient_wrt_atomic_parameters(site=True).packed()
          else:
            gs_old = None
          gs_new = fmodel.target_functor()(
            compute_gradients=True).d_target_d_site_cart().as_double()
          gfd_old, gfd_new = finite_differences_site(fmodel=fmodel)
          for which,gs,gfd in [("old", gs_old, gfd_old),
                               ("new", gs_new, gfd_new)]:
            if (gs is None): continue
            cc = flex.linear_correlation(gs, gfd).coefficient()
            if (0 or verbose):
              print which
              print "ana:", list(gs)
              print "fin:", list(gfd)
              print target, "corr:", cc
              print
            diff = gs - gfd
            tolerance = 1.e-6
            assert approx_equal(abs(flex.min(diff) ), 0.0, tolerance)
            assert approx_equal(abs(flex.mean(diff)), 0.0, tolerance)
            assert approx_equal(abs(flex.max(diff) ), 0.0, tolerance)
            assert approx_equal(cc, 1.0, tolerance)
          fmodel.model_error_ml()

def run_call_back(flags, space_group_info):
  exercise(
    space_group_info=space_group_info,
    verbose=flags.Verbose,
    quick=flags.quick)

def run():
  debug_utils.parse_options_loop_space_groups(
    argv=sys.argv[1:], call_back=run_call_back, keywords=("quick",))

if (__name__ == "__main__"):
  run()
