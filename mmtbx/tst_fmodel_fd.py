from cctbx.array_family import flex
import mmtbx.f_model
import mmtbx.targets
from cctbx.regression.tst_miller import generate_random_hl
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import adptbx
from libtbx.test_utils import approx_equal
import random
import sys
from cctbx import xray
from mmtbx import masks

random.seed(0)
flex.set_random_seed(0)

def finite_differences_site(target_functor, eps=1.e-5):
  fmodel = target_functor.manager
  structure = fmodel.xray_structure
  unit_cell = structure.unit_cell()
  gs = flex.double()
  for i_scatterer in xrange(structure.scatterers().size()):
    sc = structure.scatterers()[i_scatterer]
    site_orig = sc.site
    d_target_d_site = []
    for ix in xrange(3):
      ts = []
      for signed_eps in [eps, -eps]:
        site_cart = list(unit_cell.orthogonalize(site_orig))
        site_cart[ix] += signed_eps
        sc.site = unit_cell.fractionalize(site_cart)
        fmodel.update_xray_structure(update_f_calc=True)
        ts.append(target_functor().target_work())
      gs.append((ts[0]-ts[1])/(2*eps))
    sc.site = site_orig
  return gs

sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
sfg_params.algorithm = "direct"
sfg_params.cos_sin_table = False

def exercise(space_group_info,
             n_elements       = 10,
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
      scatterers = xray_structure.scatterers()
      if (anomalous_flag):
        assert scatterers.size() >= 7
        for i in [1,7]:
          scatterers[i].fp = -0.2
          scatterers[i].fdp = 5
        have_non_zero_fdp = True
      else:
        for i in [1,7]:
          scatterers[i].fp = 0
          scatterers[i].fdp = 0
        have_non_zero_fdp = False
      f_obs = abs(xray_structure.structure_factors(
                               d_min          = d_min,
                               anomalous_flag = anomalous_flag,
                               cos_sin_table  = sfg_params.cos_sin_table,
                               algorithm      = sfg_params.algorithm).f_calc())
      f_obs_comp = f_obs.structure_factors_from_scatterers(
                            xray_structure = xray_structure,
                            algorithm      = sfg_params.algorithm,
                            cos_sin_table  = sfg_params.cos_sin_table).f_calc()
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
            if (have_non_zero_fdp): continue # XXX gradients not correct!
            experimental_phases = generate_random_hl(miller_set=f_obs)
          else:
            experimental_phases = None
          if (target == "ml_sad"
                and (not anomalous_flag or mmtbx.targets.phaser is None)):
            continue
          print "  ",target
          xray.set_scatterer_grad_flags(
                                   scatterers = xrs.scatterers(),
                                   site       = True)
          fmodel = mmtbx.f_model.manager(
            xray_structure               = xrs,
            f_obs                        = f_obs,
            r_free_flags                 = flags,
            target_name                  = target,
            abcd                         = experimental_phases,
            sf_and_grads_accuracy_params = sfg_params,
            k_sol                        = k_sol,
            b_sol                        = b_sol,
            b_cart                       = b_cart,
            mask_params                  = masks.mask_master_params.extract())
          fmodel.update_xray_structure(
            xray_structure=xrs,
            update_f_calc=True,
            update_f_mask=True)
          xray.set_scatterer_grad_flags(
            scatterers=fmodel.xray_structure.scatterers(),
            site=True)
          fmodel.update()
          t_f = fmodel.target_functor()
          t_f.prepare_for_minimization()
          gs = t_f(compute_gradients=True).d_target_d_site_cart().as_double()
          gfd = finite_differences_site(target_functor=t_f)
          cc = flex.linear_correlation(gs, gfd).coefficient()
          if (0 or verbose):
            print "ana:", list(gs)
            print "fin:", list(gfd)
            print "rat:", [f/a for a,f in zip(gs,gfd)]
            print target, "corr:", cc, space_group_info
            print
          diff = gs - gfd
          diff /= max(1, flex.max(flex.abs(gfd)))
          tolerance = 1.2e-5
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
