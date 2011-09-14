from cctbx.array_family import flex
import mmtbx.f_model
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import sgtbx
import random
import sys
from mmtbx.scaling import outlier_rejection
from cctbx.xray import observation_types
from cctbx.development import debug_utils

if (1):
  random.seed(0)
  flex.set_random_seed(value=0)

def exercise(d_min            = 3.5,
             k_sol            = 0.3,
             b_sol            = 60.0,
             b_cart           = [0,0,0,0,0,0],
             anomalous_flag   = False,
             scattering_table = "it1992",
             space_group_info = None):
  space_groups = [ str(space_group_info) ]
  for sg in space_groups:
      ### get random structure
      xray_structure = random_structure.xray_structure(
                          space_group_info       = sgtbx.space_group_info(sg),
                          elements               = (("O","C","N")*200),
                          volume_per_atom        = 100,
                          min_distance           = 1.5,
                          general_positions_only = True,
                          random_u_iso           = True)
      xray_structure.scattering_type_registry(table = scattering_table)
      ### Get FOBS
      for scale in [0.0001, 1.0, 1000.0]:
          dummy = abs(xray_structure.structure_factors(
                                   d_min          = d_min,
                                   anomalous_flag = anomalous_flag).f_calc())
          flags = dummy.generate_r_free_flags(fraction = 0.1,
                                              max_free = 99999999)
          fmodel = mmtbx.f_model.manager(xray_structure   = xray_structure,
                                         r_free_flags     = flags,
                                         target_name      = "ls_wunit_k1",
                                         f_obs            = dummy,
                                         k_sol            = k_sol,
                                         b_sol            = b_sol,
                                         b_cart           = b_cart)

          fmodel.update_xray_structure(xray_structure = xray_structure,
                                       update_f_calc = True,
                                       update_f_mask = True)
          f_obs = abs(fmodel.f_model())
          f_obs = f_obs.array(data = f_obs.data()*scale)
          f_obs.set_observation_type(observation_type = observation_types.amplitude())
          ### look at non-model based outliers detection
          om = outlier_rejection.outlier_manager(miller_obs   = f_obs,
                                                 r_free_flags = flags,
                                                 out          = "silent")
          tmp1 = om.basic_wilson_outliers()
          tmp2 = om.extreme_wilson_outliers()
          tmp3 = om.beamstop_shadow_outliers()

          # start loop over distorted models
          for error in [0.0,  0.8]:
              for fraction in [0.0,0.5]:
                # get distorted model
                xrs_dc = xray_structure.deep_copy_scatterers()
                sel = xrs_dc.random_remove_sites_selection(fraction = fraction)
                xrs_dc = xrs_dc.select(sel)
                xrs_dc.shake_sites_in_place(rms_difference=error)
                xrs_dc.scattering_type_registry(table = scattering_table)

                fmodel = mmtbx.f_model.manager(xray_structure    = xrs_dc,
                                               r_free_flags      = flags,
                                               target_name       = "ls_wunit_k1",
                                               f_obs             = f_obs,
                                               b_cart            = b_cart)
                for k_sol in [0.10, 0.30, 0.50]:
                  for b_sol in [60, 80]:
                     fmodel.update(k_sols = [k_sol], b_sol = b_sol)
                     a,b = fmodel.alpha_beta()
                     o_sel =  om.model_based_outliers(f_model = fmodel.f_model())
                     n_out = o_sel.data().count(False)
                     assert (n_out < 50)

def run_call_back(flags, space_group_info):
  exercise(space_group_info=space_group_info)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)


if (__name__ == "__main__"):
  run()
