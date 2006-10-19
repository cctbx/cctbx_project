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
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx.scaling import outlier_rejection
from cctbx.xray import observation_types
from cctbx.development import debug_utils


def exercise(d_min            = 3.5,
             k_sol            = 0.3,
             b_sol            = 60.0,
             b_cart           = [0,0,0,0,0,0],
             sf_algorithm     = "fft",
             sf_cos_sin_table = False,
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
                          random_u_iso           = True,
                          #u_iso                  = adptbx.b_as_u(30.0)
                          )
      xray_structure.scattering_type_registry(table = scattering_table)
      ### Get FOBS
      for scale in [0.0001, 1.0, 1000.0]:
          dummy = abs(xray_structure.structure_factors(
                                   d_min          = d_min,
                                   anomalous_flag = anomalous_flag,
                                   cos_sin_table  = sf_cos_sin_table,
                                   algorithm      = sf_algorithm).f_calc())
          flags = dummy.generate_r_free_flags(fraction = 0.1,
                                              max_free = 99999999)
          fmodel = mmtbx.f_model.manager(xray_structure   = xray_structure,
                                         sf_algorithm     = sf_algorithm,
                                         sf_cos_sin_table = sf_cos_sin_table,
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
          #print "basic_wilson_outliers    = ", tmp1.data().count(True), tmp1.data().count(False)
          #print "extreme_wilson_outliers  = ", tmp2.data().count(True), tmp2.data().count(False)
          #print "beamstop_shadow_outliers = ", tmp3.data().count(True), tmp3.data().count(False)

          # start loop over distorted models
          for error in [0.0,  0.8]:
              for fraction in [0.0,0.5]:
                # get distorted model
                xrs_dc = xray_structure.deep_copy_scatterers()
                sel = xrs_dc.random_remove_sites_selection(fraction = fraction)
                xrs_dc = xrs_dc.select(sel)
                xrs_dc.shake_sites(mean_error = error)
                xrs_dc.scattering_type_registry(table = scattering_table)

                fmodel = mmtbx.f_model.manager(xray_structure    = xrs_dc,
                                               sf_algorithm      = sf_algorithm,
                                               sf_cos_sin_table  = sf_cos_sin_table,
                                               r_free_flags      = flags,
                                               target_name       = "ls_wunit_k1",
                                               f_obs             = f_obs,
                                               b_cart            = b_cart)
                for k_sol in [0.10, 0.30, 0.50]:
                  for b_sol in [60, 80]:
                     fmodel.update(k_sol = k_sol, b_sol = b_sol)
                     #print "   scale = %12.6f mean_error = %3.1f deleted = %3.1f r_work = %6.4f "% \
                     #   (scale, error, fraction, fmodel.r_work()), k_sol, b_sol
                     a,b = fmodel.alpha_beta()
                     o_sel =  om.model_based_outliers(fmodel_manager = fmodel)
                     n_out = o_sel.data().count(False)
                     assert (n_out < 50) # some would be outliers are always expected, but not too much though
                     #print "model_based_outliers = ", o_sel.data().count(True), o_sel.data().count(False)


def run_call_back(flags, space_group_info):
  exercise(space_group_info=space_group_info)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)


if (__name__ == "__main__"):
  run()
