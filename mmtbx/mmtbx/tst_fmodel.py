from cctbx.array_family import flex
import mmtbx.f_model
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import sgtbx
from libtbx.test_utils import approx_equal
import random
import sys, math

random.seed(0)
flex.set_random_seed(0)

def test_1(xray_structure):
  # exercise almost all without dealing with particular values
  for d_min in [2.0, 2.5]:
      for sf_algorithm in ["direct", "fft"]:
          for anomalous_flag in [True, False]:
              for sf_cos_sin_table in [True, False]:
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
                  for (xrs,fc) in ((xray_structure,None),(None,f_obs_comp)):
                      ###
                      ### instantiate fmodel only
                      ###
                      fmodel = mmtbx.f_model.manager(
                                              xray_structure    = xrs,
                                              f_calc            = fc,
                                              f_obs             = f_obs,
                                              r_free_flags      = flags,
                                              target_name       = "ls_wunit_k1",
                                              sf_cos_sin_table  = sf_cos_sin_table,
                                              sf_algorithm      = sf_algorithm)
                      assert fmodel.f_obs.data().all_eq(f_obs.data())
                      assert abs(fmodel.f_calc).data().all_eq(f_obs.data())
                      assert abs(fmodel.f_model()).data().all_eq(f_obs.data())
                      assert fmodel.r_work() == 0.0
                      assert fmodel.r_free() == 0.0
                      assert fmodel.fu_aniso().all_eq(1.0)
                      assert fmodel.fu_aniso_w().all_eq(1.0)
                      assert fmodel.fu_aniso_t().all_eq(1.0)
                      assert abs(fmodel.target_w()) < 1.e-9
                      assert abs(fmodel.target_t()) < 1.e-9
                      assert fmodel.k_sol_b_sol() == (0.0,0.0)
                      assert fmodel.u_aniso == [0,0,0,0,0,0]
                      assert fmodel.u_iso() == 0.0
                      assert fmodel.f_obs_w().data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert fmodel.f_obs_t().data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_calc_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_calc_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_model_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_model_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_bulk_w()).data().all_eq(0)
                      assert abs(fmodel.f_bulk_t()).data().all_eq(0)
                      assert abs(fmodel.f_mask_w()).data().all_eq(0)
                      assert abs(fmodel.f_mask_t()).data().all_eq(0)
                      assert abs(fmodel.f_mask).data().all_eq(0)
                      assert approx_equal(fmodel.scale_k1()  , 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k1_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k1_t(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k2_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k2_t(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k3_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k3_t(), 1.0, 1.e-9)
                      assert fmodel.figures_of_merit().all_approx_equal(1.0, 1.e-6)
                      assert fmodel.phase_errors().all_approx_equal(0.0, 1.e-6)
                      assert fmodel.phase_errors_work().all_approx_equal(0.0, 1.e-6)
                      assert fmodel.phase_errors_test().all_approx_equal(0.0, 1.e-6)
                      a,b = fmodel.alpha_beta()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      a,b = fmodel.alpha_beta_w()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      a,b = fmodel.alpha_beta_t()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      assert fmodel.f_ordered_solvent.data().all_eq(0)
                      assert fmodel.f_ordered_solvent_w().data().all_eq(0)
                      assert fmodel.f_ordered_solvent_t().data().all_eq(0)
                      ###
                      ### instantiate fmodel only + update ksol & bsol
                      ###
                      fmodel = mmtbx.f_model.manager(
                                              xray_structure    = xrs,
                                              f_calc            = fc,
                                              f_obs             = f_obs,
                                              r_free_flags      = flags,
                                              target_name       = "ls_wunit_k1",
                                              sf_cos_sin_table  = sf_cos_sin_table,
                                              sf_algorithm      = sf_algorithm)
                      fmodel.update(k_sol = 0.5, b_sol = 35.0)
                      assert fmodel.f_obs.data().all_eq(f_obs.data())
                      assert abs(fmodel.f_calc).data().all_eq(f_obs.data())
                      assert abs(fmodel.f_model()).data().all_eq(f_obs.data())
                      assert fmodel.r_work() == 0.0
                      assert fmodel.r_free() == 0.0
                      assert fmodel.fu_aniso().all_eq(1.0)
                      assert fmodel.fu_aniso_w().all_eq(1.0)
                      assert fmodel.fu_aniso_t().all_eq(1.0)
                      assert abs(fmodel.target_w()) < 1.e-9
                      assert abs(fmodel.target_t()) < 1.e-9
                      assert fmodel.k_sol_b_sol() == (0.5,35.0)
                      assert fmodel.u_aniso == [0,0,0,0,0,0]
                      assert fmodel.u_iso() == 0.0
                      assert fmodel.f_obs_w().data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert fmodel.f_obs_t().data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_calc_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_calc_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_model_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_model_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_bulk_w()).data().all_eq(0)
                      assert abs(fmodel.f_bulk_t()).data().all_eq(0)
                      assert abs(fmodel.f_mask_w()).data().all_eq(0)
                      assert abs(fmodel.f_mask_t()).data().all_eq(0)
                      assert abs(fmodel.f_mask).data().all_eq(0)
                      assert approx_equal(fmodel.scale_k1()  , 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k1_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k1_t(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k2_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k2_t(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k3_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k3_t(), 1.0, 1.e-9)
                      assert fmodel.figures_of_merit().all_approx_equal(1.0, 1.e-6)
                      assert fmodel.phase_errors().all_approx_equal(0.0, 1.e-6)
                      assert fmodel.phase_errors_work().all_approx_equal(0.0, 1.e-6)
                      assert fmodel.phase_errors_test().all_approx_equal(0.0, 1.e-6)
                      a,b = fmodel.alpha_beta()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      a,b = fmodel.alpha_beta_w()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      a,b = fmodel.alpha_beta_t()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      assert fmodel.f_ordered_solvent.data().all_eq(0)
                      assert fmodel.f_ordered_solvent_w().data().all_eq(0)
                      assert fmodel.f_ordered_solvent_t().data().all_eq(0)
                      ###
                      ### instantiate fmodel only, then use deep_copy
                      ###
                      fmodel=None
                      fmodel_ = mmtbx.f_model.manager(
                                              xray_structure    = xrs,
                                              f_calc            = fc,
                                              f_obs             = f_obs,
                                              r_free_flags      = flags,
                                              target_name       = "ls_wunit_k1",
                                              sf_cos_sin_table  = sf_cos_sin_table,
                                              sf_algorithm      = sf_algorithm)
                      fmodel = fmodel_.deep_copy()
                      assert fmodel.f_obs.data().all_eq(f_obs.data())
                      assert abs(fmodel.f_calc).data().all_eq(f_obs.data())
                      assert abs(fmodel.f_model()).data().all_eq(f_obs.data())
                      assert fmodel.r_work() == 0.0
                      assert fmodel.r_free() == 0.0
                      assert fmodel.fu_aniso().all_eq(1.0)
                      assert fmodel.fu_aniso_w().all_eq(1.0)
                      assert fmodel.fu_aniso_t().all_eq(1.0)
                      assert abs(fmodel.target_w()) < 1.e-9
                      assert abs(fmodel.target_t()) < 1.e-9
                      assert fmodel.k_sol_b_sol() == (0.0,0.0)
                      assert fmodel.u_aniso == [0,0,0,0,0,0]
                      assert fmodel.u_iso() == 0.0
                      assert fmodel.f_obs_w().data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert fmodel.f_obs_t().data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_calc_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_calc_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_model_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_model_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_bulk_w()).data().all_eq(0)
                      assert abs(fmodel.f_bulk_t()).data().all_eq(0)
                      assert abs(fmodel.f_mask_w()).data().all_eq(0)
                      assert abs(fmodel.f_mask_t()).data().all_eq(0)
                      assert abs(fmodel.f_mask).data().all_eq(0)
                      assert approx_equal(fmodel.scale_k1()  , 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k1_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k1_t(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k2_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k2_t(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k3_w(), 1.0, 1.e-9)
                      assert approx_equal(fmodel.scale_k3_t(), 1.0, 1.e-9)
                      assert fmodel.figures_of_merit().all_approx_equal(1.0, 1.e-6)
                      assert fmodel.phase_errors().all_approx_equal(0.0, 1.e-6)
                      assert fmodel.phase_errors_work().all_approx_equal(0.0, 1.e-6)
                      assert fmodel.phase_errors_test().all_approx_equal(0.0, 1.e-6)
                      a,b = fmodel.alpha_beta()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      a,b = fmodel.alpha_beta_w()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      a,b = fmodel.alpha_beta_t()
                      assert a.data().all_approx_equal(1.0, 1.e-9)
                      assert b.data().all_approx_equal(0.0, 1.e-9)
                      assert fmodel.f_ordered_solvent.data().all_eq(0)
                      assert fmodel.f_ordered_solvent_w().data().all_eq(0)
                      assert fmodel.f_ordered_solvent_t().data().all_eq(0)
                      ###
                      ### instantiate fmodel only, then use resolution_filter
                      ###
                      fmodel=None
                      d_max_ = 3.0
                      d_min_ = 2.7
                      fmodel_ = mmtbx.f_model.manager(
                                              xray_structure    = xrs,
                                              f_calc            = fc,
                                              f_obs             = f_obs,
                                              r_free_flags      = flags,
                                              target_name       = "ls_wunit_k1",
                                              sf_cos_sin_table  = sf_cos_sin_table,
                                              sf_algorithm      = sf_algorithm)
                      fmodel_1 = fmodel_.resolution_filter(d_max = d_max_, d_min = d_min_)
                      if(fc is not None):
                         fc_ = fc.resolution_filter(d_max = d_max_, d_min = d_min_)
                      else:
                         fc_ = None
                      fmodel_2 = mmtbx.f_model.manager(
                                xray_structure    = xrs,
                                f_calc            = fc_,
                                f_obs             = f_obs.resolution_filter(d_max = d_max_, d_min = d_min_),
                                r_free_flags      = flags.resolution_filter(d_max = d_max_, d_min = d_min_),
                                target_name       = "ls_wunit_k1",
                                sf_cos_sin_table  = sf_cos_sin_table,
                                sf_algorithm      = sf_algorithm)
                      assert fmodel_1.f_obs.data().all_eq(fmodel_2.f_obs.data())
                      assert fmodel_1.r_free_flags.data().all_eq(fmodel_2.r_free_flags.data())
                      assert abs(fmodel_1.f_calc).data().all_eq(abs(fmodel_2.f_calc).data())
                      assert abs(fmodel_1.f_model()).data().all_eq(abs(fmodel_2.f_model()).data())
                      assert fmodel_1.fu_aniso().all_eq(fmodel_2.fu_aniso())
                      assert fmodel_1.fu_aniso_w().all_eq(fmodel_2.fu_aniso_w())
                      assert fmodel_1.fu_aniso_t().all_eq(fmodel_2.fu_aniso_t())
                      assert fmodel_1.f_obs_w().data().all_eq(
                                                     fmodel_2.f_obs_w().data())
                      assert fmodel_1.f_obs_t().data().all_eq(
                                                     fmodel_2.f_obs_t().data())
                      assert abs(fmodel_1.f_calc_w()).data().all_eq(
                                               abs(fmodel_2.f_calc_w()).data())
                      assert abs(fmodel_1.f_calc_t()).data().all_eq(
                                               abs(fmodel_2.f_calc_t()).data())
                      assert abs(fmodel_1.f_model_w()).data().all_eq(
                                              abs(fmodel_2.f_model_w()).data())
                      assert abs(fmodel_1.f_model_t()).data().all_eq(
                                              abs(fmodel_2.f_model_t()).data())
                      assert abs(fmodel_1.f_bulk_w()).data().all_eq(abs(fmodel_2.f_bulk_w()).data())
                      assert abs(fmodel_1.f_bulk_t()).data().all_eq(abs(fmodel_2.f_bulk_t()).data())
                      assert abs(fmodel_1.f_mask_w()).data().all_eq(abs(fmodel_2.f_mask_w()).data())
                      assert abs(fmodel_1.f_mask_t()).data().all_eq(abs(fmodel_2.f_mask_t()).data())
                      assert abs(fmodel_1.f_mask).data().all_eq(abs(fmodel_1.f_mask).data())
                      assert fmodel_1.figures_of_merit() .all_approx_equal(fmodel_2.figures_of_merit() )
                      assert fmodel_1.phase_errors()     .all_approx_equal(fmodel_2.phase_errors()     )
                      assert fmodel_1.phase_errors_work().all_approx_equal(fmodel_2.phase_errors_work())
                      assert fmodel_1.phase_errors_test().all_approx_equal(fmodel_2.phase_errors_test())
                      a1,b1 = fmodel_1.alpha_beta()
                      a2,b2 = fmodel_2.alpha_beta()
                      assert a1.data().all_approx_equal(a2.data())
                      assert b1.data().all_approx_equal(b2.data())
                      a1,b1 = fmodel_1.alpha_beta_w()
                      a2,b2 = fmodel_2.alpha_beta_w()
                      assert a1.data().all_approx_equal(a2.data())
                      assert b1.data().all_approx_equal(b2.data())
                      a1,b1 = fmodel_1.alpha_beta_t()
                      a2,b2 = fmodel_2.alpha_beta_t()
                      assert a1.data().all_approx_equal(a2.data())
                      assert b1.data().all_approx_equal(b2.data())
                      assert fmodel_1.f_ordered_solvent.data()    .all_eq(fmodel_2.f_ordered_solvent.data()    )
                      assert fmodel_1.f_ordered_solvent_w().data().all_eq(fmodel_2.f_ordered_solvent_w().data())
                      assert fmodel_1.f_ordered_solvent_t().data().all_eq(fmodel_2.f_ordered_solvent_t().data())

def run():
  n_elements = 70
  sgs = ["P 1", "P 4", "C 1 2/c 1"]
  for sg in sgs:
      print sg
      xray_structure = random_structure.xray_structure(
             space_group_info       = sgtbx.space_group_info(sg),
             elements               =(("O","N","C")*(n_elements/3+1))[:n_elements],
             volume_per_atom        = 100,
             min_distance           = 1.5,
             general_positions_only = True,
             anisotropic_flag       = False,
             random_u_iso           = False,
             random_occupancy       = False)
      xray_structure.scattering_type_registry(table="wk1995")
      test_1(xray_structure)

if (__name__ == "__main__"):
  run()
  print "OK"





#def run_one(space_group_info):
#  n_elements = 70
#  xray_structure = random_structure.xray_structure(
#         space_group_info       = space_group_info,
#         elements               =(("O","N","C")*(n_elements/3+1))[:n_elements],
#         volume_per_atom        = 100,
#         min_distance           = 1.5,
#         general_positions_only = True,
#         anisotropic_flag       = False,
#         random_u_iso           = False,
#         random_occupancy       = False)
#  xray_structure.scattering_type_registry(table="wk1995")
#  test_1(xray_structure)
#
#def run_call_back(flags, space_group_info):
#  run_one(space_group_info=space_group_info)
#
#def run():
#  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
#  print "OK"
#
#if (__name__ == "__main__"):
#  run()
