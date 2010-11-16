from cctbx.array_family import flex
import mmtbx.f_model
from cctbx.development import random_structure
from cctbx import sgtbx
from libtbx.test_utils import approx_equal, show_diff
from libtbx.utils import format_cpu_times
import pickle
from cStringIO import StringIO
import random
import mmtbx

random.seed(0)
flex.set_random_seed(0)

def test_1(xray_structure):
  # exercise almost all without dealing with particular values
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  for d_min in [2.0, 2.5]:
      for algorithm in ["direct", "fft"]:
          sfg_params.algorithm = algorithm
          for anomalous_flag in [True, False]:
              for cos_sin_table in [True, False]:
                  sfg_params.cos_sin_table = cos_sin_table
                  f_obs = abs(xray_structure.structure_factors(
                                       d_min          = d_min,
                                       anomalous_flag = anomalous_flag,
                                       cos_sin_table  = cos_sin_table,
                                       algorithm      = algorithm).f_calc())
                  f_obs_comp = f_obs.structure_factors_from_scatterers(
                              xray_structure = xray_structure,
                              algorithm      = algorithm,
                              cos_sin_table  = cos_sin_table).f_calc()
                  f_obs = abs(f_obs_comp)
                  flags = f_obs.generate_r_free_flags(fraction = 0.1,
                                                      max_free = 99999999)
                  f_obs.set_observation_type_xray_amplitude()
                  for (xrs,fc) in ((xray_structure,None),(None,f_obs_comp)):
                      ###
                      ### instantiate fmodel only
                      ###
                      if(xrs is None):
                         zero = flex.complex_double(f_obs.data().size(), 0.0)
                         f_mask = f_obs.array(data = zero)
                      else:
                         f_mask = None
                      fmodel = mmtbx.f_model.manager(
                                              xray_structure    = xrs,
                                              f_calc            = fc,
                                              f_mask            = f_mask,
                                              f_obs             = f_obs,
                                              r_free_flags      = flags,
                                              target_name       = "ls_wunit_k1",
                                              sf_and_grads_accuracy_params = sfg_params)
                      fmodel_info = fmodel.info()
                      assert fmodel.f_obs.data().all_eq(f_obs.data())
                      assert abs(fmodel.f_calc()).data().all_eq(f_obs.data())
                      assert abs(fmodel.f_model()).data().all_eq(f_obs.data())
                      assert approx_equal(fmodel.r_work(), 0, 1.e-9)
                      assert approx_equal(fmodel.r_free(), 0, 1.e-9)
                      assert approx_equal(fmodel.r_all(), 0, 1.e-9)
                      assert fmodel.fb_cart().all_eq(1.0)
                      assert fmodel.fb_cart_work().all_eq(1.0)
                      assert fmodel.fb_cart_t().all_eq(1.0)
                      assert abs(fmodel.target_w()) < 1.e-9
                      assert abs(fmodel.target_t()) < 1.e-9
                      assert fmodel.k_sol_b_sol() == (0.0,0.0)
                      assert approx_equal(fmodel.b_cart(),[0,0,0,0,0,0])
                      assert fmodel.b_iso() == 0.0
                      assert fmodel.f_obs_w.data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert fmodel.f_obs_t.data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_calc_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_calc_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_model_work()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_model_free()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_bulk_w()).data().all_eq(0)
                      assert abs(fmodel.f_bulk_t()).data().all_eq(0)
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
                      fmodel.model_error_ml()
                      ###
                      ### instantiate fmodel only + update ksol & bsol
                      ###
                      if(xrs is not None):
                         fc_ = None
                      else:
                         fc_ = fc
                         zero = flex.complex_double(f_obs.data().size(), 0.0)
                         f_mask = f_obs.array(data = zero)
                      fmodel = mmtbx.f_model.manager(
                                              xray_structure    = xrs,
                                              f_calc            = fc_,
                                              f_obs             = f_obs,
                                              f_mask            = f_mask,
                                              r_free_flags      = flags,
                                              target_name       = "ls_wunit_k1",
                                              sf_and_grads_accuracy_params = sfg_params)
                      fmodel_info = fmodel.info()
                      fmodel.update(k_sol = 0.5, b_sol = 35.0)
                      assert fmodel.f_obs.data().all_eq(f_obs.data())
                      assert abs(fmodel.f_calc()).data().all_eq(f_obs.data())
                      assert fmodel.fb_cart().all_eq(1.0)
                      assert fmodel.fb_cart_work().all_eq(1.0)
                      assert fmodel.fb_cart_t().all_eq(1.0)
                      assert fmodel.k_sol_b_sol() == (0.5,35.0)
                      assert approx_equal(fmodel.b_cart(),[0,0,0,0,0,0])
                      assert fmodel.b_iso() == 0.0
                      assert fmodel.f_obs_w.data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert fmodel.f_obs_t.data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_calc_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_calc_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      fmodel.model_error_ml()
                      #
                      p = pickle.dumps(fmodel, 1)
                      l = pickle.loads(p)
                      s = StringIO()
                      sl = StringIO()
                      fmodel.info().show_all(out=s)
                      l.info().show_all(out=sl)
                      assert not show_diff(sl.getvalue(), s.getvalue())
                      ###
                      ### instantiate fmodel only, then use deep_copy
                      ###
                      fmodel=None
                      fmodel_ = mmtbx.f_model.manager(
                                              xray_structure    = xrs,
                                              f_calc            = fc_,
                                              f_mask            = f_mask,
                                              f_obs             = f_obs,
                                              r_free_flags      = flags,
                                              target_name       = "ls_wunit_k1",
                                              sf_and_grads_accuracy_params = sfg_params)
                      fmodel = fmodel_.deep_copy()
                      fmodel_info = fmodel.info()
                      assert fmodel.f_obs.data().all_eq(f_obs.data())
                      assert abs(fmodel.f_calc()).data().all_eq(f_obs.data())
                      assert abs(fmodel.f_model()).data().all_eq(f_obs.data())
                      assert approx_equal(fmodel.r_work(), 0, 1.e-9)
                      assert approx_equal(fmodel.r_free(), 0, 1.e-9)
                      assert approx_equal(fmodel.r_all(), 0, 1.e-9)
                      assert fmodel.fb_cart().all_eq(1.0)
                      assert fmodel.fb_cart_work().all_eq(1.0)
                      assert fmodel.fb_cart_t().all_eq(1.0)
                      assert abs(fmodel.target_w()) < 1.e-9
                      assert abs(fmodel.target_t()) < 1.e-9
                      assert fmodel.k_sol_b_sol() == (0.0,0.0)
                      assert approx_equal(fmodel.b_cart(),[0,0,0,0,0,0])
                      assert fmodel.b_iso() == 0.0
                      assert fmodel.f_obs_w.data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert fmodel.f_obs_t.data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_calc_w()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_calc_t()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_model_work()).data().all_eq(
                                                f_obs.select(~flags.data()).data())
                      assert abs(fmodel.f_model_free()).data().all_eq(
                                                 f_obs.select(flags.data()).data())
                      assert abs(fmodel.f_bulk_w()).data().all_eq(0)
                      assert abs(fmodel.f_bulk_t()).data().all_eq(0)
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
                      fmodel.model_error_ml()
                      ###
                      ### instantiate fmodel only, then use resolution_filter
                      ###
                      fmodel=None
                      d_max_ = 3.0
                      d_min_ = 2.7
                      if(xrs is not None):
                         fc_ = None
                      else:
                         fc_ = fc
                         zero = flex.complex_double(f_obs.data().size(), 0.0)
                         f_mask = f_obs.array(data = zero)
                      fmodel_ = mmtbx.f_model.manager(
                                              xray_structure    = xrs,
                                              f_calc            = fc_,
                                              f_mask            = f_mask,
                                              f_obs             = f_obs,
                                              r_free_flags      = flags,
                                              target_name       = "ls_wunit_k1",
                                              sf_and_grads_accuracy_params = sfg_params)
                      fmodel_1 = fmodel_.resolution_filter(d_max = d_max_, d_min = d_min_, update_xray_structure=True)
                      fmodel_info = fmodel_1.info()
                      fmodel_info = fmodel_.info()
                      if(fc is not None):
                         fc_ = fc.resolution_filter(d_max = d_max_, d_min = d_min_)
                         fm_ = f_mask.resolution_filter(d_max = d_max_, d_min = d_min_)
                      else:
                         fc_ = None
                         fm_ = None
                      fmodel_2 = mmtbx.f_model.manager(
                                xray_structure    = xrs,
                                f_calc            = fc_,
                                f_mask            = fm_,
                                f_obs             = f_obs.resolution_filter(d_max = d_max_, d_min = d_min_),
                                r_free_flags      = flags.resolution_filter(d_max = d_max_, d_min = d_min_),
                                target_name       = "ls_wunit_k1",
                                sf_and_grads_accuracy_params = sfg_params)
                      assert fmodel_1.f_obs.data().all_eq(fmodel_2.f_obs.data())
                      assert fmodel_1.r_free_flags.data().all_eq(fmodel_2.r_free_flags.data())
                      assert abs(fmodel_1.f_calc()).data().all_eq(abs(fmodel_2.f_calc()).data())
                      assert abs(fmodel_1.f_model()).data().all_eq(abs(fmodel_2.f_model()).data())
                      assert fmodel_1.fb_cart().all_eq(fmodel_2.fb_cart())
                      assert fmodel_1.fb_cart_work().all_eq(fmodel_2.fb_cart_work())
                      assert fmodel_1.fb_cart_t().all_eq(fmodel_2.fb_cart_t())
                      assert fmodel_1.f_obs_w.data().all_eq(
                                                     fmodel_2.f_obs_w.data())
                      assert fmodel_1.f_obs_t.data().all_eq(
                                                     fmodel_2.f_obs_t.data())
                      assert abs(fmodel_1.f_calc_w()).data().all_eq(
                                               abs(fmodel_2.f_calc_w()).data())
                      assert abs(fmodel_1.f_calc_t()).data().all_eq(
                                               abs(fmodel_2.f_calc_t()).data())
                      assert abs(fmodel_1.f_model_work()).data().all_eq(
                                              abs(fmodel_2.f_model_work()).data())
                      assert abs(fmodel_1.f_model_free()).data().all_eq(
                                              abs(fmodel_2.f_model_free()).data())
                      assert abs(fmodel_1.f_bulk_w()).data().all_eq(abs(fmodel_2.f_bulk_w()).data())
                      assert abs(fmodel_1.f_bulk_t()).data().all_eq(abs(fmodel_2.f_bulk_t()).data())
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
                      fmodel_1.model_error_ml()
                      fmodel_2.model_error_ml()

def exercise_1():
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
             random_u_iso           = False,
             random_occupancy       = False)
      xray_structure.scattering_type_registry(table="wk1995")
      test_1(xray_structure)

def exercise_2():
  xray_structure = random_structure.xray_structure(
             space_group_info       = sgtbx.space_group_info("C 1 2/c 1"),
             elements               =("O","N","C")*50,
             volume_per_atom        = 100,
             min_distance           = 1.5,
             general_positions_only = True,
             random_u_iso           = True,
             random_occupancy       = True)
  xray_structure.scattering_type_registry(table="wk1995")
  f_obs = abs(xray_structure.structure_factors(d_min = 2.0).f_calc())
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  for algorithm in ["fft", "direct"]:
      sfg_params.algorithm=algorithm
      sfg_params.cos_sin_table = True
      sfg_params.grid_resolution_factor=1/7.
      sfg_params.wing_cutoff = 1.e-8
      flags = f_obs.generate_r_free_flags(fraction = 0.1,max_free = 99999999)
      fmodel = mmtbx.f_model.manager(xray_structure    = xray_structure,
                                     f_obs             = f_obs,
                                     r_free_flags      = flags,
                                     sf_and_grads_accuracy_params = sfg_params)
      f_calc_1 = abs(fmodel.f_calc()).data()
      f_calc_2 = abs(f_obs.structure_factors_from_scatterers(
        xray_structure = xray_structure,
        algorithm                    = sfg_params.algorithm,
        cos_sin_table                = sfg_params.cos_sin_table,
        grid_resolution_factor       = sfg_params.grid_resolution_factor,
        quality_factor               = sfg_params.quality_factor,
        u_base                       = sfg_params.u_base,
        b_base                       = sfg_params.b_base,
        wing_cutoff                  = sfg_params.wing_cutoff,
        exp_table_one_over_step_size = sfg_params.exp_table_one_over_step_size
                                                             ).f_calc()).data()
      delta = flex.abs(f_calc_1-f_calc_2)
      assert approx_equal(flex.sum(delta), 0.0)

def run():
  exercise_1()
  exercise_2()

if (__name__ == "__main__"):
  run()
  print format_cpu_times()
