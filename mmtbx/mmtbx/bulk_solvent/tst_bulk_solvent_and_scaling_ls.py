from iotbx import pdb
from cctbx.array_family import flex
from cctbx import xray
import time, math,os,sys
from libtbx.test_utils import approx_equal
from cctbx import miller
from libtbx.utils import format_cpu_times
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from mmtbx import masks
from mmtbx import bulk_solvent
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from cctbx import crystal
import random
from libtbx import adopt_init_args
import mmtbx.f_model
from libtbx import introspection
from cctbx import sgtbx
import cctbx.sgtbx.bravais_types
from cctbx import adptbx
from cctbx.development import random_structure
import libtbx.load_env

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

class alpha_beta_parameters(object):
  def __init__(self, test_ref_in_bin  = 200,
                     n_atoms_absent   = 325,
                     n_atoms_included = 0,
                     bf_atoms_absent  = 25.0,
                     final_error      = 0.0,
                     absent_atom_type = "O",
                     method           = "est",
                     verbose          = -1,
                     expert_mode      = False,
                     interpolation    = False,
                     fix_scale_for_calc_option = 1.0):
    adopt_init_args(self, locals())

class mask_parameters(object):
  def __init__(self, solvent_radius = 1.0,
                     shrink_truncation_radius = 1.0,
                     grid_step_factor = 4.0,
                     verbose = -1):
    adopt_init_args(self, locals())


#TEST-2-------------------------------------------------LS target
def data_prep(f_calc, f_mask, ss,
              k_sol = 0.33,
              b_sol = 55.0,
              u_aniso = [10.0, 20.0, -3,0, 40.0, 0.0]):

    f_bulk_data = f_mask.data() * flex.exp(-ss * b_sol) * k_sol
    fu = bulk_solvent.fu_aniso(u_aniso,
                               f_mask.indices(),
                               f_mask.unit_cell())
    f_model_data = (f_calc.data() + f_bulk_data)*fu
    f_obs = miller.array(miller_set = f_mask,
                         data = flex.abs(f_model_data))
    return k_sol, b_sol, u_aniso, f_obs

def assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags, tk,tb,tu):
    assert approx_equal(u_aniso, fmodel.u_aniso, tu)
    assert approx_equal(k_sol,   fmodel.k_sol,   tk)
    assert approx_equal(b_sol,   fmodel.b_sol,   tb)
    r_work = bulk_solvent.r_factor(
                              f_obs.select(~fmodel.r_free_flags.data()).data(),
                              fmodel.f_model_w().data())
    r_test = bulk_solvent.r_factor(
                               f_obs.select(fmodel.r_free_flags.data()).data(),
                               fmodel.f_model_t().data())
    assert approx_equal(r_work, fmodel.r_work())
    assert approx_equal(r_test, fmodel.r_free())
    assert fmodel.f_obs_w().data().all_eq(
                                     f_obs.select(~r_free_flags.data()).data())
    assert fmodel.f_obs_t().data().all_eq(
                                     f_obs.select(r_free_flags.data()).data())
    assert fmodel.f_obs is f_obs

def exercise_0(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs)

    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = True
    params.anisotropic_scaling                      = True
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = True
    params.minimization_k_sol_b_sol                 = True
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.05
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 3
    params.number_of_minimization_macro_cycles      = 5
    params.number_of_cycles_for_anisotropic_scaling = 2
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-5, tb = 1.e-2, tu = 1.e-3)
    print "OK: LS min.&grid s.: ",format_cpu_times()

def exercise_1(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = True
    params.anisotropic_scaling                      = True
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = False
    params.minimization_k_sol_b_sol                 = True
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.025
    params.b_sol_step                               = 1.0
    params.number_of_macro_cycles                   = 4
    params.number_of_minimization_macro_cycles      = 3
    params.number_of_cycles_for_anisotropic_scaling = 5
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-2, tb = 5.0, tu = 1.e-1)
    print "OK: LS minimization: ",format_cpu_times()

def exercise_2(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = True
    params.anisotropic_scaling                      = True
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = False
    params.minimization_k_sol_b_sol                 = False
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.1
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 1
    params.number_of_minimization_macro_cycles      = 5
    params.number_of_cycles_for_anisotropic_scaling = 5
    params.fix_k_sol                                = 0.33
    params.fix_b_sol                                = 55.0
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-3, tb = 1.e-3, tu = 1.e-3)
    print "OK: LS fix_ksolbsol: ",format_cpu_times()

def exercise_3(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs,
                  k_sol             = k_sol,
                  b_sol             = b_sol)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = False
    params.anisotropic_scaling                      = True
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = False
    params.minimization_k_sol_b_sol                 = False
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.1
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 1
    params.number_of_minimization_macro_cycles      = 5
    params.number_of_cycles_for_anisotropic_scaling = 5
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-3, tb = 1.e-3, tu = 1.e-3)
    print "OK: LS fix_ksolbsol: ",format_cpu_times()

def exercise_4(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs,
                  k_sol             = k_sol,
                  b_sol             = b_sol,
                  u_aniso           = u_aniso)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = False
    params.anisotropic_scaling                      = False
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = False
    params.minimization_k_sol_b_sol                 = False
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.1
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 2
    params.number_of_minimization_macro_cycles      = 2
    params.number_of_cycles_for_anisotropic_scaling = 2
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-6, tb = 1.e-6, tu = 1.e-6)
    print "OK: LS fix_all 1:    ",format_cpu_times()

def exercise_5(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = True
    params.anisotropic_scaling                      = True
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = False
    params.minimization_k_sol_b_sol                 = False
    params.minimization_u_aniso                     = False
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.1
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 2
    params.number_of_minimization_macro_cycles      = 2
    params.number_of_cycles_for_anisotropic_scaling = 2
    params.fix_k_sol                                = k_sol
    params.fix_b_sol                                = b_sol
    params.fix_u_aniso                              = u_aniso
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-6, tb = 1.e-6, tu = 1.e-6)
    print "OK: LS fix_all 2:    ",format_cpu_times()

def exercise_51(fmodel,
                f_mask,
                f_calc,
                ss,
                r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs,
                  k_sol             = k_sol,
                  b_sol             = b_sol,
                  u_aniso           = u_aniso)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = False
    params.anisotropic_scaling                      = False
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = True
    params.minimization_k_sol_b_sol                 = True
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.1
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 2
    params.number_of_minimization_macro_cycles      = 2
    params.number_of_cycles_for_anisotropic_scaling = 2
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False
    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-6, tb = 1.e-6, tu = 1.e-6)
    print "OK: LS fix_all 3:    ",format_cpu_times()


def exercise_6(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs,
                  k_sol             = k_sol,
                  b_sol             = b_sol,
                  u_aniso           = u_aniso)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = True
    params.anisotropic_scaling                      = True
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = True
    params.minimization_k_sol_b_sol                 = True
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.1
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 3
    params.number_of_minimization_macro_cycles      = 3
    params.number_of_cycles_for_anisotropic_scaling = 3
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-5, tb = 1.e-6, tu = 1.e-3)
    print "OK: LS fix_all 4:    ",format_cpu_times()

def exercise_7(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs,
                  k_sol             = k_sol,
                  b_sol             = b_sol)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = False
    params.anisotropic_scaling                      = True
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = True
    params.minimization_k_sol_b_sol                 = True
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.1
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 3
    params.number_of_minimization_macro_cycles      = 3
    params.number_of_cycles_for_anisotropic_scaling = 3
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-3, tb = 1.e-3, tu = 1.e-3)
    print "OK: exercise_7:      ",format_cpu_times()

def exercise_8(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss, k_sol=0, b_sol=0)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs,
                  k_sol             = k_sol,
                  b_sol             = b_sol)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = False
    params.anisotropic_scaling                      = True
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = True
    params.minimization_k_sol_b_sol                 = True
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.1
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 3
    params.number_of_minimization_macro_cycles      = 3
    params.number_of_cycles_for_anisotropic_scaling = 3
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-3, tb = 1.e-3, tu = 1.e-3)
    print "OK: exercise_8:      ",format_cpu_times()

def exercise_9(fmodel,
               f_mask,
               f_calc,
               ss,
               r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss)

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs,
                  u_aniso           = u_aniso)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = True
    params.anisotropic_scaling                      = False
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = True
    params.minimization_k_sol_b_sol                 = True
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.025
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 5
    params.number_of_minimization_macro_cycles      = 5
    params.number_of_cycles_for_anisotropic_scaling = 5
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-3, tb = 1.e-3, tu = 1.e-6)
    print "OK: exercise_9:      ",format_cpu_times()

def exercise_10(fmodel,
                f_mask,
                f_calc,
                ss,
                r_free_flags):
    k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss, u_aniso=[0,0,0,0,0,0])

    fmodel.update(alpha_beta_params = alpha_beta_parameters(),
                  f_obs             = f_obs,
                  u_aniso           = u_aniso)
    params = bss.solvent_and_scale_params()
    params.bulk_solvent_correction                  = True
    params.anisotropic_scaling                      = False
    params.statistical_solvent_model                = False
    params.k_sol_b_sol_grid_search                  = True
    params.minimization_k_sol_b_sol                 = True
    params.minimization_u_aniso                     = True
    params.target                                   = "ls_wunit_k1"
    params.symmetry_constraints_on_u_aniso          = True
    params.k_sol_max                                = 0.6
    params.k_sol_min                                = 0.1
    params.b_sol_max                                = 80.0
    params.b_sol_min                                = 10.0
    params.k_sol_step                               = 0.05
    params.b_sol_step                               = 5.0
    params.number_of_macro_cycles                   = 5
    params.number_of_minimization_macro_cycles      = 5
    params.number_of_cycles_for_anisotropic_scaling = 5
    params.fix_k_sol                                = None
    params.fix_b_sol                                = None
    params.fix_u_aniso                              = None
    params.start_minimization_from_k_sol            = 0.35
    params.start_minimization_from_b_sol            = 46.0
    params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
    params.apply_back_trace_of_u_aniso              = False

    fmodel.update_solvent_and_scale(params = params)
    assert_result(fmodel, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                  tk = 1.e-6, tb = 1.e-2, tu = 1.e-6)
    print "OK: exercise_10:     ",format_cpu_times()

def exercise_11(fmodel,
                f_mask,
                f_calc,
                ss,
                r_free_flags):
    k_sol_range = (0.25,0.48,0.80,0.72,0.31,0.63)
    b_sol_range = (10.4,18.9,45.2,77.1,16.0,70.45)
    u_aniso = [4,5,-9,0,40,0]
    for k_sol_,b_sol_ in zip(k_sol_range,b_sol_range):

        k_sol, b_sol, u_aniso, f_obs = data_prep(f_calc, f_mask, ss,
                                                 k_sol = k_sol_,
                                                 b_sol = b_sol_,
                                                 u_aniso = u_aniso)
        fmodel_copy = fmodel.deep_copy()
        fmodel_copy.update(alpha_beta_params = alpha_beta_parameters(),
                      f_obs             = f_obs)
        params = bss.solvent_and_scale_params()
        params.bulk_solvent_correction                  = True
        params.anisotropic_scaling                      = True
        params.statistical_solvent_model                = True
        params.k_sol_b_sol_grid_search                  = True
        params.minimization_k_sol_b_sol                 = True
        params.minimization_u_aniso                     = True
        params.target                                   = "ls_wunit_k1"
        params.symmetry_constraints_on_u_aniso          = True
        params.k_sol_max                                = 0.8
        params.k_sol_min                                = 0.1
        params.b_sol_max                                = 80.0
        params.b_sol_min                                = 10.0
        params.k_sol_step                               = 0.05
        params.b_sol_step                               = 5.0
        params.number_of_macro_cycles                   = 2
        params.number_of_minimization_macro_cycles      = 5
        params.number_of_cycles_for_anisotropic_scaling = 5
        params.fix_k_sol                                = None
        params.fix_b_sol                                = None
        params.fix_u_aniso                              = None
        params.start_minimization_from_k_sol            = 0.35
        params.start_minimization_from_b_sol            = 46.0
        params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
        params.apply_back_trace_of_u_aniso              = False

        fmodel_copy.update_solvent_and_scale(params = params)
        assert_result(fmodel_copy, k_sol, b_sol, u_aniso, f_obs, r_free_flags,
                      tk = 1.e-2, tb = 1.e-2, tu = 1.e-2)
    print "OK: closest to real: ",format_cpu_times()

def exercise_12(fmodel):
  for symbol in sgtbx.bravais_types.acentric + sgtbx.bravais_types.centric:
      space_group_info = sgtbx.space_group_info(symbol = symbol)
      random_xray_structure = random_structure.xray_structure(
                                       space_group_info  = space_group_info,
                                       elements          = ["N"]*20,
                                       volume_per_atom   = 50.0,
                                       anisotropic_flag  = False,
                                       random_u_iso      = False,
                                       u_iso             = adptbx.b_as_u(20.0))
      sg = random_xray_structure.space_group()
      uc = random_xray_structure.unit_cell()
      f_calc = random_xray_structure.structure_factors(
                                           d_min          = 2.0,
                                           anomalous_flag = False,
                                           algorithm      = "fft").f_calc()
      # to avoid known problems, re-calculate f_calc
      f_calc = f_calc.structure_factors_from_scatterers(
                                xray_structure = random_xray_structure).f_calc()
      r_free_flags = miller.array(miller_set = f_calc,
                                  data = flex.bool(f_calc.data().size(), False))

      fmodel = mmtbx.f_model.manager(
        sf_algorithm   = "fft",
        r_free_flags   = r_free_flags,
        f_obs          = f_calc.array(data=flex.double(f_calc.data().size(),1.0)),
        xray_structure = random_xray_structure,
        target_name    = "ls_wunit_k1",
        mask_params    = mask_parameters())
      fmodel.update_xray_structure(xray_structure = random_xray_structure,
                                   update_f_calc  = False,
                                   update_f_mask  = True,
                                   update_f_ordered_solvent = False)
      u_cart_p1 = adptbx.random_u_cart(u_scale=5, u_min=5)
      u_star_p1 = adptbx.u_cart_as_u_star(uc, u_cart_p1)
      u_aniso_1 = adptbx.u_star_as_u_cart(uc, u_star_p1)
      u_aniso_2 = adptbx.u_star_as_u_cart(uc, sg.average_u_star(u_star = u_star_p1))
      for u_aniso in (u_aniso_1, u_aniso_2):
          fu = bulk_solvent.fu_aniso(u_aniso,
                                     f_calc.indices(),
                                     uc)
          f_obs = miller.array(miller_set = f_calc,
                               data       = flex.abs(f_calc.data())*fu)
          fmodel_copy = fmodel.deep_copy()
          fmodel_copy.update(alpha_beta_params = alpha_beta_parameters(),
                             f_obs             = f_obs,
                             k_sol             = 0,
                             b_sol             = 0)
          for flag in (True, False):
              params = bss.solvent_and_scale_params()
              params.bulk_solvent_correction                  = False
              params.anisotropic_scaling                      = True
              params.statistical_solvent_model                = False
              params.k_sol_b_sol_grid_search                  = False
              params.minimization_k_sol_b_sol                 = False
              params.minimization_u_aniso                     = True
              params.target                                   = "ls_wunit_k1"
              params.symmetry_constraints_on_u_aniso          = flag
              params.k_sol_max                                = 0.8
              params.k_sol_min                                = 0.1
              params.b_sol_max                                = 80.0
              params.b_sol_min                                = 10.0
              params.k_sol_step                               = 0.05
              params.b_sol_step                               = 5.0
              params.number_of_macro_cycles                   = 1
              params.number_of_minimization_macro_cycles      = 0
              params.number_of_cycles_for_anisotropic_scaling = 30
              params.fix_k_sol                                = None
              params.fix_b_sol                                = None
              params.fix_u_aniso                              = None
              params.start_minimization_from_k_sol            = 0.35
              params.start_minimization_from_b_sol            = 46.0
              params.start_minimization_from_u_aniso          = [0,0,0,0,0,0]
              params.apply_back_trace_of_u_aniso              = False
              fmodel_copy.update_solvent_and_scale(params = params)
              #print
              #print flag, sg.type().number()
              #print "u        = %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f" % \
              #  (u_aniso[0],u_aniso[1],u_aniso[2],u_aniso[3],u_aniso[4],u_aniso[5])
              #print "fmodel.u = %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f" % \
              #  (fmodel_copy.u_aniso[0],fmodel_copy.u_aniso[1],fmodel_copy.u_aniso[2],
              #   fmodel_copy.u_aniso[3],fmodel_copy.u_aniso[4],fmodel_copy.u_aniso[5])
              tolerance = 1.e-4
              if(flag == False and approx_equal(u_aniso, u_aniso_1, out=None)):
                 assert approx_equal(fmodel_copy.u_aniso, u_aniso, tolerance)
              if(flag == True and approx_equal(u_aniso, u_aniso_2, out=None)):
                 assert approx_equal(fmodel_copy.u_aniso, u_aniso, tolerance)
              if(flag == False and approx_equal(u_aniso, u_aniso_2, out=None)):
                 assert approx_equal(fmodel_copy.u_aniso, u_aniso, tolerance)
              if(flag == True and approx_equal(u_aniso, u_aniso_1, out=None)):
                 for u2, ufm in zip(u_aniso_2, fmodel_copy.u_aniso):
                   if(abs(u2) < 1.e-6):
                      assert approx_equal(ufm, 0.0, tolerance)
  print "OK: uaniso constr.:  ",format_cpu_times()

def scale_from_ls(d_min = 2.2,
                  flag = True,
                test_1 = True,
                test_2 = True,
                test_3 = True,
                test_4 = True,
                test_5 = True,
                test_6 = True,
                test_7 = True,
                test_8 = True):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="regression/pdb/2ERL.pdb", test=os.path.isfile)
  if (pdb_file is None):
    print "Skipping scale_from_ml(): input file not available"
    return
  xray_structure = pdb.as_xray_structure(pdb_file)
  f_calc = xray_structure.structure_factors(d_min = d_min,
                                       anomalous_flag = False).f_calc()
  ss = 1./flex.pow2(f_calc.d_spacings().data()) / 4.
  r_free_flags = f_calc.array(
    data=flex.size_t(xrange(1,f_calc.data().size()+1)) % 10 == 0)
  #
  # create fmodel with dummy f_obs
  #
  fmodel = mmtbx.f_model.manager(
        sf_algorithm   = "fft",
        r_free_flags   = r_free_flags,
        f_obs          = f_calc.array(data=flex.double(f_calc.data().size(),1.0)),
        xray_structure = xray_structure,
        target_name    = "ls_wunit_k1",
        mask_params    = mask_parameters())
  fmodel.update_xray_structure(xray_structure = xray_structure,
                               update_f_calc = True,
                               update_f_mask = True,
                               update_f_ordered_solvent = False)
  f_mask = masks.bulk_solvent(
    xray_structure=xray_structure,
    grid_step=f_calc.d_min()/4,
    shrink_truncation_radius=mask_parameters().shrink_truncation_radius,
    solvent_radius=mask_parameters().solvent_radius).structure_factors(
      miller_set=f_calc)

###
###>> METHOD: minimization & grid search
###
  exercise_0(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)

###
###>> METHOD: minimization
###
  exercise_1(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)
###
###>> METHOD: fix_k_sol_b_sol
###
  exercise_2(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)

  exercise_3(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)
###
###>> METHOD: fix_all
###
  exercise_4(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)

  exercise_5(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)

  exercise_51(fmodel       = fmodel.deep_copy(),
              f_mask       = f_mask,
              f_calc       = f_calc,
              ss           = ss,
              r_free_flags = r_free_flags)
###
###>> METHOD: make sure nothing is changed if start with exact values
###           and ask to refine these values
###
  exercise_6(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)
###
###>> Given non-zero k_sol & b_sol fixed, find u_aniso;
###>> keep unchanged k_sol & b_sol too
###
  exercise_7(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)
###
###>> Given zero k_sol & b_sol fixed, find u_aniso;
###>> keep unchanged k_sol & b_sol too
###
  exercise_8(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)
###
###>> Given non-zero u_aniso fixed, find k_sol & b_sol;
###>> keep unchanged u_aniso too
###
  exercise_9(fmodel       = fmodel.deep_copy(),
             f_mask       = f_mask,
             f_calc       = f_calc,
             ss           = ss,
             r_free_flags = r_free_flags)
###
###>> Given zero u_aniso fixed, find k_sol & b_sol;
###>> keep unchanged u_aniso too
###
  exercise_10(fmodel       = fmodel.deep_copy(),
              f_mask       = f_mask,
              f_calc       = f_calc,
              ss           = ss,
              r_free_flags = r_free_flags)
###
###>> Find k_sol,b_sol,u_aniso starting from zero.
###>> Very real case. Not always precise!
###
  exercise_11(fmodel       = fmodel.deep_copy(),
              f_mask       = f_mask,
              f_calc       = f_calc,
              ss           = ss,
              r_free_flags = r_free_flags)
###
###>> Test symmetry constraints for U_aniso.
###
  exercise_12(fmodel = fmodel.deep_copy())


def run():
  scale_from_ls()
  print "OK: scale_from_ls(): ",format_cpu_times()

if (__name__ == "__main__"):
  run()
