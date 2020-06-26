from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.f_model
from cctbx.development import random_structure
from cctbx import sgtbx
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import random
import mmtbx
from cctbx import adptbx
import mmtbx.masks
import mmtbx.f_model
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss

random.seed(0)
flex.set_random_seed(0)

def t_1(xray_structure, d_min=3.5):
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  for deep_copy in [True, False]:
    for hl in [True, False]:
      for target_name in ["ls_wunit_k1", "ml"]:
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
              fc = f_obs.structure_factors_from_scatterers(
                xray_structure = xray_structure,
                algorithm      = algorithm,
                cos_sin_table  = cos_sin_table).f_calc()
              f_obs = abs(fc)
              flags = f_obs.generate_r_free_flags(fraction = 0.1,
                                                  max_free = 99999999)
              f_obs.set_observation_type_xray_amplitude()
              fm = f_obs.array(data = flex.complex_double(f_obs.data().size(), 0))
              for xff in [(xray_structure, None, None),(None, fc, fm)]:
                xrs, f_calc, f_mask = xff
                ###
                ### default states
                ###
                hl_coeffs = None
                if(hl):
                  hl_coeffs = f_obs.customized_copy(data = flex.hendrickson_lattman(f_obs.size()))
                fmodel = mmtbx.f_model.manager(
                  xray_structure    = xrs,
                  abcd              = hl_coeffs,
                  f_calc            = f_calc,
                  f_mask            = f_mask,
                  f_obs             = f_obs,
                  r_free_flags      = flags,
                  target_name       = target_name,
                  sf_and_grads_accuracy_params = sfg_params)
                if(deep_copy): fmodel = fmodel.deep_copy()
                #
                fi = fmodel.info()
                assert approx_equal(fi.alpha_work_max        , 1.0)
                assert approx_equal(fi.alpha_work_mean       , 1.0)
                assert approx_equal(fi.alpha_work_min        , 1.0)
                assert approx_equal(fi.beta_work_max         , 0.0)
                assert approx_equal(fi.beta_work_mean        , 0.0)
                assert approx_equal(fi.beta_work_min         , 0.0)
                assert approx_equal(fi.completeness_6_inf    , 1.0)
                assert approx_equal(fi.completeness_d_min_inf, 1.0)
                assert approx_equal(fi.completeness_in_range , 1.0)
                assert approx_equal(fi.fom_work_max          , 1.0)
                assert approx_equal(fi.fom_work_mean         , 1.0)
                assert approx_equal(fi.fom_work_min          , 1.0)
                assert approx_equal(fi.ml_coordinate_error   , 0.0)
                assert approx_equal(fi.ml_phase_error        , 0.0)
                assert approx_equal(fi.number_of_reflections , f_obs.size())
                assert fi.number_of_reflections_merged== f_obs.deep_copy().average_bijvoet_mates().data().size()
                assert fi.number_of_test_reflections  == flags.data().count(True)
                assert fi.number_of_work_reflections  == flags.data().count(False)
                assert approx_equal(fi.overall_scale_k1, 1.0)
                assert approx_equal(fi.pher_free_max   , 0.0)
                assert approx_equal(fi.pher_free_mean  , 0.0)
                assert approx_equal(fi.pher_free_min   , 0.0)
                assert approx_equal(fi.pher_work_max   , 0.0)
                assert approx_equal(fi.pher_work_mean  , 0.0)
                assert approx_equal(fi.pher_work_min   , 0.0)
                assert approx_equal(fi.r_all           , 0.0)
                assert approx_equal(fi.r_free          , 0.0)
                assert approx_equal(fi.r_work          , 0.0)
                assert fi.sf_algorithm == algorithm
                assert fi.target_name  == target_name
                assert approx_equal(fi.target_free, 0.0)
                assert approx_equal(fi.target_work, 0.0)
                assert fi.twin_fraction is None
                assert fi.twin_law      is None
                #
                if(hl): fmodel.hl_coeffs() is not None
                else:   fmodel.hl_coeffs() is None
                assert fmodel.f_model_free().data().size() == flags.data().count(True)
                assert fmodel.f_model_work().data().size() == flags.data().count(False)
                assert fmodel.f_obs().data().all_eq(f_obs.data())
                assert abs(fmodel.f_calc()).data().all_eq(f_obs.data())
                assert abs(fmodel.f_model()).data().all_eq(f_obs.data())
                assert approx_equal(fmodel.r_work(), 0, 1.e-9)
                assert approx_equal(fmodel.r_free(), 0, 1.e-9)
                assert approx_equal(fmodel.r_all(), 0, 1.e-9)
                assert fmodel.k_anisotropic().all_eq(1.0)
                assert fmodel.k_anisotropic_work().all_eq(1.0)
                assert fmodel.k_anisotropic_test().all_eq(1.0)
                assert abs(fmodel.target_w()) < 1.e-9
                assert abs(fmodel.target_t()) < 1.e-9
                assert len(fmodel.k_masks())==1
                assert fmodel.k_masks()[0].all_eq(0)
                assert fmodel.k_isotropic().all_eq(1.0)
                assert fmodel.f_obs_work().data().all_eq(f_obs.select(~flags.data()).data())
                assert fmodel.f_obs_free().data().all_eq(f_obs.select(flags.data()).data())
                assert abs(fmodel.f_calc_w()).data().all_eq(f_obs.select(~flags.data()).data())
                assert abs(fmodel.f_calc_t()).data().all_eq(f_obs.select(flags.data()).data())
                assert abs(fmodel.f_model_work()).data().all_eq(f_obs.select(~flags.data()).data())
                assert abs(fmodel.f_model_free()).data().all_eq(f_obs.select(flags.data()).data())
                assert abs(fmodel.f_bulk_w()).data().all_eq(0)
                assert abs(fmodel.f_bulk_t()).data().all_eq(0)
                assert fmodel.f_model().data().all_eq(fmodel.f_calc().data())
                assert fmodel.f_model_scaled_with_k1().data().all_eq(fmodel.f_calc().data())
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
                assert approx_equal(fmodel.model_error_ml(), 0)
                #
                # update k_isotropic, k_anisotropic, resolution_filter
                #
                fmodel.update_core(
                  k_isotropic   = flex.double(f_obs.data().size(), 1./2),
                  k_anisotropic = flex.double(f_obs.data().size(), 1./3))
                assert approx_equal(fmodel.r_work(), 0)
                assert approx_equal(fmodel.scale_k1(), 6.0)
                fmodel = fmodel.deep_copy()
                assert approx_equal(fmodel.r_work(), 0)
                assert approx_equal(fmodel.scale_k1(), 6.0)
                #
                fmodel_ = fmodel.resolution_filter(d_min=fmodel.f_obs().d_min()+0.5)
                assert fmodel_.f_obs().size() < fmodel.f_obs().size()
                fmodel_ = fmodel.resolution_filter(d_max=fmodel.f_obs().d_max_min()[0]-1.)
                assert fmodel_.f_obs().size() < fmodel.f_obs().size()

def exercise_4_f_hydrogens():
  for d_min in [1,2,3]:
    for q in [0, 0.9, 1]:
      random.seed(0)
      flex.set_random_seed(0)
      x = random_structure.xray_structure(
        space_group_info       = sgtbx.space_group_info("P 4"),
        elements               =(("O","N","C")*5 + ("H",)*95),
        volume_per_atom        = 200,
        min_distance           = 1.5,
        general_positions_only = True,
        random_u_iso           = True,
        random_occupancy       = False)
      hd_sel = x.hd_selection()
      b_isos = x.select(~hd_sel).extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.)
      mmm = b_isos.min_max_mean().as_tuple()
      b_mean = int(mmm[2])
      x = x.set_b_iso(value=b_mean, selection = hd_sel)
      x.scattering_type_registry(table="wk1995")
      x.set_occupancies(value=q, selection = hd_sel)
      fc = x.structure_factors(d_min = d_min, algorithm="direct").f_calc()
      f_obs = abs(fc)
      x = x.deep_copy_scatterers()
      x.set_occupancies(value=0.0, selection = x.hd_selection())
      sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
      sfg_params.algorithm = "direct"
      r_free_flags = f_obs.generate_r_free_flags(fraction = 0.1)
      fmodel = mmtbx.f_model.manager(
        xray_structure = x,
        f_obs          = f_obs,
        r_free_flags   = r_free_flags,
        sf_and_grads_accuracy_params = sfg_params)
      if(q==0):
        assert approx_equal(fmodel.r_work(), 0)
      else:
        assert fmodel.r_work() > 0.05, fmodel.r_work()
      params = bss.master_params.extract()
      params.bulk_solvent=False
      params.anisotropic_scaling=False
      o = fmodel.update_all_scales(fast=False, params=params)
      assert approx_equal(o.k_sol[0],0)
      assert approx_equal(o.b_sol[0],0)
      assert approx_equal(o.b_cart,[0,0,0,0,0,0])
      assert approx_equal(o.k_h, q)
      assert approx_equal(fmodel.r_work(), 0)

def exercise_1():
  n_elements = 70
  sgs = ["P 1", "P 4", "C 1 2/c 1"]
  for sg in sgs:
    print(sg)
    xray_structure = random_structure.xray_structure(
           space_group_info       = sgtbx.space_group_info(sg),
           elements               =(("O","N","C")*(n_elements//3+1))[:n_elements],
           volume_per_atom        = 100,
           min_distance           = 1.5,
           general_positions_only = True,
           random_u_iso           = True,
           random_occupancy       = False)
    xray_structure.scattering_type_registry(table="wk1995")
    t_1(xray_structure = xray_structure)

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
      flags = f_obs.generate_r_free_flags(fraction = 0.1,max_free = 99999999)
      fmodel = mmtbx.f_model.manager(xray_structure    = xray_structure,
                                     f_obs             = f_obs,
                                     r_free_flags      = flags,
                                     sf_and_grads_accuracy_params = sfg_params)
      f_calc_1 = abs(fmodel.f_calc()).data()
      f_calc_2 = abs(f_obs.structure_factors_from_scatterers(
        xray_structure = xray_structure,
        algorithm = sfg_params.algorithm).f_calc()).data()
      delta = flex.abs(f_calc_1-f_calc_2)
      assert approx_equal(flex.sum(delta), 0.0)

def exercise_5_bulk_sol_and_scaling(d_min, symbol = "C 2", k_sol = 0.37,
                                    b_sol = 64.0):
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(symbol=symbol),
    elements               =(("O","N","C")*150),
    volume_per_atom        = 200,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  x.scattering_type_registry(table="wk1995")
  f_calc = x.structure_factors(d_min = d_min, algorithm="direct").f_calc()
  mask_manager = mmtbx.masks.manager(miller_array = f_calc)
  f_mask = mask_manager.shell_f_masks(xray_structure = x)[0]
  assert flex.mean(abs(f_mask).data()) > 0
  b_cart=[-15,-5,20, 0,6,0]
  u_star = adptbx.u_cart_as_u_star(x.unit_cell(), adptbx.b_as_u(b_cart))
  k_anisotropic = mmtbx.f_model.ext.k_anisotropic(f_calc.indices(), u_star)
  ss = 1./flex.pow2(f_calc.d_spacings().data()) / 4.
  k_mask = mmtbx.f_model.ext.k_mask(ss, k_sol, b_sol)
  scale = 17.
  k_isotropic = flex.double(f_calc.data().size(), scale)
  f_model_data = scale*k_anisotropic*(f_calc.data()+k_mask*f_mask.data())
  f_model = f_calc.customized_copy(data = f_model_data)
  f_obs = abs(f_model)
  r_free_flags = f_obs.generate_r_free_flags(use_lattice_symmetry=False)
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg_params.algorithm = "direct"
  bin_selections = []
  f_calc.setup_binner(reflections_per_bin=100)
  for i_bin in f_calc.binner().range_used():
    sel = f_calc.binner().selection(i_bin)
    bin_selections.append(sel)
  # test 1
  for k_isotropic_ in [None, k_isotropic]:
    fmodel = mmtbx.f_model.manager(
      xray_structure = x,
      f_obs          = f_obs,
      k_mask         = k_mask,
      k_anisotropic  = k_anisotropic,
      k_isotropic    = k_isotropic,
      r_free_flags   = r_free_flags,
      bin_selections = bin_selections,
      sf_and_grads_accuracy_params = sfg_params)
    assert approx_equal(fmodel.r_work(), 0)
    assert approx_equal(fmodel.r_free(), 0)
    assert approx_equal(fmodel.target_w(), 0)
    assert approx_equal(fmodel.k_masks()[0], k_mask)
    assert approx_equal(fmodel.k_anisotropic(), k_anisotropic)
    assert approx_equal(fmodel.f_model_scaled_with_k1().data(), f_model_data)
    assert approx_equal(fmodel.f_model().data(), f_model_data)
  # test 2
  fmodel = mmtbx.f_model.manager(
    f_calc        = f_calc,
    f_mask        = f_mask,
    f_obs         = f_obs,
    k_mask        = k_mask,
    k_anisotropic = k_anisotropic,
    k_isotropic   = k_isotropic,
    r_free_flags  = r_free_flags,
    bin_selections = bin_selections,
    sf_and_grads_accuracy_params = sfg_params)
  assert approx_equal(fmodel.r_work(), 0)
  assert approx_equal(fmodel.r_free(), 0)
  assert approx_equal(fmodel.target_w(), 0)
  assert approx_equal(fmodel.k_masks()[0], k_mask)
  assert approx_equal(fmodel.k_anisotropic(), k_anisotropic)
  assert approx_equal(fmodel.f_model_scaled_with_k1().data(), f_model_data)
  assert approx_equal(fmodel.f_model().data(), f_model_data)
  assert fmodel.f_calc().data().all_eq(f_calc.data())
  assert fmodel.f_masks()[0].data().all_eq(f_mask.data())
  # test 3
  params = bss.master_params.extract()
  params.number_of_macro_cycles=5
  fmodel = mmtbx.f_model.manager(
    f_calc       = f_calc,
    f_mask       = f_mask,
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    bin_selections = bin_selections,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.3
  o = fmodel.update_all_scales(fast=False, params=params, remove_outliers=False)
  assert approx_equal(o.k_sol[0], k_sol,  0.01 )
  assert approx_equal(o.b_sol[0], b_sol,  0.1)
  assert approx_equal(o.b_cart,  b_cart,  1.e-1)
  assert approx_equal(fmodel.r_work(), 0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0, 1.e-3)
  # test 4 - part 1
  fmodel = mmtbx.f_model.manager(
    f_calc       = f_calc,
    f_mask       = f_mask,
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    bin_selections = bin_selections,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.3
  fmodel.update_all_scales(fast=True, remove_outliers=False)
  assert fmodel.r_work() < 0.045, [fmodel.r_work(), d_min]
  # part 2 of test 4
  d=fmodel.k_isotropic()*fmodel.k_anisotropic()*(
    f_calc.data()+fmodel.k_masks()[0]*f_mask.data())
  fmodel = mmtbx.f_model.manager(
    f_calc       = f_calc,
    f_mask       = f_mask,
    f_obs        = abs(f_calc.customized_copy(data = d)),
    r_free_flags = r_free_flags,
    bin_selections = bin_selections,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.3
  fmodel.update_all_scales(fast=True)
  assert fmodel.r_work() < 0.052, [fmodel.r_work(), d_min]

def exercise_top_largest_f_obs_f_model_differences(threshold_percent=10,
      symbol = "C 2"):
  random.seed(0)
  flex.set_random_seed(0)
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(symbol=symbol),
    elements               =(("O","N","C")*5+("H",)*10),
    volume_per_atom        = 200,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  f_obs = abs(x.structure_factors(d_min = 1.5, algorithm="fft").f_calc())
  x.shake_sites_in_place(mean_distance=1)
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs)
  fmodel.update_all_scales(update_f_part1=False)
  v, d = fmodel.top_largest_f_obs_f_model_differences(
    threshold_percent=threshold_percent)
  n = (d>v).count(True)*100./d.size()
  assert approx_equal(threshold_percent, n, 0.5)

def exercise_f_model_no_scales(symbol = "C 2"):
  random.seed(0)
  flex.set_random_seed(0)
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(symbol=symbol),
    elements               =(("O","N","C")*5+("H",)*10),
    volume_per_atom        = 200,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  f_obs = abs(x.structure_factors(d_min = 1.5, algorithm="fft").f_calc())
  x.shake_sites_in_place(mean_distance=1)
  k_iso = flex.double(f_obs.data().size(), 2)
  k_aniso = flex.double(f_obs.data().size(), 3)
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    k_isotropic    = k_iso,
    k_anisotropic  = k_aniso,
    f_obs          = f_obs)
  fc = abs(fmodel.f_calc()).data()
  fm = abs(fmodel.f_model()).data()
  fmns = abs(fmodel.f_model_no_scales()).data()
  assert approx_equal(flex.mean(fm/fc), 6)
  assert approx_equal(flex.mean(fmns/fc), 1)

def exercise_6_instantiate_consistency(symbol = "C 2"):
  random.seed(0)
  flex.set_random_seed(0)
  for scale in [1.e-4, 1.0, 1.e+4]:
    for k_sol in [0, 0.3]:
      for b_sol in [0, 50]:
        for set_h_occ_to_zero in [True, False]:
          for update_f_part1 in [True, False]:
            for apply_scale_to in ["f_obs", "f_model"]:
              # Simulate Fobs START
              x = random_structure.xray_structure(
                space_group_info       = sgtbx.space_group_info(symbol=symbol),
                elements               =(("O","N","C")*3+("H",)*10),
                volume_per_atom        = 50,
                min_distance           = 3,
                general_positions_only = True,
                random_u_iso           = True,
                random_occupancy       = False)
              x.scattering_type_registry(table="wk1995")
              x.set_occupancies(value=0.8, selection = x.hd_selection())
              f_calc = x.structure_factors(d_min = 2.0).f_calc()
              mask_manager = mmtbx.masks.manager(miller_array = f_calc)
              f_mask = mask_manager.shell_f_masks(xray_structure = x)[0]
              assert flex.mean(abs(f_mask).data()) > 0
              b_cart=adptbx.random_traceless_symmetry_constrained_b_cart(
                crystal_symmetry=x.crystal_symmetry())
              u_star = adptbx.u_cart_as_u_star(x.unit_cell(), adptbx.b_as_u(b_cart))
              k_anisotropic = mmtbx.f_model.ext.k_anisotropic(f_calc.indices(), u_star)
              ss = 1./flex.pow2(f_calc.d_spacings().data()) / 4.
              k_mask = mmtbx.f_model.ext.k_mask(ss, k_sol, b_sol)
              if(apply_scale_to=="f_model"):
                k_isotropic = flex.double(f_calc.data().size(), scale)
              else:
                k_isotropic = flex.double(f_calc.data().size(), 1)
              f_model_data = scale*k_anisotropic*(f_calc.data()+k_mask*f_mask.data())
              f_model = f_calc.customized_copy(data = f_model_data)
              f_obs = abs(f_model)
              if(apply_scale_to=="f_obs"):
                f_obs = f_obs.customized_copy(data = f_obs.data()*scale)
              r_free_flags = f_obs.generate_r_free_flags()
              # Simulate Fobs END
              if(set_h_occ_to_zero):
                x.set_occupancies(value=0.0, selection = x.hd_selection())
              x.shake_sites_in_place(mean_distance=5)
              sel = x.random_remove_sites_selection(fraction=0.3)
              x = x.select(sel)
              fmodel = mmtbx.f_model.manager(
                xray_structure = x,
                f_obs          = f_obs,
                r_free_flags   = r_free_flags)
              fmodel.update_all_scales(fast=True, show=False,
                update_f_part1=update_f_part1)
              f_part1_data = fmodel.f_calc().data()*flex.random_double(
                fmodel.f_calc().data().size())
              f_part1 = fmodel.f_calc().customized_copy(data = f_part1_data)
              fmodel.update(f_part1 = f_part1)
              r1 = fmodel.r_work()
              #
              zero=fmodel.f_calc().customized_copy(data=fmodel.f_calc().data()*0)
              fmodel_dc  = mmtbx.f_model.manager(
                f_obs         = fmodel.f_obs(),
                r_free_flags  = fmodel.r_free_flags(),
                k_isotropic   = fmodel.k_isotropic(),
                k_anisotropic = fmodel.k_anisotropic(),
                f_calc        = fmodel.f_model_no_scales(),
                f_part1       = fmodel.f_part1(),
                f_part2       = fmodel.f_part2(),
                f_mask        = zero)
              r2 = fmodel_dc.r_work()
              if(0):
                print("r1=%8.6f r2=%8.6f fp1=%6.3f fp2=%6.3f fc=%6.3f"%(r1, r2,
                  flex.mean(abs(fmodel.f_part1()).data()), \
                  flex.mean(abs(fmodel.f_part2()).data()), \
                  flex.mean(abs(fmodel.f_calc()).data())), \
                  "set_h_occ_to_zero=", set_h_occ_to_zero,\
                  "update_f_part1=", update_f_part1)
              assert approx_equal(r1, r2), [r1, r2]

def run():
  for d_min in [2, 4]:
    exercise_5_bulk_sol_and_scaling(d_min = d_min)
  exercise_6_instantiate_consistency()
  exercise_f_model_no_scales()
  exercise_top_largest_f_obs_f_model_differences()
  exercise_1()
  exercise_2()
  exercise_4_f_hydrogens()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
