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
                if(xrs is not None): assert (fmodel.estimate_f000() > 0.0)
                #
                # update k_isotropic, k_anisotorpic, resolution_filter
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

def t_exercise_3_1(x, f_obs, fc1, fc2, r_free_flags, sfg_params):
  #
  assert not fc1.data().all_eq(fc2.data())
  #
  for case in [1,2,3]:
    print "case:",case
    fmodel = mmtbx.f_model.manager(
      xray_structure = x,
      f_obs          = f_obs,
      r_free_flags   = r_free_flags,
      sf_and_grads_accuracy_params = sfg_params)
    r1 = fmodel.r_work()
    #
    if(case==1):
      fmodel2 = mmtbx.f_model.manager(
        xray_structure = x,
        f_obs          = f_obs,
        f_part1        = fc1,
        r_free_flags   = r_free_flags,
        sf_and_grads_accuracy_params = sfg_params)
    elif(case==2):
      fmodel2 = fmodel.deep_copy()
      fmodel2.update_core(f_part1=fc1)
    elif(case==3):
      fmodel2 = fmodel.deep_copy()
      fmodel2.update_core(f_part1=fc1)
      fmodel2 = fmodel2.deep_copy()
      fmodel2.update_xray_structure(xray_structure=x, update_f_calc=True)
    assert fmodel2.f_part1().data().all_eq(fc1.data())
    r2 = fmodel2.r_work()
    #
    if(case==1):
      fmodel3 = mmtbx.f_model.manager(
        xray_structure = x,
        f_obs          = f_obs,
        f_part2        = fc1,
        r_free_flags   = r_free_flags,
        sf_and_grads_accuracy_params = sfg_params)
    elif(case==2):
      fmodel3 = fmodel.deep_copy()
      fmodel3.update_core(f_part2=fc1)
    elif(case==3):
      fmodel3 = fmodel.deep_copy()
      fmodel3.update_core(f_part2=fc1)
      fmodel3 = fmodel3.deep_copy()
      fmodel3.update_xray_structure(update_f_calc=True)
    assert fmodel3.f_part2().data().all_eq(fc1.data())
    r3 = fmodel3.r_work()
    #
    if(case==1):
      fmodel4 = mmtbx.f_model.manager(
        xray_structure = x,
        f_obs          = f_obs,
        f_part1        = fc1,
        f_part2        = fc2,
        r_free_flags   = r_free_flags,
        sf_and_grads_accuracy_params = sfg_params)
    elif(case==2):
      fmodel4 = fmodel.deep_copy()
      fmodel4.update_core(f_part1=fc1, f_part2=fc2)
    elif(case==3):
      fmodel4 = fmodel.deep_copy()
      fmodel4.update_core(f_part1=fc1, f_part2=fc2)
      fmodel4 = fmodel4.deep_copy()
      fmodel2.update_xray_structure(xray_structure=x, update_f_calc=True)
    assert fmodel4.f_part1().data().all_eq(fc1.data())
    assert fmodel4.f_part2().data().all_eq(fc2.data())
    r4 = fmodel4.r_work()
    #
    assert r1>0.3
    assert approx_equal(r2, r3)
    assert r2 < r1
    assert approx_equal(r4, 0)

def exercise_3_f_part1_and_f_part2():
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P 4"),
    elements               =(("O","N","C")*10),
    volume_per_atom        = 200,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  x.scattering_type_registry(table="wk1995")
  fc = x.structure_factors(d_min = 2.0, algorithm="direct").f_calc()
  #
  x1 = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P 4"),
    unit_cell              = x.unit_cell(),
    elements               =(("C")*10),
    volume_per_atom        = 200,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  x1.scattering_type_registry(table="wk1995")
  fc1 = x1.structure_factors(d_min = 2.0, algorithm="direct").f_calc()
  #
  x2 = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P 4"),
    unit_cell              = x.unit_cell(),
    elements               =(("C")*10),
    volume_per_atom        = 200,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  x2.scattering_type_registry(table="wk1995")
  fc2 = x2.structure_factors(d_min = 2.0, algorithm="direct").f_calc()
  #
  fc_all = fc.customized_copy(data = fc.data()+fc1.data()+fc2.data())
  f_obs = abs(fc_all.deep_copy())
  r_free_flags = f_obs.generate_r_free_flags(fraction = 0.1)
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg_params.algorithm = "direct"
  #
  t_exercise_3_1(x=x, f_obs=f_obs, fc1=fc1, fc2=fc2, r_free_flags=r_free_flags,
    sfg_params=sfg_params)

def exercise_4_f_hydrogens():
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P 4"),
    elements               =(("O","N","C")*5 + ("H",)*95),
    volume_per_atom        = 200,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  x.scattering_type_registry(table="wk1995")
  x.set_occupancies(value=0.9, selection = x.hd_selection())
  fc = x.structure_factors(d_min = 2.5, algorithm="direct").f_calc()
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
  assert fmodel.r_work() > 0.15
  fmodel.update_f_hydrogens()
  assert approx_equal(fmodel.k_h, 0.9)
  assert approx_equal(fmodel.b_h, 0)
  assert approx_equal(fmodel.r_work(), 0)

def exercise_1():
  n_elements = 70
  sgs = ["P 1", "P 4", "C 1 2/c 1"]
  for sg in sgs:
    print sg
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

def get_random_and_traceless_b_cart(crystal_symmetry):
  symbol = crystal_symmetry.space_group().type().lookup_symbol()
  point_group = sgtbx.space_group_info(
    symbol=symbol).group().build_derived_point_group()
  adp_constraints = sgtbx.tensor_rank_2_constraints(
    space_group=point_group,
    reciprocal_space=True)
  u_star = adptbx.u_cart_as_u_star(crystal_symmetry.unit_cell(),
    adptbx.random_u_cart(u_scale=1,u_min=0.1))
  u_indep = adp_constraints.independent_params(all_params=u_star)
  u_star = adp_constraints.all_params(independent_params=u_indep)
  b_cart = adptbx.u_as_b(adptbx.u_star_as_u_cart(
    crystal_symmetry.unit_cell(), u_star))
  tr = (b_cart[0]+b_cart[1]+b_cart[2])/3
  b_cart = [b_cart[0]-tr, b_cart[1]-tr, b_cart[2]-tr,
           b_cart[3],b_cart[4],b_cart[5]]
  return b_cart

def exercise_5_bulk_sol_and_scaling(symbol = "C 2"):
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(symbol=symbol),
    elements               =(("O","N","C")*150),
    volume_per_atom        = 200,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  x.scattering_type_registry(table="wk1995")
  f_calc = x.structure_factors(d_min = 2.5, algorithm="direct").f_calc()
  mask_manager = mmtbx.masks.manager(miller_array = f_calc)
  f_mask = mask_manager.shell_f_masks(xray_structure = x)[0]
  assert flex.mean(abs(f_mask).data()) > 0
  b_cart=[-15,-5,20, 0,6,0]
  #b_cart=get_random_and_traceless_b_cart(crystal_symmetry=x.crystal_symmetry())
  u_star = adptbx.u_cart_as_u_star(x.unit_cell(), adptbx.b_as_u(b_cart))
  k_anisotropic = mmtbx.f_model.ext.k_anisotropic(f_calc.indices(), u_star)
  ss = 1./flex.pow2(f_calc.d_spacings().data()) / 4.
  k_mask = mmtbx.f_model.ext.k_mask(ss, 0.37, 64.0)
  scale = 17.
  k_isotropic = flex.double(f_calc.data().size(), scale)
  f_model_data = scale*k_anisotropic*(f_calc.data()+k_mask*f_mask.data())
  f_model = f_calc.customized_copy(data = f_model_data)
  f_obs = abs(f_model)
  r_free_flags = f_obs.generate_r_free_flags()
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg_params.algorithm = "direct"
  # test 1
  for k_isotropic_ in [None, k_isotropic]:
    fmodel = mmtbx.f_model.manager(
      xray_structure = x,
      f_obs          = f_obs,
      k_mask         = k_mask,
      k_anisotropic  = k_anisotropic,
      k_isotropic    = k_isotropic,
      r_free_flags   = r_free_flags,
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
  fmodel = mmtbx.f_model.manager(
    f_calc       = f_calc,
    f_mask       = f_mask,
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.3
  fmodel.update_solvent_and_scale(fast=False)
  assert approx_equal(fmodel.r_work(), 0)
  assert approx_equal(fmodel.r_free(), 0)
  assert approx_equal(fmodel.k_masks()[0], k_mask, 1.e-4)
  assert approx_equal(fmodel.k_anisotropic(), k_anisotropic, 1.e-4)
  # test 4 - part 1
  fmodel = mmtbx.f_model.manager(
    f_calc       = f_calc,
    f_mask       = f_mask,
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.3
  fmodel.update_solvent_and_scale(fast=True)
  assert fmodel.r_work() < 0.015
  # part 2 of test 4
  d=fmodel.k_isotropic()*fmodel.k_anisotropic()*(
    f_calc.data()+fmodel.k_masks()[0]*f_mask.data())
  fmodel = mmtbx.f_model.manager(
    f_calc       = f_calc,
    f_mask       = f_mask,
    f_obs        = abs(f_calc.customized_copy(data = d)),
    r_free_flags = r_free_flags,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.3
  fmodel.update_solvent_and_scale(fast=True)
  assert fmodel.r_work() < 0.005

def exercise_5_bulk_sol_and_scaling_and_H(symbol = "C 2"):
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
  x.scattering_type_registry(table="wk1995")
  x.set_occupancies(value=0.9, selection = x.hd_selection())
  f_calc = x.structure_factors(d_min = 1.5, algorithm="direct").f_calc()
  mask_manager = mmtbx.masks.manager(miller_array = f_calc)
  f_mask = mask_manager.shell_f_masks(xray_structure = x)[0]
  assert flex.mean(abs(f_mask).data()) > 0
  b_cart=get_random_and_traceless_b_cart(crystal_symmetry=x.crystal_symmetry())
  u_star = adptbx.u_cart_as_u_star(x.unit_cell(), adptbx.b_as_u(b_cart))
  k_anisotropic = mmtbx.f_model.ext.k_anisotropic(f_calc.indices(), u_star)
  ss = 1./flex.pow2(f_calc.d_spacings().data()) / 4.
  k_mask = mmtbx.f_model.ext.k_mask(ss, 0.37, 64.0)
  scale = 17.
  k_isotropic = flex.double(f_calc.data().size(), scale)
  f_model_data = scale*k_anisotropic*(f_calc.data()+k_mask*f_mask.data())
  f_model = f_calc.customized_copy(data = f_model_data)
  f_obs = abs(f_model)
  r_free_flags = f_obs.generate_r_free_flags()
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg_params.algorithm = "direct"
  x.set_occupancies(value=0.0, selection = x.hd_selection())
  # test 1
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs,
    k_mask         = k_mask,
    k_anisotropic  = k_anisotropic,
    k_isotropic    = k_isotropic,
    r_free_flags   = r_free_flags,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.06, fmodel.r_work()
  fmodel.update_f_hydrogens()
  assert approx_equal(fmodel.k_h, 0.9)
  assert approx_equal(fmodel.b_h, 0)
  assert approx_equal(fmodel.r_work(), 0)
  # test 2
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.25
  fmodel.update_all_scales(cycles=6, fast=False)
  assert approx_equal(fmodel.k_h, 0.9)
  assert approx_equal(fmodel.b_h, 0)
  assert approx_equal(fmodel.r_work(), 0)
  # test 3
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    sf_and_grads_accuracy_params = sfg_params)
  assert fmodel.r_work() > 0.25
  fmodel.update_all_scales(cycles=6, fast=True, show=False)
  assert approx_equal(fmodel.k_h, 0.9)
  assert approx_equal(fmodel.b_h, 0)
  assert fmodel.r_work() < 0.01

def run(args):
  assert len(args) == 0
  exercise_5_bulk_sol_and_scaling()
  exercise_5_bulk_sol_and_scaling_and_H()
  exercise_1()
  exercise_2()
  exercise_3_f_part1_and_f_part2()
  exercise_4_f_hydrogens()
  print format_cpu_times()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
