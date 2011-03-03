from cctbx import miller
from cctbx.array_family import flex
import mmtbx.f_model
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")

def fill_missing_f_obs(fmodel, fill_mode):
  assert fill_mode in ["fobs_mean_mixed_with_dfmodel",
                       "random",
                       "fobs_mean",
                       "dfmodel"]
  from mmtbx import masks
  from mmtbx.max_lik import maxlik
  use_f_part = fmodel.k_part() > 0
  bss_params = bss.master_params.extract()
  bss_params.k_sol_b_sol_grid_search = False
  bss_params.b_sol_max = 150.0
  bss_params.number_of_macro_cycles = 1
  f_model = fmodel.f_model()
  n_refl_orig = f_model.data().size()
  complete_set = f_model.complete_set(d_min = f_model.d_min(), d_max=None)
  f_calc_atoms = complete_set.structure_factors_from_scatterers(
    xray_structure = fmodel.xray_structure).f_calc()
  f_calc_atoms_lone = f_calc_atoms.lone_set(other = f_model)
  n_refl_lone = f_calc_atoms_lone.data().size()
  if(fmodel.k_sol()>0):
    f_mask = masks.manager(
      miller_array   = f_calc_atoms,
      mask_params    = fmodel.mask_params,
      xray_structure = fmodel.xray_structure).f_mask()
  else:
    f_mask = f_calc_atoms.array(data = f_calc_atoms.data()*0)
  f_mask_lone = f_mask.lone_set(other = f_model)
  ss = 1./flex.pow2(f_mask_lone.d_spacings().data())/4.
  r_free_flags_lone = f_mask_lone.array(
    data = flex.bool(f_mask_lone.size(), False))
  f_model_core = ext.core(
    f_calc      = f_calc_atoms_lone.data(),
    f_mask      = f_mask_lone.data(),
    k_sol       = fmodel.k_sol(),
    b_sol       = fmodel.b_sol(),
    f_part_base = f_calc_atoms_lone.data()*0,
    k_part      = 0,
    b_part      = 0,
    u_star      = fmodel.u_star(),
    hkl         = f_calc_atoms_lone.indices(),
    uc          = f_mask_lone.unit_cell(),
    ss          = ss)
  f_obs_orig = fmodel.f_obs().deep_copy()
  r_free_flags_orig = fmodel.r_free_flags()
  # compose new fileld fmodel
  filled_f_obs_selection = flex.bool(n_refl_orig, False)
  f_model_lone = abs(miller.array(
    miller_set = f_mask_lone,
    data       = f_model_core.f_model * fmodel.scale_k1()))
  new_f_obs = fmodel.f_obs().concatenate(other = f_model_lone)
  new_r_free_flags = fmodel.r_free_flags().concatenate(
    other = r_free_flags_lone)
  filled_f_obs_selection = filled_f_obs_selection.concatenate(
    flex.bool(n_refl_lone, True))
  new_abcd = None
  if(fmodel.hl_coeffs() is not None):
    new_abcd = fmodel.hl_coeffs().customized_copy(
      indices = new_f_obs.indices(),
      data = fmodel.hl_coeffs().data().concatenate(
        flex.hendrickson_lattman(n_refl_lone, [0,0,0,0])))
  fmodel = mmtbx.f_model.manager(
    xray_structure = fmodel.xray_structure,
    r_free_flags   = new_r_free_flags,
    target_name    = "ls_wunit_k1",
    f_obs          = new_f_obs,
    abcd           = new_abcd,
    mask_params    = fmodel.mask_params,
    k_sol          = fmodel.k_sol(),
    b_sol          = fmodel.b_sol(),
    b_cart         = fmodel.b_cart())
  fmodel.update_solvent_and_scale(params = bss_params, optimize_mask=False)
  if(use_f_part): fmodel.update_f_part()
  # replace 'F_obs' -> alpha * 'F_obs' for filled F_obs
  alpha, beta = maxlik.alpha_beta_est_manager(
    f_obs                    = fmodel.f_obs(),
    f_calc                   = fmodel.f_model_scaled_with_k1(),
    free_reflections_per_bin = 100,
    flags                    = fmodel.r_free_flags().data(),
    interpolation            = True).alpha_beta()
  apply_alpha_sel = flex.bool(n_refl_orig, False).concatenate(
    flex.bool(n_refl_lone, True)) # assume order did not change
  assert apply_alpha_sel.size() == fmodel.f_obs().data().size()
  alpha = alpha.select(apply_alpha_sel)
  # compose new fileld fmodel
  f_model_lone = abs(miller.array(
    miller_set = f_mask_lone,
    data       = f_model_core.f_model * fmodel.scale_k1()*alpha.data()))
  new_f_obs = f_obs_orig.concatenate(other = f_model_lone)
  new_r_free_flags = r_free_flags_orig.concatenate(
    other = r_free_flags_lone)
  # XXX implement and use fmodel.customized_copy() instead of creating a new one
  new_f_obs.set_observation_type_xray_amplitude()
  assert new_f_obs.data().size() == filled_f_obs_selection.size()
  #
  if(fill_mode == "fobs_mean"):
    sel = new_f_obs.sort_permutation(by_value = "resolution")
    new_f_obs = new_f_obs.select(selection = sel)
    new_r_free_flags = new_r_free_flags.select(selection = sel)
    new_data = flex.double()
    i_max = new_f_obs.data().size()
    d_spacings = new_f_obs.d_spacings().data()
    new_f_obs_data = new_f_obs.data()
    for i_seq, fo in enumerate(new_f_obs_data):
      if(filled_f_obs_selection[i_seq]):
        #
        x = flex.double()
        y = flex.double()
        i = i_seq
        counter = 0
        while True:
          i +=1
          if i > i_max-1: break
          if(not filled_f_obs_selection[i]):
            x.append(d_spacings[i])
            y.append(new_f_obs_data[i])
            counter += 1
          if(counter == 5): break
        #
        i = i_seq
        counter = 0
        while True:
          i -=1
          if i < 0: break
          if(not filled_f_obs_selection[i]):
            x.append(d_spacings[i])
            y.append(new_f_obs_data[i])
            counter += 1
          if(counter == 5): break
        #
        assert y.size() > 0 and x.size() == y.size()
        assert x.size() <= 10, x.size()
        new_data.append( flex.mean(y) )
      else:
        new_data.append(fo)
    new_f_obs._data = new_data
  if(fill_mode == "random"):
    new_data = flex.double()
    for i_seq, fo in enumerate(new_f_obs.data()):
      if(filled_f_obs_selection[i_seq]):
        if(fo > 1.):
          new_data.append( random.randrange(int(fo-fo/2),int(fo+fo/2)) )
        else: new_data.append( abs(random.random()) )
      else:
        new_data.append(fo)
    new_f_obs._data = new_data
  if(fill_mode in ["fobs_interpolated_mixed_with_dfmodel",
                   "fobs_mean_mixed_with_dfmodel"]):
    sel = new_f_obs.sort_permutation(by_value = "resolution")
    new_f_obs = new_f_obs.select(selection = sel)
    new_r_free_flags = new_r_free_flags.select(selection = sel)
    new_data = flex.double()
    i_max = new_f_obs.data().size()
    d_spacings = new_f_obs.d_spacings().data()
    new_f_obs_data = new_f_obs.data()
    for i_seq, fo in enumerate(new_f_obs_data):
      if(filled_f_obs_selection[i_seq]):
        #
        x = flex.double()
        y = flex.double()
        i = i_seq
        counter = 0
        d_i_seq = d_spacings[i_seq]
        while True:
          i +=1
          if i > i_max-1: break
          if(not filled_f_obs_selection[i]):
            x.append(d_spacings[i])
            y.append(new_f_obs_data[i])
            counter += 1
          if(counter == 5): break
          if(abs(d_i_seq-d_spacings[i]) > 0.1): break
        #
        i = i_seq
        counter = 0
        while True:
          i -=1
          if i < 0: break
          if(not filled_f_obs_selection[i]):
            x.append(d_spacings[i])
            y.append(new_f_obs_data[i])
            counter += 1
          if(counter == 5): break
          if(abs(d_i_seq-d_spacings[i]) > 0.1): break
        #
        assert x.size() <= 10, x.size()
        #
        if(x.size() < 10):
          i = i_seq
          j = i_seq
          while True:
            i +=1
            j -=1
            if(x.size() >= 10): break
            if((i <= i_max-1 and i >= 0) and filled_f_obs_selection[i] and
               abs(d_i_seq-d_spacings[i]) < 0.1):
              x.append(d_spacings[i])
              y.append(new_f_obs_data[i])
            if(x.size() >= 10): break
            if((j <= i_max-1 and j >= 0) and filled_f_obs_selection[j] and
               abs(d_i_seq-d_spacings[j]) < 0.1):
              x.append(d_spacings[j])
              y.append(new_f_obs_data[j])
        #
        assert y.size() > 0 and x.size() == y.size()
        assert x.size() == 10, x.size()
        new_data.append( flex.mean(y) )
      else:
        new_data.append(fo)
    new_f_obs._data = new_data
  #
  fmodel_result = mmtbx.f_model.manager(
    xray_structure         = fmodel.xray_structure,
    r_free_flags           = new_r_free_flags,
    target_name            = fmodel.target_name,
    f_obs                  = new_f_obs,
    abcd                   = new_abcd,
    k_sol                  = fmodel.k_sol(),
    b_sol                  = fmodel.b_sol(),
    b_cart                 = fmodel.b_cart(),
    mask_params            = fmodel.mask_params,
    filled_f_obs_selection = filled_f_obs_selection)
  fmodel_result.update_solvent_and_scale(params = bss_params,
    optimize_mask=False)
  if(use_f_part): fmodel.update_f_part()
  return fmodel_result
