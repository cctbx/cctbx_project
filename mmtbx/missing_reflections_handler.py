from cctbx import miller
from cctbx.array_family import flex
import mmtbx.f_model
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")
from libtbx import group_args

def compute_new_f_obs(fmodel):
  from mmtbx.max_lik import maxlik
  csc = fmodel.complete_set_core()
  f_obs = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags()
  f_model_scaled_with_k1 = csc.f_model.array(
    data = csc.f_model.data()*fmodel.scale_k1())
  f_obs_lone = abs(f_model_scaled_with_k1.lone_set(other = f_obs))
  new_f_obs = f_obs.concatenate(other = f_obs_lone)
  new_f_obs, f_model_scaled_with_k1 = new_f_obs.common_sets(f_model_scaled_with_k1)
  r_free_flags, dummy = r_free_flags.common_sets(f_model_scaled_with_k1)
  new_r_free_flags = r_free_flags.data().concatenate(
    flex.bool(f_obs_lone.size(), False))
  alpha, beta = maxlik.alpha_beta_est_manager(
    f_obs                    = new_f_obs,
    f_calc                   = f_model_scaled_with_k1,
    free_reflections_per_bin = 100,
    flags                    = new_r_free_flags,
    interpolation            = True).alpha_beta()
  apply_alpha_sel = flex.bool(f_obs.size(), False).concatenate(
    flex.bool(f_obs_lone.size(), True))
  alpha = alpha.select(apply_alpha_sel)
  f_obs_lone = f_obs_lone.array(data=f_obs_lone.data()*alpha.data())
  new_f_obs = f_obs.concatenate(other = f_obs_lone)
  new_abcd = None
  if(fmodel.hl_coeffs() is not None):
    new_abcd = fmodel.hl_coeffs().customized_copy(
      indices = new_f_obs.indices(),
      data = fmodel.hl_coeffs().data().concatenate(
        flex.hendrickson_lattman(f_obs_lone.size(), [0,0,0,0])))
  return group_args(
    f_obs = new_f_obs,
    r_free_flags = new_f_obs.array(data=new_r_free_flags),
    filled_f_obs_selection = apply_alpha_sel,
    hl_coeffs = new_abcd)

def fill_missing_f_obs(fmodel, fill_mode):
  assert fill_mode in ["fobs_mean_mixed_with_dfmodel",
                       "random",
                       "fobs_mean",
                       "dfmodel"]
  from mmtbx import masks
  use_f_part = fmodel.k_part() > 0
  bss_params = bss.master_params.extract()
  bss_params.k_sol_b_sol_grid_search = False
  bss_params.b_sol_max = 150.0
  bss_params.number_of_macro_cycles = 1

  new_f_obs_obj = compute_new_f_obs(fmodel)
  new_f_obs = new_f_obs_obj.f_obs
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
    r_free_flags           = new_f_obs_obj.r_free_flags,
    target_name            = fmodel.target_name,
    f_obs                  = new_f_obs,
    abcd                   = new_f_obs_obj.hl_coeffs,
    k_sol                  = fmodel.k_sol(),
    b_sol                  = fmodel.b_sol(),
    b_cart                 = fmodel.b_cart(),
    mask_params            = fmodel.mask_params,
    filled_f_obs_selection = new_f_obs_obj.filled_f_obs_selection)
  fmodel_result.update_solvent_and_scale(params = bss_params,
    optimize_mask=False)
  if(use_f_part): fmodel_result.update_f_part()
  return fmodel_result
