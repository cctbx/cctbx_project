from cctbx.array_family import flex
import mmtbx.f_model
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")
from libtbx import group_args
from cctbx import miller
from cctbx import maptbx

def compute_new_f_obs(fmodel):
  from mmtbx.max_lik import maxlik
  csc = fmodel.complete_set_core()
  f_obs = fmodel.f_obs()
  abcd = fmodel.hl_coeffs()
  r_free_flags = fmodel.r_free_flags()
  f_model_scaled_with_k1 = csc.f_model.array(
    data = csc.f_model.data()*fmodel.scale_k1())
  f_model_scaled_with_k1_lone, f_obs_lone = f_model_scaled_with_k1.lone_sets(f_obs)
  ###
  #C = abs(f_model_scaled_with_k1.deep_copy())
  #Cd, Ci = C.data(), C.indices()
  #O = f_obs.deep_copy()
  #for mi, d in zip(O.indices(), O.data()):
  #  for i, mic in enumerate(Ci):
  #    if(mi == mic):
  #      Cd[i]=d
  #      break
  #new_f_obs = miller.set(
  #  crystal_symmetry=C.crystal_symmetry(),
  #  indices = C.indices(),
  #  anomalous_flag=False).array(data=Cd)
  ###
  f_obs_lone = abs(f_model_scaled_with_k1_lone)
  new_f_obs = f_obs.concatenate(other = f_obs_lone)
  new_f_obs, f_model_scaled_with_k1 = new_f_obs.common_sets(f_model_scaled_with_k1)
  # XXX strangely, this doen't work
  #f_model_scaled_with_k1, new_f_obs = f_model_scaled_with_k1.common_sets(new_f_obs)
  r_free_flags, dummy = r_free_flags.common_sets(new_f_obs)
  new_r_free_flags = r_free_flags.data().concatenate(
    flex.bool(f_obs_lone.size(), False))
  #alpha, beta = maxlik.alpha_beta_est_manager(                  #C
  #  f_obs                    = new_f_obs,                       #C
  #  f_calc                   = f_model_scaled_with_k1,          #C
  #  free_reflections_per_bin = 100,                             #C
  #  flags                    = new_r_free_flags,                #C
  #  interpolation            = True).alpha_beta()               #C
  apply_alpha_sel = flex.bool(f_obs.size(), False).concatenate(
    flex.bool(f_obs_lone.size(), True))
  #alpha = alpha.select(apply_alpha_sel)                                     #C
  #f_obs_lone = f_obs_lone.array(data=f_obs_lone.data()*alpha.data())        #C
  #new_f_obs = f_obs.concatenate(other = f_obs_lone)                         #C
  #f_obs_lone = f_obs_lone.array(data=f_obs_lone.data()*85)        #C
  #new_f_obs = f_obs.concatenate(other = f_obs_lone)                         #C
  #
  new_abcd = None
  if(fmodel.hl_coeffs() is not None):
    abcd_lone = f_obs_lone.customized_copy(
      data = flex.hendrickson_lattman(f_obs_lone.size(), [0,0,0,0]))
    new_abcd = abcd.concatenate(other = abcd_lone)
    new_abcd, new_f_obs = new_abcd.common_sets(new_f_obs)
  #####
  rff_lone = f_obs_lone.customized_copy(data=flex.bool(f_obs_lone.size(),False))
  new_rff = r_free_flags.concatenate(other = rff_lone)
  new_rff, new_f_obs = new_rff.common_sets(new_f_obs)
  #####
  return group_args(
    f_obs = new_f_obs,
    r_free_flags = new_rff,
    filled_f_obs_selection = apply_alpha_sel,
    hl_coeffs = new_abcd)

def select_by_map_cc(fmodel):
  residue_selections = flex.bool(
    fmodel.xray_structure.scatterers().size(),True).iselection()
  #
  f_calc_1 = fmodel.f_model_scaled_with_k1()
  if (f_calc_1.anomalous_flag()) :
    f_calc_1 = f_calc_1.average_bijvoet_mates()
  fft_map_1 = f_calc_1.fft_map(resolution_factor=0.25)
  fft_map_1.apply_volume_scaling()
  map_1 = fft_map_1.real_map()#_unpadded()
  #
  coeffs = fmodel.electron_density_map(
    fill_missing_f_obs = False).map_coefficients(map_type = "2mFo-DFc")
  if (coeffs.anomalous_flag()) :
    coeffs = coeffs.average_bijvoet_mates()
  fft_map_2 = miller.fft_map(
    crystal_gridding = fft_map_1,
    fourier_coefficients = coeffs)
  fft_map_2.apply_volume_scaling()
  map_2 = fft_map_2.real_map()#_unpadded()
  #
  if(fmodel.f_obs().d_min()<=2.5): rad = 1.5
  else: rad = 2.0
  assert map_1.size() == map_2.size()
  assert map_1.focus() == map_2.focus()
  assert map_1.all() == map_2.all()
  sites_cart = fmodel.xray_structure.sites_cart()
  result = flex.double()
  for rsel in residue_selections:
    rsel = flex.size_t([rsel])
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = fmodel.xray_structure.unit_cell(),
      fft_n_real = map_1.focus(),
      fft_m_real = map_1.all(),
      sites_cart = sites_cart.select(rsel),
      site_radii = flex.double(rsel.size(), rad))
    m1 = map_1.select(sel)
    m2 = map_2.select(sel)
    assert m1.size() == m2.size()
    corr = flex.linear_correlation(x = m1, y = m2).coefficient()
    result.append(corr)
  selection = result >= 0.9
  #print "Atoms:", fmodel.xray_structure.scatterers().size()
  #print "Atoms kept:", selection.count(True)
  keep = 100.*selection.count(True)/selection.size()
  if(keep>30.):
    xrs = fmodel.xray_structure.select(selection = selection)
    fmodel.update_xray_structure(xray_structure = xrs, update_f_calc = True,
      update_f_mask = True)
  return fmodel


def fill_missing_f_obs(fmodel, fill_mode, update_scaling=True):
  assert fill_mode in ["fobs_mean_mixed_with_dfmodel",
                       "random",
                       "fobs_mean",
                       "dfmodel"]
  use_f_part = flex.max(abs(fmodel.f_part1()).data()) > 0
  bss_params = bss.master_params.extract()
  bss_params.k_sol_b_sol_grid_search = False
  bss_params.b_sol_max = 150.0
  bss_params.number_of_macro_cycles = 1
  new_f_obs_obj = compute_new_f_obs(
    fmodel = select_by_map_cc(fmodel = fmodel.deep_copy()))
  new_f_obs = new_f_obs_obj.f_obs
  new_r_free_flags = new_f_obs_obj.r_free_flags
  filled_f_obs_selection = new_f_obs_obj.filled_f_obs_selection
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
  assert new_f_obs.indices().all_eq(new_r_free_flags.indices())
  fmodel_result = mmtbx.f_model.manager(
    xray_structure         = fmodel.xray_structure,
    r_free_flags           = new_r_free_flags,
    target_name            = fmodel.target_name,
    f_obs                  = new_f_obs,
    abcd                   = new_f_obs_obj.hl_coeffs, # ??? This is the problem!
    k_sol                  = fmodel.k_sols(),
    b_sol                  = fmodel.b_sol(),
    b_cart                 = fmodel.b_cart(),
    mask_params            = fmodel.mask_params)
  if(update_scaling):
    fmodel_result.update_solvent_and_scale(params = bss_params,
      optimize_mask=False)
  if(use_f_part): fmodel_result.update_f_part1()
  return fmodel_result
