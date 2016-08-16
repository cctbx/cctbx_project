from __future__ import division
from cctbx.array_family import flex
import mmtbx.f_model
from cctbx import maptbx
from libtbx.test_utils import approx_equal
import mmtbx.masks
import mmtbx.bulk_solvent

def multi_mask_bulk_solvent(fmodel, log=None):
  # commonly used objects
  f_calc       = fmodel.f_calc()
  f_obs        = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags()
  ss = 1./flex.pow2(f_obs.d_spacings().data()) / 4.
  # crystal_gridding
  resolution_factor = 1./5
  grid_step = f_obs.d_min()*resolution_factor
  if(grid_step < 0.15): grid_step = 0.15
  grid_step = min(0.8, grid_step)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell = f_obs.unit_cell(),
    space_group_info = f_obs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = grid_step)
  n_real = crystal_gridding.n_real()
  # compute ASU mask
  mask_data_asu = mmtbx.masks.mask_from_xray_structure(
    xray_structure        = fmodel.xray_structure,
    p1                    = True,
    for_structure_factors = True,
    n_real                = n_real,
    in_asu                = True).mask_data
  # mask connectivity analysis
  co = maptbx.connectivity(map_data=mask_data_asu, threshold=0.9)
  conn = co.result()
  conn.reshape(mask_data_asu.accessor()) #XXX should be part of connectivity?
  z = zip(co.regions(),range(0,co.regions().size()))
  sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
  cntr=0
  f_masks = []
  accumulate = False
  m_final = 0
  for p in sorted_by_volume:
    v, i = p
    if(i == 0): continue # exclude 1 region
    s = conn==i
    #XXX this is 3 loops, may be slow. move to C++ if slow.
    mask_data_asu_i = mask_data_asu.set_selected(s, 1).set_selected(~s, 0)
    fr = mask_data_asu_i.count(1)*100./mask_data_asu_i.size()
    if(log is not None):
      print >> log, "region: %5d volume: %5d fraction: %6.2f  counter: %5d"%(
        i, v, fr, cntr)
      log.flush()
    if(fr>1.):
      m_final = m_final + mask_data_asu_i
      accumulate = True
    else:
      if(accumulate):
        f_mask_i = f_calc.structure_factors_from_asu_map(
          asu_map_data = m_final, n_real = n_real)
        accumulate = False
      else:
        f_mask_i = f_calc.structure_factors_from_asu_map(
          asu_map_data = mask_data_asu_i, n_real = n_real)
      f_masks.append(f_mask_i)
    cntr += 1 # MUST BE LAST
  if(len(f_masks)==0): # just one region
    f_mask_i = f_calc.structure_factors_from_asu_map(
      asu_map_data = m_final, n_real = n_real)
    f_masks.append(f_mask_i)
  # filter sub-masks
  # account for biggest lake first...
  fmodel = mmtbx.f_model.manager(
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    f_calc       = f_calc,
    f_mask       = [f_masks[0]])
  fmodel.update_all_scales(remove_outliers=False, update_f_part1=False)
  k_anisotropic = fmodel.k_isotropic()*fmodel.k_anisotropic()
  sc = fmodel.scale_k1()
  f_mask_biggest = fmodel.f_bulk().data()
  one = flex.double(ss.size(),1.)
  ks = flex.double([k/100 for k in range(0,105,1)])
  bs = flex.double([0,])
  #bs = flex.double(range(0,50,5))
  # ... then score others
  for i, f in enumerate(f_masks):
    if(i==0):
      fm = [f_masks[0]]
      r = 1.e+9
    else:
      res = mmtbx.bulk_solvent.ksol_bsol_grid_search(
        f_obs.data(),
        f_calc.data() + f_mask_biggest,
        f.data(),
        ks,
        bs,
        ss,
        sc,
        one,
        k_anisotropic,
        r)
      k_sol, b_sol, r = res
      f_mask_biggest += f.data()*mmtbx.f_model.ext.k_mask(ss, k_sol, b_sol)
      if(log is not None):
        print >> log, "mask: %3d k_sol: %4.2f b_sol: %6.2f r: %8.6f"%(
          i+1, k_sol, b_sol, r)
        log.flush()
      if(k_sol>0.): fm.append(f)
  del f_masks
  # accumulated remaining (selected) masks
  for i, f in enumerate(fm):
    if(i==0):
      F = f_calc.deep_copy()
      FB = 0
    else:
      FB = FB+fmodel.f_bulk().data()
      F = F.array(data = f_calc.data() + FB)
    fmodel = mmtbx.f_model.manager(
      f_obs        = f_obs,
      r_free_flags = r_free_flags,
      f_calc       = F,
      f_mask       = [f])
    fmodel.update_all_scales(remove_outliers=False, update_f_part1=False)
    if(log is not None):
      print >> log, "mask: %3d r_work: %6.4f r_free: %6.4f"%(i, fmodel.r_work(),
        fmodel.r_free()), fmodel.r_work_low()
      log.flush()
  # consistency check and reconstruct result fmodel object
  r_work_start = fmodel.r_work()
  r_free_start = fmodel.r_free()
  fmodel = mmtbx.f_model.manager(
    f_obs         = f_obs,
    r_free_flags  = r_free_flags,
    f_calc        = f_calc,
    k_isotropic   = fmodel.k_isotropic(),
    k_anisotropic = fmodel.k_anisotropic(),
    k_mask        = one,
    f_mask        = [f_obs.array(data = FB + fmodel.f_bulk().data())])
  assert approx_equal(r_work_start, fmodel.r_work())
  assert approx_equal(r_free_start, fmodel.r_free())
  #
  return fmodel
