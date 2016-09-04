from __future__ import division
from cctbx.array_family import flex
import mmtbx.f_model
from cctbx import maptbx
from libtbx.test_utils import approx_equal
import mmtbx.masks
import mmtbx.bulk_solvent
import boost.python
asu_map_ext = boost.python.import_ext("cctbx_asymmetric_map_ext")

def ccp4_map(cg, file_name, mc=None, map_data=None):
  assert [mc, map_data].count(None)==1
  if(map_data is None):
    map_data = get_map(mc=mc, cg=cg)
  from iotbx import ccp4_map
  ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=cg.unit_cell(),
      space_group=cg.space_group(),
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))


def helper_1(
      f_obs,
      r_free_flags,
      f_calc,
      f_masks,
      xrs,
      ss,
      log):
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
  return fmodel

def helper_2(
      f_obs,
      r_free_flags,
      f_calc,
      f_masks,
      xrs,
      ss,
      log):
  from mmtbx.bulk_solvent import kbu_refinery
  #TGCO = kbu_refinery.tgc(f_obs=f_obs, f_calc=f_calc, f_masks=[f_masks[0]], ss=ss,
  #  k_sols=[0.35], b_sols=[50], u_star=[0,0,0,0,0,0])
  ##TGCO.minimize_kbu()
  #TGCO.minimize_kbu_sequential(use_curvatures_options=[False,True], n_cycles=5)
  #TGCO.show_kbu()
  ##print TGCO.k_sols, TGCO.b_sols, "LOOK (u_star_start)"
  #u_start_start = TGCO.kbu.u_star()
  #print
  fmodel = mmtbx.f_model.manager(
    f_obs        = f_obs,
    r_free_flags = r_free_flags,
    f_calc       = f_calc,
    f_mask       = [f_masks[0]])
  fmodel.update_all_scales(remove_outliers=False, update_f_part1=False)
  k_anisotropic = fmodel.k_isotropic()*fmodel.k_anisotropic() * fmodel.scale_k1()



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
  f_masks = fm





  n_masks = len(f_masks)
  bin_selections = f_obs.log_binning(n_reflections_in_lowest_resolution_bin=n_masks*10)
  f_bulk_data = flex.complex_double(f_obs.data().size(), 0)

  #sc = fmodel.scale_k1()
  #k_anisotropic = mmtbx.f_model.ext.k_anisotropic(f_obs.indices(), u_start_start)
  f_obs_ = f_obs.array(data = f_obs.data()/(k_anisotropic))

  for sel in bin_selections:
    TGCO = kbu_refinery.tgc(
      f_obs   = f_obs_.select(sel),
      f_calc  = f_calc.select(sel),
      f_masks = [fm.select(sel) for fm in f_masks],
      ss      = ss.select(sel),
      k_sols  = [0.35]*n_masks,
      b_sols  = [0]*n_masks,
      u_star  = [0,0,0, 0,0,0])
    for it in [1,2]:
      TGCO.minimize_k_once(use_curvatures=False)
      TGCO.minimize_k_once(use_curvatures=True)
    TGCO.show_k_sols()
    fbd = 0
    for k, fm in zip(TGCO.kbu.k_sols(), f_masks):
      #if(k>0.):
        fbd += k * fm.data().select(sel)
    f_bulk_data = f_bulk_data.set_selected(sel.iselection(), fbd)
  f_bulk = f_obs.array(data = f_bulk_data)
  one = flex.double(ss.size(),1.)
  fmodel = mmtbx.f_model.manager(
    f_obs         = f_obs,
    r_free_flags  = r_free_flags,
    f_calc        = f_calc,
    #k_isotropic   = fmodel.k_isotropic(),
    #k_anisotropic = fmodel.k_anisotropic(),
    k_mask        = one,
    f_mask        = f_bulk)
  fmodel.update_all_scales(remove_outliers=False, update_f_part1=False)
  return fmodel

################################################################################

def helper_3(
      f_obs,
      r_free_flags,
      f_calc,
      f_masks,
      xrs,
      ss,
      log):
  from mmtbx.bulk_solvent import kbu_refinery
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.grid_step_factor=5.
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    mask_params    = mask_params,
    xray_structure = xrs)
  fmodel.update_all_scales(remove_outliers=True, update_f_part1=False)
  k_anisotropic = fmodel.k_isotropic()*fmodel.k_anisotropic()*fmodel.scale_k1()
  #
  f_obs        = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags()
  f_calc       = fmodel.f_calc()
  ss           = fmodel.ss
  for i in xrange(len(f_masks)):
    f_masks[i] = f_masks[i].common_set(f_obs)
  #
  f_masks = fmodel.f_masks() + f_masks[0:]
  f_bulk_data = flex.complex_double(f_obs.data().size(), 0)
  f_obs_ = f_obs.array(data = f_obs.data()/(k_anisotropic))
  for sel in fmodel.bin_selections:
    k_masks = []
    for i, f in enumerate(f_masks):
      if(i==0):
        F = f_calc.deep_copy().data()
        FB = 0
      else:
        FB = FB + FB_
        F = f_calc.data() + FB
      obj = mmtbx.bulk_solvent.bulk_solvent_scale_coefficients_analytical(
        f_obs     = f_obs_.data(),
        f_calc    = F,
        f_mask    = f.data(),
        selection = sel)
      k_mask = obj.x_best
      k_masks.append(k_mask)
      r = obj.r_best
      if(abs(k_mask) < 1.e-6):
        one = flex.double(ss.select(sel).size(),1.)
        ks = flex.double([k/100 for k in range(-105,105,1)])
        bs = flex.double([0,])
        res = mmtbx.bulk_solvent.ksol_bsol_grid_search(
          f_obs_.data().select(sel),
          F.select(sel),
          f.data().select(sel),
          ks,
          bs,
          ss.select(sel),
          1,
          one,
          one,
          r)
        k_mask, b_sol, r = res
      FB_ = k_mask * f.data()
      k_masks.append(k_mask)
    FB = FB + k_mask * f.data()
    f_bulk_data = f_bulk_data.set_selected(sel.iselection(), FB.select(sel))
  f_bulk = f_obs.array(data = f_bulk_data)
  one = flex.double(ss.size(),1.)
  fmodel = mmtbx.f_model.manager(
    f_obs         = f_obs,
    r_free_flags  = r_free_flags,
    f_calc        = f_calc,
    k_mask        = one,
    f_mask        = f_bulk)
  fmodel.update_all_scales(remove_outliers=False, update_f_part1=False)
  return fmodel
################################################################################

def multi_mask_bulk_solvent(fmodel, log=None):
  # Commonly used objects
  f_calc       = fmodel.f_calc()
  f_obs        = fmodel.f_obs()
  r_free_flags = fmodel.r_free_flags()
  ss = 1./flex.pow2(f_obs.d_spacings().data()) / 4.
  xrs = fmodel.xray_structure
  sgt = xrs.space_group().type()
  # Crystal_gridding
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
  # Compute mask in P1
  mask_data_p1 = mmtbx.masks.mask_from_xray_structure(
    xray_structure        = fmodel.xray_structure,
    p1                    = True,
    for_structure_factors = True,
    n_real                = n_real,
    in_asu                = False).mask_data
  maptbx.unpad_in_place(map=mask_data_p1)
  # Mask connectivity analysis
  co = maptbx.connectivity(map_data=mask_data_p1, threshold=0.9)
  conn = co.result().as_double()
  # Convert result of connectivity analysis from P1 to ASU (in-place)
  conn = asu_map_ext.asymmetric_map(sgt, conn).data()
  # Find unique indices and regions in reduced (P1->ASU) conn
  region_indices = flex.double()
  region_volumes = flex.double()
  for i in conn:
    if not i in region_indices: region_indices.append(i)
  for l in region_indices:
    szl = conn.count(l)*100./conn.size()
    region_volumes.append(szl)
  s = flex.sort_permutation(region_volumes, reverse=True)
  region_volumes = region_volumes.select(s)
  region_indices = region_indices.select(s)
  # Convert P1 mask into ASU
  mask_data_asu = asu_map_ext.asymmetric_map(sgt, mask_data_p1).data()
  conn.reshape(mask_data_asu.accessor()) #XXX still need it?
  f_masks = []
  all_zero_found = False
  print >> log, "Number of regions:", len(region_indices)
  for ii, i in enumerate(region_indices):
    s = conn==i
    if(not all_zero_found and mask_data_asu.select(s.iselection()).count(0.)>0):
      all_zero_found = True
      continue
    #XXX this is 4 loops, may be slow. move to C++ if slow.
    mask_data_asu_i = mask_data_asu.deep_copy()
    mask_data_asu_i = mask_data_asu_i.set_selected(s, 1).set_selected(~s, 0)
    if(log is not None):
      print >> log, "region: %5d fraction: %8.4f"%(ii, region_volumes[ii])
      log.flush()
    f_mask_i = f_calc.structure_factors_from_asu_map(
      asu_map_data = mask_data_asu_i, n_real = n_real)
    f_masks.append(f_mask_i)
  #
  return helper_3(
      f_obs        = f_obs,
      r_free_flags = r_free_flags,
      f_calc       = f_calc,
      f_masks      = f_masks,
      xrs          = xrs,
      ss           = ss,
      log          = log)
