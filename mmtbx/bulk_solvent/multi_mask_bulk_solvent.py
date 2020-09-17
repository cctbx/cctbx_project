from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import mmtbx.f_model
from cctbx import maptbx
import mmtbx.masks
import mmtbx.bulk_solvent
import boost_adaptbx.boost.python as bp
from six.moves import zip
from six.moves import range
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
from mmtbx import map_tools
import cctbx.miller
from libtbx.test_utils import approx_equal

def ccp4_map(cg, file_name, mc=None, map_data=None):
  assert [mc, map_data].count(None)==1
  if(map_data is None):
    map_data = get_map(mc=mc, cg=cg)
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
      file_name=file_name,
      unit_cell=cg.unit_cell(),
      space_group=cg.space_group(),
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))

def get_k_mask(method, f_obs, f_calc, f_mask, ss, sel):
  obj = method(
    f_obs     = f_obs,
    f_calc    = f_calc,
    f_mask    = f_mask,
    selection = sel)
  k_mask = obj.x_best
  r = obj.r_best
  if(abs(k_mask) < 1.e-6):
    one = flex.double(ss.select(sel).size(),1.)
    ks = flex.double([k/100 for k in range(-105,105,1)])
    #ks = flex.double(
    #  [k_mask]+[k/100. for k in range(int((k_mask-k_mask)*100.), int((k_mask+k_mask)*100.))])
    bs = flex.double([0,])
    res = mmtbx.bulk_solvent.ksol_bsol_grid_search(
      f_obs.select(sel),
      f_calc.select(sel),
      f_mask.select(sel),
      ks,
      bs,
      ss.select(sel),
      1,
      one,
      one,
      r)
    k_mask, b_sol, r = res
    #print k_mask
  return k_mask, r

def loop_work(fmodel, f_masks, method):
  #
  sel_free =  fmodel.r_free_flags().data()
  sel_work = ~sel_free
  f_bulk_data_all = flex.complex_double(fmodel.f_obs().data().size(), 0)

  k_anisotropic = fmodel.k_isotropic()*fmodel.k_anisotropic()*fmodel.scale_k1()

  f_obs_w          = fmodel.f_obs().select(sel_work)
  f_calc_w         = fmodel.f_calc().select(sel_work)
  f_masks_w        = [f.select(sel_work) for f in f_masks]
  k_anisotropic_w  = k_anisotropic.select(sel_work)
  bin_selections_w = [b.select(sel_work) for b in fmodel.bin_selections]
  ss_w             = fmodel.ss.select(sel_work)

  f_obs_f          = fmodel.f_obs().select(sel_free)
  f_calc_f         = fmodel.f_calc().select(sel_free)
  f_masks_f        = [f.select(sel_free) for f in f_masks]
  k_anisotropic_f  = k_anisotropic.select(sel_free)
  bin_selections_f = [b.select(sel_free) for b in fmodel.bin_selections]
  ss_f             = fmodel.ss.select(sel_free)
  #

  f_obs_w = f_obs_w.array(data = f_obs_w.data()/(k_anisotropic_w))
  f_obs_f = f_obs_f.array(data = f_obs_f.data()/(k_anisotropic_f))
  #
  F   = f_calc_w.deep_copy().data()
  F_f = f_calc_f.deep_copy().data()
  FB   = 0
  FB_f = 0
  r_best = 999
  for i, f in enumerate(f_masks_w):
    f_f = f_masks_f[i]
    f_bulk_data_w = flex.complex_double(f_obs_w.data().size(), 0)
    f_bulk_data_f = flex.complex_double(f_obs_f.data().size(), 0)
    k_masks = flex.double()
    for sel_w, sel_f in zip(bin_selections_w, bin_selections_f):
      k_mask, r = get_k_mask(
        method = method,
        f_obs  = f_obs_w.data(),
        f_calc = F,
        f_mask = f.data(),
        ss     = ss_w,
        sel    = sel_w)
      k_masks.append(k_mask)
    #if(k_masks[0]<0 or abs(k_masks[0])<1.e-3): k_masks *= 0.
    for i, k_mask in enumerate(k_masks):
      #if k_mask<0: k_mask=0.
      sel_w, sel_f = bin_selections_w[i], bin_selections_f[i]
      f_bulk_data_w = f_bulk_data_w.set_selected(sel_w, k_mask*f.data().select(sel_w))
      f_bulk_data_f = f_bulk_data_f.set_selected(sel_f, k_mask*f_f.data().select(sel_f))
    #
    F    = F    + f_bulk_data_w
    F_f  = F_f  + f_bulk_data_f
    FB   = FB   + f_bulk_data_w
    FB_f = FB_f + f_bulk_data_f
    #
    rw1 = mmtbx.bulk_solvent.r_factor(f_obs_w.data(),   F)
    rf1 = mmtbx.bulk_solvent.r_factor(f_obs_f.data(), F_f)
    if(rf1<r_best):
      r_best = rf1
    else:
      FB   = FB   - f_bulk_data_w
      FB_f = FB_f - f_bulk_data_f
      F   = F   - f_bulk_data_w
      F_f = F_f - f_bulk_data_f
    rw2 = mmtbx.bulk_solvent.r_factor(f_obs_w.data(),   F)
    rf2 = mmtbx.bulk_solvent.r_factor(f_obs_f.data(), F_f)
    #print "%7.5f %7.5f <> %7.5f %7.5f %7.5f"%(rw1, rf1, rw2, rf2,  rf2-rw2)
  #
  f_bulk_data_all = f_bulk_data_all.set_selected(sel_work, FB)
  f_bulk_data_all = f_bulk_data_all.set_selected(sel_free, FB_f)
  f_bulk = fmodel.f_obs().array(data = f_bulk_data_all)
  return f_bulk

def run_method(fmodel, f_masks, one, method):
  f_bulk = loop_work(fmodel=fmodel, f_masks=f_masks, method = method)
  fmodel = mmtbx.f_model.manager(
    f_obs         = fmodel.f_obs(),
    r_free_flags  = fmodel.r_free_flags(),
    f_calc        = fmodel.f_calc(),
    k_mask        = one,
    f_mask        = f_bulk)
  fmodel.update_all_scales(remove_outliers=False, update_f_part1=False)
  return fmodel

def helper_3(
      fmodel,
      f_masks,
      log):
  #f_masks = fmodel.f_masks() #+ f_masks[0:]
  for i in range(len(f_masks)):
    f_masks[i] = f_masks[i].common_set(fmodel.f_obs())
  one = flex.double(fmodel.ss.size(), 1.)
  method_1 = \
    mmtbx.bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical
  method_2 = mmtbx.bulk_solvent.bulk_solvent_scale_coefficients_analytical
  # Run methods
  fmodel_1 = run_method(fmodel=fmodel, f_masks=f_masks, one=one, method=method_1)
  fmodel_2 = run_method(fmodel=fmodel, f_masks=f_masks, one=one, method=method_2)
  #
  #print "="*79
  #fmodel_1.show()
  #print fmodel_1.r_work(), fmodel_1.r_free(), fmodel_1.r_work_low()
  #print "="*79
  #fmodel_2.show()
  #print fmodel_2.r_work(), fmodel_2.r_free(), fmodel_2.r_work_low()
  #
  if(fmodel_1.r_free()<fmodel_2.r_free()):
    if(fmodel_1.r_work_low()<fmodel_2.r_work_low()):
      return fmodel_1, 1
    else:
      return fmodel_2, 2
  else:
    if(fmodel_2.r_work_low()<fmodel_1.r_work_low()):
      return fmodel_2, 2
    else:
      return fmodel_1, 1
################################################################################

def compute_map(fmodel, crystal_gridding, map_type):
  map_coefficients = map_tools.electron_density_map(
    fmodel = fmodel).map_coefficients(
      map_type         = map_type,
      isotropize       = True,
      fill_missing     = False)
  fft_map = cctbx.miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = map_coefficients)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def get_fmodel_and_grid_step(f_obs, r_free_flags, xrs):
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.grid_step_factor=5.
  fmodel_1 = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    mask_params    = mask_params,
    xray_structure = xrs)
  fmodel_1.update_all_scales(remove_outliers=True, update_f_part1=False)
  rw1, rf1, rwl1 = fmodel_1.r_work(), fmodel_1.r_free(), fmodel_1.r_work_low()
  #
  fmodel_2 = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xrs)
  fmodel_2.update_all_scales(remove_outliers=True, update_f_part1=False)
  rw2, rf2, rwl2 = fmodel_2.r_work(), fmodel_2.r_free(), fmodel_2.r_work_low()
  if(rf1<rf2 and rwl1<rwl2): return fmodel_1, 5.
  else:                      return fmodel_2, 4.

def get_mask_1(fmodel, grid_step_factor):
  grid_step = fmodel.f_obs().d_min()*(1./grid_step_factor)
  if(grid_step < 0.15): grid_step = 0.15
  grid_step = min(0.8, grid_step)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell = fmodel.xray_structure.unit_cell(),
    space_group_info = fmodel.xray_structure.space_group_info(),
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
  return mask_data_p1, n_real, crystal_gridding

def get_mask_2(fmodel, grid_step_factor):
  sgt = fmodel.xray_structure.space_group().type()
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.grid_step_factor = grid_step_factor
  asu_mask_obj = mmtbx.masks.asu_mask(
    xray_structure = fmodel.xray_structure,
    d_min          = fmodel.f_obs().d_min(),
    mask_params    = mask_params).asu_mask
  mask_data_p1 = asu_mask_obj.mask_data_whole_uc()
  maptbx.unpad_in_place(map=mask_data_p1)
  n_real = mask_data_p1.all()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell             = fmodel.xray_structure.unit_cell(),
    space_group_info      = fmodel.xray_structure.space_group_info(),
    symmetry_flags        = maptbx.use_space_group_symmetry,
    pre_determined_n_real = n_real)

  #n_real = mask_data_p1.all()
  #mask_data_p1 = asu_map_ext.asymmetric_map(sgt, mask_data_p1, n_real).symmetry_expanded_map()
  #maptbx.unpad_in_place(map=mask_data_p1)

  return mask_data_p1, n_real, crystal_gridding

class multi_mask_bulk_solvent(object):
  def __init__(self, fmodel, log=None):
    # Commonly used objects
    xrs = fmodel.xray_structure
    sgt = xrs.space_group().type()
    # Compute default fmodel and decide on grid step
    fmodel, self.grid_step_factor = get_fmodel_and_grid_step(
      f_obs        = fmodel.f_obs(),
      r_free_flags = fmodel.r_free_flags(),
      xrs          = xrs)
    #fmodel.show()
    #print fmodel.r_work(), fmodel.r_free()
    ###
    mask_data_p1, n_real, crystal_gridding = get_mask_1(fmodel=fmodel,
      grid_step_factor=self.grid_step_factor)
    #ccp4_map(cg=crystal_gridding, file_name="m1.ccp4", map_data=mask_data_p1_)
    #xxx1 = fmodel.f_obs().structure_factors_from_map(map=mask_data_p1,
    #   use_scale = True, anomalous_flag = False, use_sg = False)

    #mask_data_p1, n_real, crystal_gridding = get_mask_2(fmodel=fmodel,
    #  grid_step_factor=self.grid_step_factor)
    #print n_real
    #STOP()
    #xxx2 = fmodel.f_obs().structure_factors_from_map(map=mask_data_p1,
    #   use_scale = True, anomalous_flag = False, use_sg = False)
    #
    #assert approx_equal(xxx1.data(), xxx2.data())

    #print mask_data_p1.all(), mask_data_p1.focus(), mask_data_p1.origin()
    #print mask_data_p1_2.all(), mask_data_p1_2.focus(), mask_data_p1_2.origin()
    #print mask_data_p1_1.count(0), mask_data_p1_2.count(0)
    #assert approx_equal(mask_data_p1_1, mask_data_p1_2)
    #STOP()
    #####
    # Mask connectivity analysis
    co = maptbx.connectivity(map_data=mask_data_p1, threshold=0.01)
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
    if(log is not None): print("Number of regions:", len(region_indices), file=log)
    mi,ma,me,diff_map_asu = None,None,None,None
    for ii, i in enumerate(region_indices):
      s = conn==i
      si = s.iselection()
      if(not all_zero_found and mask_data_asu.select(si).count(0.)>0):
        all_zero_found = True
        continue
      # DIFF MAP START
      if(region_volumes[ii]<1 and diff_map_asu is None):#(ii == 2):
        fmodel_tmp = mmtbx.f_model.manager(
          f_obs          = fmodel.f_obs(),
          r_free_flags   = fmodel.r_free_flags(),
          f_calc         = fmodel.f_calc(),
          f_mask         = f_masks[len(f_masks)-1])
        fmodel_tmp.update_all_scales(remove_outliers=False, update_f_part1=False)
        diff_map_p1 = compute_map(
          fmodel           = fmodel_tmp,
          crystal_gridding = crystal_gridding,
          map_type         = "mFo-DFc")
        diff_map_asu = asu_map_ext.asymmetric_map(sgt, diff_map_p1).data()
      if(diff_map_asu is not None):
        mi,ma,me = diff_map_asu.select(si).min_max_mean().as_tuple()
        if(ma<0. or me<0.):
          continue
      # DIFF MAP END

      #XXX this is 4 loops, may be slow. move to C++ if slow.
      mask_data_asu_i = mask_data_asu.deep_copy()
      #mask_data_asu_i = mask_data_asu_i.set_selected(s, 1).set_selected(~s, 0)
      mask_data_asu_i = mask_data_asu_i.set_selected(~s, 0)

      #if(mi is None):
      #  print "region: %5d fraction: %8.4f"%(ii, region_volumes[ii]), len(region_volumes)
      #else:
      #  print "region: %5d fraction: %8.4f"%(ii, region_volumes[ii]), len(region_volumes), "%7.3f %7.3f %7.3f"%(mi,ma,me)

      if(log is not None):
        print("region: %5d fraction: %8.4f"%(ii, region_volumes[ii]), file=log)
        log.flush()
      f_mask_i = fmodel.f_obs().structure_factors_from_asu_map(
        asu_map_data = mask_data_asu_i, n_real = n_real)
      if(len(f_masks)>0 and region_volumes[ii]>1):
        f_masks[len(f_masks)-1] = f_masks[len(f_masks)-1].array(data = f_masks[len(f_masks)-1].data()+
          f_mask_i.data())
      else:
        f_masks.append(f_mask_i)
    #
    self.fmodel_result, self.method = helper_3(
        fmodel  = fmodel,
        f_masks = f_masks,
        log     = log)
    #self.fmodel_result.show()
    #
    self.n_regions = len(region_volumes[1:])
    self.region_volumes = " ".join(["%8.4f"%(v) for v in region_volumes[1:][:10]]) # top 10
