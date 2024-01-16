"""
Utility script to process bulk (model,data) pairs
"""

from __future__ import absolute_import, division, print_function
import os, sys
import iotbx.pdb
from libtbx.utils import null_out
from scitbx.array_family import flex
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from libtbx.utils import null_out
import mmtbx.f_model
from libtbx import easy_mp
from mmtbx.bulk_solvent import mosaic
from libtbx import group_args
from libtbx.test_utils import approx_equal
from libtbx import easy_pickle
import traceback
from iotbx import extract_xtal_data

pdb_files = "/net/cci/pdb_mirror/pdb/"
hkl_files = "/net/cci-filer2/raid1/share/pdbmtz/mtz_files/"

def get_files_sorted(pdb_files, hkl_files):
  ifn_p = open("/".join([pdb_files,"INDEX"]),"r")
  ifn_r = os.listdir(hkl_files)
  pdbs  = flex.std_string()
  mtzs  = flex.std_string()
  codes = flex.std_string()
  sizes = flex.double()
  cntr=0
  for lp in ifn_p.readlines():
    lp = lp.strip()
    pdb_file_name = "/".join([pdb_files,lp])
    assert os.path.isfile(pdb_file_name)
    pdb_code = lp[-11:-7]
    #
    # FOR DEBUGGING
    #if pdb_code != "4w71": continue
    #
    #
    lr = lp.replace("pdb","r")
    hkl_file_name = "/".join([hkl_files,"%s.mtz"%pdb_code])
    if(os.path.isfile(hkl_file_name)):
      cntr+=1
      s = os.path.getsize(pdb_file_name) + os.path.getsize(hkl_file_name)
      pdbs  .append(pdb_file_name)
      mtzs  .append(hkl_file_name)
      codes .append(pdb_code)
      sizes .append(s)
      #if codes.size()==100: break
  print("Total:", cntr)
  sel = flex.sort_permutation(sizes)
  pdbs  = pdbs .select(sel)
  mtzs  = mtzs .select(sel)
  codes = codes.select(sel)
  sizes = sizes.select(sel)
  return pdbs, mtzs, codes, sizes

def get_data(pdbf, mtzf):
  pdb_inp = iotbx.pdb.input(file_name=pdbf)
  xrs = pdb_inp.xray_structure_simple()
  #
  selection = xrs.scatterers().extract_occupancies() > 0
  xrs = xrs.select(selection)
  selection = ~xrs.hd_selection()
  xrs = xrs.select(selection)
  #
  #xrs.switch_to_neutron_scattering_dictionary()
  #
  reflection_file = reflection_file_reader.any_reflection_file(
    file_name=mtzf, ensure_read_access=False)
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry=xrs.crystal_symmetry(),
    force_symmetry=True,
    reflection_files=[reflection_file],
    err=null_out())
  determine_data_and_flags_result = extract_xtal_data.run(
    reflection_file_server  = rfs,
    keep_going              = True,
    extract_r_free_flags    = True,
    force_non_anomalous     = True,
    allow_mismatch_flags    = True)
  f_obs        = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xrs)
  fmodel.update_all_scales(
    remove_outliers         = True,
    apply_scale_k1_to_f_obs = True
    )
  # Skip
  rw = fmodel.r_work()
  rf = fmodel.r_free()
  if(fmodel.f_obs().d_min()>3):           return None
  if(fmodel.f_obs().completeness()<0.95): return None
  if(rw>0.25):                            return None
  if(rf<=rw):                             return None
  if(abs(rf-rw)*100<2):                   return None
  if(fmodel.f_obs().resolution_filter(d_min=8).completeness()<0.95): return None
  #
  def f_obs():        return fmodel.f_obs()
  def r_free_flags(): return fmodel.r_free_flags()
  def f_calc():       return fmodel.f_calc()
  def flags():        return fmodel.r_free_flags()
  return group_args(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = fmodel.xray_structure,
    f_calc         = f_calc,
    flags          = flags)

def get_fmodel(o, f_mask, remove_outliers, log):
  fo, fm = o.f_obs().common_sets(f_mask)
  fc, fm = o.f_calc().common_sets(f_mask)

  bin_selections = fo.log_binning(
        n_reflections_in_lowest_resolution_bin = 100)

  fmodel = mmtbx.f_model.manager(
    bin_selections = bin_selections,
    f_obs        = fo,
    r_free_flags = o.r_free_flags(),
    f_calc       = fc,
    f_mask       = fm)
  fmodel.update_all_scales(
    remove_outliers         = remove_outliers,
    apply_scale_k1_to_f_obs = False
    )
  fmodel.show(show_header=False, show_approx=False, log = log)
  print(fmodel.r_factors(prefix="  "), file=log)
  mc = fmodel.electron_density_map().map_coefficients(
    map_type   = "mFobs-DFmodel",
    isotropize = True,
    exclude_free_r_reflections = False)
  return group_args(fmodel = fmodel, mc = mc)

class compute(object):
  def __init__(self, D, log):
    self.log = log
    #
    print("-"*79, file=log)
    print("A-2013, all defaults", file=log)
    step_default = min(0.8, D.f_obs().d_min()/4)
    if(step_default < 0.15): step_default = 0.15
    f_mask = mosaic.get_f_mask(
      xrs  = D.xray_structure,
      ma   = D.f_obs(),
      step = step_default)
    self.fmodel_2013 = get_fmodel(
      o = D, f_mask = f_mask, remove_outliers = False, log = self.log).fmodel
    #
    # Compute masks and F_masks (Mosaic)
    print("-"*79, file=log)
    print("Mosaic", file=log)
    self.mm = mosaic.f_masks(
      xray_structure          = D.xray_structure,
      step                    = 0.6,
      volume_cutoff           = 50,
      mean_diff_map_threshold = 0.5,
      r_sol                   = 1.1, # whole
      r_shrink                = 0.9, # whole
      f_obs                   = self.fmodel_2013.f_obs(),
      f_calc                  = self.fmodel_2013.f_calc(),
      log                     = log)
    #
    if(self.mm.do_mosaic):
      ######
      print("-"*79, file=log)
      print("A-2013, opt step/radii", file=log)
      f_mask_opt = mosaic.get_f_mask(
        xrs      = D.xray_structure,
        ma       = self.fmodel_2013.f_obs(),
        step     = 0.6,
        r_shrink = 0.9,
        r_sol    = 1.1)
      o = get_fmodel(o = self.fmodel_2013, f_mask = f_mask_opt,
        remove_outliers = False, log = self.log)
      self.fmodel_2013_opt = o.fmodel
      self.mc_whole_mask   = o.mc
      ######
      print("-"*79, file=log)
      print("A-2013, opt step/radii, using largest mask only", file=log)
      o = get_fmodel(o = self.fmodel_2013, f_mask = self.mm.f_mask_0,
        remove_outliers = False, log = self.log)
      self.fmodel_0 = o.fmodel
      self.mc_0     = o.mc

  def do_mosaic(self, alg):
    print("-"*79, file=self.log)
    print("Refine k_masks", file=self.log)
    assert approx_equal(self.fmodel_2013_opt.f_obs().data(),
      self.fmodel_0.f_obs().data())
    assert approx_equal(self.fmodel_2013_opt.f_calc().data(),
      self.fmodel_0.f_calc().data())

    if(self.fmodel_0.r_two_bins_lowest() < self.fmodel_2013_opt.r_two_bins_lowest()):
      fmodel = self.fmodel_0
    else:
      fmodel = self.fmodel_2013_opt

    result = mosaic.refinery(
      fmodel  = fmodel,
      fv      = self.mm.FV,
      alg     = alg,
      log     = self.log)
    print("", file=self.log)
    result.fmodel.show(show_header=False, show_approx=False, log = self.log)
    print(result.fmodel.r_factors(prefix="  "), file=self.log)
    return result

def get_map(mc, cg):
  fft_map = mc.fft_map(crystal_gridding = cg)
  #fft_map.apply_sigma_scaling()
  md = fft_map.real_map_unpadded()
  sd = md.sample_standard_deviation()
  #return fft_map.real_map_unpadded()
  return md/sd, sd

def map_stat(m, conn, i, map_sd):
  selection = conn==i
  blob = m.select(selection.iselection())
  mi,ma,me = flex.min(blob), flex.max(blob), flex.mean(blob)
  sd = blob.sample_standard_deviation()
  return group_args(mi=mi, ma=ma, me=me, sd=sd, map_sd=map_sd)

def run_one(args):
  pdbf, mtzf, code, alg = args
  # FILTER OUT NON-P1 and non X-ray/neutron
  pdb_inp = iotbx.pdb.input(file_name=pdbf)
  cs = pdb_inp.crystal_symmetry()
  #
  if(cs.space_group_number() != 1): return
  exp_type = pdb_inp.get_experiment_type()
  if not (exp_type.is_xray() or exp_type.is_neutron()): return
  #
  #log = sys.stdout
  log = open("%s.log"%code,"w")
  try:
    D = get_data(pdbf, mtzf)
    if(D is None):
      log.close()
      if(os.path.isfile("%s.log"%code)): os.remove("%s.log"%code)
      return
    # main
    o = compute(D=D, log=log)
    log.flush()
    ### SKIP
    if(not o.mm.do_mosaic or len(o.mm.regions.values())<1):
      log.close()
      if(os.path.isfile("%s.log"%code)): os.remove("%s.log"%code)
      return
    ###
    mbs = o.do_mosaic(alg=alg)
    # write maps
    if(o.mm.mFoDFc_0 is not None):
      mtz_dataset = o.mc_0.as_mtz_dataset(column_root_label='FistMask')
      mtz_dataset.add_miller_array(
        miller_array=o.mc_whole_mask, column_root_label="WholeMask")
      mtz_dataset.add_miller_array(
        miller_array=mbs.mc, column_root_label="Mosaic")
      mtz_object = mtz_dataset.mtz_object()
      mtz_object.write(file_name = "%s_mc.mtz"%code)
      # map stats
      map_0        , sd_map_0         = get_map(mc=o.mm.mFoDFc_0,   cg=o.mm.crystal_gridding)
      map_WholeMask, sd_map_WholeMask = get_map(mc=o.mc_whole_mask, cg=o.mm.crystal_gridding)
      map_Mosaic   , sd_map_Mosaic    = get_map(mc=mbs.mc,          cg=o.mm.crystal_gridding)
      ###
      #write_map_file(cg=o.mm.crystal_gridding, mc=o.mm.mc,         file_name="first.ccp4")
      #write_map_file(cg=o.mm.crystal_gridding, mc=o.mc_whole_mask, file_name="whole.ccp4")
      #write_map_file(cg=o.mm.crystal_gridding, mc=mbs.mc,          file_name="mosaic.ccp4")
      ###
      for region in o.mm.regions.values():
        region.m_0         = map_stat(m=map_0,         conn = o.mm.conn, i=region.id, map_sd=sd_map_0        )
        region.m_WholeMask = map_stat(m=map_WholeMask, conn = o.mm.conn, i=region.id, map_sd=sd_map_WholeMask)
        region.m_Mosaic    = map_stat(m=map_Mosaic,    conn = o.mm.conn, i=region.id, map_sd=sd_map_Mosaic   )
    #
    result = group_args(
      code            = code,
      d_min           = o.fmodel_2013.f_obs().d_min(),
      r_2013          = o.fmodel_2013.r_factors(as_string=False),
      r_2013_opt      = o.fmodel_2013_opt.r_factors(as_string=False),
      r_mosaic        = mbs.fmodel.r_factors(as_string=False),
      r_largest_mask  = o.fmodel_0.r_factors(as_string=False),
      solvent_content = o.mm.solvent_content,
      regions         = o.mm.regions)
    easy_pickle.dump("%s.pkl"%code, result)
    log.close()
  except: # intentional
    print("FAILED:", file=log)
    traceback.print_exc(file=log)
    log.close()

def write_map_file(cg, mc, file_name):
  from iotbx import mrcfile
  fft_map = mc.fft_map(crystal_gridding=cg)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  mrcfile.write_ccp4_map(
    file_name   = file_name,
    unit_cell   = cg.unit_cell(),
    space_group = cg.space_group(),
    map_data    = map_data,
    labels      = flex.std_string([""]))

def run(cmdargs):
  if(len(cmdargs)==1):
    alg = cmdargs[0]
    assert alg in ["alg0", "alg2", "alg4", "alg4a"]
    if alg=="None": alg=None
    NPROC=70
    pdbs, mtzs, codes, sizes = get_files_sorted(pdb_files, hkl_files)
    argss = []
    for pdb, mtz, code in zip(pdbs, mtzs, codes):
      if(os.path.isfile("%s.log"%code) and
         os.path.isfile("%s.pkl"%code) and
         os.path.isfile("%s_mc.mtz"%code)): continue
      argss.append([pdb, mtz, code, alg])
    if(NPROC>1):
      stdout_and_results = easy_mp.pool_map(
        processes    = NPROC,
        fixed_func   = run_one,
        args         = argss,
        func_wrapper = "buffer_stdout_stderr")
    else:
      for args in argss:
        run_one(args)
  else:
    assert len(cmdargs) == 3
    # Usage: python example.py 4qnn.pdb 4qnn.mtz alg4
    pdb, mtz, alg = cmdargs
    assert alg in ["alg0", "alg2", "alg4", "alg4a"]
    if alg=="None": alg=None
    assert os.path.isfile(pdb)
    assert os.path.isfile(mtz)
    code = os.path.abspath(pdb)[:-4]
    run_one([pdb, mtz, code, alg])

if (__name__ == "__main__"):
  run(sys.argv[1:])
