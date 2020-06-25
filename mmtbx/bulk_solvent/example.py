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
import mmtbx.utils
import mmtbx.f_model
from libtbx import easy_mp
from mmtbx.bulk_solvent import mosaic
from libtbx import group_args
from libtbx.test_utils import approx_equal
from libtbx import easy_pickle
import traceback

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

def get_data(mtzf, xrs):
  reflection_file = reflection_file_reader.any_reflection_file(
    file_name=mtzf, ensure_read_access=False)
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry=xrs.crystal_symmetry(),
    force_symmetry=True,
    reflection_files=[reflection_file],
    err=null_out())
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server  = rfs,
    keep_going              = True,
    extract_r_free_flags    = False,
    log                     = null_out())
  f_obs = determine_data_and_flags_result.f_obs
  return f_obs, determine_data_and_flags_result.r_free_flags

def get_fmodel_init(pdbf, mtzf, log):
  pdb_inp = iotbx.pdb.input(file_name=pdbf)
  xrs = pdb_inp.xray_structure_simple()
  f_obs, r_free_flags = get_data(mtzf=mtzf, xrs=xrs)
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xrs)
  fmodel.update_all_scales()
  print("-"*79, file=log)
  print("Default (A-2013):", file=log)
  fmodel.show(show_header=False, show_approx=False, log = log)
  print(fmodel.r_factors(), file=log)
  return fmodel, pdb_inp.construct_hierarchy()

class compute(object):
  def __init__(self, pdbf, mtzf, log, volume_cutoff=30):
    self.log = log
    self.fmodel_2013, self.ph = get_fmodel_init(pdbf=pdbf, mtzf=mtzf, log=log)
    step = min(0.4, self.fmodel_2013.f_calc().d_min()/4)
    #####
    print("-"*79, file=log)
    print("Mosaic", file=log)
    self.mm = mosaic.mosaic_f_mask(
      xray_structure = self.fmodel_2013.xray_structure,
      miller_array   = self.fmodel_2013.f_calc(),
      step           = step,
      volume_cutoff  = volume_cutoff,
      f_obs          = self.fmodel_2013.f_obs(),
      r_free_flags   = self.fmodel_2013.r_free_flags(),
      f_calc         = self.fmodel_2013.f_calc(),
      log            = log)
    #####
    print("-"*79, file=log)
    print("A-2013 with step=0.4A", file=log)
    self.fmodel_2013_04 = mmtbx.f_model.manager(
      f_obs        = self.fmodel_2013.f_obs(),
      r_free_flags = self.fmodel_2013.r_free_flags(),
      f_calc       = self.fmodel_2013.f_calc(),
      f_mask       = self.mm.f_mask)
    self.fmodel_2013_04.update_all_scales(remove_outliers=False)
    self.fmodel_2013_04.show(show_header=False, show_approx=False, log = log)
    print(self.fmodel_2013_04.r_factors(prefix="  "), file=log)
    self.mc_whole_mask = \
      self.mm.fmodel_largest_mask.electron_density_map().map_coefficients(
        map_type   = "mFobs-DFmodel",
        isotropize = True,
        exclude_free_r_reflections = False)
    #####
    print("-"*79, file=log)
    print("A-2013 with step=0.4A, using largest mask only", file=log)
    self.mm.fmodel_largest_mask.show(show_header=False, show_approx=False, log = log)
    print(self.mm.fmodel_largest_mask.r_factors(prefix="  "), file=log)
    self.mc_largest_mask = \
      self.mm.fmodel_largest_mask.electron_density_map().map_coefficients(
        map_type   = "mFobs-DFmodel",
        isotropize = True,
        exclude_free_r_reflections = False)
    #####

  def do_mosaic(self, alg):
    print("-"*79, file=self.log)
    print("Refine k_masks", file=self.log)
    # fmodel just to do mosaic
    bin_selections = self.fmodel_2013.f_obs().log_binning(
      n_reflections_in_lowest_resolution_bin = 500)
    fmodel = mmtbx.f_model.manager(
      f_obs          = self.fmodel_2013.f_obs(),
      r_free_flags   = self.fmodel_2013.r_free_flags(),
      f_calc         = self.fmodel_2013.f_calc(),
      f_mask         = self.mm.f_mask,
      bin_selections = bin_selections)
    fmodel.update_all_scales(remove_outliers=False)
    #
    result = mosaic.refinery(fmodel=fmodel, fv=self.mm.FV, alg=alg,
      log=self.log)
    print("", file=self.log)
    #result.fmodel.show(show_header=False, show_approx=False, log = self.log)
    result.fmodel.show_short(show_k_mask=False, prefix="  ", log = self.log)
    print(result.fmodel.r_factors(prefix="  "), file=self.log)
    return result

def get_map(mc, cg):
  fft_map = mc.fft_map(crystal_gridding = cg)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def map_stat(m, conn, i):
  selection = conn==i
  blob = m.select(selection.iselection())
  mi,ma,me = flex.min(blob), flex.max(blob), flex.mean(blob)
  sd = blob.sample_standard_deviation()
  return group_args(mi=mi, ma=ma, me=me, sd=sd)

def run_one(args):
  pdbf, mtzf, code, alg = args
  # FILER OUT NON-P1 and non X-ray/neutron
  pdb_inp = iotbx.pdb.input(file_name=pdbf)
  cs = pdb_inp.crystal_symmetry()
  if(cs.space_group_number() != 1): return
  if(not pdb_inp.get_experiment_type() in
     ["X-RAY DIFFRACTION", "NEUTRON DIFFRACTION"]): return
  #
  #log = sys.stdout
  log = open("%s.log"%code,"w")
  try:
    # main
    o = compute(pdbf=pdbf, mtzf=mtzf, log=log)
    log.flush()
    mbs = o.do_mosaic(alg=alg)
    log.flush()
    ### SKIP
    if(o.fmodel_2013.r_work()>0.45 or
       len(o.mm.regions.values())<1 or
       not o.mm.do_mosaic):
      log.close()
      os.remove("%s.log"%code)
      return
    ###
    # write maps
    if(o.mm.mc is not None):
      assert approx_equal(o.mm.mc.data(), o.mc_largest_mask.data())
      mtz_dataset = o.mm.mc.as_mtz_dataset(column_root_label='FistMask')
      mtz_dataset.add_miller_array(
        miller_array=o.mc_whole_mask, column_root_label="WholeMask")
      mtz_dataset.add_miller_array(
        miller_array=mbs.mc, column_root_label="Mosaic")
      mtz_object = mtz_dataset.mtz_object()
      mtz_object.write(file_name = "%s_mc.mtz"%code)
      # map stats
      map_FirstMask = get_map(mc=o.mm.mc,         cg=o.mm.crystal_gridding)
      map_WholeMask = get_map(mc=o.mc_whole_mask, cg=o.mm.crystal_gridding)
      map_Mosaic    = get_map(mc=mbs.mc,          cg=o.mm.crystal_gridding)
      cntr = 0
      for region in o.mm.regions.values():
        region.m_FistMask  = map_stat(m=map_FirstMask, conn = o.mm.conn, i=region.id)
        region.m_WholeMask = map_stat(m=map_WholeMask, conn = o.mm.conn, i=region.id)
        region.m_Mosaic    = map_stat(m=map_Mosaic,    conn = o.mm.conn, i=region.id)
        if(region.diff_map.me is not None):
          cntr += 1
          assert approx_equal(region.diff_map.me, region.m_FistMask.me)
          assert approx_equal(region.diff_map.sd, region.m_FistMask.sd)
          assert approx_equal(region.diff_map.ma, region.m_FistMask.ma)
      assert cntr>0
    #
    result = group_args(
      code            = code,
      d_min           = o.fmodel_2013.f_obs().d_min(),
      r_2013          = o.fmodel_2013.r_factors(as_string=False),
      r_2013_04       = o.fmodel_2013_04.r_factors(as_string=False),
      r_mosaic        = mbs.fmodel.r_factors(as_string=False),
      solvent_content = o.mm.solvent_content,
      regions         = o.mm.regions)
    easy_pickle.dump("%s.pkl"%code, result)
    log.close()
  except: # intentional
    print("FAILED:", file=log)
    traceback.print_exc(file=log)
    log.close()

def run(cmdargs):
  if(len(cmdargs)==0):
    NPROC=50
    pdbs, mtzs, codes, sizes = get_files_sorted(pdb_files, hkl_files)
    argss = []
    for pdb, mtz, code in zip(pdbs, mtzs, codes):
      argss.append([pdb, mtz, code, "alg4"])
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
    # Usage: python example.py 4qnn.pdb 4qnn.mtz alg4
    pdb, mtz, alg = cmdargs
    code = os.path.abspath(pdb)[:-4]
    run_one([pdb, mtz, code, alg])

if (__name__ == "__main__"):
  run(sys.argv[1:])
