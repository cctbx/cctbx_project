"""
Utility script to process bulk (model,data) pairs
"""

from __future__ import absolute_import, division, print_function
import os, sys
import iotbx.pdb
from scitbx.array_family import flex
import mmtbx.f_model
from libtbx import easy_mp
from mmtbx.bulk_solvent import mosaic
from mmtbx.bulk_solvent import example
from libtbx import group_args
from libtbx import easy_pickle
import traceback
from cctbx import maptbx
from mmtbx import masks
import boost_adaptbx.boost.python as bp
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")


from cctbx.masks import vdw_radii_from_xray_structure
ext = bp.import_ext("mmtbx_masks_ext")

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

def get_fmodel(o, f_mask, remove_outliers):
  fo, fm = o.f_obs().common_sets(f_mask)
  fc, fm = o.f_calc().common_sets(f_mask)
  fmodel = mmtbx.f_model.manager(
    f_obs  = fo,
    f_calc = fc,
    f_mask = fm)
  fmodel.update_all_scales(
    remove_outliers         = remove_outliers,
    apply_scale_k1_to_f_obs = False)
  return fmodel

class compute(object):
  def __init__(self, pdbf, mtzf):
    #
    D = example.get_data(pdbf, mtzf)
    #step = 0.4
    step = D.f_obs().d_min()/4
    #
    f_mask = mosaic.get_f_mask(
      xrs  = D.xray_structure,
      ma   = D.f_obs(),
      step = step)
    fmodel_2013 = get_fmodel(
      o = D, f_mask = f_mask, remove_outliers = True)
    fmodel_2013.show(show_header=False, show_approx=False)
    #print(fmodel_2013.r_work(d_min=11.188))
    #print(fmodel_2013.resolution_filter(d_min=11.188).r_work())
    #STOP()
    #
    rr = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    #
    r_sol_best = None
    r_best = 99999
    for r_solvent in rr:
      fmodel = try_triplet(
        xrs         = D.xray_structure,
        fmodel      = fmodel_2013,
        step        = 0.4,
        r_solvent   = 0.8,
        r_shrink    = 0.9,
        filter_mask = True)
      fmodel.show(show_header=False, show_approx=False)
      STOP()
      r = fmodel.resolution_filter(d_min=10).r_work()
      if(r<r_best):
        r_sol_best = r_solvent
        r_best = r
      print(r, r_solvent)
    print(r_sol_best)



def run_one(args):
  pdbf, mtzf, code = args
  # FILTER OUT NON-P1 and non X-ray/neutron
  pdb_inp = iotbx.pdb.input(file_name=pdbf)
  cs = pdb_inp.crystal_symmetry()
  #
  if(cs.space_group_number() != 1): return
  exp_type = pdb_inp.get_experiment_type()
  if not (exp_type.is_xray() or exp_type.is_neutron()): return

  #log = sys.stdout
  #log = open("%s.log"%code,"w")
  if 1:#try:
    # main
    o = compute(pdbf=pdbf, mtzf=mtzf)
#    #log.flush()
#    ### SKIP
#    if(not o.mm.do_mosaic or o.fmodel_2013.r_work()>0.30 or
#       len(o.mm.regions.values())<1):
#      log.close()
#      if(os.path.isfile("%s.log"%code)): os.remove("%s.log"%code)
#      return
#    ###
#    result = group_args(
#      code            = code,
#      d_min           = o.fmodel_2013.f_obs().d_min(),
#      r_2013          = o.fmodel_2013.r_factors(as_string=False),
#      r_2013_04       = o.fmodel_2013_04.r_factors(as_string=False),
#      r_filtered      = o.fmodel_filtered.r_factors(as_string=False),
#      r_mosaic        = mbs.fmodel.r_factors(as_string=False),
#      r_largest_mask  = r_largest_mask,
#      solvent_content = o.mm.solvent_content,
#      regions         = o.mm.regions)
#    easy_pickle.dump("%s.pkl"%code, result)
#    log.close()
#  except: # intentional
#    print("FAILED:", file=log)
#    traceback.print_exc(file=log)
#    log.close()

def get_mask_p1(xrs, step, r_shrink, r_solvent):
  #crystal_gridding = maptbx.crystal_gridding(
  #  unit_cell        = xrs.unit_cell(),
  #  space_group_info = xrs.space_group_info(),
  #  symmetry_flags   = maptbx.use_space_group_symmetry,
  #  step             = step)
  #n_real = crystal_gridding.n_real()
  #mask_p1 = mmtbx.masks.mask_from_xray_structure(
  #  xray_structure           = xrs,
  #  p1                       = True,
  #  for_structure_factors    = False,
  #  n_real                   = n_real,
  #  solvent_radius           = r_solvent,
  #  shrink_truncation_radius = r_shrink,
  #  in_asu                   = False).mask_data
  #maptbx.unpad_in_place(map=mask_p1)

  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = xrs.unit_cell(),
    space_group_info = xrs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = step)
  n_real = crystal_gridding.n_real()
  atom_radii = vdw_radii_from_xray_structure(xray_structure = xrs)
  mask_params = masks.mask_master_params.extract()
  grid_step_factor = 1.962/step
  mask_params.solvent_radius           = r_solvent
  mask_params.shrink_truncation_radius = r_shrink

  o = mmtbx.masks.bulk_solvent(
      xray_structure              = xrs,
      ignore_zero_occupancy_atoms = False,
      solvent_radius              = mask_params.solvent_radius,
      shrink_truncation_radius    = mask_params.shrink_truncation_radius,
      ignore_hydrogen_atoms       = False,
      gridding_n_real             = n_real,
      atom_radii                  = atom_radii)
  mask_p1 = o.data.as_double()
  print(mask_p1)

  return mask_p1

def try_triplet(xrs, fmodel, step, r_solvent, r_shrink, filter_mask):
  mask_p1 = get_mask_p1(
    xrs       = xrs,
    step      = step,
    r_shrink  = r_shrink,
    r_solvent = r_solvent)
  if(filter_mask):
    mask_p1 = mosaic.filter_mask(
      mask_p1          = mask_p1,
      volume_cutoff    = 100,
      crystal_symmetry = xrs.crystal_symmetry())
  mask = asu_map_ext.asymmetric_map(
    xrs.crystal_symmetry().space_group().type(), mask_p1).data()
  f_mask = fmodel.f_obs().structure_factors_from_asu_map(
    asu_map_data = mask_p1, n_real = mask_p1.accessor().all())
  fmodel = mmtbx.f_model.manager(
    f_obs          = fmodel.f_obs(),
    r_free_flags   = fmodel.r_free_flags(),
    f_calc         = fmodel.f_calc(),
    bin_selections = fmodel.bin_selections,
    f_mask         = f_mask)
  fmodel.update_all_scales(remove_outliers=False, apply_scale_k1_to_f_obs=False)
  return fmodel

def run(cmdargs):
  if(len(cmdargs)==1):
    NPROC=70
    pdbs, mtzs, codes, sizes = get_files_sorted(pdb_files, hkl_files)
    argss = []
    for pdb, mtz, code in zip(pdbs, mtzs, codes):
      argss.append([pdb, mtz, code])
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
    assert len(cmdargs) == 2
    # Usage: python example.py 4qnn.pdb 4qnn.mtz
    pdb, mtz = cmdargs
    assert os.path.isfile(pdb)
    assert os.path.isfile(mtz)
    code = os.path.abspath(pdb)[:-4]
    run_one([pdb, mtz, code])

if (__name__ == "__main__"):
  run(sys.argv[1:])
