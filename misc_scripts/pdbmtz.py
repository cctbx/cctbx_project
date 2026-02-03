"""
Utility script to process bulk (model,data) pairs
"""

from __future__ import absolute_import, division, print_function
import os, sys
from scitbx.array_family import flex
from iotbx import reflection_file_utils
from libtbx.utils import null_out
import mmtbx.f_model
import mmtbx.model
from libtbx import easy_mp
from collections import OrderedDict
from mmtbx.command_line import cif_as_mtz
import iotbx.cif
from iotbx import crystal_symmetry_from_any
from iotbx import reflection_file_utils
from cctbx import french_wilson
import traceback
import iotbx.pdb
import os

cif_files = os.getenv('PDB_MIRROR_MMCIF')
hkl_files = os.getenv('PDB_MIRROR_STRUCTURE_FACTORS')

def get_model_files_dict(path):
  ifn = open("/".join([path,"INDEX"]),"r")
  result = {}
  for l in ifn.readlines():
    l = l.strip()
    file_name = "/".join([path,l])
    #assert os.path.isfile(file_name) # Terribly runtime expensive
    code = l[-11:-7]
    result[code] = file_name
  return result

def get_hkl_files_dict(path):
  ifn = open("/".join([path,"INDEX"]),"r")
  result     = OrderedDict()
  codes      = flex.std_string()
  file_names = flex.std_string()
  sizes      = flex.double()
  cntr = 0
  for l in ifn.readlines():
    l = l.strip()
    if not l.endswith(".gz"): continue
    cntr += 1
    code = l[-13:-9]
    #
    # FOR DEBUGGING
    #
    #if code != "6a1o": continue
    #if code != "6a1q": continue
    #if(cntr==1000): break
    #
    file_name = "/".join([path,l])
    #assert os.path.isfile(file_name) # Terribly runtime expensive
    codes.append(code)
    file_names.append(file_name)
    sizes.append(os.path.getsize(file_name))
  sel = flex.sort_permutation(sizes)
  codes      = codes     .select(sel)
  file_names = file_names.select(sel)
  sizes      = sizes     .select(sel)
  for c,f,s in zip(codes, file_names, sizes):
    result[c] = [f,s]
  return result

def apply_filter(ma, pub):
  info = ma.info()
  label_str = info.label_string()
  if("FOM" in label_str): return None
  if("FC" in label_str): return None
  x = ma.deep_copy()
  x = x.resolution_filter(d_max = pub.low, d_min = pub.high)
  x = x.resolution_filter(d_min = pub.resolution)
  x = x.sigma_filter(cutoff_factor = pub.sigma)
  x.set_info(info)
  return x

def regularize_cv_flags(ma):
  info = ma.info()
  fv = reflection_file_utils.guess_r_free_flag_value(ma)
  if(fv is None): return None
  fraction = ma.data().count(fv)*100./ma.data().size()
  if(fraction > 15): return None
  ma = ma.array(data = ma.data()==fv)
  ma.set_info(info)
  return ma

def i_to_f(d):
  if not d.is_xray_intensity_array(): return d
  result = french_wilson.french_wilson_scale(
    miller_array                   = d,
    params                         = None, #self.parameters.french_wilson,
    sigma_iobs_rejection_criterion = None,
    log                            = null_out())
  return result

def equalize_anom_flags(r_free_flags, f_obs):
  if r_free_flags.anomalous_flag()==f_obs.anomalous_flag():
    return r_free_flags
  if f_obs.anomalous_flag(): return r_free_flags.as_anomalous_array()
  else:                      return r_free_flags.as_non_anomalous_array()

def run_one(args):
  cif, hkl, code = args
  #
  try:
    inp = iotbx.pdb.input(file_name=cif)
    pub = inp.get_r_rfree_sigma()
    if(not inp.get_experiment_type().is_xray()): return None
    cs1 = crystal_symmetry_from_any.extract_from(cif)
    cs2 = crystal_symmetry_from_any.extract_from(hkl)
    if([cs1,cs2].count(None)==0):
      assert cs1.is_similar_symmetry(cs2)
    o = cif_as_mtz.extract(
      file_name                       = hkl,
      crystal_symmetry                = cs1,
      wavelength_id                   = None,
      crystal_id                      = None,
      show_details_if_error           = True,
      output_r_free_label             = "R-free-flags",
      merge_non_unique_under_symmetry = True,
      map_to_asu                      = True,
      remove_systematic_absences      = True,
      all_miller_arrays               = None,
      incompatible_flags_to_work_set  = False,
      ignore_bad_sigmas               = True,
      extend_flags                    = False,
      return_as_miller_arrays         = False,
      log                             = null_out())
    #
    # PRE-SELECT AND PRE-FILTER DATA AND FLAGS
    #
    flags = []
    data  = []
    for ma in o.as_miller_arrays():
      ma = apply_filter(ma=ma, pub=pub)
      if ma is None: continue
      if("R-free-flags" in ma.info().label_string()):
        ma = regularize_cv_flags(ma)
        if(ma is None): continue
        flags.append(ma.deep_copy())
      if(ma.is_xray_data_array()):
        ma = i_to_f(ma)
        if ma is None: continue
        data.append(ma.deep_copy())
    #
    # STOP HERE IF NO FLAGS FOUND
    #
    if len(flags)==0: return None
    #
    # THIS IS X-RAY MODEL
    #
    model = mmtbx.model.manager(
      model_input         = inp,
      stop_for_unknowns   = False,
      skip_ss_annotations = True,
      process_biomt       = False)
    xrs = model.get_xray_structure()
    #
    # TRY ALL-VS-ALL AND KEEP PLAUSIBLE PAIRS
    #
    rw_best = pub.r_work
    if rw_best is None: rw_best = 999
    best_f_obs = None
    best_flags = None
    for f_obs in data:
      for r_free_flags in flags:
        r_free_flags_dc = equalize_anom_flags(
          r_free_flags = r_free_flags,
          f_obs        = f_obs)
        f_obs_dc, r_free_flags_dc = f_obs.common_sets(r_free_flags_dc)
        fmodel = mmtbx.f_model.manager(
          f_obs          = f_obs_dc.deep_copy(),
          r_free_flags   = r_free_flags_dc.deep_copy(),
          xray_structure = xrs)
        fmodel.update_all_scales()
        rw, rf = fmodel.r_work(), fmodel.r_free()
        if rw < rf and rw < rw_best+0.03:
          rw_best = rw
          best_f_obs = f_obs_dc.deep_copy()
          best_flags = r_free_flags_dc.deep_copy()
    #
    # DUMP RESULT INTO FINAL MTZ
    #
    if best_f_obs is not None:
      assert best_flags is not None
      mtz_dataset = best_f_obs.as_mtz_dataset(column_root_label = "F-obs")
      mtz_dataset.add_miller_array(
        miller_array      = best_flags,
        column_root_label = "R-free-flags")
      mtz_object = mtz_dataset.mtz_object()
      mtz_object.write(file_name = "%s.mtz"%code)
  #
  except Exception as e:
    of = open("%s.log"%code, "w")
    traceback.print_exc(file=of)
    of.close()

def run(cmdargs, NPROC=120):
  #
  processed = []
  for f in os.listdir("."):
    if(f.endswith(".mtz")):
      processed.append(f[:f.index(".")])
  #
  cifs = get_model_files_dict(path=cif_files)
  hkls = get_hkl_files_dict(  path=hkl_files)
  print("Number of cifs:", len(cifs.keys()))
  print("Number of hkls:", len(hkls.keys()))
  argss = []
  for code in hkls.keys():
    if(code in processed): continue
    hkl,size = hkls[code]
    cif = cifs[code]
    argss.append([cif, hkl, code])
  print("NEW FILES TO PROCESS:", len(argss))
  if(NPROC>1):
    stdout_and_results = easy_mp.pool_map(
      processes    = NPROC,
      fixed_func   = run_one,
      args         = argss,
      func_wrapper = "buffer_stdout_stderr")
  else:
    for args in argss:
      run_one(args)

if (__name__ == "__main__"):
  run(sys.argv[1:])
