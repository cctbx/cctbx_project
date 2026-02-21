"""
Utility script to process bulk (model,data) pairs
"""

from __future__ import absolute_import, division, print_function
import os, sys
from scitbx.array_family import flex
import mmtbx.f_model
import mmtbx.model
from libtbx import easy_mp
from collections import OrderedDict
import iotbx.cif
import traceback
import iotbx.pdb
import os
from iotbx import reflection_file_reader
from mmtbx import map_tools

cif_files = os.getenv('PDB_MIRROR_MMCIF')
pdbmtz = os.getenv('PDBMTZ')

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

def get_mtz_dict(path):
  result     = OrderedDict()
  codes      = flex.std_string()
  file_names = flex.std_string()
  sizes      = flex.double()
  cntr = 0
  for l in os.listdir(path):
    l = l.strip()
    if not l.endswith(".mtz"): continue
    cntr += 1
    code = l[:-4]
    #
    # FOR DEBUGGING
    #
    #if code != "6a1o": continue
    #if code != "6a1q": continue
    if(cntr==10): break
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

def get_map(fmodel, map_type):
  return map_tools.electron_density_map(
    fmodel = fmodel).map_coefficients(
      map_type         = map_type,
      isotropize       = True,
      fill_missing     = False)

def run_one(args):
  cif, hkl, code = args
  #
  try:
    inp = iotbx.pdb.input(file_name=cif)
    model = mmtbx.model.manager(
      model_input         = inp,
      stop_for_unknowns   = False,
      skip_ss_annotations = True,
      process_biomt       = False)
    xrs = model.get_xray_structure()
    #
    f_obs, r_free_flags = None, None
    miller_arrays = reflection_file_reader.any_reflection_file(file_name =
      hkl).as_miller_arrays()
    for ma in miller_arrays:
      if("F-obs" in ma.info().label_string()): f_obs = ma
      if("R-free-flags" in ma.info().label_string()): r_free_flags = ma
    assert [f_obs, r_free_flags].count(None)==0
    r_free_flags = r_free_flags.array(data = r_free_flags.data()==1)
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs,
      r_free_flags   = r_free_flags,
      xray_structure = xrs)
    fmodel.update_all_scales()
    print(fmodel.r_work(), fmodel.r_free())
    mc1 = get_map(fmodel=fmodel, map_type="2mFo-DFc")
    mc2 = get_map(fmodel=fmodel, map_type="mFo-DFc")
    mtz_dataset = mc1.as_mtz_dataset(column_root_label = "2FOFCWT")
    mtz_dataset.add_miller_array(
      miller_array=mc2,
      column_root_label="FOFCWT")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = "%s_mc.mtz"%code)
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
  hkls = get_mtz_dict(path=pdbmtz)
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
