"""
Utility script to process bulk (model,data) pairs
"""

from __future__ import absolute_import, division, print_function
import os, sys
from scitbx.array_family import flex
import mmtbx.model
from libtbx import easy_mp
from collections import OrderedDict
import iotbx.cif
import traceback
import os
from iotbx import reflection_file_reader
from libtbx import easy_run

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
    #if(cntr==10): break
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

def run_one(args):
  cif, hkl, code = args
  #
  try:
    miller_arrays = reflection_file_reader.any_reflection_file(file_name =
      hkl).as_miller_arrays()
    d_min = miller_arrays[0].d_min()
    if d_min < 3.0: return None
    if d_min > 4.5: return None
    cmd = "phenix.find_reference %s use_abs=true use_master=true include_models=csm output.prefix=%s --json >& %s.zlog"%(cif, code, code)
    easy_run.call(cmd)
  #
  except Exception as e:
    of = open("%s.log"%code, "w")
    traceback.print_exc(file=of)
    of.close()

def run(cmdargs, NPROC=128):
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
    try:
      cif = cifs[code]
      argss.append([cif, hkl, code])
    except: pass # intentional
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
