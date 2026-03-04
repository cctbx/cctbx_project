"""
Utility script to process bulk (model,data) pairs
"""

from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_mp
import traceback
from libtbx import easy_run

cif_files = os.getenv('PDB_MIRROR_MMCIF')

def get_model_files_dict(path):
  ifn = open("/".join([path,"INDEX"]),"r")
  result = {}
  for l in ifn.readlines():
    l = l.strip()
    file_name = "/".join([path,l])
    #assert os.path.isfile(file_name) # Terribly runtime expensive
    code = l[-11:-7]
    #
    #
    #if code != "13mj": continue # XXX DEBUGGING
    #
    #
    result[code] = file_name
  return result

def fetch_sorry(fn, code, cmd):
  sorry = None
  with open(fn, "r") as fo:
    for l in fo.readlines():
      if l.startswith("Sorry:"):
        sorry = l
        break
  if sorry is not None:
    with open("%s.error.log"%code, "w") as fo:
      print(cmd, file=fo)
      print(sorry, file=fo)
    return True
  return False

def run_one(args):
  cif, code = args
  #
  ofn = "%s_H.cif"%code
  zlg = "%s.zlog"%code
  cmd = "mmtbx.reduce2 %s output.filename=%s >& %s"%(cif, ofn, zlg)
  try:
    rc = easy_run.call(cmd)
    if rc!=0:
      found = fetch_sorry(fn=zlg, code=code, cmd=cmd)
      if not found:
        raise RuntimeError("Error code is failure, but Sorry not found.")
  except Exception as e:
    of = open("%s.error.log"%code, "w")
    print(cmd, file=of)
    print(file=of)
    traceback.print_exc(file=of)
    of.close()

def run(NPROC=128):
  #
  processed = []
  for f in os.listdir("."):
    if(f.endswith(".pdb")):
      processed.append(f[:f.index(".")])
  #
  cifs = get_model_files_dict(path=cif_files)
  print("Number of cifs:", len(cifs.keys()))
  argss = []
  for code in cifs.keys():
    if(code in processed): continue
    cif = cifs[code]
    argss.append([cif, code])
    #if len(argss)==20: break
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
  run()
