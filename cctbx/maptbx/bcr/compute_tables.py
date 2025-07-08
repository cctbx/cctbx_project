from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from scitbx.array_family import flex
from cctbx.eltbx import xray_scattering
from libtbx import group_args, easy_mp
import sys

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_maptbx_bcr_bcr_ext")

std_labels = xray_scattering.standard_labels_list()

import json

def add_entry(results, d_min, dist, err, table, R, B, C):
  results[str(d_min)] = {
      "dist": dist,
      "err": err,
      "scat_factor": table,
      "R": R,
      "B": B,
      "C": C
  }

def rfactor(a,b):
  n = flex.sum(flex.abs(a-b))
  d = flex.sum(flex.abs(a+b))
  return n/d*100*2, flex.max(flex.abs(a-b))

def run_one(args):
  e, d_min, table = args
  dist = 18.0
  o = maptbx.atom_curves(scattering_type=e, scattering_table=table)
  b = o.bcr_approx(
    d_min       = d_min,
    radius_max  = dist,
    radius_step = 0.01,
    kprot=2,
    kpres=0,
    mxp=1000,
    epsc=0.001,
    )
  r, err = rfactor(b.image_values, b.bcr_approx_values)
  #print(r, err)
  #for R, B, C in zip(b.R, b.B, b.C):
  #  print("%20.16f %20.16f %20.16f"%(R, B, C) )
  return group_args(
   d_min = d_min,
   dist  = dist,
   err   = err,
   table = table,
   R     = b.R,
   B     = b.B,
   C     = b.C)

if (__name__ == "__main__"):
  NPROC=3
  #
  for e in ["H", "O", "C", "N", "S"]:
    print(e)
    sys.stdout.flush()
    results = {}
    argss = []
    for d_min in [round(0.3+i*0.01, 2) for i in range(int((6-0.3)/0.01)+1)]:
      argss.append([e, d_min, "wk1995"])
    #
    if(NPROC>1):
      stdout_and_results = easy_mp.pool_map(
        processes    = NPROC,
        fixed_func   = run_one,
        args         = argss,
        func_wrapper = "buffer_stdout_stderr")
      for it in stdout_and_results:
        o = it[1]
        add_entry(
          results = results,
          d_min   = o.d_min,
          dist    = o.dist,
          err     = o.err,
          table   = o.table,
          R       = o.R,
          B       = o.B,
          C       = o.C)
    else:
      for args in argss:
        o = run_one(args)
        add_entry(
          results = results,
          d_min   = o.d_min,
          dist    = o.dist,
          err     = o.err,
          table   = o.table,
          R       = o.R,
          B       = o.B,
          C       = o.C)
    #
    with open("%s.json"%e, "w") as f:
      json.dump(results, f, indent=2)
