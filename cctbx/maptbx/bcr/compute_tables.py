from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from scitbx.array_family import flex
from cctbx.eltbx import xray_scattering
from libtbx import group_args, easy_mp
import sys
import traceback

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
  dist = 14.0
  o = maptbx.atom_curves(scattering_type=e, scattering_table=table)
  #
  err1, err2 = None, None
  #
  try:
    b1 = o.bcr_approx(
      d_min       = d_min,
      radius_max  = dist,
      radius_step = 0.01,
      mxp   = 1000,
      epsc  = 0.001,
      epsp  = 0.000,
      edist = 1.0E-13,
      kpres = 1,
      kprot = 111
      )
    r1, err1 = rfactor(b1.image_values, b1.bcr_approx_values)
  except Exception as e:
    of = open("%s_%s.111.error.log"%(e, str(d_min)), "w")
    traceback.print_exc(file=of)
    of.close()
  #
  if err1 > 0.01:
    try:
      b2 = o.bcr_approx(
        d_min       = d_min,
        radius_max  = dist,
        radius_step = 0.01,
        mxp   = 1000,
        epsc  = 0.001,
        epsp  = 0.000,
        edist = 1.0E-13,
        kpres = 1,
        kprot = 112
        )
      r2, err2 = rfactor(b2.image_values, b2.bcr_approx_values)
    except Exception as e:
      of = open("%s_%s.112.error.log"%(e, str(d_min)), "w")
      traceback.print_exc(file=of)
      of.close()
  #
  if not None in [err1, err2]:
    if err1 < err2: b, err = b1, err1
    else:           b, err = b2, err2
  elif [err1, err2].count(None) == 2: return None
  else:
    if err1 is not None: b, err = b1, err1
    else:                b, err = b2, err2
  #
  return group_args(
   d_min = d_min,
   dist  = dist,
   err   = err,
   table = table,
   R     = b.R,
   B     = b.B,
   C     = b.C)

if (__name__ == "__main__"):
  #
  NPROC=120
  #
  for table in ["electron", "wk1995"]:
    for e in ["H", "O", "C", "N", "S"]:
      print(e)
      sys.stdout.flush()
      results = {}
      argss = []
      for d_min in [round(0.3+i*0.01, 2) for i in range(int((6-0.3)/0.01)+1)]:
        argss.append([e, d_min, table])
      #
      if(NPROC>1):
        stdout_and_results = easy_mp.pool_map(
          processes    = NPROC,
          fixed_func   = run_one,
          args         = argss,
          func_wrapper = "buffer_stdout_stderr")
        for it in stdout_and_results:
          o = it[1]
          if o is None:
            print("FAILED"*5)
            continue
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
      with open("%s_%s.json"%(e, table), "w") as f:
        json.dump(results, f, indent=2)
