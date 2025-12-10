from __future__ import absolute_import, division, print_function
from cctbx.eltbx import xray_scattering
from libtbx import easy_mp

from cctbx.maptbx.bcr import bcr

std_labels = xray_scattering.standard_labels_list()

def run_one(args):
  bcr.compute_tables(
      MinResolution = 4.0,
      MaxResolution = 4.0,
      scattering_table = "wk1995",
      TypesAtoms = args)


if (__name__ == "__main__"):
  #
  NPROC=120
  #
  argss = []
  for e in [["H","C"], ["O","N"], ["P","S"]]:
    argss.append(e)
  print(argss)
  #
  if(NPROC>1):
    stdout_and_results = easy_mp.pool_map(
      processes    = NPROC,
      fixed_func   = run_one,
      args         = argss,
      func_wrapper = "buffer_stdout_stderr")
  else:
    for args in argss:
      o = run_one(args)
