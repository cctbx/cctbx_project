from __future__ import absolute_import, division, print_function
from cctbx.eltbx import xray_scattering
from libtbx import easy_mp
from cctbx import crystal
from cctbx.array_family import flex
from cctbx import xray
from cctbx.maptbx.bcr import bcr

std_labels = xray_scattering.standard_labels_list()


def chunk_list(input_list, N=3):
    """
    Splits the input list into a list of lists, where each sub-list has at
    least N elements.
    Returns: A list of sub-lists, where each sub-list contains at least N
             elements, except possibly the last one if not enough elements
            remain.
    """
    if N <= 0: raise ValueError("N must be a positive integer.")
    result = []
    # Iterate through the list while ensuring we create sub-lists of at least N elements.
    for i in range(0, len(input_list), N):
      sub_list = input_list[i:i + N]
      # Adjust if the last sub-list is smaller than N but not empty
      if len(sub_list) < N and len(result) > 0:
        result[-1].extend(sub_list)  # Combine with the last sub-list
      else:
        result.append(sub_list)
    return result

def make_xrs(s):
  cs1 = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 1")
  sp1 = crystal.special_position_settings(cs1)
  scatterers1 = flex.xray_scatterer((
    xray.scatterer(s, (0.5, 0, 0)),
    xray.scatterer(s, (0, 0, 0))))
  return xray.structure(sp1, scatterers1)

def run_one_x(args):
  bcr.compute_tables(
    MinResolution    = 0.996,
    MaxResolution    = 10.1,
    DistMax          = 11.0,
    Ngrid            = 1101,
    scattering_table = "wk1995",
    TypesAtoms       = args)

def run_one_e(args):
  bcr.compute_tables(
    MinResolution    = 0.996,
    MaxResolution    = 10.1,
    DistMax          = 11.0,
    Ngrid            = 1101,
    scattering_table = "electron",
    TypesAtoms       = args)

def run_one():
   bcr.compute_tables(
     MinResolution   = 0.996,
     MaxResolution    = 1.1,
     DistMax          = 11.0,
     Ngrid            = 1101,
     scattering_table = "wk1995",
     TypesAtoms       = ["N", "S", "C"])

def run_all(xray=True, electron=True):
  #
  e_list = []
  for l in std_labels:
    xrs = make_xrs(s=l)
    r = xrs.scattering_type_registry(table="electron")
    if len(list(r.unassigned_types()))==0:
      e_list.append(l)
  #
  e_list = chunk_list(input_list=e_list)
  x_list = chunk_list(input_list=std_labels)
  #
  NPROC=120
  #
  if(NPROC>1):
    if xray:
      stdout_and_results = easy_mp.pool_map(
        processes    = NPROC,
        fixed_func   = run_one_x,
        args         = x_list,
        func_wrapper = "buffer_stdout_stderr")
    if electron:
      stdout_and_results = easy_mp.pool_map(
        processes    = NPROC,
        fixed_func   = run_one_e,
        args         = e_list,
        func_wrapper = "buffer_stdout_stderr")
  else: # XXX BROKEN
    for args in argss:
      o = run_one(args)

if (__name__ == "__main__"):
  if True: run_all()
  else:    run_one()
