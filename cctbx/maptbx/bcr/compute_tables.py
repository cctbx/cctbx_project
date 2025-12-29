from __future__ import absolute_import, division, print_function
from cctbx.eltbx import xray_scattering
from libtbx import easy_mp
from cctbx import crystal
from cctbx.array_family import flex
from cctbx import xray
from cctbx.maptbx.bcr import bcr


std_labels = xray_scattering.standard_labels_list()

def chunk_list(original_list):
    L = len(original_list)

    # Check if the list is smaller than the minimum required length for meaningful chunks
    if L < 3:
        raise ValueError("The length of the list must be at least 3.")

    chunks = []
    start_index = 0

    # Loop until we reach the end of the original list
    while start_index < L:
        # Ensure there's enough remaining elements for at least a size of 3
        end_index = start_index + 3  # Default size of each chunk

        # Adjust the end_index if there aren't enough elements left
        if end_index > L:
            end_index = L

        chunks.append(original_list[start_index:end_index])

        # Move the start index forward by the length of the current chunk
        start_index += 3  # Ensure each chunk has at least 3 items

    return chunks

def make_xrs(s):
  cs1 = crystal.symmetry((10, 20, 30, 90, 90, 90), "P 1")
  sp1 = crystal.special_position_settings(cs1)
  scatterers1 = flex.xray_scatterer((
    xray.scatterer(s, (0.5, 0, 0)),
    xray.scatterer(s, (0, 0, 0))))
  return xray.structure(sp1, scatterers1)

def run_one_x(args):
  bcr.compute_tables(
    MinResolution = 0.996,
    MaxResolution = 10.1,
    scattering_table = "wk1995",
    TypesAtoms = args)

def run_one_e(args):
  bcr.compute_tables(
    MinResolution = 0.996,
    MaxResolution = 10.1,
    scattering_table = "electron",
    TypesAtoms = args)


if (__name__ == "__main__"):
  #
  e_list = []
  for l in std_labels:
    xrs = make_xrs(s=l)
    r = xrs.scattering_type_registry(table="electron")
    if len(list(r.unassigned_types()))==0:
      e_list.append(l)
    #xrs.structure_factors(d_min=2)
  #
  e_list = chunk_list(original_list=e_list)
  x_list = chunk_list(original_list=std_labels)
  #
  NPROC=120
  #
  #argss = []
  #for e in [["H","C"], ["O","N"], ["P","S"], ["Fe","Fe2+","Fe3+"], ["Zn","Zn2+"]]:
  #  argss.append(e)
  #print(e_list)
  #
  if(NPROC>1):
    stdout_and_results = easy_mp.pool_map(
      processes    = NPROC,
      fixed_func   = run_one_x,
      args         = x_list,
      func_wrapper = "buffer_stdout_stderr")
    stdout_and_results = easy_mp.pool_map(
      processes    = NPROC,
      fixed_func   = run_one_e,
      args         = e_list,
      func_wrapper = "buffer_stdout_stderr")
  else:
    for args in argss:
      o = run_one(args)
