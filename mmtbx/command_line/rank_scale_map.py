from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.rank_scale_map

from cctbx import maptbx
import iotbx.mrcfile
from cctbx import crystal
from libtbx.utils import Sorry
from scitbx.array_family import flex
import sys, time

def show_overall_statistics(m, header):
  s = maptbx.more_statistics(m)
  print(header)
  print("  min/max/mean: %6.4f %6.4f %6.4f"%(s.min(), s.max(), s.mean()))
  print("  kurtosis    : %6.4f" % s.kurtosis())
  print("  skewness    : %6.4f" % s.skewness())
  print("  sigma       : %6.4f" % s.sigma())

def show_citation():
  print("-"*79)
  msg = """Compute rank-scaled (histogram equalized) map.

Input: CCP4 formatted map file.

Example:
  phenix.rank_scale_map map.ccp4

Citation:
  Map comparison and statistics. For details see:
  Acta Cryst. (2014). D70, 2593-2606
  Metrics for comparison of crystallographic maps
  A. Urzhumtsev, P. V. Afonine, V. Y. Lunin, T. C. Terwilliger and P. D. Adams"""
  print(msg)
  print("-"*79)

def run(args):
  show_citation()
  if(len(args)!=1): raise Sorry("Need to provide CCP4 formatted map file.")
  # map
  try:
    ccp4_map = iotbx.mrcfile.map_reader(file_name=args[0])
  except Exception: # XXX should probably be RuntimeError?
    raise Sorry("Not a valid file (provide CCP4 formatted map file).")
  cs = crystal.symmetry(ccp4_map.unit_cell().parameters(),
    ccp4_map.space_group_number)
  m = ccp4_map.data.as_double()
  # show general statistics
  show_overall_statistics(m=m, header="Map basic info (%s):"%args[0])
  # HE
  m_he = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
  show_overall_statistics(m=m_he, header="Rank-scaled (HE) map info:")
  #
  iotbx.mrcfile.write_ccp4_map(
    file_name=args[0]+"_rank_scaled.ccp4",
    unit_cell=cs.unit_cell(),
    space_group=cs.space_group(),
    #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
    #gridding_last=n_real,  # This causes a bug (map gets shifted)
    map_data=m_he,
    labels=flex.std_string([""]))

if (__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print("Time: %-8.3f"%(time.time()-t0))
  print("All done.")
