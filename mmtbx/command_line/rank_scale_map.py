"""Compute rank-scaled (histogram equalized) map"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.rank_scale_map

from cctbx import maptbx
from libtbx.utils import Sorry
import sys, time
from iotbx.data_manager import DataManager

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
  dm = DataManager()
  dm.set_overwrite(True)
  map_manager = dm.get_real_map(args[0])
  map_manager.shift_origin()

  cs = map_manager.crystal_symmetry()
  m = map_manager.map_data().as_double()
  # show general statistics
  show_overall_statistics(m=m, header="Map basic info (%s):"%args[0])
  # HE
  m_he = maptbx.volume_scale(map = m,  n_bins = 10000).map_data()
  show_overall_statistics(m=m_he, header="Rank-scaled (HE) map info:")
  #
  file_name=args[0]+"_rank_scaled.ccp4"
  he_map_manager = map_manager.customized_copy(map_data = m_he)
  he_map_manager.add_label("Histogram-equalized map")
  dm.write_real_map_file(he_map_manager, file_name)

if (__name__ == "__main__"):
  t0 = time.time()
  run(sys.argv[1:])
  print("Time: %-8.3f"%(time.time()-t0))
  print("All done.")

