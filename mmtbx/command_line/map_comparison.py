from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_comparison

from cctbx import maptbx
import iotbx.ccp4_map
from cctbx import crystal
from libtbx.utils import Sorry
from scitbx.array_family import flex
import sys

def show_overall_statistics(m, header):
  s = maptbx.more_statistics(m)
  print header
  print "  min/max/mean: %6.4f %6.4f %6.4f"%(s.min(), s.max(), s.mean())
  print "  kurtosis    : %6.4f" % s.kurtosis()
  print "  skewness    : %6.4f" % s.skewness()
  print "  sigma       : %6.4f" % s.sigma()

def show_citation():
  print "-"*79
  msg = """Map comparison and statistics. For details see:
  Acta Cryst. (2014). D70, 2593-2606
  Metrics for comparison of crystallographic maps
  A. Urzhumtsev, P. V. Afonine, V. Y. Lunin, T. C. Terwilliger and P. D. Adams"""
  print msg
  print "-"*79

def run(args):
  show_citation()
  if(len(args)!=2): raise Sorry("Need to provide two CCP4 formatted map files.")
  # map 1
  try:
    ccp4_map_1 = iotbx.ccp4_map.map_reader(file_name=args[0])
  except Exception: # XXX should probably be RuntimeError?
    raise Sorry("Not a valid file (provide CCP4 formatted map file).")
  cs_1 = crystal.symmetry(ccp4_map_1.unit_cell().parameters(),
    ccp4_map_1.space_group_number)
  m1 = ccp4_map_1.data.as_double()
  # map 2
  try:
    ccp4_map_2 = iotbx.ccp4_map.map_reader(file_name=args[1])
  except Exception: # XXX see above
    raise Sorry("Not a valid file (provide CCP4 formatted map file).")
  cs_2 = crystal.symmetry(ccp4_map_2.unit_cell().parameters(),
    ccp4_map_2.space_group_number)
  m2 = ccp4_map_2.data.as_double()
  # sanity checks
  assert cs_1.is_similar_symmetry(cs_2)
  assert m1.accessor().all() == m2.accessor().all()
  assert m1.accessor().focus() == m2.accessor().focus()
  assert m1.accessor().origin() == m2.accessor().origin()
  # show general statistics
  show_overall_statistics(m=m1, header="Map 1 (%s):"%args[0])
  show_overall_statistics(m=m2, header="Map 2 (%s):"%args[1])
  print "CC, input maps: %6.4f" % flex.linear_correlation(
    x = m1.as_1d(), y = m2.as_1d()).coefficient()
  # compute CCpeak
  m1_he = maptbx.volume_scale(map = m1,  n_bins = 10000).map_data()
  m2_he = maptbx.volume_scale(map = m2,  n_bins = 10000).map_data()
  print "CC, quantile rank-scaled (histogram equalized) maps: %6.4f" % \
    flex.linear_correlation(x = m1_he.as_1d(), y = m2_he.as_1d()).coefficient()
  print "cutoff  CCpeak"
  for cutoff in [i/100. for i in range(0,100,5)]+[0.99, 1.0]:
    print "%3.2f   %7.4f" % (cutoff,
      maptbx.cc_peak(map_1=m1_he, map_2=m2_he, cutoff=cutoff))
  # actual calcs
  h1 = maptbx.histogram(map=m1, n_bins=10000)
  h2 = maptbx.histogram(map=m2, n_bins=10000)
  print "Map 1 (%s)     Map 2 (%s)"%(args[0],args[1])
  print "(map_vale,cdf,frequency) <> (map_vale,cdf,frequency)"
  for a1,c1,v1, a2,c2,v2 in zip(h1.arguments(), h1.c_values(), h1.values(),
                                h2.arguments(), h2.c_values(), h2.values()):
    print "(%9.5f %9.5f %9.5f) (%9.5f %9.5f %9.5f)"%(a1,c1,v1, a2,c2,v2)

if (__name__ == "__main__"):
  run(sys.argv[1:])
