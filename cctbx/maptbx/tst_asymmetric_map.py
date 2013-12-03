from __future__ import division
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost.python
ext = boost.python.import_ext("cctbx_asymmetric_map_ext")
from cctbx_asymmetric_map_ext import *
from cctbx.array_family import flex

def run():
  group = space_group_info("P21");
  print "Group: ", group.type().lookup_symbol()
  elements = ('C', 'N', 'O', 'H')*11
  struc = random_structure.xray_structure(
    space_group_info = group,
    volume_per_atom = 25.,
    general_positions_only = False,
    elements = elements,
    min_distance = 1.0)
  print "Structure: ", struc
  struc.show_summary()
  d_min = 2.
  fc = struc.structure_factors(d_min=d_min).f_calc()
  print "FC: ", fc
  fc.show_summary()
  fftmap = fc.fft_map()
  print "Map: ", fftmap
  fftmap.statistics().show_summary()
  amap = asymmetric_map(struc.space_group().type(), fftmap.real_map())
  print "Asu map: ", amap
  afc = amap.structure_factors(fc.indices())
  print "ASU FC: ", afc
  afftmap = amap.map_for_fft()
  print "ASU map for fft: ", afftmap
  adata = amap.data()
  print "ASU data: ", adata
  df = flex.abs(afc - fc.data())
  r1 = flex.sum(df) / flex.sum(flex.abs(fc.data()))
  print "R1: ", r1
  assert r1 < 1.e-4

if (__name__ == "__main__"):
  run()
