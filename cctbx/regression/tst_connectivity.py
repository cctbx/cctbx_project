from __future__ import division
from scitbx.array_family import flex
import iotbx.pdb
from cctbx import maptbx
from libtbx.test_utils import approx_equal

def exercise():
  pdb_str="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1       5.000   7.000  -9.000  1.00 20.00           C
END
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  fc = xrs.structure_factors(d_min = 1., algorithm = "direct").f_calc()
  fft_map = fc.fft_map(resolution_factor = 0.2)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  co = maptbx.connectivity(map_data = map_data, threshold=1)
  map_result = co.result()
  print map_result.as_1d().min_max_mean().as_tuple()
  print (map_result==1).count(True), (map_result==1).count(False)
  print map_result.size()

if __name__ == "__main__" :
  exercise()
