from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
from cctbx import maptbx
from libtbx.test_utils import approx_equal
import sys

def exercise1():
  """
  Exercise to see how much values on generated map are different on
  different platforms
  """
  pdb_str="""
CRYST1   10.000  10.000   10.000  90.00  90.00  90.00 P 4
HETATM    1  C    C      1       2.000   2.000   2.000  1.00 20.00           C
HETATM    2  C    C      2       4.000   4.000   4.000  1.00 20.00           C
END
"""

  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  # xrs.show_summary()
  d_min = 1.
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  symmetry_flags = maptbx.use_space_group_symmetry
  fftmap = fc.fft_map(symmetry_flags = symmetry_flags)
  rmup = fftmap.real_map_unpadded()
  nx, ny, nz = rmup.accessor().all() # 32,32,32
  if not approx_equal((nx,ny,nz), (32,32,32)):
    print('Even map dimensions are different???:', file=sys.stderr)
    print('  ', (nx, ny, nz), file=sys.stderr)
    print('  ', (32,32,32), file=sys.stderr)

  # print('dimensions:', nx,ny,nz)
  # min, max and mean values on the map
  mmm = rmup.as_1d().min_max_mean().as_tuple()

  # Mac results:
  mac_results = [[(6,6,7),   1411.8],
                 [(6,6,8),   1026.0],
                 [(0,0,0),    -23.3],
                 [(31,31,31), -23.9],
                ]

  if not approx_equal(mmm, (-25.5, 1531.7, 0.), 1):
    print('approx_equal failed:', file=sys.stderr)
    print('  ', mmm, file=sys.stderr)
    print('  ', (-25.5, 1531.7, 0.), file=sys.stderr)
  for coord, value in mac_results:
    if not approx_equal(rmup[coord], value, 1):
      print('approx_equal failed:', file=sys.stderr)
      print('  ', rmup[coord], file=sys.stderr)
      print('  ', value, file=sys.stderr)

if __name__ == "__main__":
  t0 = time.time()
  exercise1()
  print("OK time =%8.3f"%(time.time() - t0))
