from __future__ import absolute_import, division, print_function

from iotbx.file_reader import any_file
from cctbx.maptbx import shift_origin_if_needed
import sys

def run(args):
  map_fname = args[0]
  af = any_file(map_fname)
  assert af.file_type == "ccp4_map"
  ccp4_map = af.file_content
  print("origin:", ccp4_map.origin)
  # see how access this info in cctbx_project/iotbx/ccp4_map/__init__.py: def show_summary
  print("summary:", ccp4_map.show_summary())
  xc = yc = zc = 1 # real coordinates
  # fractional coordinates
  xf,yf,zf = ccp4_map.unit_cell().fractionalize([xc,yc,zc])
  print( "map value:", ccp4_map.map_data().eight_point_interpolation([xf, yf, zf]))
  shifted_map_data = shift_origin_if_needed(ccp4_map.map_data()).map_data
  print( "map value on shifted map:", shifted_map_data.eight_point_interpolation([xf, yf, zf]))
  print( "shifted origin:", shifted_map_data.origin())

if __name__ == '__main__':
  run(sys.argv[1:])
