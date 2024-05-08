from __future__ import absolute_import, division, print_function
import iotbx.mrcfile
from cctbx.array_family import flex
from cctbx import crystal
from libtbx.utils import Sorry
import math

def run():
  m = flex.double(flex.grid(10,20,30), 1)
  m[1,2,3] = float('nan')
  assert math.isnan(m[1,2,3])
  cs = crystal.symmetry((10,20,30, 90,90,90), "P1")
  iotbx.mrcfile.write_ccp4_map(
    file_name   = "map_with_nan.mrc",
    unit_cell   = cs.unit_cell(),
    space_group = cs.space_group(),
    map_data    = m,
    labels      = flex.std_string(["Some text"]))
  exception = None
  try: o = iotbx.mrcfile.map_reader(file_name="map_with_nan.mrc")
  except Sorry as e: exception = str(e)
  assert exception == "Input map file contains 'nan':  not a valid map file"

if(__name__ == "__main__"):
  run()
