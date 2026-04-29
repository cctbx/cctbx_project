from __future__ import absolute_import, division, print_function
import iotbx.mrcfile
from cctbx.array_family import flex
from cctbx import crystal
from iotbx.data_manager import DataManager
from libtbx.utils import Sorry
import math
import mrcfile
import numpy as np

def test_nan_in_mrc():
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

def test_float16():
  ''' Verify that we can read and write float16 MRC files '''
  m = flex.double(flex.grid(10,20,30), 1)
  cs = crystal.symmetry((10,20,30, 90,90,90), "P1")
  iotbx.mrcfile.write_ccp4_map(
    file_name   = "map_with_nan.mrc",
    unit_cell   = cs.unit_cell(),
    space_group = cs.space_group(),
    map_data    = m,
    labels      = flex.std_string(["Some text"]))

  with mrcfile.open('map_with_nan.mrc') as src:
    d = src.data.copy()
    h = src.header.copy()
    eh = src.extended_header.copy()

  with mrcfile.new('test_float16.mrc', overwrite=True) as dst:
    dst.set_data(d.astype(np.float16, copy=False))
    dst.header.cella = h.cella
    dst.set_extended_header(eh)

  dm = DataManager()
  dm.process_real_map_file('test_float16.mrc')

def run():
  test_nan_in_mrc()
  test_float16()

if(__name__ == "__main__"):
  run()
