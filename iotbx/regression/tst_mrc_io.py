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

def test_foreign_byte_order():
  ''' A map written with the opposite endianness must read back with the
  same values. On big-endian s390x the little-endian regression maps read
  as byte-swapped garbage, which usually contains NaNs and so triggers the
  "contains 'nan'" Sorry (cctbx_project issue #1216). '''
  import sys
  native = '<' if sys.byteorder == 'little' else '>'
  foreign = '>' if sys.byteorder == 'little' else '<'
  d = (np.arange(10*20*30, dtype=np.float32) / 7.0).reshape(30, 20, 10)
  for byte_order, file_name in [(native, 'native_byte_order.mrc'),
                                (foreign, 'foreign_byte_order.mrc')]:
    with mrcfile.new(file_name, overwrite=True) as dst:
      dst.set_data(d.astype(np.dtype(byte_order + 'f4')))
      dst.header.cella = (10, 20, 30)
    with mrcfile.open(file_name) as src:
      assert np.dtype(src.data.dtype).newbyteorder('=') == np.dtype('f4')
      assert src.data.dtype.byteorder in (byte_order, '=')
      assert np.array_equal(src.data, d)

  m_native = iotbx.mrcfile.map_reader(file_name='native_byte_order.mrc')
  m_foreign = iotbx.mrcfile.map_reader(file_name='foreign_byte_order.mrc')
  assert m_native.map_data().all() == m_foreign.map_data().all()
  assert m_native.map_data().all() == (10, 20, 30), m_native.map_data().all()
  assert math.isclose(flex.max(m_native.map_data()), float(d.max()),
    rel_tol=1e-6)
  assert flex.sum(flex.abs(m_native.map_data() - m_foreign.map_data())) == 0

  dm = DataManager()
  dm.process_real_map_file('foreign_byte_order.mrc')
  mm = dm.get_real_map('foreign_byte_order.mrc')
  assert flex.sum(flex.abs(mm.map_data() - m_native.map_data())) == 0

def run():
  test_nan_in_mrc()
  test_float16()
  test_foreign_byte_order()

if(__name__ == "__main__"):
  run()
