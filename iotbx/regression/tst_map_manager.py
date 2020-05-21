from __future__ import absolute_import, division, print_function
import os
from iotbx.data_manager import DataManager

import libtbx.load_env
import libtbx.phil

def test_01():

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_ccp4 = os.path.join(data_dir, 'data',
                          'non_zero_origin_map.ccp4')
  data_pdb = os.path.join(data_dir, 'data',
                          'non_zero_origin_map.ccp4')

  dm = DataManager(['miller_array','real_map', 'phil'])
  dm.set_overwrite(True)
  dm.process_real_map_file(data_ccp4)

  # test writing and reading file
  mm = dm.get_real_map()
  mm.shift_origin()
  mm.show_summary()
  dm.write_map_with_map_manager(mm, filename='test.ccp4', overwrite=True)


  dm.process_real_map_file('test.ccp4')
  new_mm=dm.get_real_map('test.ccp4')
  new_mm.show_summary()
  assert (not new_mm.is_similar(mm))
  new_mm.shift_origin()
  new_mm.show_summary()
  assert new_mm.is_similar(mm)
  os.remove('test.ccp4')

  # Convert to map coeffs, write out, read back, convert back to map

  map_coeffs = mm.map_as_fourier_coefficients(high_resolution = 3)
  mtz_dataset = map_coeffs.as_mtz_dataset(column_root_label='F')
  mtz_object=mtz_dataset.mtz_object()
  dm.write_miller_array_file(mtz_object, filename="map_coeffs.mtz")
  # Note these Fourier coeffs correspond to working map (not original position)

  array_labels=dm.get_miller_array_labels("map_coeffs.mtz")
  labels=array_labels[0]
  dm.get_reflection_file_server(filenames=["map_coeffs.mtz"],labels=[labels])
  miller_arrays=dm.get_miller_arrays()
  new_map_coeffs=miller_arrays[0]
  map_data_from_map_coeffs=mm.fourier_coefficients_as_map(
      map_coeffs=new_map_coeffs)

  mm_from_map_coeffs=mm.customized_copy(map_data=map_data_from_map_coeffs)
  assert mm_from_map_coeffs.is_similar(mm)


if (__name__ == '__main__'):
  test_01()
  print ("OK")
