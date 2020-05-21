from __future__ import absolute_import, division, print_function
import os
from iotbx.data_manager import DataManager
from iotbx.map_manager import map_manager
from libtbx.test_utils import approx_equal

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

  # get map_data
  map_data=mm.map_data()
  assert approx_equal(map_data[15,10,19], 0.38,eps=0.01)

  # get crystal_symmetry
  cs=mm.crystal_symmetry()
  assert approx_equal(cs.unit_cell().parameters()[0] ,22.41,eps=0.01)

  # and full cell symmetry
  full_cs=mm.unit_cell_crystal_symmetry()
  assert approx_equal(full_cs.unit_cell().parameters()[0] ,149.4066,eps=0.01)

  # get shift_tracker:
  st=mm.map_shift_tracker()
  assert st.is_similar(mm)


  # write map directly:
  mm.write_map('test.ccp4')

  # read back directly
  new_mm=map_manager('test.ccp4')
  assert (not new_mm.is_similar(mm))

  new_mm.shift_origin()
  assert mm.is_similar(new_mm)

  # deep_copy
  new_mm=mm.deep_copy()
  assert new_mm.is_similar(mm)

  # customized_copy
  new_mm=mm.customized_copy(map_data=mm.map_data().deep_copy())
  assert new_mm.is_similar(mm)


  # Initialize with existing map_manager
  mm_with_existing=map_manager(map_manager_object=mm)
  assert mm_with_existing.is_similar(mm)

  # Initialize with existing map_manager and map_data
  mm_with_existing=map_manager(map_manager_object=mm,map_data=mm.map_data().deep_copy())
  assert mm_with_existing.is_similar(mm)

  # Initialize with nothing:
  mm_empty=map_manager()
  assert (not mm_empty.is_similar(mm))


  # Initialize with parameters
  mm_para=map_manager(
     unit_cell_grid= mm.unit_cell_grid,
     unit_cell_parameters= mm.unit_cell_parameters,
     space_group_number= mm.space_group_number,
     origin_grid_units= mm.origin_shift_grid_units,
     working_map_n_xyz= mm.working_map_n_xyz,
     map_data=mm._map_data,
     high_resolution=mm.high_resolution)
  assert mm_para.is_similar(mm)

  # Adjust origin and gridding:
  mm_read=map_manager()
  mm_read.read_map(data_ccp4)
  mm_read.set_origin_and_gridding((10,10,10),gridding=(100,100,100))
  assert (not mm_read.is_similar(mm))
  assert (not mm_read.already_shifted())

  # Adjust origin and gridding should fail if origin already shifted:
  mm_read=map_manager()
  mm_read.read_map(data_ccp4)
  mm_read.shift_origin()
  mm_read.set_origin_and_gridding((10,10,10),gridding=(100,100,100))
  assert (mm_read.is_similar(mm))  # not shifted as it failed
  assert (mm_read.already_shifted())

  # Read a map directly
  mm_read=map_manager()
  mm_read.read_map(data_ccp4)
  mm_read.shift_origin()
  assert mm_read.is_similar(mm)


  # Add map_data
  mm_read.add_map_data(map_data=mm.map_data().deep_copy())
  assert mm_read.is_similar(mm)



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

  # Read mtz data directly
  mm_read=mm.map_shift_tracker()
  mm_read.read_mtz("map_coeffs.mtz")
  assert mm_read.is_similar(mm)


  # Generate some map data
  empty_mm=map_manager()
  empty_mm.generate_map()
  assert approx_equal(mm.map_data()[15,10,19],0.38,eps=0.01)

if (__name__ == '__main__'):
  test_01()
  print ("OK")
