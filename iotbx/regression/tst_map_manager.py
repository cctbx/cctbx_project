from __future__ import absolute_import, division, print_function
import os
from iotbx.data_manager import DataManager
from iotbx.map_manager import map_manager
from libtbx.test_utils import approx_equal

def test_01():

  # Source data (map and model)

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_ccp4 = os.path.join(data_dir, 'data',
                          'non_zero_origin_map.ccp4')
  data_pdb = os.path.join(data_dir, 'data',
                          'non_zero_origin_model.pdb')


  # Read in map data with data_manager
  dm = DataManager(['real_map'])
  dm.set_overwrite(True)

  # Next step uses map_manager to do the actual reading
  dm.process_real_map_file(data_ccp4)
  mm = dm.get_real_map()

  # Shift the origin of the map; starts at (100,100,100)
  print (mm.map_data().origin())
  assert mm.map_data().origin() == (100,100,100)
  assert mm.origin_shift_grid_units == (0,0,0)
  assert not mm.shifted()
  mm.shift_origin()
  assert mm.map_data().origin() == (0,0,0)
  assert mm.origin_shift_grid_units == (100,100,100)
  assert mm.shifted()
  mm.show_summary()

  # test cc_to_other_map
  assert mm.cc_to_other_map_manager(mm) == 1

  # test writing and reading file
  dm.write_real_map_file(mm, filename = 'test_map_manager.ccp4',
     overwrite = True)
  dm.process_real_map_file('test_map_manager.ccp4')
  new_mm=dm.get_real_map('test_map_manager.ccp4')
  os.remove('test_map_manager.ccp4')
  new_mm.shift_origin()
  # Check whether gridding and crystal_symmetry are similar for mm, new_mm
  assert new_mm.is_similar(mm)
  assert approx_equal(new_mm.map_data()[3125],mm.map_data()[3125])


  # test writing and reading file without shifting origin
  dm = DataManager(['real_map'])
  dm.set_overwrite(True)
  dm.process_real_map_file(data_ccp4)
  mm = dm.get_real_map()
  mm.show_summary()
  dm.write_real_map_file(mm, filename = 'test_map_manager.ccp4',
     overwrite = True)
  new_mm = map_manager('test_map_manager.ccp4')
  assert (new_mm.is_similar(mm))
  new_mm.shift_origin()
  assert (not new_mm.is_similar(mm))

  # get map_data
  dm = DataManager(['real_map'])
  dm.set_overwrite(True)
  dm.process_real_map_file(data_ccp4)
  mm = dm.get_real_map()
  mm.shift_origin()
  map_data = mm.map_data()
  assert approx_equal(map_data[15, 10, 19], 0.38, eps = 0.01)

  # get crystal_symmetry
  cs = mm.crystal_symmetry()
  assert approx_equal(cs.unit_cell().parameters()[0] , 22.41, eps = 0.01)

  # and full cell symmetry
  full_cs = mm.unit_cell_crystal_symmetry()
  assert approx_equal(full_cs.unit_cell().parameters()[0] ,
       149.4066, eps = 0.01)

  # write map directly:
  mm.write_map('test_direct.ccp4')

  # read back directly
  new_mm = map_manager('test_direct.ccp4')
  assert (not new_mm.is_similar(mm))

  new_mm.shift_origin()
  assert mm.is_similar(new_mm)
  assert approx_equal(new_mm.map_data()[3125],mm.map_data()[3125])

  # deep_copy
  new_mm = mm.deep_copy()
  assert new_mm.is_similar(mm)
  assert approx_equal(new_mm.map_data()[3125],mm.map_data()[3125])

  # deep_copy a map without shifting origin
  # Make a DataManager that can write a map coeffs file too
  dm = DataManager(['miller_array', 'real_map'])
  dm.set_overwrite(True)
  dm.process_real_map_file(data_ccp4)
  omm = dm.get_real_map()
  omm.show_summary()
  new_omm = omm.deep_copy()
  assert new_omm.is_similar(omm)
  assert (not new_omm.is_similar(mm))

  # customized_copy
  new_mm = mm.customized_copy(map_data = mm.map_data().deep_copy())
  assert new_mm.is_similar(mm)


  # Initialize with parameters
  mm_para = map_manager(
     unit_cell_grid =  mm.unit_cell_grid,
     unit_cell_crystal_symmetry =  mm.unit_cell_crystal_symmetry(),
     origin_shift_grid_units =  mm.origin_shift_grid_units,
     map_data = mm.map_data(),
     wrapping = False)
  assert mm_para.is_similar(mm)

  # Adjust origin and gridding:
  mm_read = map_manager(data_ccp4)
  mm_read.shift_origin()
  mm.show_summary()
  mm_read.show_summary()
  mm_read.set_original_origin_and_gridding((10, 10, 10),
       gridding = (100, 100, 100))
  mm_read.show_summary()
  assert (not mm_read.is_similar(mm))
  assert (mm_read.origin_is_zero())

  # Set program name
  mm_read.set_program_name('test program')
  assert mm_read.program_name == 'test program'

  # Set limitation
  mm_read.add_limitation('map_is_sharpened')
  assert mm_read.limitations == ['map_is_sharpened']

  # Add a label
  mm_read.add_label('TEST LABEL')
  assert mm_read.labels[0] == 'TEST LABEL'
  mm_read.write_map('map_with_labels.mrc')
  new_mm = map_manager('map_with_labels.mrc')
  assert 'TEST LABEL' in new_mm.labels
  assert new_mm.is_in_limitations('map_is_sharpened')
  assert new_mm.labels[0].find('test program')>-1

  # Remove backslashes
  mm_read.add_label('TEST \\LABEL')
  assert mm_read.labels[0] == 'TEST \\LABEL'
  mm_read.write_map('map_with_labels.mrc')
  new_mm = map_manager('map_with_labels.mrc')
  assert 'TEST /LABEL' in new_mm.labels
  assert new_mm.is_in_limitations('map_is_sharpened')
  assert new_mm.labels[0].find('test program')>-1

  # change the cell dimensions
  mm_read = map_manager(data_ccp4)
  mm_read.shift_origin()
  assert mm_read.is_similar(mm)
  assert approx_equal(mm_read.pixel_sizes(),
        (0.7470, 0.7231, 0.7374) , eps = 0.001)
  from cctbx import crystal
  new_uc_params = list(
    mm_read.unit_cell_crystal_symmetry().unit_cell().parameters())
  new_uc_params[0]+= 10
  new_cs = crystal.symmetry(new_uc_params, 1)
  mm_read.set_unit_cell_crystal_symmetry(new_cs)
  assert not mm_read.crystal_symmetry().is_similar_symmetry(
      mm.crystal_symmetry())
  assert not mm_read.is_similar(mm)
  mm_read.show_summary()
  assert approx_equal(mm_read.pixel_sizes(),
       (0.7970, 0.7231, 0.7374) , eps = 0.001)

  # Read a map directly
  mm_read = map_manager(data_ccp4)
  mm_read.shift_origin()
  assert mm_read.is_similar(mm)

  # Set log
  import sys
  mm.set_log(sys.stdout)

  # Add map_data
  new_mm = mm_read.customized_copy(map_data = mm.map_data().deep_copy())
  assert new_mm.is_similar(mm)

  # replace data
  new_mm.set_map_data(map_data = mm.map_data().deep_copy())
  assert new_mm.is_similar(mm)

  # create a full-sized map from this one
  mm_full_size = mm_read.deep_copy().as_full_size_map()
  assert not mm_full_size.is_similar(mm_read)
  print (mm_full_size.map_data().origin(), mm_read.map_data().origin())
  print (mm_full_size.map_data().all(), mm_read.map_data().all())

  # Apply a mask to edges of a map
  assert approx_equal(new_mm.map_data().as_1d().min_max_mean().max,
                          mm.map_data().as_1d().min_max_mean().max)
  assert approx_equal( (new_mm.map_data()[0], mm.map_data()[0]),
         (0.0, 0.0))
  new_mm.create_mask_around_edges(boundary_radius = 3)
  new_mm.soft_mask(soft_mask_radius = 3)
  assert approx_equal(new_mm.map_data().as_1d().as_1d().min_max_mean().max,
                          mm.map_data().as_1d().as_1d().min_max_mean().max)
  new_mm.apply_mask(set_outside_to_mean_inside = True)
  assert approx_equal( (new_mm.map_data()[0], mm.map_data()[0]),
         (0.0116108613152, 0.0))


  dm.process_real_map_file('test_map_manager.ccp4')
  new_mm = dm.get_real_map('test_map_manager.ccp4')
  new_mm.show_summary()
  assert (not new_mm.is_similar(mm))
  new_mm.shift_origin()
  new_mm.show_summary()
  assert new_mm.is_similar(mm)
  os.remove('test_map_manager.ccp4')

  # Check origin_shifts
  print (new_mm.origin_shift_grid_units)
  print (new_mm.shift_cart())
  assert approx_equal(new_mm.origin_shift_grid_units, (100, 100, 100))
  assert approx_equal(new_mm.shift_cart(),
       (-74.70333099365234, -72.30750274658205, -73.7437515258789))

  # Convert to map coeffs, write out, read back, convert back to map

  map_coeffs = mm.map_as_fourier_coefficients(d_min = 3)
  mtz_dataset = map_coeffs.as_mtz_dataset(column_root_label = 'F')
  mtz_object = mtz_dataset.mtz_object()
  dm.write_miller_array_file(mtz_object, filename = "map_coeffs.mtz")
  # Note these Fourier coeffs correspond to working map (not original position)

  array_labels = dm.get_miller_array_labels("map_coeffs.mtz")
  labels = array_labels[0]
  dm.get_reflection_file_server(filenames = ["map_coeffs.mtz"],
       labels = [labels])
  miller_arrays = dm.get_miller_arrays()
  new_map_coeffs = miller_arrays[0]
  mm_from_map_coeffs = mm.fourier_coefficients_as_map_manager(
      map_coeffs = new_map_coeffs)

  assert mm_from_map_coeffs.is_similar(mm)

  # Find map symmetry in a map
  data_d7 = os.path.join(data_dir, 'data',
                          'D7.ccp4')
  dm = DataManager(['real_map',  'model'])
  dm.process_real_map_file(data_d7)
  dm.process_model_file(data_pdb)
  mm = dm.get_real_map(data_d7)
  model = dm.get_model(data_pdb)
  mm.shift_origin()
  mm.set_original_origin_and_gridding(original_origin=(0,0,0))

  # Box it so it is not so easy to find symmetry
  from cctbx.maptbx.box import with_bounds
  box=with_bounds(mm,lower_bounds=(2,2,2),
   upper_bounds=(43,43,43))
  new_mm=box.map_manager()
  new_mm.find_map_symmetry(symmetry='d7',min_ncs_cc=0.8,
    include_helical_symmetry=False)
  ncs_obj=new_mm.ncs_object()
  assert ncs_obj is not None
  print ("NCS: ",new_mm.ncs_object().as_ncs_spec_string())
  another_mm = map_manager(
     unit_cell_grid =  new_mm.unit_cell_grid,
     unit_cell_crystal_symmetry =  new_mm.unit_cell_crystal_symmetry(),
     origin_shift_grid_units =  new_mm.origin_shift_grid_units,
     map_data = new_mm.map_data(),
     ncs_object = ncs_obj,
     wrapping = False)
  assert another_mm.is_similar(new_mm)
  assert ncs_obj.is_similar_ncs_object(another_mm.ncs_object())
  assert new_mm.is_similar(another_mm)

  # Resample with different gridding
  resampled_mm = another_mm.resample_on_different_grid(
    target_grid_spacing = .3)
  assert tuple(resampled_mm.unit_cell_grid) == (67, 67, 67)
  assert tuple(resampled_mm.origin_shift_grid_units) == (3 , 3, 3)

  # Reset crystal_symmetry of map with ncs object
  another_mm.set_unit_cell_crystal_symmetry(
     unit_cell_crystal_symmetry=another_mm.unit_cell_crystal_symmetry())
  assert not new_cs.is_similar_symmetry(
    another_mm.unit_cell_crystal_symmetry())
  shift_cart =  another_mm.shift_cart()
  another_mm.set_unit_cell_crystal_symmetry(unit_cell_crystal_symmetry=new_cs)
  assert shift_cart !=  another_mm.shift_cart()
  assert approx_equal(another_mm.ncs_object().shift_cart(), another_mm.shift_cart())

  # Get resolution
  assert approx_equal(new_mm.resolution(force=True, method='d99') ,
    3.73886)
  assert approx_equal(new_mm.resolution(force=True, method='d_min') ,
    0.888888888889)
  assert approx_equal(new_mm.resolution(force=True, method='d9') ,
    0.888888888889)
  assert approx_equal(new_mm.resolution(force=True, method='d99') ,
    3.73886)
  assert approx_equal(new_mm.resolution() ,
    3.73886)

  # Adjust model and ncs symmetry to match this map
  assert model.shift_cart()  is None
  assert not model.shifted()
  new_mm.set_model_symmetries_and_shift_cart_to_match_map(model)
  assert model.shifted()
  assert approx_equal (model.shift_cart() ,
       (-0.888888888888889, -0.8888888888888891, -0.888888888888889))
  sc = model.get_sites_cart()
  from scitbx.matrix import col
  sc_absolute = sc - col(model.shift_cart())
  assert approx_equal(new_mm.sites_cart_to_sites_cart_absolute(sc),
     sc_absolute)
  assert approx_equal(new_mm.sites_cart_absolute_to_sites_cart(sc_absolute),
     sc)
  assert new_mm.is_compatible_ncs_object(ncs_obj)
  ncs_obj.set_shift_cart((0,0,0))
  assert not new_mm.is_compatible_ncs_object(ncs_obj)

  new_mm.set_ncs_object_shift_cart_to_match_map(ncs_obj)
  new_mm.set_ncs_object(ncs_obj)
  assert new_mm.is_compatible_ncs_object(new_mm.ncs_object())
  new_mm.show_summary()

  new_mm.shift_origin(desired_origin=(11,1,1))
  print (new_mm.shift_cart(),new_mm.ncs_object().shift_cart())
  assert new_mm.is_compatible_ncs_object(new_mm.ncs_object())
  new_mm.shift_origin()
  assert new_mm.is_compatible_ncs_object(new_mm.ncs_object())

  # filter a map
  dm = DataManager()
  mm = dm.get_real_map(data_d7)

  low_pass_filtered=mm.deep_copy()
  low_pass_filtered.resolution_filter(d_min=2.5)

  high_pass_filtered=mm.deep_copy()
  high_pass_filtered.resolution_filter(d_max=2.5)

  gaussian=mm.deep_copy()
  gaussian.gaussian_filter(smoothing_radius=1)

  randomize=mm.deep_copy()
  randomize.randomize(random_seed=1)
  assert approx_equal(randomize.map_map_cc(mm),0.7,.10)

  binary=mm.deep_copy()
  binary.binary_filter(threshold=0.5)

  scaled=mm.deep_copy()
  scaled.scale_map(10)

  assert approx_equal(
     (mm.map_data().as_1d()[1073],low_pass_filtered.map_data().as_1d()[1073],
       high_pass_filtered.map_data().as_1d()[1073],
       gaussian.map_data().as_1d()[1073],binary.map_data().as_1d()[1073],
       scaled.map_data().as_1d()[1073]),
      (0.0171344596893,0.0227163900537,-0.0072717454565,0.0149086679298,0.0,
        0.171344596893))

  info=mm.get_density_along_line((5,5,5),(10,10,10))
  assert approx_equal([info.along_density_values[4]]+list(info.along_sites[4]),
    [-0.562231123447 , 8.0, 8.0, 8.0])
  from iotbx.map_model_manager import map_model_manager
  extra_map_manager_id_list = ["low_pass_filtered",
    # "high_pass_filtered",
    "gaussian",
    "binary"]

  expected_cc= [ 1,
     # -0.038380049155,
     0.975273714961,
     0.866785173385]
  mam=map_model_manager(
    map_manager=mm,
    extra_map_manager_list =  [low_pass_filtered,
     # high_pass_filtered,
     gaussian,
     binary],
    extra_map_manager_id_list = extra_map_manager_id_list,)
  for other_id,cc in zip(extra_map_manager_id_list,expected_cc):
   assert approx_equal(cc,
      mam.map_map_cc(map_id='map_manager',other_map_id=other_id), eps = 0.05 )

# this test requires the solve_resolve module
def test_02():

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_d7 = os.path.join(data_dir, 'data', 'D7.ccp4')

  # find_separated atoms in a map
  dm = DataManager()
  mm = dm.get_real_map(data_d7)
  sites_cart = mm.trace_atoms_in_map(dist_min=1,n_atoms=10)
  assert sites_cart.size() == 10 # Note: zero if not available

if (__name__  ==  '__main__'):
  test_01()
  test_02()
  print ("OK")
