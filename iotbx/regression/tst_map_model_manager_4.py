from __future__ import absolute_import, division, print_function
import sys, os
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal


data_dir = os.path.dirname(os.path.abspath(__file__))
random_map_a = os.path.join(data_dir, 'data', "mmma.ccp4")
random_map_b = os.path.join(data_dir, 'data', "mmmb.ccp4")

def get_map_model_managers():
  # Set up source data
  from iotbx.map_model_manager import map_model_manager
  mmm = map_model_manager()
  mmm.generate_map(wrapping=True, d_min=3)
  second_model=mmm.model().deep_copy()
  mmm.box_all_maps_with_bounds_and_shift_origin(
    lower_bounds=[-1,-1,0], upper_bounds=[30, 40, 32])
  mmm.map_manager().set_wrapping(False)
  mmm.set_log(sys.stdout)

  #now get a second one

  from scitbx import matrix
  r=matrix.sqr((-0.8090,0.5878,0.0000,
   -0.5878,-0.8090,-0.0000,
   0.0000,-0.0000,1.0000,))
  t = matrix.col((100,0,0))
  new_sites_cart = r.elems*mmm.model().get_sites_cart() + t.elems
  second_model.set_sites_cart(new_sites_cart)
  from cctbx.maptbx.box import shift_and_box_model
  second_model = shift_and_box_model(second_model)

  second_mmm = map_model_manager(model=second_model)
  second_mmm.generate_map(model=second_model,wrapping=True, d_min=3)
  second_mmm.box_all_maps_around_model_and_shift_origin(box_cushion=10)
  second_mmm.map_manager().set_wrapping(False)
  second_mmm.set_log(sys.stdout)

  print(mmm.model().get_sites_cart()[0])
  print(second_mmm.model().get_sites_cart()[0])
  return mmm, second_mmm

def exercise( out = sys.stdout):


  # test shift_aware_rt

  mmm1, mmm2 = get_map_model_managers()
  initial_shift_aware_rt_info= mmm1.shift_aware_rt_to_superpose_other(mmm2)
  initial_rt_info=initial_shift_aware_rt_info.working_rt_info(from_obj=mmm2,to_obj=mmm1)

  model_2=mmm2.model().apply_selection_string("resseq 222:235")
  mmm2.set_model(model_2)
  shift_aware_rt_info= mmm1.shift_aware_rt_to_superpose_other(mmm2)
  rt_info=shift_aware_rt_info.working_rt_info(from_obj=mmm2,to_obj=mmm1)
  assert shift_aware_rt_info.is_similar(initial_shift_aware_rt_info,tol=0.002)

  shift_aware_rt = mmm1.shift_aware_rt(working_rt_info=rt_info,
     from_obj = mmm2, to_obj = mmm1)

  shift_aware_rt = mmm1.map_manager().shift_aware_rt(working_rt_info=rt_info,
     from_obj = mmm2, to_obj = mmm1)
  print (mmm1, mmm2)
  sites_cart_2 = mmm2.model().get_sites_cart()
  mapped_sites_cart = shift_aware_rt.apply_rt(sites_cart=sites_cart_2,
    from_obj=mmm2, to_obj=mmm1)
  assert approx_equal(mapped_sites_cart,mmm1.model().apply_selection_string("resseq 222:235").get_sites_cart(), eps=0.01)
  working_rt_info = shift_aware_rt.working_rt_info(from_obj=mmm2, to_obj=mmm1)
  mapped_sites_cart =working_rt_info.r.elems * mmm2.model().get_sites_cart() + working_rt_info.t.elems
  assert approx_equal(mapped_sites_cart,mmm1.model().apply_selection_string("resseq 222:235").get_sites_cart(), eps=0.01)

  inverse_shift_aware_rt = shift_aware_rt.inverse()
  mapped_sites_cart = inverse_shift_aware_rt.apply_rt(sites_cart=mmm1.model().apply_selection_string("resseq 222:235").get_sites_cart(),from_obj=mmm1,to_obj=mmm2)
  assert approx_equal(mapped_sites_cart,mmm2.model().get_sites_cart(), eps=0.01)

  mmm1, mmm2 = get_map_model_managers()

  # get r,t to map mmm2 model on mmm1 model
  shift_aware_rt_info= mmm1.shift_aware_rt_to_superpose_other(mmm2)
  rt_info=shift_aware_rt_info.working_rt_info(from_obj=mmm2,to_obj=mmm1)

  # get mmm2 map superimposed on mmm1 map (in region where it is defined, zero
  #   outside that region)

  new_mm = mmm1.superposed_map_manager_from_other(other=mmm2)
  new_mm.write_map('super.ccp4')
  mmm1.write_map('orig.ccp4')
  mmm1.write_model('orig.pdb')

  new_mm = mmm1.superposed_map_manager_from_other(other=mmm2,
    selection_string="resseq 221:225")
  assert approx_equal(new_mm.map_map_cc(mmm1.map_manager()),0.994645868918,eps=0.01)
  new_mm.write_map('super_221-225.ccp4')

  new_mm = mmm1.superposed_map_manager_from_other(other=mmm2,
     working_rt_info = rt_info)
  assert approx_equal(new_mm.map_map_cc(mmm1.map_manager()),0.994645868918,eps=0.01)
  new_mm.write_map('super_221-225.ccp4')

  # get a local resolution map (this one should look pretty constant!)
  mmm1.set_resolution(3)
  mmma = mmm1.deep_copy()
  model = mmm1.model()
  mmma.remove_model_by_id('model')
  mmmb = mmma.deep_copy()

  mmma.map_manager().randomize(random_seed=23412,d_min=3,high_resolution_fourier_noise_fraction=10,low_resolution_noise_cutoff=5)
  mmmb.map_manager().randomize(random_seed=887241,d_min=3,high_resolution_fourier_noise_fraction=10,low_resolution_noise_cutoff=5)

  #  We are going to read in these maps so that we have a constant value
  from iotbx.data_manager import DataManager
  dm = DataManager()
  mmma_map = dm.get_real_map(random_map_a)
  mmmb_map = dm.get_real_map(random_map_b)
  mmma_map.shift_origin()
  mmmb_map.shift_origin()
  mmma.add_map_manager_by_id(map_id="map_manager", map_manager=mmma_map)
  mmmb.add_map_manager_by_id(map_id="map_manager", map_manager=mmmb_map)

  assert approx_equal(mmma.map_manager().map_map_cc(mmmb.map_manager()),
    0.16,0.10)
  from iotbx.map_model_manager import map_model_manager
  model.set_b_iso(flex.double(model.get_sites_cart().size(),0))
  local_mmm=map_model_manager(map_manager_1=mmma.map_manager(),
    map_manager_2=mmmb.map_manager(), model = model)
  local_mmm.set_resolution(3)
  local_mmm.local_fsc()

  from iotbx.data_manager import DataManager
  dm = DataManager()
  dm.set_overwrite(True)

  cc_before = local_mmm.map_model_cc()
  print ("Working with randomized maps cc = ",cc_before)
  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  cc_before = dc.map_model_cc()
  dc.half_map_sharpen(n_bins=15)
  cc_after = dc.map_model_cc(map_id='map_manager_scaled')
  print("CC before, after half map sharpen: ",cc_before,cc_after)
  assert approx_equal((cc_before,cc_after), (0.80, 0.80), eps=0.10)

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  cc_before = dc.map_model_cc()
  dc.model_sharpen(n_bins=15, local_sharpen=False, anisotropic_sharpen=False,
    optimize_b_eff=False)
  cc_after = dc.map_model_cc(map_id='map_manager_scaled')
  print("CC before, after std model sharpen: ",cc_before,cc_after)
  assert approx_equal ((cc_before,cc_after), (0.80,0.90), eps =0.10)
  model_sharpened_mm = dc.get_map_manager_by_id(map_id='map_manager_scaled')

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  cc_before = dc.map_model_cc()
  dc.model_sharpen(local_sharpen=True,n_boxes=1,n_bins=15)
  cc_after = dc.map_model_cc(map_id='map_manager_scaled')
  print("CC before, after local model sharpen n_boxes=1: ",cc_before,cc_after)
  assert approx_equal ((cc_before,cc_after), (0.80,0.90), eps =0.10)
  model_sharpened_mm = dc.get_map_manager_by_id(map_id='map_manager_scaled')

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  dc.add_map_manager_by_id(model_sharpened_mm,'external_map')
  cc_before = dc.map_map_cc(map_id='map_manager',other_map_id='external_map')
  dc.external_sharpen(n_bins=15,map_id_external_map='external_map')
  cc_after = dc.map_map_cc(map_id='map_manager_scaled',other_map_id='external_map')
  print("CC before, after external sharpen n_boxes=1: ",cc_before,cc_after)
  assert approx_equal ((cc_before,cc_after), (0.7,0.95),eps=0.10)

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  dc.add_map_manager_by_id(model_sharpened_mm,'external_map')
  cc_before = dc.map_map_cc(map_id='map_manager',other_map_id='external_map')
  dc.external_sharpen(local_sharpen=True,n_boxes=1,
     n_bins=15,map_id_external_map='external_map')
  cc_after = dc.map_map_cc(map_id='map_manager_scaled',other_map_id='external_map')
  print("CC before, after external sharpen local n_boxes=1: ",cc_before,cc_after)
  assert approx_equal ((cc_before,cc_after), (0.70,0.95),eps=0.10)


  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  dc._local_sharpen(map_id_scaled_list=['map_manager_scaled'], map_id_to_be_scaled_list=['map_manager'],
     n_bins=15, n_boxes=1)
  cc = dc.map_model_cc()
  assert approx_equal (cc, 0.80, eps=0.1)

  # create a mask around density
  dc.create_mask_around_density(soft_mask=False)
  count = dc.get_map_manager_by_id('mask').map_data().count(1)
  print(count)
  assert 8000 < count < 14000
  dc.expand_mask(buffer_radius = 2)
  count = dc.get_map_manager_by_id('mask').map_data().count(1)
  print(count)
  assert count == 1

  # Test mask and map info functions
  mmm1, mmm2 = get_map_model_managers()
  mmm1.create_mask_around_density(soft_mask=False)
  mask_info = mmm1.mask_info(quiet=True)
  map_info = mmm1.map_info()
  mask_info_by_id = mmm1.mask_info(mask_id = 'mask',quiet=True)
  map_info_by_id = mmm1.map_info(map_id = 'map_manager')
  assert mask_info() == mask_info_by_id()
  assert map_info() == map_info_by_id()
  assert approx_equal(mask_info.fraction_marked, 0.207070707071, eps=0.01)

  assert approx_equal(map_info.fraction_above_sigma_cutoff, 0.0577876984127)

  # create a spherical mask around a point
  print("Spherical masks", )
  dc = mmm1.deep_copy()
  dc.mask_info(quiet=True)
  print (dc.mask_info(quiet=True).marked_points )
  assert dc.mask_info(quiet=True).marked_points in  [9184, 9318, 8862]
  dc.create_spherical_mask()
  dc.mask_info(quiet=True)
  print (dc.mask_info(quiet=True).marked_points)
  assert dc.mask_info(quiet=True).marked_points in [1311, 1286, 1354]
  dc.create_spherical_mask(soft_mask_radius=1)
  dc.mask_info(quiet=True)
  print (dc.mask_info(quiet=True).marked_points)
  assert dc.mask_info(quiet=True).marked_points in [8990,]
  dc.create_spherical_mask(soft_mask=False)
  dc.mask_info(quiet=True)
  print (dc.mask_info(quiet=True).marked_points)
  assert dc.mask_info(quiet=True).marked_points in [1566, 1458, 1812]
  dc.create_spherical_mask(mask_radius = 4)
  dc.mask_info(quiet=True)
  print (dc.mask_info(quiet=True).marked_points)
  assert dc.mask_info(quiet=True).marked_points in [886, 914, 820]
  dc.create_spherical_mask(soft_mask=False, mask_radius = 4)
  dc.mask_info(quiet=True)
  print (dc.mask_info(quiet=True).marked_points)
  assert dc.mask_info(quiet=True).marked_points == 654



if __name__ == "__main__":
  try:
    from phenix.command_line import superpose_pdbs  # special import
    exercise()
  except ImportError:
    print('Skipping tst_map_model_manager_4 tests becuase Phenix is not available.')
  print ("OK")
