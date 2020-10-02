from __future__ import absolute_import, division, print_function
import sys
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal


def get_map_model_managers():
  # Set up source data
  from iotbx.map_model_manager import map_model_manager
  mmm = map_model_manager()
  mmm.generate_map(wrapping=True)
  second_model=mmm.model().deep_copy()
  mmm.box_all_maps_around_model_and_shift_origin()
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
  second_mmm.generate_map(model=second_model,wrapping=True)
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
  print (rt_info)

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
  mmma = mmm1.deep_copy()
  model = mmm1.model()
  mmma.remove_model_by_id('model')
  mmmb = mmma.deep_copy()

  mmma.map_manager().randomize(random_seed=23412)
  mmmb.map_manager().randomize(random_seed=887241)
  assert approx_equal(mmma.map_manager().map_map_cc(mmmb.map_manager()),
    0.40,0.10)
  from iotbx.map_model_manager import map_model_manager
  local_mmm=map_model_manager(map_manager_1=mmma.map_manager(),
    map_manager_2=mmmb.map_manager())
  local_mm = local_mmm.local_fsc()
  cc = local_mm.map_map_cc(local_mmm.map_manager())

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  dc.set_model(model)
  cc_before = dc.map_model_cc()
  dc.model_sharpen(n_bins=15)
  cc = dc.map_model_cc()
  assert approx_equal ((cc_before,cc), (0.93,0.86), eps =0.01)
  model_sharpened_mm = dc.map_manager()

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  dc.set_model(model)
  cc_before = dc.map_model_cc()
  dc.model_sharpen(local_sharpen=True,n_boxes=1,n_bins=15)
  cc = dc.map_model_cc()
  assert approx_equal ((cc_before,cc), (0.93,0.86), eps =0.01)
  model_sharpened_mm = dc.map_manager()

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  dc.add_map_manager_by_id(model_sharpened_mm,'external_map')
  cc_before = dc.map_map_cc(map_id='map_manager',other_map_id='external_map')
  dc.external_sharpen(n_bins=15,map_id_external_map='external_map')
  cc = dc.map_map_cc(map_id='map_manager',other_map_id='external_map')
  assert approx_equal ((cc_before,cc), (0.98,0.99),eps=0.01)

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  dc.add_map_manager_by_id(model_sharpened_mm,'external_map')
  cc_before = dc.map_map_cc(map_id='map_manager',other_map_id='external_map')
  dc.external_sharpen(local_sharpen=True,n_boxes=1,
     n_bins=15,map_id_external_map='external_map')
  cc = dc.map_map_cc(map_id='map_manager',other_map_id='external_map')
  assert approx_equal ((cc_before,cc), (0.98,0.99),eps=0.01)

  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  model.set_b_iso(flex.double(model.get_sites_cart().size(),0))
  dc.set_model(model)
  cc_before = dc.map_model_cc()
  dc.half_map_sharpen(n_bins=15)
  cc_after = dc.map_model_cc()
  assert approx_equal((cc_before,cc_after), (0.9, 0.9), eps=0.10)


  dc = local_mmm.deep_copy()
  dc.set_log(sys.stdout)
  dc.set_model(model)
  dc.local_sharpen(n_bins=15, n_boxes=1)
  cc = dc.map_model_cc()
  assert approx_equal (cc, 0.90, eps=0.1)

  # create a mask around density
  dc.create_mask_around_density()
  count = dc.get_map_manager_by_id('mask').map_data().count(1)
  assert 9000 < count < 10000
  dc.expand_mask(buffer_radius = 2)
  count = dc.get_map_manager_by_id('mask').map_data().count(1)
  assert 20000 < count < 25000

if __name__ == "__main__":
  exercise()
  print ("OK")
