from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import os, sys
from libtbx.utils import Sorry
from libtbx.test_utils import approx_equal
from mmtbx.model import manager as model_manager
from iotbx.data_manager import DataManager


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
  rt_info= mmm1.rt_to_superpose_other(mmm2)

  shift_aware_rt = mmm1.map_manager().shift_aware_rt(working_rt_info=rt_info,
     from_obj = mmm2, to_obj = mmm1)
  print (mmm1, mmm2)
  sites_cart_2 = mmm2.model().get_sites_cart()
  mapped_sites_cart = shift_aware_rt.apply_rt(sites_cart=sites_cart_2,
    from_obj=mmm2, to_obj=mmm1)
  assert approx_equal(mapped_sites_cart,mmm1.model().get_sites_cart(), eps=0.01)
  working_rt_info = shift_aware_rt.working_rt_info(from_obj=mmm2, to_obj=mmm1)
  mapped_sites_cart =working_rt_info.r.elems * mmm2.model().get_sites_cart() + working_rt_info.t.elems
  assert approx_equal(mapped_sites_cart,mmm1.model().get_sites_cart(), eps=0.01)

  inverse_shift_aware_rt = shift_aware_rt.inverse()
  mapped_sites_cart = inverse_shift_aware_rt.apply_rt(sites_cart=mmm1.model().get_sites_cart(),from_obj=mmm1,to_obj=mmm2)
  assert approx_equal(mapped_sites_cart,mmm2.model().get_sites_cart(), eps=0.01)


  return # ZZ


  mmm1, mmm2 = get_map_model_managers()

  # get r,t to map mmm2 model on mmm1 model
  rt_info= mmm1.rt_to_superpose_other(mmm2)
  print (rt_info)

  # get mmm2 map superimposed on mmm1 map (in region where it is defined, zero
  #   outside that region)

  new_mm = mmm1.superposed_map_manager_from_other(other=mmm2)
  new_mm.write_map('super.ccp4')
  mmm1.write_map('orig.ccp4')
  mmm1.write_model('orig.pdb')

if __name__ == "__main__":
  exercise()
  print ("OK")
