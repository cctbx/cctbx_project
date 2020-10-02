from __future__ import absolute_import, division, print_function
from iotbx.map_model_manager import map_model_manager
from cctbx.development.create_models_or_maps import *
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import random

def tst_01():

  print (" test utilities in create_models_or_maps")
  flex.set_random_seed(0)
  random.seed(0)

  model=generate_model()
  assert approx_equal (
    model.get_sites_cart()[0],
    (14.476, 10.57, 8.342))

  s_model=shake_model(model,shake=2)
  assert approx_equal (
    s_model.get_sites_cart()[0],
    (14.162085804614943, 11.403509966523153, 6.450881839681677))
  map_coeffs=generate_map_coefficients(model=model)
  assert approx_equal(map_coeffs.data()[0],
     (3.70494534745-0.185333495539j))
  map_manager=generate_map(map_coeffs=map_coeffs)

  mm_2= generate_map(map_coeffs=map_coeffs,
      d_min=3.5,
      gridding=map_manager.map_data().all(),
      wrapping=True,
      origin_shift_grid_units=(100,0,0),
      low_resolution_fourier_noise_fraction=1,
      high_resolution_fourier_noise_fraction=1,
      low_resolution_real_space_noise_fraction=1,
      high_resolution_real_space_noise_fraction=1,
      )
  assert approx_equal (mm_2.map_data()[323], -0.0784650534896)
  mm_2.shift_origin()
  model.set_shift_cart(mm_2.shift_cart())
  mam=map_model_manager(map_manager=mm_2,model=model)
  mam.write_map('map.mrc')
  mam.write_model('model.pdb')

  new_mam=read_map_and_model('map.mrc','model.pdb')
  new_mam_2=read_map_and_model('model.pdb','map.mrc')
  assert new_mam.map_manager().cc_to_other_map_manager(new_mam_2.map_manager())==1




if (__name__ == "__main__"):
  import libtbx.load_env
  if libtbx.env.find_in_repositories(relative_path="chem_data") is None:
    print("Skipping exercise(): chem_data directory not available")
  else:
    tst_01()
  print ("OK")
