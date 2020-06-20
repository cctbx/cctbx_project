from __future__ import absolute_import, division, print_function
import os
import libtbx.load_env
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import match_map_model_ncs, map_model_manager 
from libtbx.test_utils import approx_equal

def test_01():

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_ccp4 = os.path.join(data_dir, 'data',
                          'non_zero_origin_map.ccp4')
  data_pdb = os.path.join(data_dir, 'data',
                          'non_zero_origin_model.pdb')
  data_ncs_spec = os.path.join(data_dir, 'data',
                          'non_zero_origin_ncs_spec.ncs_spec')

  dm = DataManager(['ncs_spec','model', 'real_map', 'phil'])
  dm.set_overwrite(True)

  # Read in map and model

  map_file=data_ccp4
  dm.process_real_map_file(map_file)
  mm = dm.get_real_map(map_file)

  model_file=data_pdb
  dm.process_model_file(model_file)
  model = dm.get_model(model_file)

  ncs_file=data_ncs_spec
  dm.process_ncs_spec_file(ncs_file)
  ncs = dm.get_ncs_spec(ncs_file)

  mmmn = match_map_model_ncs()
  mmmn.add_map_manager(mm)
  mmmn.add_model(model)
  mmmn.add_ncs_object(ncs)

  original_ncs=mmmn.ncs_object()
  assert approx_equal((24.0528, 11.5833, 20.0004),
     tuple(original_ncs.ncs_groups()[0].translations_orth()[-1]),
     eps=0.1)

  assert tuple(mmmn._map_manager.origin_shift_grid_units) == (0,0,0)

  # Shift origin to (0,0,0)
  mmmn.shift_origin()
  assert tuple(mmmn._map_manager.origin_shift_grid_units) == (100,100,100)

  mmmn.write_model('s.pdb')
  mmmn.write_map('s.mrc')

  shifted_ncs=mmmn.ncs_object()
  assert approx_equal((-153.758, -74.044, -127.487),
      tuple(shifted_ncs.ncs_groups()[0].translations_orth()[-1]),eps=0.1)


  # Shift a model and shift it back

  model=mmmn.model()
  shifted_model=mmmn.shift_model_to_match_working_map(model=model)
  model_in_original_position=mmmn.shift_model_to_match_original_map(
      model=shifted_model)
  assert (approx_equal(model.get_sites_cart(), # not a copy
                      shifted_model.get_sites_cart()))
  assert approx_equal(model.get_sites_cart(),
                      model_in_original_position.get_sites_cart())


  # Generate a map and model

  mmm=map_model_manager()
  mmm.generate_map()
  model=mmm.model()
  mm=mmm.map_manager()
  assert approx_equal(
     model.get_sites_cart()[0], (14.476, 10.57, 8.34) ,eps=0.01)
  assert approx_equal(mm.map_data()[10,10,10],-0.0195,eps=0.001)

# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01()
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")
