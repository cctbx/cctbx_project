from __future__ import absolute_import, division, print_function
import os, sys
import libtbx.load_env
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager

def test_01(method = 'model_sharpen'):

  # Source data

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_ccp4 = os.path.join(data_dir, 'data',
                          'non_zero_origin_map.ccp4')
  data_pdb = os.path.join(data_dir, 'data',
                          'non_zero_origin_model.pdb')
  data_ncs_spec = os.path.join(data_dir, 'data',
                          'non_zero_origin_ncs_spec.ncs_spec')

  # Read in data

  dm = DataManager(['ncs_spec','model', 'real_map', 'phil'])
  dm.set_overwrite(True)

  map_file=data_ccp4
  dm.process_real_map_file(map_file)
  mm = dm.get_real_map(map_file)

  model_file=data_pdb
  dm.process_model_file(model_file)
  model = dm.get_model(model_file)

  ncs_file=data_ncs_spec
  dm.process_ncs_spec_file(ncs_file)
  ncs = dm.get_ncs_spec(ncs_file)

  mmm=map_model_manager(
    model = model,
    map_manager_1 = mm.deep_copy(),
    map_manager_2 = mm.deep_copy(),
    ncs_object = ncs,
    wrapping = False)
  mmm.add_map_manager_by_id(
     map_id='external_map',map_manager=mmm.map_manager().deep_copy())
  mmm.set_resolution(3)
  mmm.set_log(sys.stdout)

  dc = mmm.deep_copy()

  sharpening_method = getattr(mmm,method)

  # Model sharpening
  mmm = dc.deep_copy()
  sharpening_method(anisotropic_sharpen = False, n_bins=10)
  sharpening_method(anisotropic_sharpen = False, n_bins=10,
     local_sharpen = True)
  sharpening_method(anisotropic_sharpen = True, n_bins=10)
  sharpening_method(anisotropic_sharpen = True, n_bins=10,
     local_sharpen = True, n_boxes = 1)


# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01(method = 'model_sharpening')
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")

