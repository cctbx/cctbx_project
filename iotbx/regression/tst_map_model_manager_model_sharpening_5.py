from __future__ import absolute_import, division, print_function
import os, sys
import libtbx.load_env
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager

def test_01(method = 'model_sharpen',
   expected_results=None):

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
  if method == 'local_resolution_map':
    mmm.generate_map(map_id = 'map_manager_2')
    dd = mmm.resolution()
  else:
    mmm.set_resolution(3)
  mmm.set_log(sys.stdout)

  dc = mmm.deep_copy()

  sharpen_method = getattr(mmm,method)

  # sharpen by method (can be model_sharpen, half_map_sharpen or
  #     external_sharpen and also local_resolution_map)

  if method == 'local_resolution_map':
    mm = sharpen_method()
    x = mm.map_data().as_1d().min_max_mean()
    from libtbx.test_utils import approx_equal
    assert approx_equal((x.min,x.max,x.mean),
      (1.9835316316768663, 2.3146054110530336, 2.1600110833032353),
       eps = 0.01)
    mm = mmm.local_resolution_map(map_id_1 = None, map_id_2 = None,
     model_id = 'model', map_id = 'map_manager_1')
    xx = mm.map_data().as_1d().min_max_mean()
    from libtbx.test_utils import approx_equal
    assert approx_equal((x.mean, xx.mean),
      (2.1600110833032353, 2.16),
       eps = 0.01)

  else: # usual
    sharpen_method(anisotropic_sharpen = False, n_bins=10)
    assert mmm.map_model_cc() > 0.9
    sharpen_method(anisotropic_sharpen = False, n_bins=10,
       local_sharpen = True)
    assert mmm.map_model_cc() > 0.9
    sharpen_method(anisotropic_sharpen = True, n_bins=10)
    assert mmm.map_model_cc() > 0.9
    sharpen_method(anisotropic_sharpen = True, n_bins=10,
       local_sharpen = True, n_boxes = 1)
    assert mmm.map_model_cc() > 0.9


def test_02():
  """Summary of anisotropic scaling when anisotropy is not available"""

  from libtbx.test_utils import approx_equal

  # Cheap map_model_manager: no data files, no chem_data
  mmm = map_model_manager()
  mmm.set_log(None)
  mmm.generate_map(d_min = 3)
  mmm.add_map_manager_by_id(map_id = 'previous_map',
     map_manager = mmm.map_manager().deep_copy())

  # Supply the anisotropy of each map ourselves, keyed on map_id so that
  #  the two calls made by _get_aniso_before_and_after are distinguishable
  aniso_by_map_id = {}
  def get_aniso_of_map(d_min = None, map_id = None):
    return aniso_by_map_id[map_id]
  mmm._get_aniso_of_map = get_aniso_of_map

  prev_b_cart = (10., 20., 30., 1., 2., 3.)
  new_b_cart = (4., 5., 6., 0., 1., 2.)

  # Both available: b_sharpen is previous minus new, element by element
  aniso_by_map_id['previous_map'] = prev_b_cart
  aniso_by_map_id['map_manager'] = new_b_cart
  info = mmm._get_aniso_before_and_after(d_min = 3,
     map_id = 'map_manager', previous_map_id = 'previous_map')
  assert approx_equal(tuple(info.b_sharpen), (6., 15., 24., 1., 1., 1.))
  assert info.text.find('B-cart') > -1
  assert info.text.find('B-sharpen') > -1
  assert info.text.find('not available') < 0

  # Not available for one map or for both: return normally, b_sharpen is
  #  None and the text says the information is not available
  for prev_value, new_value in [
      (None, new_b_cart),
      (prev_b_cart, None),
      (None, None)]:
    aniso_by_map_id['previous_map'] = prev_value
    aniso_by_map_id['map_manager'] = new_value
    info = mmm._get_aniso_before_and_after(d_min = 3,
       map_id = 'map_manager', previous_map_id = 'previous_map')
    assert info.b_sharpen is None
    assert info.text.find('not available') > -1


# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  test_02()
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01(method = 'model_sharpen')
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")

