from __future__ import absolute_import, division, print_function
import os, sys
import libtbx.load_env
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager
from libtbx.test_utils import approx_equal

def test_01():

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


  # Model sharpening
  mmm = dc.deep_copy()
  tls_info=mmm.tls_from_map(n_bins=10,
    model_id = 'model',
    map_id = 'map_manager',
    iterations = 1,
   )
  tlso = tls_info.tlso_list[0]
  print ("t:", tlso.t)
  print ("l:", tlso.l)
  print ("s:", tlso.s)
  print ("origin:", tlso.origin)
  assert approx_equal(tlso.t,
 (0.7920199173476214, 0.742281514408794, 0.7103008342583756, -0.05199687072329786, 0.04301326317889638, -0.0032498605215769945))
  assert approx_equal(tlso.l,
  (-0.0004103603606922303, -0.00042929108338180655, 0.0001519656028327732, -2.8489076942333132e-06, -6.23198622708519e-05, 3.694504506269563e-05))
  assert approx_equal(tlso.s,
   (5.562000065270741e-09, -5.108813278707348e-10, -4.132731326301999e-09, 7.925572910527853e-10, -2.7018794222798323e-09, 1.0742616538181975e-09, 7.501166703517675e-10, -1.0661760561475498e-09, -2.860120638319051e-09))
  assert approx_equal(tlso.origin,
  (-64.703319312974, -62.30575036040377, -63.74368724016462))

  print("TLS: ",tlso.t,tlso.l,tlso.s,tlso.origin)


# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01()
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")

