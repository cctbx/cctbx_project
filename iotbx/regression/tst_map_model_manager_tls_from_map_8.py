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
(1.180418902258779, 1.1747521845606608, 1.178996799712174, -0.08474662674769494, -0.022609295693646402, 0.0649209491344932)
 ,eps=.001)
  assert approx_equal(tlso.l,
(-0.002159404807991249, -0.002107964765763024, 0.0008301439376854558, -5.973347993775719e-05, -0.000134276871934738, -9.05515898670584e-05)
 ,eps=.001)
  assert approx_equal(tlso.s,
(2.9348223335616302e-08, 5.52441087256425e-09, -5.382459681103171e-09, 4.3530347434547015e-09, -2.3559464233595e-08, 4.217968590464982e-09, -4.380707049750269e-09, 1.9232725033868253e-09, -5.788759082799497e-09)
 ,eps=.001)
  assert approx_equal(tlso.origin,
(-64.70331931297399, -62.30573551948903, -63.743687240164604)
 ,eps=.001)


  print("TLS: ",tlso.t,tlso.l,tlso.s,tlso.origin)


# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01()
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")

