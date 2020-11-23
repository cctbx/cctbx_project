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
 (0.6518723599417712, 0.6807846368236251, 0.6161941485135081, -0.04791178588048965, -0.0014794039180157132, 0.032774655367019095) )
  assert approx_equal(tlso.l,
  (-0.00020092054127279798, -9.557441568256003e-05, -1.711699526358822e-05, -4.8893794663274e-05, -2.026091444762368e-05, -2.194054393244713e-05))
  assert approx_equal(tlso.s,
   (1.3393493443427932e-09, 3.660975429526433e-10, -6.07051859483934e-10, 4.1665859321329886e-10, -1.35856164931005e-09, 5.409502867516901e-10, -5.4005830782848e-10, -1.2017586805521653e-09, 1.9212302485258283e-11))
  assert approx_equal(tlso.origin,
   (-64.70331931297407, -62.305747062422725, -63.74368724016457))

  print("TLS: ",tlso.t,tlso.l,tlso.s,tlso.origin)


# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01()
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")

