from __future__ import absolute_import, division, print_function
import libtbx.load_env
from iotbx.regression.tst_map_model_manager_model_sharpening_5 import test_01
# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01(method = 'model_sharpen')
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")

