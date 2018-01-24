from __future__ import division, print_function

import os

import libtbx.phil

from libtbx.utils import Sorry
from libtbx.data_manager import DataManager

def test_data_manager():
  a = DataManager()

  a.add_model('a', 'b')
  a.add_model('c', 'd')
  assert (a.get_model() == 'b')
  assert (a.get_model('a') == 'b')
  assert (a.get_model('c') == 'd')

  assert (a.has_models())
  assert (a.has_models(exact_count=True, expected_n=2))
  assert (not a.has_models(expected_n=3, raise_sorry=False))

  # test exporting phil
  working_phil = a.export_phil_scope()
  assert(working_phil.extract().data_manager.model_files == ['a', 'c'])

  try:
    a.has_models(expected_n=3, raise_sorry=True)
  except Sorry:
    pass

  try:
    a.has_models(exact_count=True, raise_sorry=True)
  except Sorry:
    pass

  a.set_default_model('c')
  assert (a.get_model() == 'd')

  assert ( (a.get_model_names() == ['a', 'c']) or
           (a.get_model_names() == ['c', 'a']) )

  a.remove_model('c')
  try:
    a.get_model()
  except Sorry:
    pass
  try:
    a.get_model('missing')
  except Sorry:
    pass
  try:
    a.set_default_model('missing')
  except Sorry:
    pass

  assert (a.get_sequence_names() == [])

  # test loading from phil
  test_phil_str = '''
data_manager {
  phil_files = data_manager_test.eff
}
'''
  with open('data_manager_test.eff', 'w') as f:
    f.write(test_phil_str)

  test_phil = libtbx.phil.parse(test_phil_str)
  a.load_phil_scope(test_phil)

  assert('data_manager_test.eff' in a.get_phil_names())

  os.remove('data_manager_test.eff')

if (__name__ == '__main__'):

  test_data_manager()

  print('OK')
