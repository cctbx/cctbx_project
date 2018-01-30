from __future__ import division, print_function

import os

import libtbx.phil

from libtbx.utils import Sorry
from libtbx.data_manager import DataManager

def test_data_manager():
  a = DataManager(datatypes=['model'])

  a.add_model('a', 'b')
  a.add_model('c', 'd')
  assert (a.get_model() == 'b')
  assert (a.get_model('a') == 'b')
  assert (a.get_model('c') == 'd')

  assert (a.has_models())
  assert (a.has_models(exact_count=True, expected_n=2))
  assert (not a.has_models(expected_n=3, raise_sorry=False))

  # exporting phil
  working_phil = a.export_phil_scope()
  assert(working_phil.extract().data_manager.model_files == ['a', 'c'])

  # data tracking
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

  a = DataManager(datatypes=['sequence', 'phil'])
  assert (a.get_sequence_names() == [])
  assert (not hasattr(a, 'get_model'))

  # phil functions
  test_phil_str = '''
data_manager {
  phil_files = data_manager_test.eff
}
'''
  with open('data_manager_test.eff', 'w') as f:
    f.write(test_phil_str)

  # loading file with get function
  assert(len(a.get_phil_names()) == 0)
  p = a.get_phil('data_manager_test.eff')
  assert(type(p) == libtbx.phil.scope)
  assert('data_manager_test.eff' in a.get_phil_names())

  # loading file with phil
  a = DataManager(datatypes=['phil'])
  test_phil = libtbx.phil.parse(test_phil_str)
  a.load_phil_scope(test_phil)

  assert('data_manager_test.eff' in a.get_phil_names())
  assert(a.get_default_phil_name() == 'data_manager_test.eff')

  os.remove('data_manager_test.eff')

  # writing
  a = DataManager(datatypes=['model', 'phil', 'sequence'])
  a.add_model('a','b')
  a.add_phil('c','d')
  a.add_sequence('e','f')

  a.write_model_file('a.dat', a.get_model(), overwrite=True)
  a.write_phil_file('c.dat', a.get_phil(), overwrite=True)
  a.write_sequence_file('e.dat', a.get_sequence(), overwrite=True)

  with open('a.dat', 'r') as f:
    lines = f.readlines()
  assert(lines[0] == 'b')

  os.remove('a.dat')
  os.remove('c.dat')
  os.remove('e.dat')

if (__name__ == '__main__'):

  test_data_manager()

  print('OK')
