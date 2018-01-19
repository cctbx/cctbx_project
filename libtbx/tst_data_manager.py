from __future__ import division, print_function

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
  assert (not a.has_models(expected_n=3))

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

if (__name__ == '__main__'):

  test_data_manager()

  print('OK')
