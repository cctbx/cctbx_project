from __future__ import absolute_import, division, print_function
from mmtbx.wwpdb import rcsb_entry_request
import requests
from libtbx.test_utils import Exception_expected
from libtbx.utils import Sorry

def exercise_1():
  """
  Exercise 1, experimental:
  """
  info = rcsb_entry_request.get_info(pdb_ids=['1ucs'])
  assert len(info) == 1
  assert info[0].get_rwork() == 0.133, info[0].get_rwork()
  assert info[0].get_rama_outliers() == 0
  assert info[0].get_rota_outliers() == 1.82
  assert info[0].get_clashscore() == 19.44
  assert not info[0].is_computational()

def exercise_2():
  """
  Exercise 2, computational:
  """
  info = rcsb_entry_request.get_info(pdb_ids=['AF_AFP12102F1'])
  assert len(info) == 1
  assert info[0].is_computational()
  assert info[0].get_rwork() == None, info[0].get_rwork()

def exercise_3():
  """
  Exercise 3, non-existing:
  """
  try:
    info = rcsb_entry_request.get_info(pdb_ids=['1234567890'])
  except Sorry as e:
    assert str(e).find("There are 1 invalid pdb ids for which RCSB did not return result.") >= 0
  else:
    raise Exception_expected

def exercise_4():
  """
  Exercise 1 + 2 (mix)
  """

  info = rcsb_entry_request.get_info(pdb_ids=['1ucs', 'AF_AFP12102F1'])
  assert len(info) == 2
  assert info[0].get_rwork() == 0.133, info[0].get_rwork()
  assert info[0].get_rama_outliers() == 0
  assert info[0].get_rota_outliers() == 1.82
  assert info[0].get_clashscore() == 19.44
  assert not info[0].is_computational()

  assert info[1].is_computational()
  assert info[1].get_plddt() == 97.37, info[1].get_plddt()
  assert info[1].get_rwork() == None, info[1].get_rwork()

if (__name__ == "__main__"):
  # thorough_exercise()
  # check if internet and rcsb are available
  exception_occured = False
  try:
    r = requests.get('https://search.rcsb.org/')
  except Exception:
    print("OK but exception.")
    exception_occured = True
  if not exception_occured and r.ok and len(r.text) > 100:
    exercise_1()
    exercise_2()
    exercise_3()
    exercise_4()
    print("OK")
  else:
    print("OK but skipped.")
