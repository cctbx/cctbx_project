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
  # print(info[0].data)
  # STOP()
  assert len(info) == 1
  assert info[0].get_rwork() == 0.133, info[0].get_rwork()
  assert info[0].get_rama_outliers() == 0
  assert info[0].get_rota_outliers() == 0, info[0].get_rota_outliers()
  assert info[0].get_clashscore() == 9.72, info[0].get_clashscore()
  assert not info[0].is_computational()

def exercise_11():
  """
  Exercise 1, experimental, cryoem, actually with outliers:
  """
  info = rcsb_entry_request.get_info(pdb_ids=['9SKV'])
  # print(info[0].data)
  # STOP()
  assert len(info) == 1
  assert info[0].get_rwork() == None, info[0].get_rwork()
  assert info[0].get_rama_outliers() == 0.74, info[0].get_rama_outliers()
  assert info[0].get_rota_outliers() == 3.53, info[0].get_rota_outliers()
  assert info[0].get_clashscore() == 4.87, info[0].get_clashscore()
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
  assert info[0].get_rota_outliers() == 0
  assert info[0].get_clashscore() == 9.72
  assert not info[0].is_computational()

  assert info[1].is_computational()
  assert info[1].get_plddt() == 97.49, info[1].get_plddt()
  assert info[1].get_rwork() == None, info[1].get_rwork()

def _info(pdb_id, methods):
  """Build rcsb_entry_info offline from a minimal RCSB-like dict."""
  data = {'rcsb_id': pdb_id, 'rcsb_ma_qa_metric_global': None}
  if methods is not None:
    data['exptl'] = [{'method': m} for m in methods]
  return rcsb_entry_request.rcsb_entry_info(data)

def exercise_5():
  """
  Exercise 5, offline: is_xray() is True only for pure X-ray entries.
  Joint X-ray/neutron (either listing order), pure neutron, MicroED, EM and
  entries without exptl are not x-ray.
  """
  assert _info('1UCS', ['X-RAY DIFFRACTION']).is_xray()
  assert not _info('5MO0', ['NEUTRON DIFFRACTION']).is_xray()
  assert not _info('3OTJ', ['NEUTRON DIFFRACTION', 'X-RAY DIFFRACTION']).is_xray()
  assert not _info('XXXX', ['X-RAY DIFFRACTION', 'NEUTRON DIFFRACTION']).is_xray()
  assert not _info('5K7R', ['ELECTRON CRYSTALLOGRAPHY']).is_xray()
  assert not _info('9GSD', ['ELECTRON MICROSCOPY']).is_xray()
  assert not _info('NONE', None).is_xray()
  assert not _info('EMPT', []).is_xray()
  # all methods are available, in RCSB order
  m = _info('3OTJ', ['NEUTRON DIFFRACTION', 'X-RAY DIFFRACTION']).get_experimental_methods()
  assert m == ['NEUTRON DIFFRACTION', 'X-RAY DIFFRACTION'], m
  assert _info('NONE', None).get_experimental_methods() == []

def exercise_12():
  """
  Exercise 12, live: joint X-ray/neutron and MicroED entries are not x-ray.
  """
  info = rcsb_entry_request.get_info(pdb_ids=['3OTJ', '5K7R', '1ucs'])
  assert len(info) == 3
  m = info[0].get_experimental_methods()
  assert m == ['NEUTRON DIFFRACTION', 'X-RAY DIFFRACTION'], m
  assert not info[0].is_xray()
  assert info[1].get_experimental_methods() == ['ELECTRON CRYSTALLOGRAPHY'], \
      info[1].get_experimental_methods()
  assert not info[1].is_xray()
  assert info[2].is_xray()

if (__name__ == "__main__"):
  # thorough_exercise()
  exercise_5()
  # check if internet and rcsb are available
  exception_occured = False
  try:
    r = requests.get('https://search.rcsb.org/')
  except Exception:
    print("OK but exception.")
    exception_occured = True
  if not exception_occured and r.ok and len(r.text) > 100:
    exercise_1()
    exercise_11()
    exercise_12()
    exercise_2()
    exercise_3()
    exercise_4()
    print("OK")
  else:
    print("OK but skipped.")
