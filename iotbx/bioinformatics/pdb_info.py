from __future__ import division
from __future__ import print_function

import sys
import libtbx.utils
import xml.etree.ElementTree as ET
from libtbx import easy_pickle
from libtbx.utils import Sorry

web_urls = {"rcsb": "http://www.rcsb.org/pdb/rest/customReport.xml?"+\
    "pdbids={pdb_list}"+\
    "&customReportColumns=structureId,refinementResolution,rWork,rFree"
}


def get_experimental_pdb_info(pdbids, site="rcsb"):
  """
  returns list of tuples (pdb_id, resolution, rwork, rfree) and dict
  pdbid: (resolution, rwork, rfree)
  """
  rlist = []
  rdict = {}
  assert site in ["rcsb"]
  pdb_list = ",".join(pdbids)
  url = web_urls.get(site).format(pdb_list=pdb_list)
  data = libtbx.utils.urlopen(url)
  str_data = data.read()
  root = ET.fromstring(str_data)
  for record in root:
    pdbid = record[0].text
    resolution = None if record[1].text == 'null' else float(record[1].text)
    rwork = None if record[2].text == 'null' else float(record[2].text)
    rfree = None if record[3].text == 'null' else float(record[3].text)
    tup = (pdbid, resolution, rwork, rfree)
    rlist.append(tup)
    rdict[record[0].text] = tup[1:]
  return rlist, rdict

class pdb_info_local(object):
  def __init__(self):
    """
    Loads pickle with data. Path is temporary in current work dir.
    Should be centralized somewhere else upon going to production.
    """
    self.db_dict = easy_pickle.load("pdb_dict.pickle")
    # print self.db_dict

  def _get_info(self, pdbid, skip_none=True, raise_if_absent=False):
    info = self.db_dict.get(pdbid.upper(), None)
    if info is None and raise_if_absent:
      raise Sorry("Not in database")
    if skip_none and info is not None and info[0] is None:
      info = None
    return info

  def get_info_list(self, pdbids, skip_none=True, raise_if_absent=False):
    """
    Get info about pdbids (list of strings) in form of list of tuples
    (pdbid, resolution, rwork, rfree). Easy to sort.
    """
    result = []
    for pdbid in pdbids:
      info = self._get_info(pdbid, raise_if_absent=raise_if_absent)
      if info is not None:
        result.append( (pdbid,) + info)
    return result

  def get_info_dict(self, pdbids, skip_none=True, raise_if_absent=False):
    """
    Get info about pdbids (list of strings) in form of dict
    pdbid: (resolution, rwork, rfree). Easy to lookup.
    """
    result = {}
    for pdbid in pdbids:
      info = self._get_info(pdbid, raise_if_absent=raise_if_absent)
      if info is not None:
        result[pdbid] = info
    return result

def get_all_experimental_pdb_info_to_pkl():
  """
  Get info (resolution, rwork, rfree) for all PDB from RCSB and dump into
  pickle file:
  pdb_dict 5.1 Mb.
  Takes ~15 seconds from LBL.
  Use only xray diffraction.
  """
  get_all = "http://www.rcsb.org/pdb/rest/customReport.xml?"+\
      "pdbids=*&customReportColumns=" +\
      "structureId,refinementResolution,rWork,rFree,experimentalTechnique,rObserved"
  print(get_all)
  rdict = {}
  data = libtbx.utils.urlopen(get_all)
  str_data = data.read()
  # print str_data
  root = ET.fromstring(str_data)
  n_bad = 0
  for record in root:
    if record[4].text != "X-RAY DIFFRACTION":
      continue
    pdbid = record[0].text
    resolution = None if record[1].text == 'null' else float(record[1].text)
    rwork = None if record[2].text == 'null' else float(record[2].text)
    rfree = None if record[3].text == 'null' else float(record[3].text)
    if rwork is None:
      # put rObserved
      rwork = None if record[5].text == 'null' else float(record[5].text)
    tup = (pdbid, resolution, rwork, rfree)
    rdict[record[0].text] = tup[1:]
    # print tup
    if tup.count(None) > 0:
      print(tup)
      n_bad += 1
  print("Total bad records", n_bad)
  easy_pickle.dump(file_name='pdb_dict.pickle', obj=rdict)

def tst_pdb_info_local():
  # Enable before running.
  # get_all_experimental_pdb_info_to_pkl()
  info_local = pdb_info_local()
  ans_dict_1 = {'1yjp': (1.8, 0.181, 0.19), '1ucs': (0.62, 0.133, 0.155)}
  ans_list_1 = [('1ucs', 0.62, 0.133, 0.155), ('1yjp', 1.8, 0.181, 0.19)]
  assert info_local.get_info_dict(["1ucs", "1yjp"]) == ans_dict_1
  assert info_local.get_info_list(["1ucs", "1yjp"]) == ans_list_1
  ans_dict_2 = {'1YJP': (1.8, 0.181, 0.19), '1UCS': (0.62, 0.133, 0.155)}
  ans_list_2 = [('1UCS', 0.62, 0.133, 0.155), ('1YJP', 1.8, 0.181, 0.19)]
  rlist, rdict = get_experimental_pdb_info(["1ucs", "1yjp"])
  assert rlist == ans_list_2
  assert rdict == ans_dict_2


def run(args):
  tst_pdb_info_local()

if __name__ == '__main__':
  run(sys.argv[1:])
