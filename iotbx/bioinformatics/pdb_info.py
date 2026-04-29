"""Get experimental information summary for PDB IDS from the PDB"""

from __future__ import absolute_import, division, print_function

import sys
from libtbx import easy_pickle
from libtbx.utils import Sorry
from libtbx import smart_open
import csv
import requests

web_urls = {"rcsb": (
  "https://data.rcsb.org/graphql", """
  {{
    entries(entry_ids: {pdb_list} )
    {{
      rcsb_id
      refine {{
        ls_R_factor_R_free
        ls_R_factor_R_work
        ls_R_factor_obs
        ls_d_res_high
      }}
    }}
  }}"""
)}

def get_experimental_pdb_info(pdbids, site="rcsb"):
  """
  returns list of tuples (pdb_id, resolution, rwork, rfree) and dict
  pdbid: (resolution, rwork, rfree)


  OBSOLETED. New functionality is mmtbx.wwpdb.rcsb_entry_request.get_info
  """
  rlist = []
  rdict = {}
  assert site in ["rcsb"]
  url = web_urls[site][0]
  request = web_urls[site][1]

  pdb_list = "%s" % pdbids
  pdb_list = pdb_list.replace("'", '"')
  request = request.format(pdb_list=pdb_list)
  r = requests.post(url, json={"query":request}, timeout=10)
  r_json = r.json()
  for res in r_json["data"]["entries"]:
    pdb_id = str(res["rcsb_id"])
    resolution, rwork, rfree = (None, None, None)
    if res["refine"] is not None:
      resolution = None if res["refine"][0]["ls_d_res_high"] is None else float(res["refine"][0]["ls_d_res_high"])
      rwork = None if res["refine"][0]["ls_R_factor_R_work"] is None else float(res["refine"][0]["ls_R_factor_R_work"])
      rfree = None if res["refine"][0]["ls_R_factor_R_free"] is None else float(res["refine"][0]["ls_R_factor_R_free"])
      if rwork is None:
        rwork = None if res["refine"][0]["ls_R_factor_obs"] is None else float(res["refine"][0]["ls_R_factor_obs"])
    tup = (pdb_id, resolution, rwork, rfree)
    rlist.append(tup)
    rdict[pdb_id] = tup[1:]
  return rlist, rdict

class pdb_info_local(object):
  def __init__(self):
    """
    Loads pickle with data. Path is temporary in current work dir.
    Should be centralized somewhere else upon going to production.
    """
    db_dict = {}

    import iotbx
    from pathlib import Path
    data_dir = Path(iotbx.__file__).parent / 'bioinformatics'
    pdb_info_file = str( data_dir / 'pdb_info.csv.gz')


    csv_file = smart_open.for_reading(file_name=pdb_info_file, gzip_mode="rt")
    csv_reader = csv.reader(csv_file,delimiter=";")
    for row in csv_reader:
      db_dict[row[0]] = (row[1],row[2],row[3],row[4],row[5])
    self.db_dict = db_dict

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
  pdb_dict 5.6 Mb.
  Takes ~1 minute from home.
  Use only xray diffraction.
  """

  base_url = "https://search.rcsb.org/rcsbsearch/v1/query?json="
  q = {
  "query": {
    "type": "terminal",
    "service": "text",
    "parameters": {
      "attribute": "exptl.method",
      "operator": "exact_match",
      "value": "X-RAY DIFFRACTION"
    }
  },
  "request_options": {
    "return_all_hits": True,
  },
  "return_type": "entry"
  }

  # First get all x-ray pdb ids
  r1 = requests.post(base_url, json=q)
  r1_json = r1.json()
  print ('Total:', r1_json["total_count"])
  res_ids = []
  for res in r1_json["result_set"]:
    res_ids.append(str(res["identifier"]))
  print ('total resids', len(res_ids))
  # Now get the info:
  rlist, rdict = get_experimental_pdb_info(res_ids)
  n_bad = 0
  for tup in rlist:
    if tup.count(None) > 0:
      print(tup)
      n_bad += 1
  print("Total bad records", n_bad)
  easy_pickle.dump(file_name='pdb_dict.pickle', obj=rdict)

def tst_pdb_info_local():
  # Enable before running.
  # get_all_experimental_pdb_info_to_pkl()

  # I don't know why there are 5 values now in the table for each PDB.
  # info_local = pdb_info_local()
  # ans_dict_1 = {'1yjp': (1.8, 0.181, 0.19), '1ucs': (0.62, 0.133, 0.155)}
  # ans_list_1 = [('1ucs', 0.62, 0.133, 0.155), ('1yjp', 1.8, 0.181, 0.19)]
  # assert info_local.get_info_dict(["1ucs", "1yjp"]) == ans_dict_1,
  # assert info_local.get_info_list(["1ucs", "1yjp"]) == ans_list_1
  ans_dict_2 = {'1YJP': (1.8, 0.18086, 0.19014), '1UCS': (0.62, 0.133, 0.155)}
  ans_list_2 = [('1UCS', 0.62, 0.133, 0.155), ('1YJP', 1.8, 0.18086, 0.19014)]
  try:
    rlist, rdict = get_experimental_pdb_info(["1ucs", "1yjp"])
  except requests.exceptions.ReadTimeout:
    print("Skipped test: transient read timeout, can't run test right now")
    return
  assert rlist == ans_list_2, rlist
  assert rdict == ans_dict_2, rdict


def run(args):
  tst_pdb_info_local()

if __name__ == '__main__':
  run(sys.argv[1:])
