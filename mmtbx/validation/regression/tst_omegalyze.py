from __future__ import absolute_import, division, print_function
from mmtbx.validation import omegalyze
from libtbx.test_utils import show_diff, approx_equal
from iotbx import pdb
from iotbx.data_manager import DataManager
from libtbx.test_utils import convert_string_to_cif_long
import libtbx.load_env
import os
import json
import time

ref_omegalyze_give_text = """residues:type:omega:conformation:mc_bmax
 A  40  PHE to  A  41  PRO: Pro     : -14.27:Cis     :28.76
 A 206  VAL to  A 207  LEU: non-Pro : 123.05:Twisted :29.35
 A 504  ILE to  A 505  PRO: Pro     :  -0.31:Cis     :27.55
 A 584  ASN to  A 585  LYS: non-Pro : -12.68:Cis     :26.39
 A 603  ILE to  A 604  GLY: non-Pro :-144.84:Twisted :27.51
 B 929  ARG to  B 930  LEU: non-Pro :-136.14:Twisted :26.97
 B1331  LYS to  B1332  ASP: non-Pro : 134.92:Twisted :27.81
 B1474  GLU to  B1475  LYS: non-Pro : -19.25:Cis     :43.95
 B1593  LYS to  B1594  PRO: Pro     :   8.04:Cis     :26.47
SUMMARY: 3 cis prolines out of 77 PRO
SUMMARY: 0 twisted prolines out of 77 PRO
SUMMARY: 2 other cis residues out of 1464 nonPRO
SUMMARY: 4 other twisted residues out of 1464 nonPRO
"""

class omegalyze_test_string():
  #I wrote the regression test to use a class with a custom .write() method as a
  #  proof of principle for learning OOP and to see if I could. Possible because
  #  all my print functions accept an optional writeto= variable.
  def write(self,string):
    self.output += str(string)
  def __init__(self):
    self.output = ""

def exercise_omegalyze():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2hr0.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_omegalyze(): input pdb (2hr0.pdb) not available")
    return
  #-----
  pdb_io = pdb.input(regression_pdb)
  pdbid = os.path.basename(regression_pdb)
  hierarchy = pdb_io.construct_hierarchy()

  text_test = omegalyze_test_string()
  outliers = omegalyze.omegalyze(
    pdb_hierarchy=hierarchy,
    nontrans_only=True,
    out=text_test,
    quiet=False)
  outliers.show_old_output(out=text_test, verbose=True)

  assert not show_diff(text_test.output , ref_omegalyze_give_text)

def exercise_omegalyze_json(test_mmcif=False):
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2hr0.pdb",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_omegalyze(): input pdb (2hr0.pdb) not available")
    return
  #-----
  dm = DataManager()
  if test_mmcif:
    with open(regression_pdb) as f:
      pdb_2hr0_str = f.read()
    pdb_2hr0_str = convert_string_to_cif_long(pdb_2hr0_str, chain_addition="LONGCHAIN")
    dm.process_model_str("1", pdb_2hr0_str)
    m = dm.get_model("1")
    #print(pdb_2hr0_str)
  else:
    m = dm.get_model(regression_pdb)

  omegalyze_json = omegalyze.omegalyze(pdb_hierarchy=m.get_hierarchy(), nontrans_only=True).as_JSON()
  omjson_dict = json.loads(omegalyze_json)
  #import pprint
  #pprint.pprint(omjson_dict)
  assert len(omjson_dict['flat_results'])==9, "tst_omegalyze json output not returning correct number of nontrans residues, now: "+str(len(omjson_dict['flat_results']))
  assert approx_equal(omjson_dict['flat_results'][0]['omega'], -14.27418253081719), "tst_omegalyze json output first calculated omega dihedral angle not matching previous value, now: "+str(omjson_dict['flat_results'][0]['omega'])
  assert omjson_dict['flat_results'][0]['omega_type']=='Cis', "tst_omegalyze json output first omega_type not matching previous value, now: "+str(omjson_dict['flat_results'][0]['omega_type'])
  assert approx_equal(omjson_dict['flat_results'][0]['highest_mc_b'], 28.76), "tst_omegalyze json output first calculated highest_mc_b not matching previous value, now: "+str(omjson_dict['flat_results'][0]['highest_mc_b'])
  assert approx_equal(omjson_dict['flat_results'][8]['omega'], 8.043663329121266), "tst_omegalyze json output last calculated omega dihedral angle not matching previous value, now: "+str(omjson_dict['flat_results'][8]['omega'])
  assert omjson_dict['flat_results'][8]['omega_type']=='Cis', "tst_omegalyze json output last omega_type not matching previous value, now: "+str(omjson_dict['flat_results'][8]['omega_type'])

  from mmtbx.validation import test_utils
  #assert count(omjson_dict['hierarchical_results'], "PRO")==3, "tst_omegalyze json hierarchical output number of Pro omega outliers changed"
  assert test_utils.count_dict_values(omjson_dict['hierarchical_results'], "Cis")==5, "tst_omegalyze json hierarchical output total number of omega Cis outliers changed to: "+str(test_utils.count_dict_values(omjson_dict['hierarchical_results'], "Cis"))
  assert test_utils.count_dict_values(omjson_dict['hierarchical_results'], "Twisted")==4, "tst_omegalyze json hierarchical output total number of omega Twisted outliers changed to: "+str(test_utils.count_dict_values(omjson_dict['hierarchical_results'], "Twisted"))
  assert omjson_dict['summary_results'][""]['num_cis_proline']==3, "tst_omegalyze json summary results num cis prolines changed to: " + str(omjson_dict['summary_results'][""]['num_cis_proline'])
  return omjson_dict

def run():
  t0 = time.time()
  exercise_omegalyze()
  om_dict = exercise_omegalyze_json()
  om_dict_cif = exercise_omegalyze_json(test_mmcif=True)
  assert om_dict['summary_results'] == om_dict_cif['summary_results'], "tst_omegalyze summary results changed between pdb and cif version"
  print("OK. Time: %8.3f"%(time.time()-t0))

if (__name__ == "__main__"):
  run()
