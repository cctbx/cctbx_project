from __future__ import division
import time
import os

import remediator
import iotbx.pdb
import libtbx

def test_pdb_to_v3(old_pdb_file, v3_pdb_file):
  print("............checking conversion of "+old_pdb_file+" to new PDB format")
  remed = remediator.Remediator()
  old_pdb_model = remediator.get_model_from_file(old_pdb_file)
  #write_model_to_file(old_pdb_model.model_as_pdb(), "old_file.pdb")
  old_pdb_model_v3 = remed.remediate_model_object(old_pdb_model, True)
  v3_pdb_model = remediator.get_model_from_file(v3_pdb_file)
  #write_model_to_file(old_pdb_model_v3.model_as_pdb(), "old_file_remed.pdb")
  #write_model_to_file(v3_pdb_model.model_as_pdb(), "v3_file.pdb")
  assert old_pdb_model_v3.get_hierarchy().as_pdb_string() == v3_pdb_model.get_hierarchy().as_pdb_string(), \
    "remediation test of "+old_pdb_file+" to v3 failed"
  print("OK")
  #print(old_pdb_model_v3.get_hierarchy().is_similar_hierarchy(v3_pdb_model.get_hierarchy()))
  #print(old_pdb_model.get_hierarchy().is_similar_hierarchy(v3_pdb_model.get_hierarchy()))

def test_pdb_to_v23(pdb_file):
  remediator.remediate(pdb_file, False)

def test_pdb_remediate_cycle(pdb_file):
  t0 = time.time()
  print("testing remediation cycle of "+pdb_file)
  remediatorObj = remediator.Remediator()
  pdb_model = remediator.get_model_from_file(pdb_file)
  pdb_model = remediatorObj.remediate_model_object(pdb_model, True)
  #write_model_to_file(pdb_model.model_as_pdb(), "input_model.pdb")
  print("............checking if remediated "+pdb_file+" is v3")
  assert remediatorObj.is_model_v3(pdb_model), "remediation of "+pdb_file+ " during remediation cycle test failed"
  pdb_model = remediatorObj.remediate_model_object(pdb_model, False)
  #write_model_to_file(pdb_model.model_as_pdb(), "old_model.pdb")
  pdb_model = remediatorObj.remediate_model_object(pdb_model, True)
  #write_model_to_file(pdb_model.model_as_pdb(), "reold_model.pdb")
  print("............checking remediated-to-v2.3 and remediated "+pdb_file+" is v3")
  assert remediatorObj.is_model_v3(pdb_model), "cycling of "+pdb_file+" from v3 to v2.3 back to v3 failed"
  print("OK. "+pdb_file+" time: %8.3f"%(time.time()-t0))

def write_model_to_file(model, file_name):
  with open("tst/"+file_name, 'w+') as fh:
    fh.write(model)

def get_regression_folder_file_path(file_name):
  regression_folder = os.path.join(libtbx.env.dist_path("iotbx"), 'pdb', 'remediation', 'tst')
  file_path = os.path.join(regression_folder, file_name)
  assert os.path.isfile(file_path), file_name + " test file is missing from " + regression_folder
  return file_path

def run_tests():
  test_pdb_to_v3(get_regression_folder_file_path("1ubqFHv23.pdb"), get_regression_folder_file_path("1ubqFH.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("1ubqFH.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("1ubqFHv23.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("404dH.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("mg_tst.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("2f55FH.pdb"))
  #test_pdb_remediate_cycle("tst/4v9dFH.pdb")

if __name__=="__main__":
  run_tests()
  print("OK")
