import time

import remediator
import iotbx.pdb
import mmtbx

def test_pdb_to_v3(old_pdb_file, v3_pdb_file):
  remed = remediator.Remediator()
  old_pdb_model = remediator.get_model_from_file(old_pdb_file)
  old_pdb_model_v3 = remed.remediate_model_object(old_pdb_model, True)
  v3_pdb_model = remediator.get_model_from_file(v3_pdb_file)
  write_model_to_file(old_pdb_model_v3.model_as_pdb(), "old_file.pdb")
  write_model_to_file(v3_pdb_model.model_as_pdb(), "v3_file.pdb")
  print(old_pdb_model_v3.get_hierarchy().as_pdb_string() == v3_pdb_model.get_hierarchy().as_pdb_string())
  print(old_pdb_model_v3.get_hierarchy().is_similar_hierarchy(v3_pdb_model.get_hierarchy()))
  print(old_pdb_model.get_hierarchy().is_similar_hierarchy(v3_pdb_model.get_hierarchy()))

def test_pdb_to_v23(pdb_file):
  remediator.remediate(pdb_file, False)

def test_pdb_remediate_cycle(pdb_file):
  t0 = time.time()
  remediatorObj = remediator.Remediator()
  pdb_model = remediator.get_model_from_file(pdb_file)
  pdb_model = remediatorObj.remediate_model_object(pdb_model, True)
  #write_model_to_file(pdb_model.model_as_pdb(), "input_model.pdb")
  print("............checking if remediated "+pdb_file+" is v3")
  assert remediatorObj.is_model_v3(pdb_model)
  pdb_model = remediatorObj.remediate_model_object(pdb_model, False)
  #write_model_to_file(pdb_model.model_as_pdb(), "old_model.pdb")
  pdb_model = remediatorObj.remediate_model_object(pdb_model, True)
  #write_model_to_file(pdb_model.model_as_pdb(), "reold_model.pdb")
  print("............checking remediated-to-v2.3 and remediated "+pdb_file+" is v3")
  assert remediatorObj.is_model_v3(pdb_model)
  print("OK. "+pdb_file+" time: %8.3f"%(time.time()-t0))

def write_model_to_file(model, file_name):
  with open("tst/"+file_name, 'w+') as fh:
    fh.write(model)

def run_tests():
  test_pdb_to_v3("tst/1ubqFHv23.pdb", "tst/1ubqFH.pdb")
  test_pdb_remediate_cycle("tst/1ubqFH.pdb")
  test_pdb_remediate_cycle("tst/1ubqFHv23.pdb")
  test_pdb_remediate_cycle("tst/404dH.pdb")
  test_pdb_remediate_cycle("tst/mg_tst.pdb")
  test_pdb_remediate_cycle("tst/4v9dFH.pdb")

if __name__=="__main__":
  run_tests()
  print("OK")
