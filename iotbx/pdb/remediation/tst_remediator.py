from __future__ import division
import time
import os
import difflib

import remediator
import iotbx.pdb
import libtbx
from iotbx.data_manager import DataManager
from libtbx import easy_run

multimodel_tst_str = """MODEL        1
ATOM     21  P     A C   2      33.564  25.095  49.676  0.50 58.86           P
ATOM     22  O1P   A C   2      32.847  23.800  49.561  0.50 58.82           O
ATOM     23  O2P   A C   2      33.003  26.303  49.020  0.50 59.18           O
ATOM     24  O5*   A C   2      33.761  25.416  51.223  0.50 57.39           O
ATOM     25  C5*   A C   2      34.402  24.477  52.075  0.50 54.86           C
ATOM     26  C4*   A C   2      34.562  25.046  53.463  0.50 53.12           C
ATOM     27  O4*   A C   2      35.472  26.177  53.431  0.50 52.63           O
ATOM     28  C3*   A C   2      33.304  25.614  54.097  0.50 51.70           C
ATOM     29  O3*   A C   2      32.491  24.581  54.644  0.50 48.48           O
ATOM     30  C2*   A C   2      33.889  26.530  55.163  0.50 51.30           C
ATOM     31  O2*   A C   2      34.349  25.814  56.290  0.50 51.90           O
ATOM     32  C1*   A C   2      35.090  27.121  54.420  0.50 51.50           C
ATOM     33  N9    A C   2      34.791  28.391  53.755  0.50 50.06           N
ATOM     34  C8    A C   2      34.427  28.594  52.447  0.50 49.53           C
ATOM     35  N7    A C   2      34.216  29.853  52.145  0.50 49.19           N
ATOM     36  C5    A C   2      34.462  30.526  53.334  0.50 48.49           C
ATOM     37  C6    A C   2      34.410  31.888  53.680  0.50 47.61           C
ATOM     38  N6    A C   2      34.090  32.858  52.824  0.50 47.71           N
ATOM     39  N1    A C   2      34.705  32.225  54.953  0.50 47.40           N
ATOM     40  C2    A C   2      35.031  31.253  55.813  0.50 48.40           C
ATOM     41  N3    A C   2      35.118  29.941  55.609  0.50 48.16           N
ATOM     42  C4    A C   2      34.817  29.639  54.334  0.50 48.73           C
ATOM      0 1H5*   A C   2      35.271  24.244  51.713  0.50 54.86           H   new
ATOM      0 2H5*   A C   2      33.882  23.659  52.113  0.50 54.86           H   new
ATOM      0  H4*   A C   2      34.865  24.283  53.979  0.50 53.12           H   new
ATOM      0  H3*   A C   2      32.712  26.073  53.481  0.50 51.70           H   new
ATOM      0  H2*   A C   2      33.246  27.170  55.506  0.50 51.30           H   new
ATOM      0 2HO*   A C   2      34.347  26.318  56.962  0.50 51.90           H   new
ATOM      0  H1*   A C   2      35.790  27.300  55.067  0.50 51.50           H   new
ATOM      0  H8    A C   2      34.338  27.902  51.832  0.50 49.53           H   new
ATOM      0 1H6    A C   2      34.073  33.676  53.088  0.50 47.71           H   new
ATOM      0 2H6    A C   2      33.900  32.666  52.008  0.50 47.71           H   new
ATOM      0  H2    A C   2      35.225  31.535  56.678  0.50 48.40           H   new
ATOM    210  P     U C  11       5.031  37.445  41.706  0.50 38.57           P
ATOM    211  O1P   U C  11       4.328  37.521  40.399  0.50 38.83           O
ATOM    212  O2P   U C  11       4.502  38.190  42.876  0.50 39.51           O
ATOM    213  O5*   U C  11       5.164  35.910  42.092  0.50 39.44           O
ATOM    214  C5*   U C  11       5.615  34.962  41.129  0.50 42.64           C
ATOM    215  C4*   U C  11       5.702  33.586  41.744  0.50 44.82           C
ATOM    216  O4*   U C  11       6.788  33.536  42.710  0.50 45.43           O
ATOM    217  C3*   U C  11       4.497  33.138  42.552  0.50 45.57           C
ATOM    218  O3*   U C  11       3.432  32.687  41.721  0.50 47.09           O
ATOM    219  C2*   U C  11       5.098  32.011  43.379  0.50 45.57           C
ATOM    220  O2*   U C  11       5.261  30.823  42.630  0.50 46.33           O
ATOM    221  C1*   U C  11       6.472  32.593  43.723  0.50 45.59           C
ATOM    222  N1    U C  11       6.468  33.279  45.023  0.50 45.31           N
ATOM    223  C2    U C  11       6.566  32.494  46.157  0.50 44.65           C
ATOM    224  O2    U C  11       6.667  31.278  46.111  0.50 43.78           O
ATOM    225  N3    U C  11       6.539  33.183  47.344  0.50 44.04           N
ATOM    226  C4    U C  11       6.431  34.550  47.509  0.50 45.77           C
ATOM    227  O4    U C  11       6.433  35.025  48.646  0.50 45.26           O
ATOM    228  C5    U C  11       6.341  35.294  46.286  0.50 45.83           C
ATOM    229  C6    U C  11       6.364  34.649  45.114  0.50 45.15           C
ATOM      0 1H5*   U C  11       6.484  35.226  40.789  0.50 42.64           H   new
ATOM      0 2H5*   U C  11       5.008  34.947  40.373  0.50 42.64           H   new
ATOM      0  H4*   U C  11       5.805  33.013  40.968  0.50 44.82           H   new
ATOM      0  H3*   U C  11       4.090  33.842  43.080  0.50 45.57           H   new
ATOM      0  H2*   U C  11       4.547  31.761  44.137  0.50 45.57           H   new
ATOM      0 2HO*   U C  11       5.328  30.169  43.153  0.50 46.33           H   new
ATOM      0  H1*   U C  11       7.121  31.874  43.775  0.50 45.59           H   new
ATOM      0  H3    U C  11       6.595  32.713  48.062  0.50 44.04           H   new
ATOM      0  H5    U C  11       6.267  36.221  46.308  0.50 45.83           H   new
ATOM      0  H6    U C  11       6.308  35.146  44.330  0.50 45.15           H   new
ENDMDL
MODEL        2
ATOM     21  P     A C   2      33.564  25.095  49.676  0.50 58.86           P
ATOM     22  O1P   A C   2      32.847  23.800  49.561  0.50 58.82           O
ATOM     23  O2P   A C   2      33.003  26.303  49.020  0.50 59.18           O
ATOM     24  O5*   A C   2      33.761  25.416  51.223  0.50 57.39           O
ATOM     25  C5*   A C   2      34.402  24.477  52.075  0.50 54.86           C
ATOM     26  C4*   A C   2      34.562  25.046  53.463  0.50 53.12           C
ATOM     27  O4*   A C   2      35.472  26.177  53.431  0.50 52.63           O
ATOM     28  C3*   A C   2      33.304  25.614  54.097  0.50 51.70           C
ATOM     29  O3*   A C   2      32.491  24.581  54.644  0.50 48.48           O
ATOM     30  C2*   A C   2      33.889  26.530  55.163  0.50 51.30           C
ATOM     31  O2*   A C   2      34.349  25.814  56.290  0.50 51.90           O
ATOM     32  C1*   A C   2      35.090  27.121  54.420  0.50 51.50           C
ATOM     33  N9    A C   2      34.791  28.391  53.755  0.50 50.06           N
ATOM     34  C8    A C   2      34.427  28.594  52.447  0.50 49.53           C
ATOM     35  N7    A C   2      34.216  29.853  52.145  0.50 49.19           N
ATOM     36  C5    A C   2      34.462  30.526  53.334  0.50 48.49           C
ATOM     37  C6    A C   2      34.410  31.888  53.680  0.50 47.61           C
ATOM     38  N6    A C   2      34.090  32.858  52.824  0.50 47.71           N
ATOM     39  N1    A C   2      34.705  32.225  54.953  0.50 47.40           N
ATOM     40  C2    A C   2      35.031  31.253  55.813  0.50 48.40           C
ATOM     41  N3    A C   2      35.118  29.941  55.609  0.50 48.16           N
ATOM     42  C4    A C   2      34.817  29.639  54.334  0.50 48.73           C
ATOM      0 1H5*   A C   2      35.271  24.244  51.713  0.50 54.86           H   new
ATOM      0 2H5*   A C   2      33.882  23.659  52.113  0.50 54.86           H   new
ATOM      0  H4*   A C   2      34.865  24.283  53.979  0.50 53.12           H   new
ATOM      0  H3*   A C   2      32.712  26.073  53.481  0.50 51.70           H   new
ATOM      0  H2*   A C   2      33.246  27.170  55.506  0.50 51.30           H   new
ATOM      0 2HO*   A C   2      34.347  26.318  56.962  0.50 51.90           H   new
ATOM      0  H1*   A C   2      35.790  27.300  55.067  0.50 51.50           H   new
ATOM      0  H8    A C   2      34.338  27.902  51.832  0.50 49.53           H   new
ATOM      0 1H6    A C   2      34.073  33.676  53.088  0.50 47.71           H   new
ATOM      0 2H6    A C   2      33.900  32.666  52.008  0.50 47.71           H   new
ATOM      0  H2    A C   2      35.225  31.535  56.678  0.50 48.40           H   new
ATOM    210  P     U C  11       5.031  37.445  41.706  0.50 38.57           P
ATOM    211  O1P   U C  11       4.328  37.521  40.399  0.50 38.83           O
ATOM    212  O2P   U C  11       4.502  38.190  42.876  0.50 39.51           O
ATOM    213  O5*   U C  11       5.164  35.910  42.092  0.50 39.44           O
ATOM    214  C5*   U C  11       5.615  34.962  41.129  0.50 42.64           C
ATOM    215  C4*   U C  11       5.702  33.586  41.744  0.50 44.82           C
ATOM    216  O4*   U C  11       6.788  33.536  42.710  0.50 45.43           O
ATOM    217  C3*   U C  11       4.497  33.138  42.552  0.50 45.57           C
ATOM    218  O3*   U C  11       3.432  32.687  41.721  0.50 47.09           O
ATOM    219  C2*   U C  11       5.098  32.011  43.379  0.50 45.57           C
ATOM    220  O2*   U C  11       5.261  30.823  42.630  0.50 46.33           O
ATOM    221  C1*   U C  11       6.472  32.593  43.723  0.50 45.59           C
ATOM    222  N1    U C  11       6.468  33.279  45.023  0.50 45.31           N
ATOM    223  C2    U C  11       6.566  32.494  46.157  0.50 44.65           C
ATOM    224  O2    U C  11       6.667  31.278  46.111  0.50 43.78           O
ATOM    225  N3    U C  11       6.539  33.183  47.344  0.50 44.04           N
ATOM    226  C4    U C  11       6.431  34.550  47.509  0.50 45.77           C
ATOM    227  O4    U C  11       6.433  35.025  48.646  0.50 45.26           O
ATOM    228  C5    U C  11       6.341  35.294  46.286  0.50 45.83           C
ATOM    229  C6    U C  11       6.364  34.649  45.114  0.50 45.15           C
ATOM      0 1H5*   U C  11       6.484  35.226  40.789  0.50 42.64           H   new
ATOM      0 2H5*   U C  11       5.008  34.947  40.373  0.50 42.64           H   new
ATOM      0  H4*   U C  11       5.805  33.013  40.968  0.50 44.82           H   new
ATOM      0  H3*   U C  11       4.090  33.842  43.080  0.50 45.57           H   new
ATOM      0  H2*   U C  11       4.547  31.761  44.137  0.50 45.57           H   new
ATOM      0 2HO*   U C  11       5.328  30.169  43.153  0.50 46.33           H   new
ATOM      0  H1*   U C  11       7.121  31.874  43.775  0.50 45.59           H   new
ATOM      0  H3    U C  11       6.595  32.713  48.062  0.50 44.04           H   new
ATOM      0  H5    U C  11       6.267  36.221  46.308  0.50 45.83           H   new
ATOM      0  H6    U C  11       6.308  35.146  44.330  0.50 45.15           H   new
ENDMDL
END
"""

def test_pdb_to_v3(old_pdb_file, v3_pdb_file):
  print("............checking conversion of "+old_pdb_file+" to new PDB format")
  remed = remediator.Remediator()
  old_pdb_model = remediator.get_model_from_file(old_pdb_file)
  #write_model_to_file(old_pdb_model.model_as_pdb(), "old_file.pdb")
  old_pdb_model_v3 = remed.remediate_model_object(old_pdb_model, True)
  v3_pdb_model = remediator.get_model_from_file(v3_pdb_file)
  #write_model_to_file(old_pdb_model_v3.model_as_pdb(), "old_file_remed.pdb")
  #write_model_to_file(v3_pdb_model.model_as_pdb(), "v3_file.pdb")
  old_pdb_model_v3_string = old_pdb_model_v3.get_hierarchy().as_pdb_string()
  v3_pdb_model_string = v3_pdb_model.get_hierarchy().as_pdb_string()
  #print(old_pdb_model_v3_string.splitlines(1))
  #print(v3_pdb_model_string.splitlines(1))
  if old_pdb_model_v3_string != v3_pdb_model_string:

    diff = list(difflib.unified_diff(old_pdb_model_v3_string.splitlines(), v3_pdb_model_string.splitlines()))
    assert False, "remediation test of "+old_pdb_file+" to v3 failed with these differences:\n" + '\n'.join(diff)

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

def tst_multi_model():
  print("............checking if multi-model files remediate properly")
  dm = DataManager()
  dm.process_model_str("1",multimodel_tst_str)
  m = dm.get_model("1")
  remediatorObj = remediator.Remediator()
  pdb_model = remediatorObj.remediate_model_object(m, True)
  assert remediatorObj.is_model_v3(pdb_model), "test if multimodel file converts to v3"

def write_model_to_file(model, file_name):
  with open("tst/"+file_name, 'w+') as fh:
    fh.write(model)

def get_regression_folder_file_path(file_name):
  regression_folder = os.path.join(libtbx.env.dist_path("iotbx"), 'pdb', 'remediation', 'tst')
  file_path = os.path.join(regression_folder, file_name)
  assert os.path.isfile(file_path), file_name + " test file is missing from " + regression_folder
  return file_path

def test_full_str_convert(old_file, new_file):
  print("testing full file remediation")
  prefix = os.path.basename(__file__).replace(".py","")
  assert not easy_run.call("iotbx.pdb_remediator %s > %s.pdb"%(old_file, prefix))
  with open(prefix+".pdb") as file:
    remediated_list = [line.rstrip() for line in file]
    remediated_pdb = "\n".join(remediated_list)
  with open(new_file) as file:
    new_list = [line.rstrip() for line in file]
    new_pdb = "\n".join(new_list)
  if new_pdb != remediated_pdb:

    diff = list(difflib.unified_diff(new_pdb.splitlines(), remediated_pdb.splitlines()))
    assert False, "remediation test of "+old_file+" to v3 failed with these differences:\n" + '\n'.join(diff)

  print("OK")

def run_tests():
  test_full_str_convert(get_regression_folder_file_path("dna-rna-testv23.pdb"), get_regression_folder_file_path("dna-rna-test.pdb"))
  test_pdb_to_v3(get_regression_folder_file_path("protein-dna-testv23.pdb"), get_regression_folder_file_path("protein-dna-test.pdb"))
  test_pdb_to_v3(get_regression_folder_file_path("dna-rna-testv23.pdb"), get_regression_folder_file_path("dna-rna-test.pdb"))

  test_pdb_remediate_cycle(get_regression_folder_file_path("protein-dna-test.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("protein-dna-testv23.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("dna-rna-test.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("dna-rna-testv23.pdb"))
  test_pdb_remediate_cycle(get_regression_folder_file_path("mg_tst.pdb"))
  tst_multi_model()

if __name__=="__main__":
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    run_tests()
    print("OK")
  else:
    print('chem_data is not available for remediator tests, skipping')
