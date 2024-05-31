from __future__ import absolute_import, division, print_function
from mmtbx.suitename import suitealyze
from iotbx.data_manager import DataManager
from libtbx.test_utils import convert_string_to_cif_long, convert_pdb_to_cif_for_pdb_str

import libtbx.load_env
import time
import json
import os

def exercise_suitename_json(test_mmcif=False):
  # derived from 2goz
  # note: chain B, residue 20 of 2goz is a DNA residue (DC). As of 2023, DNA conformations are not handled by suitename,
  #   and this residue is skipped.
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb2goz_refmac_tls.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print("Skipping exercise: input pdb (pdb2goz_refmac_tls.ent) not available")
    return
  dm = DataManager()
  if test_mmcif:
    with open(regression_pdb) as f:
      pdb_2goz_str = f.read()
    pdb_2goz_str = convert_string_to_cif_long(pdb_2goz_str, hetatm_name_addition = "", chain_addition="LONGCHAIN")
    dm.process_model_str("1", pdb_2goz_str)
    m = dm.get_model("1")
    chainA = "ALONGCHAIN"
    chainB = "BLONGCHAIN"
  else:
    m = dm.get_model(regression_pdb)
    chainA = "A"
    chainB = "B"
  sz = suitealyze.suitealyze(pdb_hierarchy=m.get_hierarchy())
  sz_dict = json.loads(sz.as_JSON())
  #import pprint
  #pprint.pprint(sz_dict)
  assert len(sz_dict['flat_results']) == 62, "tst_suitename json output not returning correct number of suites, now: "+str(len(sz_dict['flat_results']))
  assert sz_dict['flat_results'][0]["cluster"] == "__", "tst_suitename json output first cluster value changed, now: "+sz_dict['flat_results'][0]["cluster"]
  from mmtbx.validation import test_utils
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "1a")==37, "tst_suitename json hierarchical output total number of 1a changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "1a"))
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "1b")==5, "tst_suitename json hierarchical output total number of 1b changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "1b"))
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "7r")==1, "tst_suitename json hierarchical output total number of 7r changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "7r"))
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "!!")==5, "tst_suitename json hierarchical output total number of !! changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "!!"))
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "__")==2, "tst_suitename json hierarchical output total number of __ changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "__"))
  assert sz_dict['summary_results'][""]["num_outliers"] == 5, "tst_suitename json summary output total number of outliers changed"
  assert sz_dict['summary_results'][""]["num_suites"] == 62, "tst_suitename json summary output total number of suites changed"
  assert sz_dict['suitestrings'][""][chainA] == "__G1aG1aA1aU1aG1aU1aA7rC0aU1aA1aC1aC1aA1cG1bC4aU1gG1aA1[U6gG9aA1aG1aU1aC1aC1aC1aA!!A!!A2aU1bA!!G1aG1aA1aC1aG&aA1aA1aA1aC1aG1aC1cC", "model 1 chain A suitestring changed"
  assert sz_dict['suitestrings'][""][chainB] == "__G1aG1aC1aG1aU1bC!!C1aU1aG1cG1a?1aA1bU4aC1bC!!A1aA1a?1aC", "model 1 chain B suitestring changed"

multimod_2goz_pdb_str = """MODEL        1
ATOM    539  P     C A  26     -19.024  25.068  -5.945  1.00 46.81           P
ATOM    540  OP1   C A  26     -19.207  26.235  -5.055  1.00 45.04           O
ATOM    541  OP2   C A  26     -17.676  24.613  -6.369  1.00 44.77           O
ATOM    542  O5'   C A  26     -19.675  23.774  -5.260  1.00 45.98           O
ATOM    543  C5'   C A  26     -20.927  23.881  -4.629  1.00 44.99           C
ATOM    544  C4'   C A  26     -21.435  22.491  -4.330  1.00 44.44           C
ATOM    545  O4'   C A  26     -21.502  21.721  -5.562  1.00 44.23           O
ATOM    546  C3'   C A  26     -20.539  21.645  -3.436  1.00 43.30           C
ATOM    547  O3'   C A  26     -20.642  22.010  -2.045  1.00 43.60           O
ATOM    548  C2'   C A  26     -21.120  20.271  -3.754  1.00 42.13           C
ATOM    549  O2'   C A  26     -22.418  20.065  -3.239  1.00 39.90           O
ATOM    550  C1'   C A  26     -21.222  20.354  -5.267  1.00 41.42           C
ATOM    551  N1    C A  26     -19.997  19.822  -5.973  1.00 39.41           N
ATOM    552  C2    C A  26     -19.788  18.429  -6.015  1.00 39.39           C
ATOM    553  O2    C A  26     -20.593  17.642  -5.478  1.00 37.57           O
ATOM    554  N3    C A  26     -18.683  17.965  -6.652  1.00 36.84           N
ATOM    555  C4    C A  26     -17.823  18.785  -7.222  1.00 38.50           C
ATOM    556  N4    C A  26     -16.768  18.220  -7.830  1.00 39.67           N
ATOM    557  C5    C A  26     -18.010  20.208  -7.197  1.00 37.68           C
ATOM    558  C6    C A  26     -19.098  20.666  -6.567  1.00 38.34           C
ATOM    559  P     A A  27     -19.357  21.872  -1.098  1.00 45.51           P
ATOM    560  OP1   A A  27     -19.721  22.274   0.274  1.00 45.89           O
ATOM    561  OP2   A A  27     -18.210  22.527  -1.770  1.00 45.31           O
ATOM    562  O5'   A A  27     -19.069  20.294  -1.096  1.00 45.85           O
ATOM    563  C5'   A A  27     -20.076  19.401  -0.618  1.00 45.44           C
ATOM    564  C4'   A A  27     -19.669  17.972  -0.877  1.00 44.25           C
ATOM    565  O4'   A A  27     -19.522  17.744  -2.296  1.00 44.42           O
ATOM    566  C3'   A A  27     -18.303  17.615  -0.313  1.00 43.99           C
ATOM    567  O3'   A A  27     -18.508  17.353   1.048  1.00 44.86           O
ATOM    568  C2'   A A  27     -17.965  16.379  -1.157  1.00 42.87           C
ATOM    569  O2'   A A  27     -18.734  15.239  -0.850  1.00 42.69           O
ATOM    570  C1'   A A  27     -18.430  16.850  -2.521  1.00 42.15           C
ATOM    571  N9    A A  27     -17.387  17.528  -3.288  1.00 40.58           N
ATOM    572  C8    A A  27     -17.180  18.869  -3.373  1.00 40.38           C
ATOM    573  N7    A A  27     -16.163  19.191  -4.139  1.00 40.04           N
ATOM    574  C5    A A  27     -15.672  17.981  -4.572  1.00 39.49           C
ATOM    575  C6    A A  27     -14.602  17.635  -5.416  1.00 42.18           C
ATOM    576  N6    A A  27     -13.805  18.531  -5.990  1.00 41.83           N
ATOM    577  N1    A A  27     -14.374  16.326  -5.659  1.00 41.96           N
ATOM    578  C2    A A  27     -15.177  15.419  -5.080  1.00 41.06           C
ATOM    579  N3    A A  27     -16.209  15.633  -4.274  1.00 39.86           N
ATOM    580  C4    A A  27     -16.409  16.943  -4.055  1.00 39.15           C
ATOM    581  P     A A  28     -17.341  17.272   2.110  1.00 45.26           P
ATOM    582  OP1   A A  28     -17.978  17.212   3.439  1.00 48.43           O
ATOM    583  OP2   A A  28     -16.346  18.351   1.847  1.00 44.90           O
ATOM    584  O5'   A A  28     -16.694  15.843   1.817  1.00 46.52           O
ATOM    585  C5'   A A  28     -17.365  14.642   2.237  1.00 49.78           C
ATOM    586  C4'   A A  28     -16.459  13.761   3.084  1.00 51.07           C
ATOM    587  O4'   A A  28     -15.433  13.241   2.193  1.00 51.00           O
ATOM    588  C3'   A A  28     -15.746  14.526   4.211  1.00 52.58           C
ATOM    589  O3'   A A  28     -16.208  14.324   5.616  1.00 56.20           O
ATOM    590  C2'   A A  28     -14.251  14.279   3.995  1.00 51.35           C
ATOM    591  O2'   A A  28     -13.582  13.793   5.135  1.00 51.93           O
ATOM    592  C1'   A A  28     -14.184  13.259   2.857  1.00 49.66           C
ATOM    593  N9    A A  28     -13.133  13.599   1.911  1.00 48.08           N
ATOM    594  C8    A A  28     -12.647  14.842   1.637  1.00 47.72           C
ATOM    595  N7    A A  28     -11.695  14.858   0.741  1.00 47.42           N
ATOM    596  C5    A A  28     -11.535  13.530   0.409  1.00 47.59           C
ATOM    597  C6    A A  28     -10.676  12.874  -0.499  1.00 48.14           C
ATOM    598  N6    A A  28      -9.771  13.521  -1.243  1.00 47.37           N
ATOM    599  N1    A A  28     -10.770  11.527  -0.604  1.00 48.11           N
ATOM    600  C2    A A  28     -11.675  10.892   0.155  1.00 48.78           C
ATOM    601  N3    A A  28     -12.529  11.409   1.043  1.00 47.86           N
ATOM    602  C4    A A  28     -12.413  12.742   1.122  1.00 47.87           C
ATOM    603  P     A A  29     -16.571  12.975   6.417  1.00 56.21           P
ATOM    604  OP1   A A  29     -15.913  13.074   7.733  1.00 55.14           O
ATOM    605  OP2   A A  29     -16.319  11.787   5.562  1.00 58.11           O
ATOM    606  O5'   A A  29     -18.142  13.153   6.601  1.00 56.19           O
ATOM    607  C5'   A A  29     -18.988  12.025   6.745  1.00 57.91           C
ATOM    608  C4'   A A  29     -20.346  12.311   6.130  1.00 58.84           C
ATOM    609  O4'   A A  29     -20.942  13.497   6.732  1.00 57.42           O
ATOM    610  C3'   A A  29     -20.327  12.540   4.614  1.00 59.78           C
ATOM    611  O3'   A A  29     -21.378  11.750   4.001  1.00 63.25           O
ATOM    612  C2'   A A  29     -20.556  14.053   4.510  1.00 59.01           C
ATOM    613  O2'   A A  29     -21.148  14.457   3.293  1.00 58.45           O
ATOM    614  C1'   A A  29     -21.523  14.248   5.685  1.00 57.26           C
ATOM    615  N9    A A  29     -21.837  15.624   6.106  1.00 56.47           N
ATOM    616  C8    A A  29     -21.000  16.652   6.455  1.00 56.26           C
ATOM    617  N7    A A  29     -21.619  17.768   6.784  1.00 55.43           N
ATOM    618  C5    A A  29     -22.960  17.456   6.637  1.00 56.14           C
ATOM    619  C6    A A  29     -24.156  18.194   6.825  1.00 56.57           C
ATOM    620  N6    A A  29     -24.202  19.469   7.225  1.00 56.29           N
ATOM    621  N1    A A  29     -25.329  17.561   6.590  1.00 56.55           N
ATOM    622  C2    A A  29     -25.312  16.285   6.185  1.00 56.37           C
ATOM    623  N3    A A  29     -24.266  15.492   5.977  1.00 55.90           N
ATOM    624  C4    A A  29     -23.109  16.140   6.222  1.00 56.39           C
ATOM    625  P     U A  30     -21.077  10.352   3.257  1.00 65.20           P
ATOM    626  OP1   U A  30     -22.343   9.743   2.804  1.00 65.24           O
ATOM    627  OP2   U A  30     -20.167   9.523   4.084  1.00 65.39           O
ATOM    628  O5'   U A  30     -20.316  10.903   1.969  1.00 66.07           O
ATOM    629  C5'   U A  30     -20.969  11.820   1.095  1.00 66.37           C
ATOM    630  C4'   U A  30     -20.234  11.844  -0.226  1.00 66.40           C
ATOM    631  O4'   U A  30     -18.809  11.940   0.024  1.00 67.11           O
ATOM    632  C3'   U A  30     -20.352  10.547  -1.009  1.00 65.69           C
ATOM    633  O3'   U A  30     -21.627  10.410  -1.672  1.00 62.51           O
ATOM    634  C2'   U A  30     -19.112  10.638  -1.908  1.00 66.65           C
ATOM    635  O2'   U A  30     -19.229  11.597  -2.942  1.00 67.13           O
ATOM    636  C1'   U A  30     -18.088  11.109  -0.879  1.00 68.23           C
ATOM    637  N1    U A  30     -17.381  10.011  -0.110  1.00 68.83           N
ATOM    638  C2    U A  30     -15.993  10.004   0.007  1.00 69.07           C
ATOM    639  O2    U A  30     -15.275  10.851  -0.501  1.00 69.49           O
ATOM    640  N3    U A  30     -15.472   8.955   0.741  1.00 68.66           N
ATOM    641  C4    U A  30     -16.202   7.927   1.346  1.00 69.29           C
ATOM    642  O4    U A  30     -15.653   7.028   1.974  1.00 69.21           O
ATOM    643  C5    U A  30     -17.628   8.003   1.183  1.00 69.98           C
ATOM    644  C6    U A  30     -18.145   9.017   0.481  1.00 69.92           C
ENDMDL
MODEL        2
ATOM    539  P     C A  26     -19.024  25.068  -5.945  1.00 46.81           P
ATOM    540  OP1   C A  26     -19.207  26.235  -5.055  1.00 45.04           O
ATOM    541  OP2   C A  26     -17.676  24.613  -6.369  1.00 44.77           O
ATOM    542  O5'   C A  26     -19.675  23.774  -5.260  1.00 45.98           O
ATOM    543  C5'   C A  26     -20.927  23.881  -4.629  1.00 44.99           C
ATOM    544  C4'   C A  26     -21.435  22.491  -4.330  1.00 44.44           C
ATOM    545  O4'   C A  26     -21.502  21.721  -5.562  1.00 44.23           O
ATOM    546  C3'   C A  26     -20.539  21.645  -3.436  1.00 43.30           C
ATOM    547  O3'   C A  26     -20.642  22.010  -2.045  1.00 43.60           O
ATOM    548  C2'   C A  26     -21.120  20.271  -3.754  1.00 42.13           C
ATOM    549  O2'   C A  26     -22.418  20.065  -3.239  1.00 39.90           O
ATOM    550  C1'   C A  26     -21.222  20.354  -5.267  1.00 41.42           C
ATOM    551  N1    C A  26     -19.997  19.822  -5.973  1.00 39.41           N
ATOM    552  C2    C A  26     -19.788  18.429  -6.015  1.00 39.39           C
ATOM    553  O2    C A  26     -20.593  17.642  -5.478  1.00 37.57           O
ATOM    554  N3    C A  26     -18.683  17.965  -6.652  1.00 36.84           N
ATOM    555  C4    C A  26     -17.823  18.785  -7.222  1.00 38.50           C
ATOM    556  N4    C A  26     -16.768  18.220  -7.830  1.00 39.67           N
ATOM    557  C5    C A  26     -18.010  20.208  -7.197  1.00 37.68           C
ATOM    558  C6    C A  26     -19.098  20.666  -6.567  1.00 38.34           C
ATOM    559  P     A A  27     -19.357  21.872  -1.098  1.00 45.51           P
ATOM    560  OP1   A A  27     -19.721  22.274   0.274  1.00 45.89           O
ATOM    561  OP2   A A  27     -18.210  22.527  -1.770  1.00 45.31           O
ATOM    562  O5'   A A  27     -19.069  20.294  -1.096  1.00 45.85           O
ATOM    563  C5'   A A  27     -20.076  19.401  -0.618  1.00 45.44           C
ATOM    564  C4'   A A  27     -19.669  17.972  -0.877  1.00 44.25           C
ATOM    565  O4'   A A  27     -19.522  17.744  -2.296  1.00 44.42           O
ATOM    566  C3'   A A  27     -18.303  17.615  -0.313  1.00 43.99           C
ATOM    567  O3'   A A  27     -18.508  17.353   1.048  1.00 44.86           O
ATOM    568  C2'   A A  27     -17.965  16.379  -1.157  1.00 42.87           C
ATOM    569  O2'   A A  27     -18.734  15.239  -0.850  1.00 42.69           O
ATOM    570  C1'   A A  27     -18.430  16.850  -2.521  1.00 42.15           C
ATOM    571  N9    A A  27     -17.387  17.528  -3.288  1.00 40.58           N
ATOM    572  C8    A A  27     -17.180  18.869  -3.373  1.00 40.38           C
ATOM    573  N7    A A  27     -16.163  19.191  -4.139  1.00 40.04           N
ATOM    574  C5    A A  27     -15.672  17.981  -4.572  1.00 39.49           C
ATOM    575  C6    A A  27     -14.602  17.635  -5.416  1.00 42.18           C
ATOM    576  N6    A A  27     -13.805  18.531  -5.990  1.00 41.83           N
ATOM    577  N1    A A  27     -14.374  16.326  -5.659  1.00 41.96           N
ATOM    578  C2    A A  27     -15.177  15.419  -5.080  1.00 41.06           C
ATOM    579  N3    A A  27     -16.209  15.633  -4.274  1.00 39.86           N
ATOM    580  C4    A A  27     -16.409  16.943  -4.055  1.00 39.15           C
ATOM    581  P     A A  28     -17.341  17.272   2.110  1.00 45.26           P
ATOM    582  OP1   A A  28     -17.978  17.212   3.439  1.00 48.43           O
ATOM    583  OP2   A A  28     -16.346  18.351   1.847  1.00 44.90           O
ATOM    584  O5'   A A  28     -16.694  15.843   1.817  1.00 46.52           O
ATOM    585  C5'   A A  28     -17.365  14.642   2.237  1.00 49.78           C
ATOM    586  C4'   A A  28     -16.459  13.761   3.084  1.00 51.07           C
ATOM    587  O4'   A A  28     -15.433  13.241   2.193  1.00 51.00           O
ATOM    588  C3'   A A  28     -15.746  14.526   4.211  1.00 52.58           C
ATOM    589  O3'   A A  28     -16.208  14.324   5.616  1.00 56.20           O
ATOM    590  C2'   A A  28     -14.251  14.279   3.995  1.00 51.35           C
ATOM    591  O2'   A A  28     -13.582  13.793   5.135  1.00 51.93           O
ATOM    592  C1'   A A  28     -14.184  13.259   2.857  1.00 49.66           C
ATOM    593  N9    A A  28     -13.133  13.599   1.911  1.00 48.08           N
ATOM    594  C8    A A  28     -12.647  14.842   1.637  1.00 47.72           C
ATOM    595  N7    A A  28     -11.695  14.858   0.741  1.00 47.42           N
ATOM    596  C5    A A  28     -11.535  13.530   0.409  1.00 47.59           C
ATOM    597  C6    A A  28     -10.676  12.874  -0.499  1.00 48.14           C
ATOM    598  N6    A A  28      -9.771  13.521  -1.243  1.00 47.37           N
ATOM    599  N1    A A  28     -10.770  11.527  -0.604  1.00 48.11           N
ATOM    600  C2    A A  28     -11.675  10.892   0.155  1.00 48.78           C
ATOM    601  N3    A A  28     -12.529  11.409   1.043  1.00 47.86           N
ATOM    602  C4    A A  28     -12.413  12.742   1.122  1.00 47.87           C
ATOM    603  P     A A  29     -16.571  12.975   6.417  1.00 56.21           P
ATOM    604  OP1   A A  29     -15.913  13.074   7.733  1.00 55.14           O
ATOM    605  OP2   A A  29     -16.319  11.787   5.562  1.00 58.11           O
ATOM    606  O5'   A A  29     -18.142  13.153   6.601  1.00 56.19           O
ATOM    607  C5'   A A  29     -18.988  12.025   6.745  1.00 57.91           C
ATOM    608  C4'   A A  29     -20.346  12.311   6.130  1.00 58.84           C
ATOM    609  O4'   A A  29     -20.942  13.497   6.732  1.00 57.42           O
ATOM    610  C3'   A A  29     -20.327  12.540   4.614  1.00 59.78           C
ATOM    611  O3'   A A  29     -21.378  11.750   4.001  1.00 63.25           O
ATOM    612  C2'   A A  29     -20.556  14.053   4.510  1.00 59.01           C
ATOM    613  O2'   A A  29     -21.148  14.457   3.293  1.00 58.45           O
ATOM    614  C1'   A A  29     -21.523  14.248   5.685  1.00 57.26           C
ATOM    615  N9    A A  29     -21.837  15.624   6.106  1.00 56.47           N
ATOM    616  C8    A A  29     -21.000  16.652   6.455  1.00 56.26           C
ATOM    617  N7    A A  29     -21.619  17.768   6.784  1.00 55.43           N
ATOM    618  C5    A A  29     -22.960  17.456   6.637  1.00 56.14           C
ATOM    619  C6    A A  29     -24.156  18.194   6.825  1.00 56.57           C
ATOM    620  N6    A A  29     -24.202  19.469   7.225  1.00 56.29           N
ATOM    621  N1    A A  29     -25.329  17.561   6.590  1.00 56.55           N
ATOM    622  C2    A A  29     -25.312  16.285   6.185  1.00 56.37           C
ATOM    623  N3    A A  29     -24.266  15.492   5.977  1.00 55.90           N
ATOM    624  C4    A A  29     -23.109  16.140   6.222  1.00 56.39           C
ATOM    625  P     U A  30     -21.077  10.352   3.257  1.00 65.20           P
ATOM    626  OP1   U A  30     -22.343   9.743   2.804  1.00 65.24           O
ATOM    627  OP2   U A  30     -20.167   9.523   4.084  1.00 65.39           O
ATOM    628  O5'   U A  30     -20.316  10.903   1.969  1.00 66.07           O
ATOM    629  C5'   U A  30     -20.969  11.820   1.095  1.00 66.37           C
ATOM    630  C4'   U A  30     -20.234  11.844  -0.226  1.00 66.40           C
ATOM    631  O4'   U A  30     -18.809  11.940   0.024  1.00 67.11           O
ATOM    632  C3'   U A  30     -20.352  10.547  -1.009  1.00 65.69           C
ATOM    633  O3'   U A  30     -21.627  10.410  -1.672  1.00 62.51           O
ATOM    634  C2'   U A  30     -19.112  10.638  -1.908  1.00 66.65           C
ATOM    635  O2'   U A  30     -19.229  11.597  -2.942  1.00 67.13           O
ATOM    636  C1'   U A  30     -18.088  11.109  -0.879  1.00 68.23           C
ATOM    637  N1    U A  30     -17.381  10.011  -0.110  1.00 68.83           N
ATOM    638  C2    U A  30     -15.993  10.004   0.007  1.00 69.07           C
ATOM    639  O2    U A  30     -15.275  10.851  -0.501  1.00 69.49           O
ATOM    640  N3    U A  30     -15.472   8.955   0.741  1.00 68.66           N
ATOM    641  C4    U A  30     -16.202   7.927   1.346  1.00 69.29           C
ATOM    642  O4    U A  30     -15.653   7.028   1.974  1.00 69.21           O
ATOM    643  C5    U A  30     -17.628   8.003   1.183  1.00 69.98           C
ATOM    644  C6    U A  30     -18.145   9.017   0.481  1.00 69.92           C
ENDMDL
END
"""

def exercise_multimodel_suitename_json(test_mmcif=False):
  # derived from 2goz
  dm = DataManager()
  dm.process_model_str("1",multimod_2goz_pdb_str)
  m = dm.get_model("1")
  sz = suitealyze.suitealyze(pdb_hierarchy=m.get_hierarchy())
  sz_dict = json.loads(sz.as_JSON())
  if test_mmcif:
    model_1 = "1"
    model_2 = "2"
  else:
    model_1 = "   1"
    model_2 = "   2"
  #import pprint
  #pprint.pprint(sz_dict)
  assert len(sz_dict['flat_results']) == 10, "tst_suitename json output not returning correct number of suites, now: "+str(len(sz_dict['flat_results']))
  assert sz_dict['flat_results'][0]["cluster"] == "__", "tst_suitename json output first cluster value changed, now: "+sz_dict['flat_results'][0]["cluster"]
  from mmtbx.validation import test_utils
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "1a")==2, "tst_suitename json hierarchical output total number of 1a changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "1a"))
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "2a")==2, "tst_suitename json hierarchical output total number of 1b changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "2a"))
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "!!")==4, "tst_suitename json hierarchical output total number of !! changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "!!"))
  assert test_utils.count_dict_values(sz_dict['hierarchical_results'], "__")==2, "tst_suitename json hierarchical output total number of __ changed, now: "+str(test_utils.count_dict_values(sz_dict['hierarchical_results'], "__"))
  assert sz_dict['summary_results'][model_1]["num_outliers"] == 2, "tst_suitename json summary output total number of outliers changed"
  assert sz_dict['summary_results'][model_1]["num_suites"] == 5, "tst_suitename json summary output total number of suites changed"
  assert sz_dict['summary_results'][model_2]["num_outliers"] == 2, "tst_suitename json summary output total number of outliers changed"
  assert sz_dict['summary_results'][model_2]["num_suites"] == 5, "tst_suitename json summary output total number of suites changed"
  return sz_dict

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_suitename_json()
  suite_json = exercise_multimodel_suitename_json()
  exercise_suitename_json(test_mmcif=True)
  convert_pdb_to_cif_for_pdb_str(locals(), chain_addition="LONGCHAIN", hetatm_name_addition = "", key_str="multimod_", print_new_string = False)
  suite_json_cif = exercise_multimodel_suitename_json(test_mmcif=True)
  assert suite_json['summary_results']['   1'] == suite_json_cif['summary_results']['1'], "tst_suitename summary results changed between pdb and cif version"
  assert suite_json['summary_results']['   2'] == suite_json_cif['summary_results']['2'], "tst_suitename summary results changed between pdb and cif version"

  print("OK. Time: %8.3f"%(time.time()-t0))
