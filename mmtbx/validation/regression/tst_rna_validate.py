
from __future__ import absolute_import, division, print_function
import mmtbx.validation.rna_validate
import iotbx.pdb
from mmtbx.validation.rna_validate import rna_validation
from libtbx.easy_pickle import loads, dumps
from libtbx.test_utils import approx_equal
import libtbx.load_env
from six.moves import cStringIO as StringIO
from iotbx.data_manager import DataManager
from libtbx.test_utils import convert_string_to_cif_long
import sys, os
import json
import time
from mmtbx.validation import test_utils


# This actually tests expected output - the remaining tests guard against
# fixed bugs.
def exercise_1(test_mmcif=False):
  # derived from 2goz
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
    pdb_2goz_str = convert_string_to_cif_long(pdb_2goz_str, chain_addition="LONGCHAIN", hetatm_name_addition = "")
    dm.process_model_str("1", pdb_2goz_str)
    m = dm.get_model("1")
  else:
    m = dm.get_model(regression_pdb)
  rv = rna_validation(m.get_hierarchy())
  assert len(rv.puckers.results) == 2, len(rv.puckers.results)
  assert len(rv.bonds.results) == 4, len(rv.bonds.results)
  assert len(rv.angles.results) == 14, len(rv.angles.results)
  assert len(rv.suites.results) == 5, len(rv.suites.results)
  #assert approx_equal(rv.suites.average_suiteness(), 0.55, eps=0.01)
  pickle_unpickle(rv)
  pdb_in = iotbx.pdb.input(regression_pdb)
  result = mmtbx.validation.rna_validate.rna_validation(
    pdb_hierarchy=pdb_in.construct_hierarchy(),
    geometry_restraints_manager=None,
    params=None)
  pickle_unpickle(result)
  bonds_json = json.loads(rv.bonds.as_JSON())
  angles_json = json.loads(rv.angles.as_JSON())
  puckers_json = json.loads(rv.puckers.as_JSON())
  suites_json = json.loads(rv.suites.as_JSON())
  #import pprint
  #pprint.pprint(puckers_json)
  assert len(bonds_json["flat_results"]) == 4, "tst_rna_validate json output not returning correct number of bond outliers, changed to: "+str(len(bonds_json["flat_results"]))
  assert approx_equal(bonds_json['flat_results'][0]["score"], 4.40664), "tst_rna_validate json output first bond outlier score changed to: "+str(bonds_json['flat_results'][0]["score"])
  assert approx_equal(bonds_json['flat_results'][0]["xyz"][0], -17.61949), "tst_rna_validate json output first bond x changed to: "+str(bonds_json['flat_results'][0]["xyz"][0])
  assert approx_equal(bonds_json['flat_results'][0]["xyz"][1], -2.9325), "tst_rna_validate json output first bond y changed to: "+str(bonds_json['flat_results'][0]["xyz"][1])
  assert approx_equal(bonds_json['flat_results'][0]["xyz"][2], -23.262), "tst_rna_validate json output first bond z changed to: "+str(bonds_json['flat_results'][0]["xyz"][2])
  assert test_utils.count_dict_values(bonds_json['hierarchical_results'], True)==4, "tst_rna_validate json hierarchical output total number of bond outliers changed to: "+str(test_utils.count_dict_values(bonds_json['hierarchical_results'], True))
  assert bonds_json['summary_results']['']["num_outliers_too_large"]==4, "tst_rna_validate json output summary result number of bond outliers too large changed to: "+str(bonds_json['summary_results']['']["num_outliers_too_large"])
  assert len(angles_json["flat_results"]) == 14, "tst_rna_validate json output not returning correct number of angle outliers, changed to: " + str(len(angles_json["flat_results"]))
  assert approx_equal(angles_json['flat_results'][0]["score"], -4.765809), "tst_rna_validate json output first angle outlier score changed to: " + str(angles_json['flat_results'][0]["score"])
  assert approx_equal(angles_json['flat_results'][0]["xyz"][0], -16.161), "tst_rna_validate json output first angle x score changed to: " + str(angles_json['flat_results'][0]["xyz"][0])
  assert approx_equal(angles_json['flat_results'][0]["xyz"][1], -6.863), "tst_rna_validate json output first angle y score changed to: " + str(angles_json['flat_results'][0]["xyz"][1])
  assert approx_equal(angles_json['flat_results'][0]["xyz"][2], -20.773), "tst_rna_validate json output first angle z score changed to: " + str(angles_json['flat_results'][0]["xyz"][2])
  assert test_utils.count_dict_values(angles_json['hierarchical_results'], [1.0, 1.0, 1.0])==14, "tst_rna_validate json hierarchical output number of angle outliers changed to: "+str(test_utils.count_dict_values(angles_json['hierarchical_results'], [1.0, 1.0, 1.0]))
  assert angles_json['summary_results']['']["num_outliers_too_large"]==6, "tst_rna_validate json output summary result number of angle outliers too large changed to: "+str(angles_json['summary_results']['']["num_outliers_too_large"])
  assert len(puckers_json["flat_results"]) == 2, "tst_rna_validate json output not returning correct number of pucker outliers, changed to: " + str(len(puckers_json["flat_results"]))
  assert approx_equal(puckers_json['flat_results'][0]["delta_angle"], 106.454213), "tst_rna_validate json output first pucker outlier delta changed to: " + str(puckers_json['flat_results'][0]["delta_angle"])
  assert test_utils.count_dict_values(puckers_json['hierarchical_results'], " ")==2, "tst_rna_validate json hierarchical output number of pucker outliers changed to: "+str(test_utils.count_dict_values(puckers_json['hierarchical_results'], " "))
  assert puckers_json['summary_results']['']["num_outliers"]==2, "tst_rna_validate json output summary result number of pucker outliers changed to: "+str(puckers_json['summary_results']['']["num_outliers"])
  assert len(suites_json["flat_results"]) == 5, "tst_rna_validate json output not returning correct number of suite outliers, changed to: " + str(len(suites_json["flat_results"]))
  assert suites_json['flat_results'][0]["suiteness"]==0, "tst_rna_validate json output first suite outlier suiteness changed to: " + str(suites_json['flat_results'][0]["suiteness"])
  assert test_utils.count_dict_values(suites_json['hierarchical_results'], "!!")==5, "tst_rna_validate json hierarchical output number of suite outliers changed to: "+str(test_utils.count_dict_values(suites_json['hierarchical_results'], "!!"))
  #assert suites_json['summary_results']['']["num_outliers"]==2, "tst_rna_validate json output summary result number of pucker outliers changed to: "+str(puckers_json['summary_results']['']["num_outliers"])
  return rv

def exercise_pdbvcif(rv, rv_cif):
  bonds_json = json.loads(rv.bonds.as_JSON())
  angles_json = json.loads(rv.angles.as_JSON())
  puckers_json = json.loads(rv.puckers.as_JSON())
  suites_json = json.loads(rv.suites.as_JSON())
  bonds_json_cif = json.loads(rv_cif.bonds.as_JSON())
  angles_json_cif = json.loads(rv_cif.angles.as_JSON())
  puckers_json_cif = json.loads(rv_cif.puckers.as_JSON())
  suites_json_cif = json.loads(rv_cif.suites.as_JSON())
  assert bonds_json['summary_results'] == bonds_json_cif['summary_results'], "tst_rna_validate summary results changed between pdb and cif version"
  assert angles_json['summary_results'] == angles_json_cif['summary_results'], "tst_rna_validate summary results changed between pdb and cif version"
  assert puckers_json['summary_results'] == puckers_json_cif['summary_results'], "tst_rna_validate summary results changed between pdb and cif version"
  assert suites_json['summary_results'] == suites_json_cif['summary_results'], "tst_rna_validate summary results changed between pdb and cif version"

def exercise_2():
  # fragment from 3g8t
  pdb_raw = """\
ATOM   7975  P     G Q 140      10.347 137.422  73.792  1.00118.69           P
ATOM   7976  OP1   G Q 140       9.348 138.439  74.195  1.00118.63           O
ATOM   7977  OP2   G Q 140      11.208 137.681  72.617  1.00118.60           O
ATOM   7978  O5'   G Q 140      11.286 137.083  75.051  1.00119.34           O
ATOM   7979  C5'   G Q 140      11.015 137.599  76.358  1.00120.43           C
ATOM   7980  C4'   G Q 140      11.617 136.733  77.454  1.00121.27           C
ATOM   7981  O4'   G Q 140      11.723 135.361  77.003  1.00121.51           O
ATOM   7982  C3'   G Q 140      13.038 137.080  77.885  1.00121.72           C
ATOM   7983  O3'   G Q 140      13.050 138.143  78.836  1.00122.07           O
ATOM   7984  C2'   G Q 140      13.517 135.770  78.511  1.00121.87           C
ATOM   7985  O2'   G Q 140      13.160 135.640  79.877  1.00121.91           O
ATOM   7986  C1'   G Q 140      12.793 134.717  77.671  1.00121.84           C
ATOM   7987  N9    G Q 140      13.651 134.055  76.684  1.00121.95           N
ATOM   7988  C8    G Q 140      14.226 134.620  75.568  1.00121.96           C
ATOM   7989  N7    G Q 140      14.944 133.787  74.870  1.00121.97           N
ATOM   7990  C5    G Q 140      14.843 132.588  75.563  1.00121.99           C
ATOM   7991  C6    G Q 140      15.412 131.320  75.281  1.00122.03           C
ATOM   7992  O6    G Q 140      16.142 131.004  74.331  1.00121.99           O
ATOM   7993  N1    G Q 140      15.064 130.361  76.235  1.00122.12           N
ATOM   7994  C2    G Q 140      14.264 130.602  77.331  1.00122.14           C
ATOM   7995  N2    G Q 140      14.033 129.557  78.139  1.00122.06           N
ATOM   7996  N3    G Q 140      13.721 131.787  77.604  1.00122.09           N
ATOM   7997  C4    G Q 140      14.050 132.734  76.684  1.00122.00           C
ATOM   7998  P     A Q 141      14.015 139.407  78.645  1.00122.32           P
ATOM   7999  OP1   A Q 141      13.141 140.571  78.383  1.00122.35           O
ATOM   8000  OP2   A Q 141      15.073 139.069  77.665  1.00122.37           O
ATOM   8001  O5'   A Q 141      14.690 139.587  80.088  1.00122.53           O
ATOM   8002  C5'   A Q 141      15.914 138.934  80.438  1.00122.90           C
ATOM   8003  C4'   A Q 141      15.673 137.727  81.331  1.00123.23           C
ATOM   8004  O4'   A Q 141      16.899 137.378  82.025  1.00123.39           O
ATOM   8005  C3'   A Q 141      14.605 137.909  82.412  1.00123.36           C
ATOM   8006  O3'   A Q 141      13.551 136.969  82.215  1.00123.31           O
ATOM   8007  C2'   A Q 141      15.347 137.686  83.732  1.00123.49           C
ATOM   8008  O2'   A Q 141      14.582 136.976  84.689  1.00123.46           O
ATOM   8009  C1'   A Q 141      16.551 136.857  83.291  1.00123.61           C
ATOM   8010  N9    A Q 141      17.714 136.912  84.186  1.00123.79           N
ATOM   8011  C8    A Q 141      17.797 137.531  85.406  1.00123.80           C
ATOM   8012  N7    A Q 141      18.964 137.415  85.992  1.00123.83           N
ATOM   8013  C5    A Q 141      19.706 136.662  85.099  1.00123.90           C
ATOM   8014  C6    A Q 141      21.037 136.187  85.130  1.00123.95           C
ATOM   8015  N6    A Q 141      21.875 136.421  86.146  1.00123.97           N
ATOM   8016  N1    A Q 141      21.476 135.460  84.076  1.00123.93           N
ATOM   8017  C2    A Q 141      20.634 135.225  83.057  1.00123.90           C
ATOM   8018  N3    A Q 141      19.364 135.620  82.917  1.00123.84           N
ATOM   8019  C4    A Q 141      18.953 136.340  83.979  1.00123.84           C
TER    8020        A Q 141
"""
  dm = DataManager()
  #print(help(dm))
  dm.process_model_str("", pdb_raw)
  rv = rna_validation(dm.get_model().get_hierarchy())
  assert len(rv.puckers.results) == 1
  pickle_unpickle(rv)

def exercise_3():
  # derived from 3bbi
  pdb_raw = """\
ATOM      1  O5'A  U A   1      39.826  29.792  61.182  0.50 82.88           O
ATOM      2  O5'B  U A   1      39.852  29.856  60.945  0.50 82.91           O
ATOM      3  C5'A  U A   1      40.022  30.855  62.130  0.50 80.10           C
ATOM      4  C5'B  U A   1      40.079  31.072  61.650  0.50 79.19           C
ATOM      5  C4'A  U A   1      38.772  31.691  62.255  0.50 78.82           C
ATOM      6  C4'B  U A   1      38.796  31.774  62.020  0.50 77.88           C
ATOM      7  O4'A  U A   1      38.856  32.884  61.429  0.50 80.02           O
ATOM      8  O4'B  U A   1      38.669  32.983  61.220  0.50 79.40           O
ATOM      9  C3'A  U A   1      37.538  30.959  61.777  0.50 78.46           C
ATOM     10  C3'B  U A   1      37.499  31.010  61.758  0.50 77.50           C
ATOM     11  O3'A  U A   1      37.037  30.176  62.849  0.50 75.37           O
ATOM     12  O3'B  U A   1      37.158  30.120  62.826  0.50 74.93           O
ATOM     13  C2'A  U A   1      36.641  32.101  61.304  0.50 79.86           C
ATOM     14  C2'B  U A   1      36.483  32.143  61.671  0.50 78.29           C
ATOM     15  O2'A  U A   1      35.961  32.782  62.341  0.50 81.49           O
ATOM     16  O2'B  U A   1      36.053  32.616  62.932  0.50 80.66           O
ATOM     17  C1'A  U A   1      37.673  33.025  60.649  0.50 79.69           C
ATOM     18  C1'B  U A   1      37.299  33.233  60.979  0.50 77.90           C
ATOM     19  N1 A  U A   1      38.048  32.567  59.307  0.50 77.53           N
ATOM     20  N1 B  U A   1      37.046  33.343  59.533  0.50 75.93           N
ATOM     21  C2 A  U A   1      37.329  32.949  58.191  0.50 76.19           C
ATOM     22  C2 B  U A   1      37.795  32.584  58.638  0.50 75.96           C
ATOM     23  O2 A  U A   1      36.364  33.685  58.191  0.50 76.18           O
ATOM     24  O2 B  U A   1      38.754  31.891  58.945  0.50 77.97           O
ATOM     25  N3 A  U A   1      37.803  32.400  57.033  0.50 74.94           N
ATOM     26  N3 B  U A   1      37.376  32.694  57.351  0.50 75.08           N
ATOM     27  C4 A  U A   1      38.865  31.551  56.889  0.50 74.13           C
ATOM     28  C4 B  U A   1      36.340  33.467  56.890  0.50 74.34           C
ATOM     29  O4 A  U A   1      39.163  31.127  55.805  0.50 72.11           O
ATOM     30  O4 B  U A   1      36.162  33.540  55.698  0.50 72.95           O
ATOM     31  C5 A  U A   1      39.560  31.215  58.099  0.50 74.56           C
ATOM     32  C5 B  U A   1      35.632  34.263  57.853  0.50 72.59           C
ATOM     33  C6 A  U A   1      39.134  31.742  59.244  0.50 76.24           C
ATOM     34  C6 B  U A   1      36.006  34.156  59.108  0.50 73.72           C
ATOM     35  P     C A   2      36.230  28.821  62.547  1.00 72.47           P
ATOM     36  OP1   C A   2      35.510  28.509  63.805  1.00 72.07           O
ATOM     37  OP2   C A   2      37.093  27.782  61.920  1.00 67.52           O
ATOM     38  O5'   C A   2      35.141  29.294  61.493  1.00 61.16           O
ATOM     39  C5'   C A   2      34.065  30.126  61.886  1.00 53.00           C
ATOM     40  C4'   C A   2      33.217  30.452  60.692  1.00 54.42           C
ATOM     41  O4'   C A   2      34.058  31.113  59.709  1.00 57.20           O
ATOM     42  C3'   C A   2      32.696  29.239  59.934  1.00 54.16           C
ATOM     43  O3'   C A   2      31.508  28.706  60.508  1.00 49.29           O
ATOM     44  C2'   C A   2      32.447  29.817  58.551  1.00 50.69           C
ATOM     45  O2'   C A   2      31.256  30.569  58.491  1.00 53.09           O
ATOM     46  C1'   C A   2      33.635  30.762  58.406  1.00 51.94           C
ATOM     47  N1    C A   2      34.752  30.152  57.665  1.00 50.45           N
ATOM     48  C2    C A   2      34.662  30.120  56.283  1.00 49.86           C
ATOM     49  O2    C A   2      33.659  30.593  55.738  1.00 50.93           O
ATOM     50  N3    C A   2      35.654  29.571  55.565  1.00 47.16           N
ATOM     51  C4    C A   2      36.709  29.052  56.179  1.00 47.93           C
ATOM     52  N4    C A   2      37.659  28.519  55.412  1.00 47.02           N
ATOM     53  C5    C A   2      36.834  29.060  57.597  1.00 46.21           C
ATOM     54  C6    C A   2      35.840  29.624  58.297  1.00 47.22           C
"""
  dm = DataManager()
  #print(help(dm))
  dm.process_model_str("", pdb_raw)
  rv = rna_validation(dm.get_model().get_hierarchy())
  pickle_unpickle(rv)

def pickle_unpickle(result):
  result2 = loads(dumps(result))
  out1 = StringIO()
  out2 = StringIO()
  result.show(out=out1)
  result2.show(out=out2)
  assert (out1.getvalue() == out2.getvalue())

def run():
  t0 = time.time()
  verbose = "--verbose" in sys.argv[1:]
  rv = exercise_1()
  rv_cif = exercise_1(test_mmcif=True)
  exercise_pdbvcif(rv, rv_cif)
  exercise_2()
  exercise_3()
  print("OK. Time: %8.3f"%(time.time()-t0))

if (__name__ == "__main__"):
  run()
