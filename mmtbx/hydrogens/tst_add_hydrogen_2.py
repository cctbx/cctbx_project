from __future__ import division, print_function
import time
import mmtbx.model
import iotbx.pdb
import mmtbx.hydrogens
from libtbx.utils import null_out


def run():
  correct_H_position_with_cdl()


def get_model(pdb_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp,
                              log         = null_out())
  model_with_h = mmtbx.hydrogens.add(model = model)
  return model_with_h


def correct_H_position_with_cdl():
  """CDL changes ideal angles --> make sure H is placed correctly"""
  fail = False
  try:
    model_with_h = get_model(pdb_str = pdb_str_1)
  except Exception as e:
    fail = True

  assert(fail == False), "Fix requires correct indentation of atom labels in hydrogenate."

  #number_h = model_with_h.get_hd_selection().count(True)


pdb_str_1 = """
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21
SCALE1      0.013843  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013887  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011496        0.00000
ATOM      1  N   PRO H  14      52.628 -74.147  33.427  1.00 20.43           N
ATOM      2  CA  PRO H  14      53.440 -73.630  34.533  1.00 20.01           C
ATOM      3  C   PRO H  14      54.482 -72.584  34.124  1.00 20.76           C
ATOM      4  O   PRO H  14      55.025 -72.627  33.021  1.00 16.34           O
ATOM      5  CB  PRO H  14      54.055 -74.895  35.134  1.00 22.06           C
ATOM      6  CG  PRO H  14      54.084 -75.862  33.972  1.00 25.16           C
ATOM      7  CD  PRO H  14      52.770 -75.608  33.294  1.00 17.36           C
ATOM      8  N   SER H  15      54.727 -71.646  35.038  1.00 21.70           N
ATOM      9  CA  SER H  15      55.670 -70.537  34.874  1.00 25.33           C
ATOM     10  C   SER H  15      55.049 -69.401  34.057  1.00 24.78           C
ATOM     11  O   SER H  15      55.581 -68.291  34.023  1.00 27.51           O
ATOM     12  CB  SER H  15      56.982 -71.005  34.219  1.00 25.20           C
ATOM     13  OG  SER H  15      56.914 -70.938  32.802  1.00 28.91           O
ATOM     14  N   GLN H  16      53.918 -69.678  33.412  1.00 24.55           N
ATOM     15  CA  GLN H  16      53.224 -68.673  32.611  1.00 29.39           C
ATOM     16  C   GLN H  16      52.340 -67.778  33.475  1.00 28.13           C
ATOM     17  O   GLN H  16      52.234 -67.987  34.681  1.00 26.35           O
ATOM     18  CB  GLN H  16      52.371 -69.346  31.533  1.00 31.67           C
ATOM     19  CG  GLN H  16      53.196 -70.112  30.524  1.00 44.80           C
ATOM     20  CD  GLN H  16      54.379 -69.303  30.030  1.00 48.55           C
ATOM     21  OE1 GLN H  16      54.213 -68.269  29.386  1.00 52.45           O
ATOM     22  NE2 GLN H  16      55.584 -69.766  30.342  1.00 55.07           N
TER
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
