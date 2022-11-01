from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
import time

"""
Test multiplication of hierarchy and SS annotations in different combinations
of MTRIX and BIOMT records presence.
Corner case where MTRIX and BIOMT are identical. We want to silently skip BIOMT multiplication
if MTRIX was done, but not skip if MTRIX was not done.
"""

# pdb id: 7a5v
pdb_lines = """
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1

MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.309017 -0.951057  0.000000      239.91182
MTRIX2   2  0.951057  0.309017  0.000000      -37.99830
MTRIX3   2  0.000000  0.000000  1.000000        0.00000
MTRIX1   3 -0.809017 -0.587785  0.000000      350.18719
MTRIX2   3  0.587785 -0.809017  0.000000      178.42929
MTRIX3   3  0.000000  0.000000  1.000000        0.00000

REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2  0.309017 -0.951057  0.000000      239.91182
REMARK 350   BIOMT2   2  0.951057  0.309017  0.000000      -37.99830
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   3 -0.809017 -0.587785  0.000000      350.18719
REMARK 350   BIOMT2   3  0.587785 -0.809017  0.000000      178.42929
REMARK 350   BIOMT3   3  0.000000  0.000000  1.000000        0.00000

ATOM     21  N   SER A  10     121.825 143.282  92.198  1.00 49.53           N
ATOM     22  CA  SER A  10     122.763 143.750  91.142  1.00 53.90           C
ATOM     23  C   SER A  10     122.515 145.240  90.852  1.00 44.64           C
ATOM     24  O   SER A  10     123.502 145.971  90.659  1.00 49.23           O
ATOM     25  CB  SER A  10     122.649 142.922  89.891  1.00 54.86           C
ATOM     26  OG  SER A  10     121.430 143.193  89.228  1.00 60.17           O
ATOM     27  N   PHE A  11     121.253 145.670  90.877  1.00 39.99           N
ATOM     28  CA  PHE A  11     120.852 147.084  90.685  1.00 46.43           C
ATOM     29  C   PHE A  11     121.347 147.931  91.873  1.00 53.53           C
ATOM     30  O   PHE A  11     121.866 149.042  91.635  1.00 45.84           O
ATOM     31  CB  PHE A  11     119.339 147.230  90.531  1.00 45.10           C
ATOM     32  CG  PHE A  11     118.884 148.665  90.434  1.00 50.94           C
ATOM     33  CD1 PHE A  11     119.380 149.498  89.437  1.00 58.07           C
ATOM     34  CD2 PHE A  11     118.032 149.212  91.385  1.00 57.16           C
ATOM     35  CE1 PHE A  11     118.988 150.831  89.362  1.00 56.76           C
ATOM     36  CE2 PHE A  11     117.629 150.541  91.299  1.00 59.88           C
ATOM     37  CZ  PHE A  11     118.121 151.352  90.299  1.00 55.92           C
ATOM     38  N   VAL A  12     121.144 147.455  93.110  1.00 43.35           N
ATOM     39  CA  VAL A  12     121.615 148.167  94.333  1.00 40.73           C
ATOM     40  C   VAL A  12     123.143 148.194  94.308  1.00 37.70           C
ATOM     41  O   VAL A  12     123.700 149.281  94.511  1.00 41.01           O
ATOM     42  CB  VAL A  12     121.073 147.560  95.640  1.00 46.83           C
ATOM     43  CG1 VAL A  12     121.686 148.244  96.866  1.00 40.11           C
ATOM     44  CG2 VAL A  12     119.551 147.638  95.701  1.00 41.22           C
ATOM     45  N   LYS A  13     123.785 147.091  93.923  1.00 36.58           N
ATOM     46  CA  LYS A  13     125.268 147.018  93.854  1.00 41.77           C
ATOM     47  C   LYS A  13     125.761 148.054  92.839  1.00 42.08           C
ATOM     48  O   LYS A  13     126.776 148.724  93.128  1.00 38.35           O
ATOM     49  CB  LYS A  13     125.757 145.612  93.490  1.00 46.55           C
ATOM     50  CG  LYS A  13     127.268 145.492  93.315  1.00 46.87           C
ATOM     51  CD  LYS A  13     127.740 144.113  92.898  1.00 58.69           C
ATOM     52  CE  LYS A  13     129.223 144.071  92.578  1.00 64.30           C
ATOM     53  NZ  LYS A  13     130.049 143.879  93.793  1.00 72.87           N
ATOM     54  N   GLU A  14     125.080 148.160  91.694  1.00 42.55           N
ATOM     55  CA  GLU A  14     125.407 149.148  90.630  1.00 47.24           C
ATOM     56  C   GLU A  14     125.269 150.559  91.210  1.00 34.67           C
ATOM     57  O   GLU A  14     126.169 151.376  90.970  1.00 40.45           O
ATOM     58  CB  GLU A  14     124.498 148.982  89.407  1.00 52.19           C
ATOM     59  CG  GLU A  14     125.009 147.963  88.402  1.00 82.30           C
ATOM     60  CD  GLU A  14     124.265 147.939  87.071  1.00104.62           C
"""

def exercise_biomt_only():
  """
  Use-case for phenix.pdb.biomt_reconstruction - we want to expand only using biomt
  """
  inp = iotbx.pdb.input(lines=pdb_lines, source_info=None)
  model = mmtbx.model.manager(
    model_input = inp,
    expand_with_mtrix = False)
  assert model.get_number_of_atoms() == 40, model.get_number_of_atoms()
  model.expand_with_BIOMT_records()
  assert model.get_number_of_atoms() == 120, model.get_number_of_atoms()

def exercise_mtrix_only():
  """
  Use-case for phenix.pdb.mtrix_reconstruction - we want to expand only using mtrix
  """
  inp = iotbx.pdb.input(lines=pdb_lines, source_info=None)
  model = mmtbx.model.manager(
    model_input = inp)
  assert model.get_number_of_atoms() == 120, model.get_number_of_atoms()

  # use-case for phenix.real_space_refine ??? when normally we want to expand by MTRIX first
  # and by BIOMT second, but not in case when BIOMT == MTRIX...
  model.expand_with_BIOMT_records()
  assert model.get_number_of_atoms() == 120, model.get_number_of_atoms()

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_biomt_only()
  exercise_mtrix_only()
  print("Total time: %8.3f"%(time.time() - t0))
  print("OK.")
