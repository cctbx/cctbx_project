from __future__ import absolute_import, division, print_function

from libtbx import easy_run

pdb = """
ATOM     59  N   ALA A   9      44.498  67.689  40.761  1.00 46.22           N
ATOM     60  CA  ALA A   9      43.052  67.784  40.923  1.00 45.87           C
ATOM     61  C   ALA A   9      42.363  67.105  39.740  1.00 37.56           C
ATOM     62  O   ALA A   9      41.727  67.767  38.913  1.00 49.01           O
ATOM     63  CB  ALA A   9      42.634  69.251  41.068  1.00 41.14           C
ATOM     64  N   PRO A  10      42.481  65.784  39.622  1.00 54.52           N
ATOM     65  CA  PRO A  10      41.870  65.108  38.475  1.00 49.99           C
ATOM     66  C   PRO A  10      40.358  65.027  38.616  1.00 54.63           C
ATOM     67  O   PRO A  10      39.825  64.775  39.697  1.00 54.09           O
ATOM     68  CB  PRO A  10      42.511  63.718  38.495  1.00 39.39           C
ATOM     69  CG  PRO A  10      42.883  63.499  39.909  1.00 41.67           C
ATOM     70  CD  PRO A  10      43.196  64.849  40.502  1.00 44.54           C
ATOM     71  N   SER A  11      39.662  65.271  37.503  1.00 61.40           N
ATOM     72  CA  SER A  11      38.219  65.068  37.450  1.00 59.89           C
ATOM     73  C   SER A  11      37.919  63.580  37.301  1.00 57.15           C
ATOM     74  O   SER A  11      38.565  62.895  36.507  1.00 53.39           O
ATOM     75  CB  SER A  11      37.599  65.840  36.289  1.00 55.80           C
ATOM     76  N   PRO A  12      37.102  63.122  37.964  1.00 59.58           N
ATOM     77  CA  PRO A  12      36.829  61.665  37.929  1.00 63.69           C
ATOM     78  C   PRO A  12      35.870  61.273  36.808  1.00 63.21           C
ATOM     79  O   PRO A  12      34.706  60.913  37.019  1.00 68.01           O
ATOM     80  CB  PRO A  12      36.236  61.415  39.313  1.00 64.16           C
ATOM     81  CG  PRO A  12      35.516  62.685  39.620  1.00 68.82           C
ATOM     82  CD  PRO A  12      36.321  63.804  39.007  1.00 64.77           C
ATOM     83  N   THR A  13      36.378  61.338  35.577  1.00 65.04           N
ATOM     84  CA  THR A  13      35.605  61.101  34.371  1.00 72.28           C
ATOM     85  C   THR A  13      36.415  60.173  33.463  1.00 76.81           C
ATOM     86  O   THR A  13      37.017  60.604  32.483  1.00 75.12           O
ATOM     87  CB  THR A  13      35.271  62.418  33.666  1.00 75.72           C
ATOM     88  OG1 THR A  13      36.491  63.024  33.216  1.00 83.26           O
ATOM     89  CG2 THR A  13      34.574  63.397  34.618  1.00 68.25           C
ATOM     90  N   GLY A  14      36.457  58.895  33.806  1.00 79.08           N
ATOM     91  CA  GLY A  14      37.045  57.902  32.929  1.00 78.29           C
ATOM     92  C   GLY A  14      38.493  57.587  33.245  1.00 76.16           C
ATOM     93  O   GLY A  14      39.013  57.869  34.331  1.00 79.41           O
"""

def run():
  f=open('cdl_esd.pdb', 'w')
  f.write(pdb)
  f.close()
  cmd = 'phenix.pdb_interpretation cdl_esd.pdb'
  rc = easy_run.call(cmd)
  print('OK')
  return rc

if __name__=="__main__":
  run()
