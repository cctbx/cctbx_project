from __future__ import absolute_import, division, print_function
from mmtbx.command_line import molprobity

pdb_txt = """\
CRYST1  104.400  104.400  124.250  90.00  90.00 120.00 H 3
SCALE1      0.009579  0.005530  0.000000        0.00000
SCALE2      0.000000  0.011060  0.000000        0.00000
SCALE3      0.000000  0.000000  0.008048        0.00000
ATOM    567  N   TYR    75      20.675  45.830  21.914  1.00 29.26      A    N
ATOM    568  CA  TYR    75      21.642  45.771  20.827  1.00 28.13      A    C
ATOM    569  C   TYR    75      21.877  44.324  20.468  1.00 31.45      A    C
ATOM    570  O   TYR    75      22.331  43.553  21.310  1.00 34.90      A    O
ATOM    571  CB  TYR    75      22.954  46.407  21.219  1.00 29.62      A    C
ATOM    572  CG  TYR    75      23.961  46.389  20.096  1.00 35.61      A    C
ATOM    573  CD1 TYR    75      23.759  47.148  18.949  1.00 36.53      A    C
ATOM    574  CD2 TYR    75      25.102  45.613  20.170  1.00 36.83      A    C
ATOM    575  CE1 TYR    75      24.669  47.150  17.921  1.00 35.23      A    C
ATOM    576  CE2 TYR    75      26.023  45.604  19.139  1.00 44.16      A    C
ATOM    577  CZ  TYR    75      25.801  46.379  18.015  1.00 43.29      A    C
ATOM    578  OH  TYR    75      26.709  46.377  16.975  1.00 52.94      A    O
ATOM    579  N   ASN    76      21.606  43.962  19.226  1.00 33.77      A    N
ATOM    580  CA  ASN    76      22.031  42.666  18.719  1.00 35.78      A    C
ATOM    581  C   ASN    76      21.481  41.532  19.576  1.00 32.53      A    C
ATOM    582  O   ASN    76      22.202  40.610  19.959  1.00 33.59      A    O
ATOM    583  CB  ASN    76      23.558  42.596  18.642  1.00 37.23      A    C
ATOM    584  CG  ASN    76      24.081  42.961  17.278  1.00 45.44      A    C
ATOM    585  OD1 ASN    76      24.351  44.126  16.999  1.00 46.64      A    O
ATOM    586  ND2 ASN    76      24.191  41.966  16.399  1.00 46.66      A    N
ATOM    587  N   GLY    77      20.191  41.620  19.897  1.00 32.36      A    N
ATOM    588  CA  GLY    77      19.529  40.601  20.686  1.00 33.39      A    C
ATOM    589  C   GLY    77      19.785  40.649  22.181  1.00 32.77      A    C
ATOM    590  O   GLY    77      19.179  39.865  22.918  1.00 34.11      A    O
END
"""

def test1(prefix="tst_mp4"):
  """
  Illustrating failure in twin refine tutorial data, discovered 2/26/2018.
  """
  with open("%s.pdb" % prefix, 'w') as f:
    f.write(pdb_txt)
  result = molprobity.run(args=['%s.pdb' % prefix])

if (__name__ == "__main__"):
  test1()
  print("OK")
