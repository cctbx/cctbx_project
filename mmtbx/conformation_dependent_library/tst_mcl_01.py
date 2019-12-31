from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run

pdb_string = '''
CRYST1  230.463  123.139  153.201  90.00 131.05  90.00 C 1 2 1
SCALE1      0.004339  0.000000  0.003779        0.00000
SCALE2      0.000000  0.008121  0.000000        0.00000
SCALE3      0.000000  0.000000  0.008655        0.00000
ATOM      1  N   CYS A 861     -12.347 -25.689  76.808  1.00127.20           N
ATOM      2  CA  CYS A 861     -13.104 -25.567  78.048  1.00133.90           C
ATOM      3  C   CYS A 861     -14.052 -26.748  78.202  1.00134.10           C
ATOM      4  O   CYS A 861     -14.691 -27.168  77.239  1.00131.28           O
ATOM      5  CB  CYS A 861     -13.889 -24.252  78.073  1.00135.25           C
ATOM      6  SG  CYS A 861     -14.037 -23.445  79.699  1.00155.43           S
ATOM      7  N   CYS A 864     -17.542 -24.913  77.978  1.00139.78           N
ATOM      8  CA  CYS A 864     -17.949 -24.175  76.780  1.00132.76           C
ATOM      9  C   CYS A 864     -17.994 -25.041  75.522  1.00125.77           C
ATOM     10  O   CYS A 864     -18.762 -24.765  74.604  1.00125.97           O
ATOM     11  CB  CYS A 864     -17.017 -22.980  76.533  1.00134.20           C
ATOM     12  SG  CYS A 864     -16.835 -21.826  77.907  1.00126.17           S
ATOM     13  N   CYS A1054     -20.150 -20.841  82.906  1.00100.48           N
ATOM     14  CA  CYS A1054     -18.688 -20.783  82.915  1.00101.74           C
ATOM     15  C   CYS A1054     -18.182 -19.501  83.572  1.00103.40           C
ATOM     16  O   CYS A1054     -18.656 -18.409  83.261  1.00107.65           O
ATOM     17  CB  CYS A1054     -18.128 -20.893  81.493  1.00111.10           C
ATOM     18  SG  CYS A1054     -16.308 -20.868  81.398  1.00 99.22           S
ATOM     19  N   HIS A1060     -13.988 -16.716  77.115  1.00102.85           N
ATOM     20  CA  HIS A1060     -14.171 -17.807  76.158  1.00107.46           C
ATOM     21  C   HIS A1060     -14.871 -17.324  74.890  1.00101.11           C
ATOM     22  O   HIS A1060     -14.674 -17.888  73.814  1.00 96.36           O
ATOM     23  CB  HIS A1060     -14.975 -18.973  76.754  1.00122.04           C
ATOM     24  CG  HIS A1060     -14.199 -19.831  77.706  1.00131.75           C
ATOM     25  ND1 HIS A1060     -14.110 -19.567  79.055  1.00134.67           N
ATOM     26  CD2 HIS A1060     -13.424 -20.921  77.487  1.00136.31           C
ATOM     27  CE1 HIS A1060     -13.342 -20.474  79.633  1.00135.85           C
ATOM     28  NE2 HIS A1060     -12.911 -21.303  78.702  1.00135.21           N
TER
HETATM   29 ZN    ZN A2083     -15.163 -21.423  79.459  1.00340.34          Zn
END
'''

def main():
  with open('tst_mcl_01.pdb', 'w') as f:
    f.write(pdb_string)
  cmd = 'phenix.pdb_interpretation tst_mcl_01.pdb write_geo=1'
  print (cmd)
  rc = easy_run.go(cmd)
  assert os.path.exists('tst_mcl_01.pdb.geo')
  return rc.return_code

if __name__ == '__main__':
  rc = main()
  assert not rc
  print('OK')
