from __future__ import absolute_import, division, print_function
import sys
from libtbx import easy_run

pdb_6mu9_enol = '''
CRYST1   35.307   42.094   48.542  97.73 105.52 113.34 P 1
SCALE1      0.028323  0.012221  0.011626        0.00000
SCALE2      0.000000  0.025874  0.007309        0.00000
SCALE3      0.000000  0.000000  0.022217        0.00000
ATOM      8  N   ALA A 295      -7.890  -8.768  10.104  1.00  4.63           N
ATOM      9  CA  ALA A 295      -8.605  -8.087  11.181  1.00  4.72           C
ATOM     10  C   ALA A 295      -8.140  -8.599  12.543  1.00  4.93           C
ATOM     11  O   ALA A 295      -7.881  -7.797  13.461  1.00  5.29           O
ATOM     12  CB  ALA A 295     -10.089  -8.264  11.003  1.00  5.17           C
ATOM     13  H   ALA A 295      -8.397  -9.173   9.540  1.00  4.64           H
ATOM     14  HA  ALA A 295      -8.410  -7.138  11.137  1.00  4.72           H
ATOM     15  HB1 ALA A 295     -10.550  -7.807  11.724  1.00  5.17           H
ATOM     16  HB2 ALA A 295     -10.352  -7.885  10.150  1.00  5.17           H
ATOM     17  HB3 ALA A 295     -10.298  -9.211  11.023  1.00  5.17           H
ATOM     18  N   GLU A 296      -8.040  -9.914  12.713  1.00  5.04           N
ATOM     19  C   GLU A 296      -6.182 -10.125  14.321  1.00  4.85           C
ATOM     20  O   GLU A 296      -5.839  -9.875  15.483  1.00  5.72           O
ATOM     21  CA AGLU A 296      -7.647 -10.459  14.018  0.60  5.09           C
ATOM     22  CB AGLU A 296      -7.875 -11.968  14.096  0.60  5.76           C
ATOM     23  CG AGLU A 296      -9.328 -12.413  14.130  0.60  7.89           C
ATOM     24  CD AGLU A 296      -9.471 -13.771  14.845  0.60 10.09           C
ATOM     25  OE1AGLU A 296      -9.314 -13.813  16.066  0.60 11.70           O
ATOM     26  OE2AGLU A 296      -9.703 -14.799  14.216  0.60 10.03           O
ATOM     27  H  AGLU A 296      -8.181 -10.501  12.107  0.60  5.04           H
ATOM     28  HA AGLU A 296      -8.192 -10.043  14.704  0.60  5.17           H
ATOM     29  HB2AGLU A 296      -7.463 -12.380  13.321  0.60  5.90           H
ATOM     30  HB3AGLU A 296      -7.451 -12.301  14.902  0.60  5.90           H
ATOM     31  HG2AGLU A 296      -9.855 -11.757  14.612  0.60  7.94           H
ATOM     32  HG3AGLU A 296      -9.657 -12.509  13.223  0.60  7.94           H
ATOM     33  CA BGLU A 296      -7.625 -10.460  14.013  0.40  5.17           C
ATOM     34  CB BGLU A 296      -7.851 -11.964  14.083  0.40  5.89           C
ATOM     35  CG BGLU A 296      -9.332 -12.273  14.088  0.40  7.94           C
ATOM     36  CD BGLU A 296      -9.643 -13.744  14.210  0.40 10.23           C
ATOM     37  OE1BGLU A 296      -8.894 -14.454  14.914  0.40 13.79           O
ATOM     38  OE2BGLU A 296     -10.651 -14.166  13.605  0.40 14.65           O
ATOM     39  H  BGLU A 296      -8.194 -10.501  12.110  0.40  5.04           H
ATOM     40  HA BGLU A 296      -8.171 -10.052  14.704  0.40  5.17           H
ATOM     41  HB2BGLU A 296      -7.451 -12.388  13.308  0.40  5.90           H
ATOM     42  HB3BGLU A 296      -7.461 -12.312  14.900  0.40  5.90           H
ATOM     43  HG2BGLU A 296      -9.744 -11.819  14.840  0.40  7.94           H
ATOM     44  HG3BGLU A 296      -9.721 -11.955  13.259  0.40  7.94           H
ATOM     45  N   ALA A 298      -4.628  -7.432  13.375  1.00  4.46           N
ATOM     46  CA  ALA A 298      -4.642  -6.018  13.758  1.00  4.55           C
ATOM     47  C   ALA A 298      -5.093  -5.864  15.208  1.00  5.06           C
ATOM     48  O   ALA A 298      -4.498  -5.097  15.981  1.00  5.62           O
ATOM     49  CB  ALA A 298      -5.539  -5.231  12.825  1.00  4.75           C
ATOM     50  H   ALA A 298      -5.188  -7.633  12.758  1.00  4.46           H
ATOM     51  HA  ALA A 298      -3.743  -5.661  13.683  1.00  4.56           H
ATOM     52  HB1 ALA A 298      -5.535  -4.299  13.095  1.00  4.75           H
ATOM     53  HB2 ALA A 298      -5.203  -5.315  11.919  1.00  4.75           H
ATOM     54  HB3 ALA A 298      -6.440  -5.587  12.879  1.00  4.75           H
ATOM     55  N   LYS A 299      -6.151  -6.592  15.571  1.00  5.45           N
ATOM     56  C   LYS A 299      -5.553  -6.995  17.929  1.00  5.57           C
ATOM     57  O   LYS A 299      -5.364  -6.373  18.988  1.00  6.24           O
ATOM     58  CA ALYS A 299      -6.650  -6.563  16.943  0.70  5.72           C
ATOM     59  CB ALYS A 299      -7.887  -7.435  17.038  0.70  6.52           C
ATOM     60  CG ALYS A 299      -8.615  -7.361  18.360  0.70  8.93           C
ATOM     61  CD ALYS A 299      -9.836  -8.261  18.356  0.70 13.49           C
ATOM     62  CE ALYS A 299     -10.613  -8.205  19.658  0.70 19.14           C
ATOM     63  NZ ALYS A 299     -11.844  -9.041  19.615  0.70 25.08           N
ATOM    105  HNO LYS A 299      -5.898  -6.487  18.858  1.00  5.37           H
ATOM     64  H  ALYS A 299      -6.591  -7.109  15.046  0.70  5.45           H
ATOM     65  HA ALYS A 299      -6.908  -5.655  17.164  0.70  5.77           H
ATOM     66  HB2ALYS A 299      -8.510  -7.165  16.345  0.70  6.80           H
ATOM     67  HB3ALYS A 299      -7.626  -8.359  16.899  0.70  6.80           H
ATOM     68  HG2ALYS A 299      -8.023  -7.652  19.071  0.70  8.92           H
ATOM     69  HG3ALYS A 299      -8.907  -6.449  18.516  0.70  8.92           H
ATOM     70  HD2ALYS A 299     -10.430  -7.984  17.641  0.70 12.97           H
ATOM     71  HD3ALYS A 299      -9.552  -9.178  18.216  0.70 12.97           H
ATOM     72  HE2ALYS A 299     -10.051  -8.530  20.379  0.70 17.52           H
ATOM     73  HE3ALYS A 299     -10.876  -7.289  19.827  0.70 17.52           H
ATOM     74  HZ1ALYS A 299     -12.445  -8.668  19.075  0.70 21.82           H
ATOM     75  HZ2ALYS A 299     -11.647  -9.853  19.309  0.70 21.82           H
ATOM     76  HZ3ALYS A 299     -12.187  -9.116  20.432  0.70 21.82           H
ATOM     77  CA BLYS A 299      -6.652  -6.577  16.942  0.30  5.76           C
ATOM     78  CB BLYS A 299      -7.867  -7.492  17.043  0.30  6.79           C
ATOM     79  CG BLYS A 299      -8.651  -7.382  18.335  0.30  8.92           C
ATOM     80  CD BLYS A 299      -9.850  -8.314  18.304  0.30 12.96           C
ATOM     81  CE BLYS A 299     -10.624  -8.309  19.611  0.30 17.52           C
ATOM     82  NZ BLYS A 299     -11.225  -6.976  19.889  0.30 21.81           N
ATOM     83  H  BLYS A 299      -6.592  -7.105  15.043  0.30  5.45           H
ATOM     84  HA BLYS A 299      -6.933  -5.676  17.168  0.30  5.77           H
ATOM     85  HB2BLYS A 299      -8.473  -7.281  16.315  0.30  6.80           H
ATOM     86  HB3BLYS A 299      -7.568  -8.411  16.960  0.30  6.80           H
ATOM     87  HG2BLYS A 299      -8.084  -7.634  19.081  0.30  8.92           H
ATOM     88  HG3BLYS A 299      -8.970  -6.472  18.443  0.30  8.92           H
ATOM     89  HD2BLYS A 299     -10.452  -8.033  17.597  0.30 12.97           H
ATOM     90  HD3BLYS A 299      -9.543  -9.219  18.139  0.30 12.97           H
ATOM     91  HE2BLYS A 299     -11.341  -8.960  19.560  0.30 17.52           H
ATOM     92  HE3BLYS A 299     -10.022  -8.528  20.340  0.30 17.52           H
ATOM     93  HZ1BLYS A 299     -11.916  -6.835  19.346  0.30 21.82           H
ATOM     94  HZ2BLYS A 299     -11.509  -6.941  20.732  0.30 21.82           H
ATOM     95  HZ3BLYS A 299     -10.621  -6.336  19.760  0.30 21.82           H
ATOM     96  N   GLN A 300      -4.803  -8.034  17.581  1.00  5.36           N
ATOM     97  CA  GLN A 300      -3.714  -8.488  18.441  1.00  5.61           C
ATOM     98  C   GLN A 300      -2.580  -7.469  18.504  1.00  5.41           C
ATOM     99  O   GLN A 300      -2.009  -7.269  19.573  1.00  6.27           O
ATOM    100  CB  GLN A 300      -3.191  -9.851  17.989  1.00  5.98           C
ATOM    101  CG  GLN A 300      -2.052 -10.353  18.864  1.00  6.43           C
ATOM    102  CD  GLN A 300      -2.484 -10.561  20.307  1.00  6.45           C
ATOM    103  OE1 GLN A 300      -1.877 -10.015  21.224  1.00  7.66           O
ATOM    104  NE2 GLN A 300      -3.535 -11.343  20.514  1.00  7.78           N
ATOM    106  HA  GLN A 300      -4.061  -8.591  19.341  1.00  5.62           H
ATOM    107  HB2 GLN A 300      -3.912 -10.498  18.031  1.00  5.99           H
ATOM    108  HB3 GLN A 300      -2.862  -9.778  17.079  1.00  5.99           H
ATOM    109  HG2 GLN A 300      -1.748 -11.207  18.520  1.00  6.43           H
ATOM    110  HG3 GLN A 300      -1.321  -9.716  18.853  1.00  6.43           H
ATOM    111 HE21 GLN A 300      -3.935 -11.701  19.854  1.00  7.79           H
ATOM    112 HE22 GLN A 300      -3.808 -11.479  21.314  1.00  7.79           H
ATOM    113  N   VAL A 302      -2.733  -4.226  18.245  1.00  6.36           N
ATOM    114  CA  VAL A 302      -3.196  -3.172  19.170  1.00  7.64           C
ATOM    115  C   VAL A 302      -3.117  -3.633  20.625  1.00  6.92           C
ATOM    116  O   VAL A 302      -2.686  -2.885  21.506  1.00  7.89           O
ATOM    117  CB  VAL A 302      -4.658  -2.770  18.820  1.00  9.66           C
ATOM    118  CG1 VAL A 302      -5.326  -1.925  19.899  1.00  9.76           C
ATOM    119  CG2 VAL A 302      -4.642  -2.018  17.523  1.00 10.64           C
ATOM    120  H   VAL A 302      -3.323  -4.479  17.676  1.00  6.36           H
ATOM    121  HA  VAL A 302      -2.633  -2.389  19.069  1.00  7.65           H
ATOM    122  HB  VAL A 302      -5.187  -3.573  18.696  1.00  9.67           H
ATOM    123 HG11 VAL A 302      -6.125  -1.518  19.529  1.00  9.77           H
ATOM    124 HG12 VAL A 302      -5.563  -2.494  20.648  1.00  9.77           H
ATOM    125 HG13 VAL A 302      -4.708  -1.235  20.186  1.00  9.77           H
ATOM    126 HG21 VAL A 302      -5.551  -1.772  17.290  1.00 10.64           H
ATOM    127 HG22 VAL A 302      -4.100  -1.221  17.627  1.00 10.64           H
ATOM    128 HG23 VAL A 302      -4.266  -2.587  16.833  1.00 10.64           H
TER
HETATM  151  O   HOH A 612     -14.077  -8.095  20.396  1.00 30.97           O
ANISOU  151  O   HOH A 612     2034   4356   5379    389    733    650       O
HETATM  152  O   HOH A 614      -9.477 -11.830  17.688  1.00 21.91           O
ANISOU  152  O   HOH A 614     2233   3500   2591   -829    560  -1423       O
HETATM  153  O   HOH A 624      -9.375  -6.307  21.625  1.00 16.97           O
ANISOU  153  O   HOH A 624     1632   3441   1376   -252    511   -189       O
HETATM  154  O   HOH A 758       0.566   3.963 -21.726  1.00 17.32           O
ANISOU  154  O   HOH A 758     1490   2309   2782   -162   -208   -138       O
HETATM  155  O   HOH A 767      -7.259  -4.934  20.452  1.00  8.81           O
ANISOU  155  O   HOH A 767     1113   1517    719     52    126   -154       O
HETATM  156  O   HOH A 772      -6.854 -11.129  17.769  1.00  9.63           O
ANISOU  156  O   HOH A 772     1307   1551    802   -418     77    261       O
HETATM  157  O   HOH A 896     -11.730 -10.705  16.670  1.00 31.68           O
ANISOU  157  O   HOH A 896     2876   5566   3594   -409    181  -1771       O
HETATM  158  O  AHOH A 926     -11.872  -5.521  20.961  0.60 13.70           O
ANISOU  158  O  AHOH A 926     1544   2160   1503    302    417    263       O
END
'''

h_atom = 'ATOM   4292  H   GLN A 300      -4.898  -8.487  16.858  1.00  5.37           H'

def run1():
  pf='tst_pdb_6mu9_enol.pdb'
  f = open(pf, 'w')
  f.write(pdb_6mu9_enol)
  f.close()
  cmd = 'phenix.geometry_minimization %s' % pf
  print(cmd)
  rc = easy_run.go(cmd)
  find = ['Changed (significantly) 1 bond restraint(s),  added 1 bond restraint(s)',
          'Changed (significantly) 0 angle/dihedral restraint(s), added 3 angle/dihedral restraint(s)',
          'New atom/bond "ATOM     51  HNO LYS A 299 .*.     H  "',
          ]
  for f in find:
    for line in rc.stdout_lines:
      if line.find(f)>-1:
        # print(line)
        break
    else:
      assert 0, 'line not found: %s' % f
  return rc

def run2():
  pf='tst_pdb_6mu9_enol_error.pdb'
  pdb_6mu9_enol_error = pdb_6mu9_enol.replace(
    'ATOM    106  HA  GLN A 300      -4.061  -8.591  19.341  1.00  5.62           H',
    'ATOM    106  HA  GLN A 300      -4.061  -8.591  19.341  1.00  5.62           H\n%s' % h_atom,
    )
  # print(pdb_6mu9_enol_error)
  f = open(pf, 'w')
  f.write(pdb_6mu9_enol_error)
  f.close()
  cmd = 'phenix.geometry_minimization %s' % pf
  print(cmd)
  rc = easy_run.go(cmd)
  find = ['Sorry: Enol-peptide should not have a "H" hydrogen on following residue :']
  for f in find:
    for line in rc.stdout_lines:
      if line.find(f)>-1:
        break
    else:
      assert 0, 'line not found: %s' % f
  return rc

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  rc = run1()
  assert rc.return_code==0
  rc = run2()
  assert rc.return_code!=0
