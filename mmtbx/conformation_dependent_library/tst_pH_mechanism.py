from __future__ import absolute_import, division, print_function
import sys
from libtbx import easy_run

gnp = '''
HETATM 2435  PG  GNP A 201      -5.193  14.551 -21.840  1.00  9.31           P
HETATM 2436  O1G GNP A 201      -6.728  14.452 -21.462  1.00 12.25           O
HETATM 2437  O2G GNP A 201      -5.090  14.640 -23.230  1.00 18.00           O
HETATM 2438  O3G GNP A 201      -4.278  15.571 -21.150  1.00 15.01           O
HETATM 2439  N3B GNP A 201      -4.611  12.945 -21.512  1.00 13.32           N
HETATM 2440  PB  GNP A 201      -3.816  12.564 -20.004  1.00  5.89           P
HETATM 2441  O1B GNP A 201      -2.466  13.219 -20.064  1.00 12.14           O
HETATM 2442  O2B GNP A 201      -4.929  13.094 -19.079  1.00 11.92           O
HETATM 2443  O3A GNP A 201      -3.652  11.002 -20.112  1.00  6.88           O
HETATM 2444  PA  GNP A 201      -4.825   9.898 -19.654  1.00  9.13           P
HETATM 2445  O1A GNP A 201      -4.707   9.944 -18.092  1.00  9.97           O
HETATM 2446  O2A GNP A 201      -6.090  10.319 -20.321  1.00  6.24           O
HETATM 2447  O5' GNP A 201      -4.205   8.569 -20.187  1.00 11.53           O
HETATM 2448  C5' GNP A 201      -3.543   8.298 -21.351  1.00  6.70           C
HETATM 2449  C4' GNP A 201      -3.319   6.842 -21.597  1.00  9.27           C
HETATM 2450  O4' GNP A 201      -2.180   6.257 -20.945  1.00 11.15           O
HETATM 2451  C3' GNP A 201      -4.555   5.972 -21.280  1.00  5.98           C
HETATM 2452  O3' GNP A 201      -4.838   5.178 -22.402  1.00 20.95           O
HETATM 2453  C2' GNP A 201      -4.170   5.413 -19.954  1.00 12.11           C
HETATM 2454  O2' GNP A 201      -4.838   4.350 -19.264  1.00 16.85           O
HETATM 2455  C1' GNP A 201      -2.581   5.457 -19.941  1.00 10.49           C
HETATM 2456  N9  GNP A 201      -1.773   5.496 -18.743  1.00 14.77           N
HETATM 2457  C8  GNP A 201      -1.890   6.470 -17.880  1.00 11.82           C
HETATM 2458  N7  GNP A 201      -1.084   6.312 -16.856  1.00 14.77           N
HETATM 2459  C5  GNP A 201      -0.458   5.227 -17.093  1.00 13.86           C
HETATM 2460  C6  GNP A 201       0.560   4.449 -16.442  1.00 16.13           C
HETATM 2461  O6  GNP A 201       0.924   4.910 -15.409  1.00 13.44           O
HETATM 2462  N1  GNP A 201       1.054   3.385 -16.896  1.00 14.95           N
HETATM 2463  C2  GNP A 201       0.591   2.888 -18.089  1.00 12.85           C
HETATM 2464  N2  GNP A 201       1.104   1.752 -18.562  1.00 17.51           N
HETATM 2465  N3  GNP A 201      -0.378   3.514 -18.819  1.00 12.28           N
HETATM 2466  C4  GNP A 201      -0.911   4.657 -18.367  1.00 18.03           C
HETATM 2467 DOG2 GNP A 201      -4.220  14.350 -23.480  1.00 21.60           D
HETATM 2468 DNB3 GNP A 201      -4.638  12.314 -22.176  1.00 15.99           D
HETATM 2469 H5'2 GNP A 201      -4.064   8.664 -22.100  1.00  8.04           H
HETATM 2470 H5'1 GNP A 201      -2.670   8.749 -21.327  1.00  8.04           H
HETATM 2471  H4' GNP A 201      -3.158   6.751 -22.563  1.00 11.13           H
HETATM 2472  H3' GNP A 201      -5.321   6.570 -21.140  1.00  7.17           H
HETATM 2473 DO3' GNP A 201      -5.612   5.410 -22.779  1.00 25.14           D
HETATM 2474  H2' GNP A 201      -4.391   6.162 -19.360  1.00 14.54           H
HETATM 2475 DO2' GNP A 201      -5.490   4.672 -18.757  1.00 20.22           D
HETATM 2476  H1' GNP A 201      -2.361   4.562 -20.283  1.00 12.59           H
HETATM 2477  H8  GNP A 201      -2.520   7.220 -17.970  1.00 14.18           H
HETATM 2478  DN1 GNP A 201       1.698   2.933 -16.423  1.00 17.94           D
HETATM 2479 DN21 GNP A 201       1.566   1.190 -18.004  1.00 21.01           D
HETATM 2480 DN22 GNP A 201       0.995   1.532 -19.446  1.00 21.01           D
'''
llp_neutron = '''
REMARK CRYST1   55.532  124.986  130.444  90.00  90.00  90.00 P 21 21 21
HETATM    1  N   LLP B 258      -0.447 -10.941  14.499  1.00 15.19           N
HETATM    2  CA  LLP B 258      -0.166 -10.435  13.173  1.00 16.94           C
HETATM    3  C   LLP B 258       1.020  -9.463  13.270  1.00 16.29           C
HETATM    4  O   LLP B 258       1.880  -9.487  12.435  1.00 16.55           O
HETATM    5  CB  LLP B 258      -1.413  -9.787  12.490  1.00 19.72           C
HETATM    6  CG  LLP B 258      -2.474 -10.860  12.217  1.00 22.40           C
HETATM    7  CD  LLP B 258      -3.792 -10.723  13.002  1.00 23.55           C
HETATM    8  CE  LLP B 258      -4.807 -11.905  12.717  1.00 24.31           C
HETATM    9  NZ  LLP B 258      -6.088 -11.491  13.288  1.00 24.08           N
HETATM   10  C2  LLP B 258      -9.937  -9.685  14.136  1.00 18.63           C
HETATM   11  C2' LLP B 258     -10.748  -8.689  13.386  1.00 18.63           C
HETATM   12  C3  LLP B 258      -8.817 -10.412  13.596  1.00 18.26           C
HETATM   13  C4  LLP B 258      -8.134 -11.344  14.456  1.00 20.66           C
HETATM   14  C4' LLP B 258      -7.015 -12.088  13.961  1.00 24.62           C
HETATM   15  C5  LLP B 258      -8.637 -11.466  15.781  1.00 19.23           C
HETATM   16  C5' LLP B 258      -7.977 -12.415  16.725  1.00 19.27           C
HETATM   17  C6  LLP B 258      -9.708 -10.734  16.205  1.00 16.33           C
HETATM   18  D   LLP B 258      -0.978 -10.290  15.036  1.00 15.65           D
HETATM   19  DB2 LLP B 258      -1.099  -9.339  11.539  1.00 19.84           D
HETATM   20  DB3 LLP B 258      -1.817  -8.976  13.109  1.00 21.03           D
HETATM   21  DD2 LLP B 258      -3.584 -10.689  14.078  1.00 25.88           D
HETATM   22  DE2 LLP B 258      -4.436 -12.823  13.186  1.00 24.67           D
HETATM   23  DE3 LLP B 258      -4.912 -12.073  11.641  1.00 23.22           D
HETATM   24  DG2 LLP B 258      -2.076 -11.872  12.359  1.00 22.69           D
HETATM   25  DG3 LLP B 258      -2.733 -10.803  11.151  1.00 21.60           D
HETATM   26  N1  LLP B 258     -10.252  -9.936  15.386  1.00 16.38           N
HETATM   27  O3  LLP B 258      -8.483 -10.155  12.287  1.00 22.39           O
HETATM   28  OP1 LLP B 258      -5.476 -15.288  16.374  1.00 16.32           O
HETATM   29  OP2 LLP B 258      -7.046 -14.828  18.431  1.00 13.48           O
HETATM   30  OP3 LLP B 258      -4.890 -13.481  18.022  1.00 15.01           O
HETATM   31  OP4 LLP B 258      -6.764 -13.026  16.528  1.00 20.08           O
HETATM   32  P   LLP B 258      -6.052 -14.262  17.370  1.00 13.48           P
HETATM   33  H6  LLP B 258     -10.075 -10.811  17.219  1.00 17.86           H
HETATM   34 H2'3 LLP B 258     -11.798  -9.003  13.358  1.00 19.77           H
HETATM   35 H4'1 LLP B 258      -6.902 -13.121  14.273  1.00 23.11           H
HETATM   36 H5'1 LLP B 258      -8.706 -13.215  16.912  1.00 19.90           H
HETATM   37 H5'2 LLP B 258      -7.879 -11.891  17.686  1.00 18.49           H
HETATM   38  D1ZALLP B 258     -11.031  -9.415  15.735  0.94 16.15           D
HETATM   39  DA ALLP B 258       0.162 -11.264  12.527  0.90 17.45           D
HETATM   40  DD3ALLP B 258      -4.256  -9.766  12.729  0.67 24.01           D
HETATM   41 D2'1ALLP B 258     -10.434  -8.534  12.352  0.24 19.01           D
HETATM   42 D2'2ALLP B 258     -10.701  -7.714  13.885  0.40 19.18           D
HETATM   43  HA BLLP B 258       0.162 -11.264  12.527  0.10 17.45           H
HETATM   44  HD3BLLP B 258      -4.256  -9.766  12.729  0.33 24.01           H
HETATM   45  H1ZBLLP B 258     -11.031  -9.415  15.735  0.06 16.15           H
HETATM   46 H2'1BLLP B 258     -10.434  -8.534  12.352  0.76 19.01           H
HETATM   47 H2'2BLLP B 258     -10.701  -7.714  13.885  0.60 19.18           H
'''

def run1():
  f = open('tst_pH_gnp.pdb', 'w')
  f.write(gnp)
  f.close()
  cmd = 'phenix.geometry_minimization tst_pH_gnp.pdb'
  print(cmd)
  rc = easy_run.go(cmd)
  find = ['Changed (significantly) 24 bond restraint(s),  added 1 bond restraint(s)',
          'Changed (significantly) 35 angle restraint(s), added 1 angle restraint(s)',
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
  f = open('tst_LLP_neutron.pdb', 'w')
  f.write(llp_neutron)
  f.close()
  cmd = 'phenix.geometry_minimization tst_LLP_neutron.pdb'
  print(cmd)
  rc = easy_run.go(cmd)
  find = ['Changed (significantly) 29 bond restraint(s),  added 2 bond restraint(s)',
          'Changed (significantly) 13 angle restraint(s), added 4 angle restraint(s)',
          ]
  for f in find:
    for line in rc.stdout_lines:
      if line.find(f)>-1:
        # print(line)
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
  assert rc.return_code==0
