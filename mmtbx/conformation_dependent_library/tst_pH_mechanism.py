from __future__ import division
from __future__ import print_function
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

def run():
  f = file('tst_pH_gnp.pdb', 'wb')
  f.write(gnp)
  f.close()
  cmd = 'phenix.geometry_minimization tst_pH_gnp.pdb'
  rc = easy_run.go(cmd)
  find = ['Changed 28 bond restraint(s),  added 1 bond restraint(s)',
          'Changed 43 angle restraint(s), added 1 angle restraint(s)',
          ]
  for f in find:
    for line in rc.stdout_lines:
      if line.find(f)>-1:
        print(line)
        break
    else:
      assert 0, 'line not found: %s' % f
  return rc

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  rc = run(*tuple(args))
  assert rc.return_code==0
