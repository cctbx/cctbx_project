from __future__ import absolute_import, division, print_function

pdb_string = '''
CRYST1   74.120  110.220   88.210  90.00 111.08  90.00 P 1 21 1
SCALE1      0.013492  0.000000  0.005201        0.00000
SCALE2      0.000000  0.009073  0.000000        0.00000
SCALE3      0.000000  0.000000  0.012150        0.00000
ATOM      1  N   MET A  18      15.215  -9.367  61.781  1.00  4.60           N
ATOM      2  CA  MET A  18      15.994  -9.639  60.573  1.00  4.33           C
ATOM      3  C   MET A  18      16.745  -8.406  60.090  1.00  4.51           C
ATOM      4  O   MET A  18      16.143  -7.509  59.515  1.00  4.82           O
ATOM      5  CB  MET A  18      15.126 -10.183  59.422  1.00  4.76           C
ATOM      6  CG  MET A  18      14.343 -11.416  59.813  1.00  5.58           C
ATOM      7  SD  MET A  18      13.572 -12.228  58.401  1.00  6.39           S
ATOM      8  CE  MET A  18      14.906 -13.295  57.886  1.00  6.54           C
ATOM      9  N   PRO A  19      18.064  -8.360  60.319  1.00  5.05           N
ATOM     10  CA  PRO A  19      18.906  -7.408  59.619  1.00  5.21           C
ATOM     11  C   PRO A  19      19.235  -7.982  58.236  1.00  5.15           C
ATOM     12  O   PRO A  19      19.074  -9.179  58.005  1.00  5.57           O
ATOM     13  CB  PRO A  19      20.167  -7.364  60.481  1.00  5.61           C
ATOM     14  CG  PRO A  19      20.273  -8.777  60.996  1.00  5.96           C
ATOM     15  CD  PRO A  19      18.864  -9.259  61.175  1.00  5.46           C
ATOM     16  N   GLY A  20      19.718  -7.162  57.313  1.00  5.38           N
ATOM     17  CA  GLY A  20      20.123  -7.718  56.021  1.00  5.56           C
ATOM     18  C   GLY A  20      21.352  -8.618  56.134  1.00  5.51           C
ATOM     19  O   GLY A  20      21.447  -9.635  55.447  1.00  6.08           O
ATOM     20  N   MET A  45      23.492 -17.747  55.646  1.00  5.92           N
ATOM     21  CA  MET A  45      22.469 -16.820  56.107  1.00  5.49           C
ATOM     22  C   MET A  45      22.107 -17.166  57.544  1.00  5.48           C
ATOM     23  O   MET A  45      21.090 -17.797  57.833  1.00  6.28           O
ATOM     24  CB  MET A  45      21.247 -16.804  55.194  1.00  5.38           C
ATOM     25  CG  MET A  45      20.401 -15.567  55.446  1.00  6.08           C
ATOM     26  SD  MET A  45      21.221 -13.984  55.142  1.00  6.66           S
ATOM     27  CE  MET A  45      21.440 -13.998  53.389  1.00 10.50           C
ATOM     28  N  AMET A  60      17.784 -13.489  62.293  0.50  4.96           N
ATOM     29  CA AMET A  60      18.366 -14.113  61.111  0.50  4.63           C
ATOM     30  C  AMET A  60      18.517 -13.083  60.009  0.50  4.62           C
ATOM     31  O  AMET A  60      17.645 -12.235  59.816  0.50  5.04           O
ATOM     32  CB AMET A  60      17.486 -15.268  60.628  0.50  5.40           C
ATOM     33  CG AMET A  60      18.030 -15.953  59.379  0.50  5.79           C
ATOM     34  SD AMET A  60      17.293 -17.563  59.033  0.50  8.88           S
ATOM     35  CE AMET A  60      15.691 -17.103  58.425  0.50  9.84           C
ATOM     36  N  BMET A  60      17.808 -13.468  62.300  0.50  4.99           N
ATOM     37  CA BMET A  60      18.372 -14.103  61.122  0.50  4.66           C
ATOM     38  C  BMET A  60      18.520 -13.079  60.013  0.50  4.61           C
ATOM     39  O  BMET A  60      17.651 -12.228  59.822  0.50  5.05           O
ATOM     40  CB BMET A  60      17.463 -15.242  60.671  0.50  5.34           C
ATOM     41  CG BMET A  60      17.993 -16.000  59.469  0.50  5.73           C
ATOM     42  SD BMET A  60      16.961 -17.415  59.051  0.50  7.74           S
ATOM     43  CE BMET A  60      17.480 -17.628  57.356  0.50  6.86           C
ATOM     44  N   GLY A  61      19.623 -13.141  59.279  1.00  4.61           N
ATOM     45  CA  GLY A  61      19.790 -12.275  58.124  1.00  4.53           C
ATOM     46  C   GLY A  61      18.763 -12.558  57.033  1.00  4.33           C
ATOM     47  O   GLY A  61      18.268 -13.677  56.885  1.00  4.88           O
ATOM     48  N   HIS A  62      18.445 -11.511  56.263  1.00  4.67           N
ATOM     49  CA  HIS A  62      17.535 -11.671  55.130  1.00  4.67           C
ATOM     50  C   HIS A  62      18.172 -11.460  53.753  1.00  4.62           C
ATOM     51  O   HIS A  62      17.509 -11.713  52.751  1.00  5.66           O
ATOM     52  CB  HIS A  62      16.245 -10.808  55.263  1.00  5.00           C
ATOM     53  CG  HIS A  62      16.453  -9.324  55.365  1.00  5.03           C
ATOM     54  ND1 HIS A  62      17.096  -8.584  54.399  1.00  5.27           N
ATOM     55  CD2 HIS A  62      15.999  -8.429  56.281  1.00  5.73           C
ATOM     56  CE1 HIS A  62      17.081  -7.304  54.759  1.00  5.40           C
ATOM     57  NE2 HIS A  62      16.409  -7.182  55.884  1.00  5.91           N
ATOM     58  N   GLY A  63      19.414 -10.999  53.682  1.00  4.54           N
ATOM     59  CA  GLY A  63      20.077 -10.720  52.424  1.00  5.11           C
ATOM     60  C   GLY A  63      19.597  -9.432  51.792  1.00  5.28           C
ATOM     61  O   GLY A  63      18.676  -8.784  52.268  1.00 10.32           O
ATOM     62  N   MET A  64      20.188  -9.032  50.690  1.00  4.88           N
ATOM     63  CA  MET A  64      19.875  -7.722  50.110  1.00  4.61           C
ATOM     64  C   MET A  64      18.802  -7.851  49.034  1.00  4.01           C
ATOM     65  O   MET A  64      18.950  -8.608  48.072  1.00  5.02           O
ATOM     66  CB  MET A  64      21.123  -7.109  49.485  1.00  5.88           C
ATOM     67  CG  MET A  64      22.189  -6.837  50.547  1.00  7.27           C
ATOM     68  SD  MET A  64      23.536  -5.778  50.026  1.00 10.02           S
ATOM     69  CE  MET A  64      24.268  -6.809  48.769  1.00  9.52           C
ATOM     70  N   SER A  68      15.441 -11.601  47.812  1.00  3.72           N
ATOM     71  CA  SER A  68      16.098 -12.378  48.866  1.00  4.05           C
ATOM     72  C   SER A  68      15.232 -12.433  50.117  1.00  4.18           C
ATOM     73  O   SER A  68      14.960 -13.517  50.652  1.00  4.38           O
ATOM     74  CB  SER A  68      17.456 -11.763  49.191  1.00  4.18           C
ATOM     75  OG  SER A  68      18.137 -12.530  50.159  1.00  5.33           O
ATOM     76  N   CYS A  69      14.801 -11.268  50.599  1.00  4.15           N
ATOM     77  CA  CYS A  69      14.051 -11.277  51.835  1.00  4.21           C
ATOM     78  C   CYS A  69      12.676 -11.926  51.689  1.00  3.95           C
ATOM     79  O   CYS A  69      12.147 -12.460  52.674  1.00  5.01           O
ATOM     80  CB  CYS A  69      13.942  -9.887  52.432  1.00  4.97           C
ATOM     81  SG  CYS A  69      12.884  -8.739  51.550  1.00  5.68           S
ATOM     82  N   TYR A  72      13.108 -15.718  52.476  1.00  4.39           N
ATOM     83  CA  TYR A  72      13.265 -16.077  53.880  1.00  4.87           C
ATOM     84  C   TYR A  72      12.019 -15.828  54.678  1.00  4.47           C
ATOM     85  O   TYR A  72      11.644 -16.649  55.510  1.00  5.33           O
ATOM     86  CB  TYR A  72      14.470 -15.353  54.494  1.00  5.22           C
ATOM     87  CG  TYR A  72      15.777 -15.839  53.953  1.00  5.47           C
ATOM     88  CD1 TYR A  72      16.251 -17.119  54.274  1.00  6.76           C
ATOM     89  CD2 TYR A  72      16.537 -15.067  53.084  1.00  6.32           C
ATOM     90  CE1 TYR A  72      17.428 -17.590  53.743  1.00  7.24           C
ATOM     91  CE2 TYR A  72      17.727 -15.532  52.569  1.00  7.02           C
ATOM     92  CZ  TYR A  72      18.173 -16.793  52.898  1.00  6.65           C
ATOM     93  OH  TYR A  72      19.366 -17.214  52.363  1.00  8.74           O
ATOM     94  N   ARG A  87      12.058  -6.379  60.531  1.00  4.03           N
ATOM     95  CA  ARG A  87      13.280  -6.007  59.840  1.00  4.44           C
ATOM     96  C   ARG A  87      13.927  -4.837  60.581  1.00  4.37           C
ATOM     97  O   ARG A  87      13.242  -3.890  60.978  1.00  4.52           O
ATOM     98  CB  ARG A  87      13.005  -5.622  58.385  1.00  4.97           C
ATOM     99  CG  ARG A  87      14.195  -4.965  57.708  1.00  5.40           C
ATOM    100  CD  ARG A  87      13.976  -4.770  56.203  1.00  5.93           C
ATOM    101  NE  ARG A  87      14.676  -3.589  55.713  1.00  6.34           N
ATOM    102  CZ  ARG A  87      15.987  -3.377  55.798  1.00  5.76           C
ATOM    103  NH1 ARG A  87      16.787  -4.315  56.242  1.00  6.21           N
ATOM    104  NH2 ARG A  87      16.487  -2.207  55.442  1.00  6.66           N
ATOM    105  N   GLU A 181      13.993  -3.488  51.085  1.00  4.28           N
ATOM    106  CA  GLU A 181      13.274  -4.704  51.468  1.00  4.20           C
ATOM    107  C   GLU A 181      11.842  -4.541  51.980  1.00  3.76           C
ATOM    108  O   GLU A 181      11.044  -5.466  51.776  1.00  4.13           O
ATOM    109  CB  GLU A 181      14.066  -5.476  52.519  1.00  4.85           C
ATOM    110  CG  GLU A 181      15.375  -6.034  51.998  1.00  5.84           C
ATOM    111  CD  GLU A 181      16.565  -5.108  52.105  1.00  5.87           C
ATOM    112  OE1 GLU A 181      17.620  -5.461  51.537  1.00  6.33           O
ATOM    113  OE2 GLU A 181      16.481  -4.053  52.751  1.00  7.97           O
ATOM    114  N   ILE A 185       8.764  -7.518  52.166  1.00  4.11           N
ATOM    115  CA  ILE A 185       8.380  -8.207  53.399  1.00  4.15           C
ATOM    116  C   ILE A 185       6.887  -7.955  53.683  1.00  3.80           C
ATOM    117  O   ILE A 185       6.158  -8.888  54.046  1.00  4.10           O
ATOM    118  CB  ILE A 185       9.244  -7.743  54.571  1.00  4.79           C
ATOM    119  CG1 ILE A 185      10.711  -8.123  54.320  1.00  5.83           C
ATOM    120  CG2 ILE A 185       8.742  -8.362  55.870  1.00  5.63           C
ATOM    121  CD1 ILE A 185      11.673  -7.661  55.373  1.00  6.69           C
TER
HETATM  122  C1 AGOL A 301      19.343  -2.253  53.379  0.25 11.18           C
HETATM  123  C2 AGOL A 301      19.099  -3.727  53.183  0.25 10.86           C
HETATM  124  C3 AGOL A 301      19.755  -4.689  54.142  0.25 11.90           C
HETATM  125  O1 AGOL A 301      18.265  -1.622  52.704  0.25  8.81           O
HETATM  126  O2 AGOL A 301      19.693  -4.136  51.959  0.25  8.60           O
HETATM  127  O3 AGOL A 301      20.239  -5.795  53.373  0.25  7.42           O
HETATM  128  C1 BGOL A 301      21.418  -2.460  54.579  0.25 14.24           C
HETATM  129  C2 BGOL A 301      21.388  -3.327  55.820  0.25 13.76           C
HETATM  130  C3 BGOL A 301      19.940  -3.430  56.209  0.25  9.95           C
HETATM  131  O1 BGOL A 301      21.378  -3.406  53.502  0.25 10.70           O
HETATM  132  O2 BGOL A 301      21.860  -4.628  55.469  0.25 11.34           O
HETATM  133  O3 BGOL A 301      19.286  -2.256  55.666  0.25  5.65           O
HETATM  134  C1 CGOL A 301      22.179  -2.132  55.772  0.50 11.93           C
HETATM  135  C2 CGOL A 301      21.418  -3.361  55.236  0.50 11.64           C
HETATM  136  C3 CGOL A 301      20.007  -3.217  55.753  0.50 12.10           C
HETATM  137  O1 CGOL A 301      23.559  -2.106  55.373  0.50 16.86           O
HETATM  138  O2 CGOL A 301      21.976  -4.600  55.688  0.50 12.26           O
HETATM  139  O3 CGOL A 301      19.681  -4.230  56.700  0.50  8.36           O
HETATM  140  O   HOH A 422      20.632 -13.536  49.381  1.00 16.15           O
END
'''

import time
from libtbx import easy_run

def main(run_mopac=None):
  try:
    run_mopac=int(run_mopac)
  except Exception:
    run_mopac=True
  preamble = 'tst_flipping_his'
  selection = 'chain A and resid 62 and resname HIS'
  f=open('%s.pdb' % preamble, 'w')
  f.write(pdb_string)
  del f
  cmd = 'phenix.ready_set %(preamble)s.pdb' % locals()
  print('\n ~> %s' % cmd)
  rc = easy_run.go(cmd)
  assert rc.return_code==0
  cmd = 'mmtbx.quantum_interface %(preamble)s.pdb' % locals()
  print('\n ~> %s' % cmd)
  rc = easy_run.go(cmd)
  assert rc.return_code
  cmd = 'mmtbx.quantum_interface %(preamble)s.updated.pdb write_qmr_phil=1 qi.selection="%(selection)s"' % locals()
  cmd += ' format=quantum_interface'
  print('\n ~> %s' % cmd)
  rc = easy_run.go(cmd)
  assert rc.return_code==0
  t0=time.time()
  cmd = 'mmtbx.quantum_interface %(preamble)s.updated.pdb run_qmr=True %(preamble)s.updated_A_62_HIS.phil' % locals()
  cmd += ' iterate_histidine=True only_i=1'
  if type(run_mopac)==type(1):
    cmd += ' qi.nproc=%d' % run_mopac
  print('\n ~> %s' % cmd)
  if run_mopac:
    rc = easy_run.go(cmd)
    print('Time : %0.fs' % (time.time()-t0))
    assert rc.return_code==0
  print('DONE')

if __name__ == '__main__':
  import sys
  main(*tuple(sys.argv[1:]))
