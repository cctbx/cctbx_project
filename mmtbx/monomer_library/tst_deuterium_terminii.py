from __future__ import division
pdbs = {'one' : ['''
CRYST1   71.765   48.519   56.920  90.00  90.00  90.00 P 1
SCALE1      0.013934  0.000000  0.000000        0.00000
SCALE2      0.000000  0.020610  0.000000        0.00000
SCALE3      0.000000  0.000000  0.017569        0.00000
ATOM      1  N   ASNDE   1      -7.656  -1.943   3.155  1.00 15.02           N
ATOM      2  CA  ASNDE   1      -6.522  -2.828   2.831  1.00 14.10           C
ATOM      3  C   ASNDE   1      -5.241  -2.329   3.427  1.00 13.13           C
ATOM      4  O   ASNDE   1      -4.978  -1.124   3.426  1.00 11.91           O
ATOM      5  CB  ASNDE   1      -6.346  -2.985   1.341  1.00 15.38           C
ATOM      6  CG  ASNDE   1      -7.584  -3.524   0.692  1.00 14.08           C
ATOM      7  OD1 ASNDE   1      -8.025  -4.639   1.016  1.00 17.46           O
ATOM      8  ND2 ASNDE   1      -8.204  -2.711  -0.169  1.00 11.72           N
ATOM      9  D   ASNDE   1      -8.117  -1.532   2.343  1.00 15.02           D
ATOM     10  DA  ASNDE   1      -6.744  -3.809   3.252  1.00 14.10           D
ATOM     11  DB2 ASNDE   1      -6.124  -2.013   0.900  1.00 15.38           D
ATOM     12  DB3 ASNDE   1      -5.528  -3.678   1.146  1.00 15.38           D
ATOM     13 DD21 ASNDE   1      -9.052  -3.021  -0.643  1.00 11.72           D
ATOM     14 DD22 ASNDE   1      -7.829  -1.780  -0.352  1.00 11.72           D
ATOM     15  N   ASNDE   2      -4.438  -3.276   3.905  1.00 12.26           N
ATOM     16  CA  ASNDE   2      -3.193  -2.962   4.589  1.00 11.74           C
ATOM     17  C   ASNDE   2      -1.955  -3.534   3.895  1.00 11.10           C
ATOM     18  O   ASNDE   2      -1.872  -4.747   3.648  1.00 10.42           O
ATOM     19  CB  ASNDE   2      -3.259  -3.488   6.042  1.00 12.15           C
ATOM     20  CG  ASNDE   2      -2.006  -3.127   6.861  1.00 12.82           C
ATOM     21  OD1 ASNDE   2      -1.702  -1.941   7.072  1.00 15.05           O
ATOM     22  ND2 ASNDE   2      -1.271  -4.151   7.306  1.00 13.48           N
ATOM     23  D   ASNDE   2      -4.626  -4.276   3.832  1.00 12.26           D
ATOM     24  DA  ASNDE   2      -3.074  -1.879   4.586  1.00 11.74           D
ATOM     25  DB2 ASNDE   2      -4.125  -3.051   6.539  1.00 12.15           D
ATOM     26  DB3 ASNDE   2      -3.349  -4.574   6.023  1.00 12.15           D
ATOM     27 DD21 ASNDE   2      -0.429  -3.974   7.854  1.00 13.48           D
ATOM     28 DD22 ASNDE   2      -1.552  -5.109   7.097  1.00 13.48           D
ATOM     29  N   GLNDE   3      -1.005  -2.638   3.598  1.00 10.29           N
ATOM     30  CA  GLNDE   3       0.384  -2.978   3.199  1.00 10.53           C
ATOM     31  C   GLNDE   3       1.435  -2.260   4.088  1.00 10.24           C
ATOM     32  O   GLNDE   3       1.547  -1.023   4.115  1.00  8.86           O
ATOM     33  CB  GLNDE   3       0.656  -2.718   1.711  1.00  9.80           C
ATOM     34  CG  GLNDE   3       1.944  -3.408   1.213  1.00 10.25           C
ATOM     35  CD  GLNDE   3       2.504  -2.822  -0.089  1.00 12.43           C
ATOM     36  OE1 GLNDE   3       2.744  -1.598  -0.190  1.00 14.62           O
ATOM     37  NE2 GLNDE   3       2.750  -3.705  -1.091  1.00  9.05           N
ATOM     38  D   GLNDE   3      -1.162  -1.630   3.622  1.00 10.29           D
ATOM     39  DA  GLNDE   3       0.505  -4.051   3.346  1.00 10.53           D
ATOM     40  DB2 GLNDE   3      -0.180  -3.100   1.125  1.00  9.80           D
ATOM     41  DB3 GLNDE   3       0.764  -1.645   1.552  1.00  9.80           D
ATOM     42  DG2 GLNDE   3       2.713  -3.307   1.979  1.00 10.25           D
ATOM     43  DG3 GLNDE   3       1.730  -4.462   1.036  1.00 10.25           D
ATOM     44 DE21 GLNDE   3       2.562  -4.698  -0.953  1.00  9.05           D
ATOM     45 DE22 GLNDE   3       3.123  -3.375  -1.981  1.00  9.05           D
ATOM     46  N   GLNDE   4       2.154  -3.045   4.871  1.00 10.38           N
ATOM     47  CA  GLNDE   4       3.270  -2.505   5.640  1.00 11.39           C
ATOM     48  C   GLNDE   4       4.594  -3.098   5.172  1.00 11.52           C
ATOM     49  O   GLNDE   4       4.768  -4.320   5.054  1.00 12.05           O
ATOM     50  CB  GLNDE   4       3.056  -2.683   7.147  1.00 11.96           C
ATOM     51  CG  GLNDE   4       1.829  -1.916   7.647  1.00 10.81           C
ATOM     52  CD  GLNDE   4       1.344  -2.452   8.954  1.00 13.10           C
ATOM     53  OE1 GLNDE   4       0.774  -3.541   9.002  1.00 10.65           O
ATOM     54  NE2 GLNDE   4       1.549  -1.679  10.039  1.00 12.30           N
ATOM     55  D   GLNDE   4       1.995  -4.045   4.996  1.00 10.38           D
ATOM     56  DA  GLNDE   4       3.328  -1.429   5.475  1.00 11.39           D
ATOM     57  DB2 GLNDE   4       2.909  -3.741   7.365  1.00 11.96           D
ATOM     58  DB3 GLNDE   4       3.931  -2.311   7.679  1.00 11.96           D
ATOM     59  DG2 GLNDE   4       2.090  -0.866   7.782  1.00 10.81           D
ATOM     60  DG3 GLNDE   4       1.024  -2.009   6.918  1.00 10.81           D
ATOM     61 DE21 GLNDE   4       2.015  -0.777   9.942  1.00 12.30           D
ATOM     62 DE22 GLNDE   4       1.238  -1.996  10.957  1.00 12.30           D
ATOM     63  N   ASNDE   5       5.514  -2.202   4.856  1.00 11.99           N
ATOM     64  CA  ASNDE   5       6.831  -2.556   4.318  1.00 12.30           C
ATOM     65  C   ASNDE   5       7.854  -2.105   5.324  1.00 13.40           C
ATOM     66  O   ASNDE   5       8.219  -0.923   5.374  1.00 13.92           O
ATOM     67  CB  ASNDE   5       7.065  -1.850   2.993  1.00 12.13           C
ATOM     68  CG  ASNDE   5       5.961  -2.131   2.003  1.00 12.77           C
ATOM     69  OD1 ASNDE   5       5.798  -3.262   1.551  1.00 14.27           O
ATOM     70  ND2 ASNDE   5       5.195  -1.119   1.679  1.00 10.07           N
ATOM     71  D   ASNDE   5       5.383  -1.196   4.960  1.00 11.99           D
ATOM     72  DA  ASNDE   5       6.915  -3.627   4.132  1.00 12.30           D
ATOM     73  DB2 ASNDE   5       7.106  -0.774   3.162  1.00 12.13           D
ATOM     74  DB3 ASNDE   5       8.005  -2.196   2.563  1.00 12.13           D
ATOM     75 DD21 ASNDE   5       4.432  -1.249   1.014  1.00 10.07           D
ATOM     76 DD22 ASNDE   5       5.361  -0.201   2.091  1.00 10.07           D
TER
ATOM     77  D2  ASNDE   1      -8.282  -2.409   3.603  1.00 15.02           D
ATOM     78  DC  ASNDE   5       8.249  -2.763   5.965  1.00 13.40           D
TER
END
''',
'''bond pdb=" N   ASNDE   1 "
     pdb=" D2  ASNDE   1 "
  ideal  model  delta    sigma   weight residual
  0.960  0.900  0.060 2.00e-02 2.50e+03 9.04e+00
bond pdb=" C   ASNDE   5 "
     pdb=" DC  ASNDE   5 "
  ideal  model  delta    sigma   weight residual
  1.000  1.000  0.000 2.00e-02 2.50e+03 1.06e-05'''
],
}
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file

def main():
  for code, item in pdbs.items():
    print(code)
    pf = 'tst_deuterium_end_%s.pdb' % code
    with open(pf, 'w') as f:
      f.write(item[0])
    cmd = 'phenix.pdb_interpretation write_geo=True %(pf)s > %(pf)s.log' % locals()
    print(cmd)
    easy_run.go(cmd)
    print(item[1])
    for line in item[1].split('\n'):
      print(line)
      assert_lines_in_file(file_name='%s.geo' % pf, lines=line)

if __name__ == '__main__':
  main()
