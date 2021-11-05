from __future__ import absolute_import, division, print_function
import sys
from libtbx import easy_run

pdbs = {'4udx_sf4.pdb' : '''
HETATM 9499 FE1  SF4 X1001      39.663   2.346  22.723  1.00  6.06          FE
HETATM 9500 FE2  SF4 X1001      42.157   2.138  23.788  1.00  5.97          FE
HETATM 9501 FE3  SF4 X1001      41.050   0.016  22.600  1.00  6.18          FE
HETATM 9502 FE4  SF4 X1001      41.899   2.128  21.073  1.00  6.14          FE
HETATM 9503  S1  SF4 X1001      41.420   3.878  22.535  1.00  6.14           S
HETATM 9504  S2  SF4 X1001      40.249   1.068  24.522  1.00  6.06           S
HETATM 9505  S3  SF4 X1001      43.249   0.748  22.278  1.00  6.36           S
HETATM 9506  S4  SF4 X1001      39.895   0.913  20.872  1.00  6.51           S
  ''',
  '4udx_sf4_cys.pdb' : '''
ATOM    683  N   CYS X  48      44.683  -0.694  27.648  1.00  5.76           N
ATOM    684  CA  CYS X  48      44.043   0.578  27.358  1.00  5.66           C
ATOM    685  C   CYS X  48      42.533   0.474  27.608  1.00  5.45           C
ATOM    686  O   CYS X  48      41.867  -0.426  27.113  1.00  6.44           O
ATOM    687  CB  CYS X  48      44.312   1.012  25.919  1.00  5.93           C
ATOM    688  SG  CYS X  48      43.598   2.633  25.504  1.00  5.98           S
ATOM    734  N   CYS X  51      38.525   4.088  25.616  1.00  5.96           N
ATOM    735  CA  CYS X  51      37.209   4.232  25.021  1.00  6.18           C
ATOM    736  C   CYS X  51      36.381   2.999  25.357  1.00  5.94           C
ATOM    737  O   CYS X  51      36.866   2.008  25.906  1.00  6.63           O
ATOM    738  CB  CYS X  51      37.245   4.495  23.514  1.00  6.07           C
ATOM    739  SG  CYS X  51      37.502   3.012  22.476  1.00  6.27           S
ATOM    801  N   CYS X  56      41.197  -2.946  25.932  1.00  6.08           N
ATOM    802  CA  CYS X  56      42.138  -3.527  25.009  1.00  6.07           C
ATOM    803  C   CYS X  56      43.396  -3.916  25.747  1.00  5.94           C
ATOM    804  O   CYS X  56      43.802  -3.243  26.704  1.00  7.17           O
ATOM    805  CB  CYS X  56      42.512  -2.569  23.889  1.00  6.66           C
ATOM    806  SG  CYS X  56      41.086  -2.294  22.740  1.00  6.50           S
ATOM   1022  N   CYS X  70      41.337   0.053  17.877  1.00  6.60           N
ATOM   1023  CA  CYS X  70      41.399   1.221  16.973  1.00  6.52           C
ATOM   1024  C   CYS X  70      42.679   1.255  16.144  1.00  6.45           C
ATOM   1025  O   CYS X  70      42.824   2.158  15.296  1.00  7.57           O
ATOM   1026  CB  CYS X  70      41.302   2.516  17.784  1.00  6.74           C
ATOM   1027  SG  CYS X  70      42.684   2.663  18.965  1.00  6.45           S
HETATM 9499 FE1  SF4 X1001      39.663   2.346  22.723  1.00  6.06          FE
HETATM 9500 FE2  SF4 X1001      42.157   2.138  23.788  1.00  5.97          FE
HETATM 9501 FE3  SF4 X1001      41.050   0.016  22.600  1.00  6.18          FE
HETATM 9502 FE4  SF4 X1001      41.899   2.128  21.073  1.00  6.14          FE
HETATM 9503  S1  SF4 X1001      41.420   3.878  22.535  1.00  6.14           S
HETATM 9504  S2  SF4 X1001      40.249   1.068  24.522  1.00  6.06           S
HETATM 9505  S3  SF4 X1001      43.249   0.748  22.278  1.00  6.36           S
HETATM 9506  S4  SF4 X1001      39.895   0.913  20.872  1.00  6.51           S
  ''',
  '4u9h_f3s_cys.pdb' : '''
CRYST1   66.530   97.950  125.850  90.00  90.00  90.00 P 21 21 21
SCALE1      0.015031  0.000000  0.000000        0.00000
SCALE2      0.000000  0.010209  0.000000        0.00000
SCALE3      0.000000  0.000000  0.007946        0.00000
ATOM      1  N   CYS S 231     -40.840  -6.337  14.851  1.00  5.51           N
ATOM      2  CA  CYS S 231     -40.013  -5.494  14.018  1.00  5.50           C
ATOM      3  C   CYS S 231     -40.789  -4.444  13.235  1.00  5.59           C
ATOM      4  O   CYS S 231     -40.377  -3.294  13.181  1.00  6.08           O
ATOM      5  CB  CYS S 231     -39.137  -6.330  13.077  1.00  5.85           C
ATOM      6  SG  CYS S 231     -38.302  -7.647  14.022  1.00  5.36           S
ATOM      7  H   CYS S 231     -40.715  -7.176  14.722  1.00  8.43           H
ATOM      8  HA  CYS S 231     -39.406  -5.011  14.590  1.00  4.25           H
ATOM      9  HB2 CYS S 231     -39.693  -6.723  12.387  1.00  8.29           H
ATOM     10  HB3 CYS S 231     -38.460  -5.757  12.678  1.00  6.29           H
ATOM     11  N   CYS S 249     -36.838 -10.889   5.704  1.00  5.29           N
ATOM     12  CA  CYS S 249     -36.582 -11.388   7.062  1.00  5.34           C
ATOM     13  C   CYS S 249     -35.867 -12.738   6.948  1.00  5.19           C
ATOM     14  O   CYS S 249     -36.217 -13.566   6.098  1.00  5.72           O
ATOM     15  CB  CYS S 249     -37.959 -11.517   7.758  1.00  5.35           C
ATOM     16  SG  CYS S 249     -37.935 -12.287   9.420  1.00  5.03           S
ATOM     17  H   CYS S 249     -37.663 -10.835   5.492  1.00  6.38           H
ATOM     18  HA  CYS S 249     -36.025 -10.764   7.557  1.00  6.71           H
ATOM     19  HB2 CYS S 249     -38.359 -10.636   7.819  1.00  7.35           H
ATOM     20  HB3 CYS S 249     -38.524 -12.071   7.205  1.00  6.05           H
ATOM     40  N   CYS S 252     -34.737 -14.293  12.518  1.00  5.11           N
ATOM     41  CA  CYS S 252     -33.673 -14.292  13.509  1.00  4.88           C
ATOM     42  C   CYS S 252     -34.097 -14.398  14.968  1.00  5.08           C
ATOM     43  O   CYS S 252     -33.221 -14.580  15.827  1.00  5.63           O
ATOM     44  CB  CYS S 252     -32.706 -13.130  13.306  1.00  5.00           C
ATOM     45  SG  CYS S 252     -33.261 -11.513  13.953  1.00  5.02           S
ATOM     46  H   CYS S 252     -34.873 -13.554  12.121  1.00  4.62           H
ATOM     47  HA  CYS S 252     -33.162 -15.093  13.343  0.99  4.06           H
ATOM     48  HB2 CYS S 252     -31.880 -13.344  13.777  1.00  5.77           H
ATOM     49  HB3 CYS S 252     -32.490 -13.045  12.358  1.00  4.19           H
TER
HETATM   50  S1  F3S S1003     -37.373  -8.610  10.393  1.00  5.80           S
HETATM   51  S2  F3S S1003     -34.483 -11.028  10.412  1.00  5.12           S
HETATM   52  S3  F3S S1003     -37.067 -11.059  12.944  1.00  5.03           S
HETATM   53  S4  F3S S1003     -34.752  -8.309  13.065  1.00  5.36           S
HETATM   54 FE1  F3S S1003     -36.744 -10.733  10.683  1.00  4.93          Fe2+
HETATM   55 FE3  F3S S1003     -36.932  -8.794  12.568  1.00  5.27          Fe2+
HETATM   56 FE4  F3S S1003     -34.843 -10.514  12.587  1.00  4.74          Fe2+
END
'''
}

results = {'4udx_sf4.pdb' : [[722, 9714, 12002, 8287],
                             [0,0,13,4],
                             ],
           }

def run():
  fn = '4udx_sf4.pdb'
  f=open(fn, 'w')
  f.write(pdbs[fn])
  f.close()
  for i, superpose in enumerate(['None', 'all']):
    cmd = 'phenix.geometry_minimization %s superpose_ideal=%s' % (fn,
                                                                 superpose,
      )
    print(cmd)
    rc = easy_run.go(cmd)
    for line in rc.stdout_lines:
      #print line
      if line.find('bond_residual_sum')>-1:
        bond_value = round(float(line.split()[-1]))
      if line.find('angle_residual_sum')>-1:
        angle_value = round(float(line.split()[-1]))
    # in the second iteration bond_value and angle_value are expected to be 0!
    # assert bond_value, "%f %f" % (bond_value, results[fn][i][0])
    # assert angle_value, angle_value
    assert bond_value==results[fn][i][0], 'not matching %s to %s' % (
          bond_value,
          results[fn][i][0],
          )
    assert angle_value==results[fn][i][1], 'not matching %s to %s' % (
      angle_value,
      results[fn][i][1],
      )

  for i, superpose in enumerate(['None', 'all']):
    cmd = 'phenix.pdb_interpretation %s superpose_ideal=%s' % (fn,
                                                               superpose,
      )
    print(cmd)
    rc = easy_run.go(cmd)
    for line in rc.stdout_lines:
      #print line
      if line.find('bond_residual_sum')>-1:
        value = round(float(line.split()[-1]))
        print('value',value)
        assert value==results[fn][i][2], 'not matching %s to %s' % (
          value,
          results[fn][i][2],
          )
      if line.find('angle_residual_sum')>-1:
        value = round(float(line.split()[-1]))
        print('value',value)
        assert value==results[fn][i][3], 'not matching %s to %s' % (
          value,
          results[fn][i][3],
          )

  from libtbx.test_utils import assert_lines_in_file
  fn = '4udx_sf4_cys.pdb'
  f=open(fn, 'w')
  f.write(pdbs[fn])
  f.close()
  lines = [
    'Iron sulfur cluster coordination',
    'SF4 X1001',
    'pdb="FE1  SF4 X1001 " - pdb=" SG  CYS X  51 "',
    'pdb="FE4  SF4 X1001 " - pdb=" SG  CYS X  70 "',
    'pdb="FE3  SF4 X1001 " - pdb=" SG  CYS X  56 "',
    'pdb="FE2  SF4 X1001 " - pdb=" SG  CYS X  48 "',
    ]
  cmd = 'phenix.pdb_interpretation %s link_all=True > %s.log' % (fn, fn)
  print(cmd)
  assert not easy_run.call(cmd)
  for line in lines:
    print(line)
    assert_lines_in_file(file_name='%s.log' % fn, lines=line)
  fn = '4u9h_f3s_cys.pdb'
  with open(fn, 'w') as f:
    f.write(pdbs[fn])
  lines = [
    'Iron sulfur cluster coordination',
    'F3S S1003',
    'pdb="FE1  F3S S1003 " - pdb=" SG  CYS S 249 "',
    'pdb="FE3  F3S S1003 " - pdb=" SG  CYS S 231 "',
    'pdb="FE4  F3S S1003 " - pdb=" SG  CYS S 252 "',
    ]
  cmd = 'phenix.pdb_interpretation %s link_all=True > %s.log' % (fn, fn)
  print(cmd)
  assert not easy_run.call(cmd)
  for line in lines:
    print(line)
    assert_lines_in_file(file_name='%s.log' % fn, lines=line)
  return 0

if __name__=="__main__":
  rc = run()#sys.argv[1])
  assert rc==0
