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
}

results = {'4udx_sf4.pdb' : [[655, 9675, 12002, 8287],
                             [0,0,13,4],
                             ],
           }

def run():
  fn = '4udx_sf4.pdb'
  f=file(fn, 'wb')
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

  fn = '4udx_sf4_cys.pdb'
  f=file(fn, 'wb')
  f.write(pdbs[fn])
  f.close()
  # lines = '''
  # SF4/F3S coordination
  #    SF4 X1001
  #      pdb="FE1  SF4 X1001 " - pdb=" SG  CYS X  51 "
  #      pdb="FE2  SF4 X1001 " - pdb=" SG  CYS X  48 "
  #      pdb="FE3  SF4 X1001 " - pdb=" SG  CYS X  56 "
  #      pdb="FE4  SF4 X1001 " - pdb=" SG  CYS X  70 "'''

  from libtbx.test_utils import assert_lines_in_file
  lines = [
    'SF4/F3S coordination',
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
  return 0

if __name__=="__main__":
  rc = run()#sys.argv[1])
  assert rc==0

