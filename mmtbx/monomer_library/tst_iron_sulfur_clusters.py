from __future__ import division
pdbs = {'FES' : ['''
CRYST1   69.637   83.783  289.575  90.00  90.00  90.00 P 21 21 21
SCALE1      0.014360  0.000000  0.000000        0.00000
SCALE2      0.000000  0.011936  0.000000        0.00000
SCALE3      0.000000  0.000000  0.003453        0.00000
ATOM      2  N   CYS B  65      27.906  13.996 110.876  1.00 38.84           N
ATOM      3  CB  CYS B  65      25.655  13.129 111.327  1.00 45.32           C
ATOM      4  SG  CYS B  65      24.970  12.918 109.662  1.00 48.54           S
ATOM      6  CB  CYS B  70      24.347  16.552 110.907  1.00 35.56           C
ATOM      7  SG  CYS B  70      25.742  16.540 109.744  1.00 38.98           S
ATOM      9  CB  CYS B  73      28.351  16.094 103.957  1.00 38.81           C
ATOM     10  SG  CYS B  73      27.104  14.870 103.769  1.00 54.50           S
ATOM     11  CB  CYS B  85      27.375  11.085 104.799  1.00 53.46           C
ATOM     12  SG  CYS B  85      25.854  11.959 105.192  1.00 43.66           S
TER
HETATM   13  S1  FES B1001      27.661  14.411 107.356  1.00 45.15           S
HETATM   14  S2  FES B1001      24.426  14.718 106.631  1.00 50.67           S
HETATM   15 FE1  FES B1001      26.274  13.994 105.655  1.00 45.74          Fe
HETATM   16 FE2  FES B1001      25.714  14.660 108.456  1.00 54.80          Fe
END
''',
'''Iron sulfur cluster coordination
      pdb=" FES B1001 "
        pdb="FE1  FES B1001 " - pdb=" SG  CYS B  85 "
        pdb="FE1  FES B1001 " - pdb=" SG  CYS B  73 "
        pdb="FE2  FES B1001 " - pdb=" SG  CYS B  65 "
        pdb="FE2  FES B1001 " - pdb=" SG  CYS B  70 "
      Number of angles added : 6'''
],
  'SF4' : ['''
CRYST1   69.637   83.783  289.575  90.00  90.00  90.00 P 21 21 21
SCALE1      0.014360  0.000000  0.000000        0.00000
SCALE2      0.000000  0.011936  0.000000        0.00000
SCALE3      0.000000  0.000000  0.003453        0.00000
ATOM      1  CB  CYS B 158      13.213  19.933 106.619  1.00 39.32           C
ATOM      2  SG  CYS B 158      14.305  21.083 105.791  1.00 40.40           S
ATOM      4  SG  CYS B 161      19.842  19.147 102.318  1.00 39.67           S
ATOM      6  CB  CYS B 164      15.403  21.627  98.258  1.00 32.61           C
ATOM      7  SG  CYS B 164      15.348  22.627  99.764  1.00 45.66           S
ATOM      8  CB  CYS B 225      14.524  16.612  99.097  1.00 42.52           C
ATOM      9  SG  CYS B 225      14.453  16.183 100.856  1.00 47.12           S
TER
HETATM   10  S1  SF4 B1002      16.060  18.164 103.702  1.00 42.49           S
HETATM   11  S2  SF4 B1002      16.763  21.423 103.093  1.00 39.97           S
HETATM   12  S3  SF4 B1002      16.716  19.174 100.702  1.00 47.37           S
HETATM   13  S4  SF4 B1002      13.852  20.036 102.235  1.00 41.04           S
HETATM   14 FE1  SF4 B1002      15.572  20.996 101.237  1.00 44.15          Fe
HETATM   15 FE2  SF4 B1002      15.000  18.186 101.701  1.00 44.87          Fe
HETATM   16 FE3  SF4 B1002      15.185  20.158 104.048  1.00 43.71          Fe
HETATM   17 FE4  SF4 B1002      17.638  19.418 102.765  1.00 45.13          Fe
END
''',
'''Iron sulfur cluster coordination
      pdb=" SF4 B1002 "
        pdb="FE3  SF4 B1002 " - pdb=" SG  CYS B 158 "
        pdb="FE1  SF4 B1002 " - pdb=" SG  CYS B 164 "
        pdb="FE2  SF4 B1002 " - pdb=" SG  CYS B 225 "
        pdb="FE4  SF4 B1002 " - pdb=" SG  CYS B 161 "
      Number of angles added : 12'''
],
  'F3S' : ['''
CRYST1   69.637   83.783  289.575  90.00  90.00  90.00 P 21 21 21
SCALE1      0.014360  0.000000  0.000000        0.00000
SCALE2      0.000000  0.011936  0.000000        0.00000
SCALE3      0.000000  0.000000  0.003453        0.00000
ATOM      1  CB  CYS B 168      15.203  26.368  93.567  1.00 47.04           C
ATOM      2  SG  CYS B 168      13.943  27.456  92.871  1.00 47.62           S
ATOM      4  CB  CYS B 215       7.764  28.960  93.921  1.00 37.98           C
ATOM      5  SG  CYS B 215       7.441  27.195  93.528  1.00 39.80           S
ATOM     10  CB  CYS B 221      12.828  21.872  95.639  1.00 37.15           C
ATOM     11  SG  CYS B 221      11.108  21.981  95.007  1.00 44.00           S
TER
HETATM   12  S1  F3S B1003      12.639  23.857  92.203  1.00 32.42           S
HETATM   13  S2  F3S B1003      10.627  27.040  91.359  1.00 41.10           S
HETATM   14  S3  F3S B1003      11.078  25.992  94.637  1.00 40.00           S
HETATM   15  S4  F3S B1003       8.840  23.963  92.856  1.00 38.57           S
HETATM   16 FE1  F3S B1003      12.222  26.021  92.644  1.00 44.92          Fe
HETATM   17 FE3  F3S B1003      10.970  23.908  93.675  1.00 48.58          Fe
HETATM   18 FE4  F3S B1003       9.433  26.146  93.042  1.00 47.61          Fe
END
''',
'''Iron sulfur cluster coordination
      pdb=" F3S B1003 "
        pdb="FE4  F3S B1003 " - pdb=" SG  CYS B 215 "
        pdb="FE1  F3S B1003 " - pdb=" SG  CYS B 168 "
        pdb="FE3  F3S B1003 " - pdb=" SG  CYS B 221 "
      Number of angles added : 9'''
],
}

from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file

def main():
  for code, item in pdbs.items():
    print(code)
    pf = 'tst_fe_s_%s.pdb' % code
    with open(pf, 'w') as f:
      f.write(item[0])
    cmd = 'phenix.pdb_interpretation %(pf)s > %(pf)s.log' % locals()
    print(cmd)
    easy_run.go(cmd)
    print(item[1])
    for line in item[1].split('\n'):
      print(line)
      assert_lines_in_file(file_name='%s.log' % pf, lines=line)

if __name__ == '__main__':
  main()
