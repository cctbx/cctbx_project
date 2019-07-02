from __future__ import absolute_import, division, print_function
from mmtbx.validation.clashscore import check_and_add_hydrogen
import libtbx.load_env
import iotbx.pdb
import unittest
import os
from six.moves import map

__author__ = 'Youval'

test_pdb_str = '''\
CRYST1   80.020   97.150   49.850  90.00  90.00  90.00 C 2 2 21      8
SCALE1      0.012497  0.000000  0.000000        0.00000
SCALE2      0.000000  0.010293  0.000000        0.00000
SCALE3      0.000000  0.000000  0.020060        0.00000
ATOM      1  N   CYS A   1      41.272   7.927 -56.483  1.00 44.87           N
ATOM      2  CA  CYS A   1      40.571   9.126 -55.942  1.00 41.71           C
ATOM      3  C   CYS A   1      41.440  10.386 -55.943  1.00 38.08           C
ATOM      4  O   CYS A   1      41.256  11.281 -55.115  1.00 36.43           O
ATOM      5  CB  CYS A   1      40.055   8.833 -54.532  1.00 39.06           C
ATOM      6  SG  CYS A   1      41.291   8.128 -53.442  1.00 42.77           S
ATOM      7  N   ASP A   2      42.366  10.465 -56.895  1.00 41.37           N
ATOM      8  CA  ASP A   2      43.253  11.617 -57.014  1.00 43.37           C
ATOM      9  C   ASP A   2      42.583  12.834 -57.635  1.00 35.58           C
ATOM     10  O   ASP A   2      43.159  13.921 -57.670  1.00 38.51           O
ATOM     11  CB  ASP A   2      44.501  11.254 -57.812  1.00 55.07           C
ATOM     12  CG  ASP A   2      45.724  11.138 -56.938  1.00 64.10           C
ATOM     13  OD1 ASP A   2      45.964  10.034 -56.402  1.00 71.17           O
ATOM     14  OD2 ASP A   2      46.428  12.159 -56.767  1.00 62.17           O
ATOM     15  N   ALA A   3      41.370  12.640 -58.137  1.00 28.02           N
ATOM     16  CA  ALA A   3      40.609  13.720 -58.744  1.00 28.01           C
ATOM     17  C   ALA A   3      39.995  14.617 -57.670  1.00 29.12           C
ATOM     18  O   ALA A   3      39.612  15.753 -57.949  1.00 39.41           O
ATOM     19  CB  ALA A   3      39.515  13.146 -59.633  1.00 23.29           C
ATOM     20  N   PHE A   4      39.918  14.106 -56.444  1.00 23.62           N
ATOM     21  CA  PHE A   4      39.331  14.849 -55.336  1.00 12.37           C
ATOM     22  C   PHE A   4      40.343  15.632 -54.507  1.00 18.90           C
ATOM     23  O   PHE A   4      39.977  16.574 -53.802  1.00 20.87           O
ATOM     24  CB  PHE A   4      38.534  13.900 -54.442  1.00 10.74           C
ATOM     25  CG  PHE A   4      37.432  13.185 -55.164  1.00 10.86           C
ATOM     26  CD1 PHE A   4      36.180  13.772 -55.303  1.00 10.18           C
ATOM     27  CD2 PHE A   4      37.655  11.939 -55.736  1.00 14.14           C
ATOM     28  CE1 PHE A   4      35.166  13.131 -56.003  1.00 13.99           C
ATOM     29  CE2 PHE A   4      36.646  11.286 -56.441  1.00  9.42           C
ATOM     30  CZ  PHE A   4      35.401  11.882 -56.576  1.00 10.19           C
ATOM     31  N   VAL A   5      41.616  15.267 -54.629  1.00 15.08           N
ATOM     32  CA  VAL A   5      42.691  15.918 -53.884  1.00 18.34           C
ATOM     33  C   VAL A   5      42.775  17.426 -54.125  1.00 20.84           C
ATOM     34  O   VAL A   5      42.713  17.886 -55.265  1.00 38.92           O
ATOM     35  CB  VAL A   5      44.055  15.258 -54.208  1.00 21.65           C
ATOM     36  CG1 VAL A   5      45.191  15.981 -53.494  1.00 17.83           C
ATOM     37  CG2 VAL A   5      44.029  13.787 -53.800  1.00 25.16           C
ATOM     38  N   GLY A   6      42.891  18.187 -53.038  1.00 26.12           N
ATOM     39  CA  GLY A   6      42.993  19.635 -53.139  1.00 29.95           C
ATOM     40  C   GLY A   6      42.216  20.401 -52.082  1.00 26.34           C
ATOM     41  O   GLY A   6      41.670  19.815 -51.150  1.00 34.19           O
ATOM     42  N   THR A   7      42.194  21.723 -52.220  1.00 28.07           N
ATOM     43  CA  THR A   7      41.473  22.601 -51.301  1.00 22.77           C
ATOM     44  C   THR A   7      40.276  23.144 -52.072  1.00 21.63           C
ATOM     45  O   THR A   7      40.424  23.581 -53.213  1.00 30.40           O
ATOM     46  CB  THR A   7      42.360  23.765 -50.841  1.00 21.18           C
ATOM     47  OG1 THR A   7      43.614  23.245 -50.383  1.00 28.49           O
ATOM     48  CG2 THR A   7      41.698  24.521 -49.703  1.00 19.00           C
ATOM     49  N   TRP A   8      39.093  23.105 -51.459  1.00 19.58           N
ATOM     50  CA  TRP A   8      37.871  23.559 -52.125  1.00 22.39           C
ATOM     51  C   TRP A   8      37.059  24.548 -51.297  1.00 28.33           C
ATOM     52  O   TRP A   8      36.944  24.404 -50.082  1.00 36.01           O
ATOM     53  CB  TRP A   8      36.977  22.356 -52.476  1.00 25.27           C
ATOM     54  CG  TRP A   8      37.714  21.209 -53.124  1.00 28.71           C
ATOM     55  CD1 TRP A   8      38.477  20.265 -52.493  1.00 28.57           C
ATOM     56  CD2 TRP A   8      37.807  20.925 -54.526  1.00 25.39           C
ATOM     57  NE1 TRP A   8      39.050  19.422 -53.415  1.00 28.39           N
ATOM     58  CE2 TRP A   8      38.656  19.803 -54.670  1.00 23.95           C
ATOM     59  CE3 TRP A   8      37.263  21.510 -55.675  1.00 27.66           C
ATOM     60  CZ2 TRP A   8      38.973  19.258 -55.915  1.00 25.82           C
ATOM     61  CZ3 TRP A   8      37.578  20.965 -56.917  1.00 24.97           C
ATOM     62  CH2 TRP A   8      38.426  19.851 -57.023  1.00 25.28           C
ATOM     63  N   LYS A   9      36.490  25.543 -51.975  1.00 30.17           N
ATOM     64  CA  LYS A   9      35.653  26.559 -51.338  1.00 32.58           C
ATOM     65  C   LYS A   9      34.199  26.310 -51.718  1.00 28.46           C
ATOM     66  O   LYS A   9      33.906  25.943 -52.854  1.00 28.06           O
ATOM     67  CB  LYS A   9      36.023  27.969 -51.820  1.00 40.47           C
ATOM     68  CG  LYS A   9      37.324  28.552 -51.298  1.00 52.85           C
ATOM     69  CD  LYS A   9      37.391  30.038 -51.647  1.00 62.78           C
ATOM     70  CE  LYS A   9      38.673  30.705 -51.157  1.00 72.08           C
ATOM     71  NZ  LYS A   9      39.851  30.410 -52.023  1.00 74.54           N
ATOM     72  N   LEU A  10      33.290  26.516 -50.770  1.00 29.40           N
ATOM     73  CA  LEU A  10      31.866  26.349 -51.036  1.00 26.79           C
ATOM     74  C   LEU A  10      31.452  27.507 -51.944  1.00 35.56           C
ATOM     75  O   LEU A  10      31.912  28.637 -51.765  1.00 38.53           O
ATOM     76  CB  LEU A  10      31.073  26.419 -49.732  1.00 25.83           C
ATOM     77  CG  LEU A  10      29.558  26.286 -49.880  1.00 23.54           C
ATOM     78  CD1 LEU A  10      29.194  24.832 -50.107  1.00 22.85           C
ATOM     79  CD2 LEU A  10      28.870  26.812 -48.635  1.00 32.69           C
HETATM  901  N   NPH A 117      23.870  15.268 -50.490  1.00 25.06           N
HETATM  902  CA  NPH A 117      23.515  14.664 -49.210  1.00 22.94           C
HETATM  903  CB  NPH A 117      23.658  15.702 -48.113  1.00 21.44           C
HETATM  904  SG  NPH A 117      25.281  16.410 -47.839  1.00 31.29           S
HETATM  905  CD  NPH A 117      25.498  16.310 -46.059  0.50 33.70           C
HETATM  906  CE  NPH A 117      26.971  16.492 -45.820  0.50 39.92           C
HETATM  907  OZ  NPH A 117      27.815  15.348 -45.806  0.50 43.26           O
HETATM  908  NZ  NPH A 117      27.411  17.709 -45.649  0.50 43.87           N
HETATM  909  C6  NPH A 117      28.830  18.015 -45.436  0.50 51.14           C
HETATM  910  C5  NPH A 117      29.800  17.009 -45.411  0.50 51.26           C
HETATM  911  C6A NPH A 117      29.195  19.344 -45.262  0.50 52.54           C
HETATM  912  C4A NPH A 117      31.137  17.332 -45.219  0.50 52.42           C
HETATM  913  C10 NPH A 117      30.543  19.684 -45.064  0.50 54.23           C
HETATM  914  C7  NPH A 117      28.201  20.331 -45.288  0.50 53.28           C
HETATM  915  C4  NPH A 117      32.087  16.314 -45.202  0.50 53.40           C
HETATM  916  C1A NPH A 117      31.522  18.669 -45.046  0.50 55.19           C
HETATM  917  N10 NPH A 117      30.892  21.027 -44.889  0.50 53.87           N
HETATM  918  C8  NPH A 117      28.549  21.670 -45.111  0.50 54.83           C
HETATM  919  C3  NPH A 117      33.439  16.599 -45.020  0.50 53.96           C
HETATM  920  N1  NPH A 117      32.879  18.963 -44.861  0.50 57.78           N
HETATM  921  C9  NPH A 117      29.889  22.011 -44.911  0.50 55.29           C
HETATM  922  C2  NPH A 117      33.832  17.924 -44.852  0.50 55.79           C
HETATM  923  C   NPH A 117      22.061  14.206 -49.245  1.00 23.82           C
HETATM  924  O   NPH A 117      21.158  15.023 -49.423  1.00 24.62           O
ATOM    925  N   VAL A 118      21.831  12.913 -49.037  1.00 28.67           N
ATOM    926  CA  VAL A 118      20.474  12.371 -49.062  1.00 33.96           C
ATOM    927  C   VAL A 118      20.069  11.729 -47.737  1.00 36.17           C
ATOM    928  O   VAL A 118      20.854  11.018 -47.109  1.00 40.12           O
ATOM    929  CB  VAL A 118      20.297  11.336 -50.209  1.00 39.91           C
ATOM    930  CG1 VAL A 118      18.831  10.915 -50.336  1.00 38.68           C
ATOM    931  CG2 VAL A 118      20.791  11.918 -51.524  1.00 43.67           C
ATOM    932  N   MET A 119      18.829  11.986 -47.328  1.00 39.11           N
ATOM    933  CA  MET A 119      18.273  11.449 -46.092  1.00 39.84           C
ATOM    934  C   MET A 119      16.761  11.361 -46.256  1.00 45.51           C
ATOM    935  O   MET A 119      16.066  12.375 -46.173  1.00 45.08           O
ATOM    936  CB  MET A 119      18.616  12.370 -44.916  1.00 46.77           C
ATOM    937  CG  MET A 119      19.314  11.688 -43.740  1.00 45.80           C
ATOM    938  SD  MET A 119      18.225  10.980 -42.492  1.00 41.41           S
ATOM    939  CE  MET A 119      17.987   9.359 -43.118  1.00 44.09           C
ATOM    940  N   LYS A 120      16.265  10.152 -46.521  1.00 52.01           N
ATOM    941  CA  LYS A 120      14.834   9.908 -46.706  1.00 55.50           C
ATOM    942  C   LYS A 120      14.220  10.827 -47.764  1.00 61.07           C
ATOM    943  O   LYS A 120      13.324  11.621 -47.467  1.00 63.20           O
ATOM    944  CB  LYS A 120      14.079  10.073 -45.381  1.00 57.57           C
ATOM    945  CG  LYS A 120      14.393   9.026 -44.324  1.00 62.80           C
ATOM    946  CD  LYS A 120      13.804   9.426 -42.977  1.00 66.61           C
ATOM    947  CE  LYS A 120      14.422  10.730 -42.479  1.00 72.33           C
ATOM    948  NZ  LYS A 120      13.835  11.219 -41.198  1.00 69.78           N
TER    1036      ALA A 131
HETATM 1037  O   HOH A1001       8.647  20.222 -39.056  1.00 45.16           O
HETATM 1038  O   HOH A1002      16.907   6.512 -42.744  1.00 33.42           O
HETATM 1039  O   HOH A1003      28.827   2.276 -52.623  1.00 38.38           O
HETATM 1040  O   HOH A1004      30.966  10.844 -61.160  1.00 40.82           O
HETATM 1041  O   HOH A1005      30.838  16.415 -40.165  1.00 29.15           O
HETATM 1042  O   HOH A1006      31.573   1.965 -31.583  1.00 39.89           O
HETATM 1043  O   HOH A1007      37.321   3.217 -36.538  1.00 30.75           O
HETATM 1044  O   HOH A1008      38.663   3.050 -42.546  1.00 40.80           O
HETATM 1045  O   HOH A1009      38.830   9.609 -38.415  1.00 40.22           O
HETATM 1046  O   HOH A1010      14.636  22.139 -33.478  1.00 37.49           O
HETATM 1047  O   HOH A1011      16.639  25.103 -47.161  1.00 43.03           O
HETATM 1048  O   HOH A1012      17.561  26.508 -50.212  1.00 36.05           O
HETATM 1049  O   HOH A1013      20.199  26.017 -48.000  1.00 30.09           O
HETATM 1050  O   HOH A1014      19.269  19.390 -54.054  1.00 32.98           O
HETATM 1051  O   HOH A1015      24.341   3.887 -33.555  1.00 48.69           O
HETATM 1052  O   HOH A1016      24.806  10.616 -58.584  1.00 53.76           O
HETATM 1053  O   HOH A1017      29.937   3.471 -49.079  1.00 11.83           O
HETATM 1054  O   HOH A1018      31.421  10.991 -50.126  1.00 18.05           O
HETATM 1055  O   HOH A1019      35.725   3.882 -45.468  1.00  5.17           O
HETATM 1056  O   HOH A1020      35.293  11.651 -35.895  1.00 12.46           O
HETATM 1057  O   HOH A1021      36.784  16.657 -63.380  1.00 28.91           O
HETATM 1058  O   HOH A1022      38.577  25.715 -64.560  1.00 46.10           O
HETATM 1059  O   HOH A1023      39.943  20.367 -62.416  1.00 29.86           O
HETATM 1060  O   HOH A1024      42.053   3.459 -49.050  1.00 24.27           O
HETATM 1061  O   HOH A1026      43.187  22.656 -54.818  1.00 33.44           O
HETATM 1062  O   HOH A1027      47.201  12.336 -45.984  1.00 35.11           O
HETATM 1063  O   HOH A1028      53.460   7.293 -47.113  1.00 59.68           O
HETATM 1064  O   HOH A2001      21.683  15.390 -34.744  1.00 34.84           O
HETATM 1065  O   HOH A2002      33.838  15.714 -65.621  1.00 51.47           O
HETATM 1066  O   HOH A2003      36.315  13.744 -63.033  1.00 50.47           O
HETATM 1067  O   HOH A2004      47.315  14.827 -51.041  1.00 38.84           O
HETATM 1068  O   HOH A3001      18.141  28.190 -41.798  1.00 76.07           O
HETATM 1069  O   HOH A3002      26.887  17.479 -42.522  1.00 65.69           O
HETATM 1070  O   HOH A3003      27.587  21.810 -41.976  1.00 83.31           O
HETATM 1071  O   HOH A3004      29.928  13.873 -41.695  1.00 54.90           O
HETATM 1072  O   HOH A3005      32.870  14.371 -39.386  1.00 43.83           O
HETATM 1073  O   HOH A3006      33.981  15.263 -42.437  1.00 51.87           O
HETATM 1074  O   HOH A3007      49.639  13.800 -48.775  1.00 78.70           O
'''


class MyTestCase(unittest.TestCase):

  def setUp(self):
    self.file_to_delete = []
    # import files used in tests
    self.file_name = 'test_pdb_file.pdb'
    open(self.file_name,'w').write(test_pdb_str)
    self.file_to_delete.append(self.file_name)

  def test_identifying_and_addition_of_hydrogen(self):
    """ test identifying and addition of hydrogen """
    has_reduce = libtbx.env.has_module(name="reduce")
    if has_reduce:
      pdb_inp = iotbx.pdb.input(file_name=self.file_name)
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      elements = pdb_hierarchy.atoms().extract_element()
      h_count_0 = elements.count(' H') + elements.count(' D')
      new_pdb_str,_ = check_and_add_hydrogen(
        file_name=self.file_name,
        verbose=False)

      pdb_inp = iotbx.pdb.input(source_info=None, lines=new_pdb_str)
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      elements = pdb_hierarchy.atoms().extract_element()
      h_count_1 = elements.count(' H') + elements.count(' D')

      self.assertEqual(h_count_0,0)
      self.assertTrue(h_count_1>0)
    else:
      # Skip test if reduce is not present
      pass

  def tearDown(self):
    """ delete files created in during testing"""
    if self.file_to_delete:
      for fn in self.file_to_delete:
        if os.path.isfile(fn): os.remove(fn)

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_something']
  suite = unittest.TestSuite(list(map(MyTestCase, tests)))
  return suite


if __name__ == '__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
