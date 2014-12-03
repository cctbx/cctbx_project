from __future__ import division
from libtbx.test_utils import approx_equal
import iotbx.ncs
import iotbx.pdb

pdb_str_1 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   GLY C  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  GLY C  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   GLY C  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   GLY C  34     189.986 271.004 173.508  1.00  0.00           O
TER
"""

pdb_str_2 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  O   HOH S   1     109.583 203.076 175.423  1.00  0.00           O
TER
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   GLY C  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  GLY C  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   GLY C  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   GLY C  34     189.986 271.004 173.508  1.00  0.00           O
TER
ATOM      9  O   TYR D   4     189.583 273.076 175.423  1.00  0.00           O
TER
"""

pdb_str_3 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  O   HOH C  34     189.583 273.076 175.423  1.00  0.00           O
ATOM     10  O   HOH C  34     188.804 273.006 174.173  1.00  0.00           O
ATOM     11  O   HOH C  34     188.920 271.622 173.510  1.00  0.00           O
ATOM     12  O   HOH C  34     189.986 271.004 173.508  1.00  0.00           O
TER
"""

pdb_str_4 = """
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY B  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   TYR C  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  TYR C  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   TYR C  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   TYR C  34     189.986 271.004 173.508  1.00  0.00           O
TER
"""

pdb_str_5 = """
ATOM      1  N   MET A   1     158.070 173.095 147.115  1.00 50.00           N
ATOM      2  CA  MET A   1     157.408 172.627 148.359  1.00 50.00           C
ATOM      3  CB  MET A   1     157.550 171.094 148.516  1.00 50.00           C
ATOM      4  CG  MET A   1     156.748 170.503 149.691  1.00 50.00           C
ATOM      5  SD  MET A   1     154.968 170.855 149.612  1.00 50.00           S
ATOM      6  CE  MET A   1     154.505 169.913 151.091  1.00 50.00           C
ATOM      7  C   MET A   1     157.958 173.331 149.563  1.00 50.00           C
ATOM      8  O   MET A   1     157.196 173.814 150.399  1.00 50.00           O
TER
ATOM      9  N   MET B   1     174.781 155.306 150.054  1.00 50.00           N
ATOM     10  CA  MET B   1     174.332 154.630 151.298  1.00 50.00           C
ATOM     11  CB  MET B   1     175.016 153.251 151.453  1.00 50.00           C
ATOM     12  CG  MET B   1     174.481 152.410 152.628  1.00 50.00           C
ATOM     13  SD  MET B   1     172.693 152.099 152.550  1.00 50.00           S
ATOM     14  CE  MET B   1     172.601 151.052 154.028  1.00 50.00           C
ATOM     15  C   MET B   1     174.594 155.484 152.502  1.00 50.00           C
ATOM     16  O   MET B   1     173.710 155.660 153.339  1.00 50.00           O
TER
ATOM     17  N   MET C   1     148.867 195.697 144.146  1.00 50.00           N
ATOM     18  CA  MET C   1     148.080 195.499 145.390  1.00 50.00           C
ATOM     19  CB  MET C   1     147.662 194.018 145.549  1.00 50.00           C
ATOM     20  CG  MET C   1     146.701 193.755 146.723  1.00 50.00           C
ATOM     21  SD  MET C   1     145.166 194.723 146.643  1.00 50.00           S
ATOM     22  CE  MET C   1     144.395 194.012 148.122  1.00 50.00           C
ATOM     23  C   MET C   1     148.846 195.960 146.594  1.00 50.00           C
ATOM     24  O   MET C   1     148.308 196.685 147.429  1.00 50.00           O
TER
ATOM    417  N   MET 1   1     274.499 237.478  69.907  1.00 50.00           N
ATOM    418  CA  MET 1   1     275.223 237.861  71.146  1.00 50.00           C
ATOM    419  CB  MET 1   1     275.281 239.400  71.298  1.00 50.00           C
ATOM    420  CG  MET 1   1     276.159 239.886  72.466  1.00 50.00           C
ATOM    421  SD  MET 1   1     277.878 239.307  72.379  1.00 50.00           S
ATOM    422  CE  MET 1   1     278.468 240.186  73.852  1.00 50.00           C
ATOM    423  C   MET 1   1     274.593 237.238  72.356  1.00 50.00           C
ATOM    424  O   MET 1   1     275.291 236.663  73.190  1.00 50.00           O
TER
END
"""

pdb_str_6 = """\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
ATOM      6  CA  LYS A 156      10.547   9.007   9.084  1.00 55.55           C
TER
ATOM      7  CA  LYS B 157      11.545   9.413   8.000  1.00 72.27           C
ATOM      8  CA  LYS B 158      12.277  10.718   8.343  1.00 75.78           C
ATOM      9  CA  LYS B 159      11.349  11.791   8.809  1.00 75.88           C
TER
ATOM      7  CA  LYS F 157       2.154   3.953  16.298  1.00 72.27           C
ATOM      8  CA  LYS F 158       2.014   3.732  17.811  1.00 75.78           C
ATOM      9  CA  LYS F 159       2.558   2.413  18.250  1.00 75.88           C
TER
ATOM      7  CA  LYS D 157       4.334  10.965  12.119  1.00 72.27           C
ATOM      8  CA  LYS D 158       4.057  11.980  13.238  1.00 75.78           C
ATOM      9  CA  LYS D 159       3.177  11.427  14.310  1.00 75.88           C
TER
ATOM      1  CA  LYS C 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS C 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS C 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS C 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS C 155       5.801   9.685  13.736  1.00 41.68           C
ATOM      6  CA  LYS C 156       4.731   9.594  12.667  1.00 55.55           C
TER
ATOM      1  CA  LYS E 151       6.987   4.106  17.432  1.00 44.22           C
ATOM      2  CA  LYS E 152       6.017   3.539  16.502  1.00 49.42           C
ATOM      3  CA  LYS E 153       6.497   3.492  15.036  1.00 46.15           C
ATOM      4  CA  LYS E 154       6.348   2.458  14.400  1.00 52.97           C
ATOM      5  CA  LYS E 155       4.647   4.221  16.634  1.00 41.68           C
ATOM      6  CA  LYS E 156       3.552   3.605  15.788  1.00 55.55           C
TER
"""

pdb_str_7 = """\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
TER
ATOM    222  CA  LEU X  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM    223  CA  ARG X  41      62.395  51.344  80.786  1.00107.25           C
ATOM    224  CA  ARG X  42      62.395  41.344  80.786  1.00107.25           C
TER
ATOM      1  CA  THR D   1       8.111  11.080  10.645  1.00 20.00           C
ATOM      2  CA  THR D   2       8.000   9.722  10.125  1.00 20.00           C
ATOM      3  CA  THR D   3       8.075   8.694  11.249  1.00 20.00           C
ATOM      4  CA  THR D   4       8.890   8.818  12.163  1.00 20.00           C
TER
ATOM      1  CA  LYS B 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS B 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS B 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS B 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS B 155       5.801   9.685  13.736  1.00 41.68           C
TER
ATOM      1  CA  LYS C 151       6.987   4.106  17.432  1.00 44.22           C
ATOM      2  CA  LYS C 152       6.017   3.539  16.502  1.00 49.42           C
ATOM      3  CA  LYS C 153       6.497   3.492  15.036  1.00 46.15           C
ATOM      4  CA  LYS C 154       6.348   2.458  14.400  1.00 52.97           C
ATOM      5  CA  LYS C 155       4.647   4.221  16.634  1.00 41.68           C
TER
ATOM    222  CA  LEU Y  40     194.618   5.253  81.582  1.00 87.10           C
ATOM    223  CA  ARG Y  41     162.395  41.344  70.786  1.00107.25           C
ATOM    224  CA  ARG Y  42     162.395  31.344  70.786  1.00107.25           C
TER
ATOM      1  CA  THR E   1       8.111 -10.645  11.080  1.00 20.00           C
ATOM      2  CA  THR E   2       8.000 -10.125   9.722  1.00 20.00           C
ATOM      3  CA  THR E   3       8.075 -11.249   8.694  1.00 20.00           C
ATOM      4  CA  THR E   4       8.890 -12.163   8.818  1.00 20.00           C
TER
"""

pdb_str_8 = """\
ATOM      1  N   THR A   1      13.014  18.419   8.520  1.00 20.00           N
ATOM      2  CA  THR A   1      12.903  17.061   8.000  1.00 20.00           C
ATOM      3  C   THR A   1      12.978  16.033   9.124  1.00 20.00           C
ATOM      4  O   THR A   1      13.793  16.157  10.038  1.00 20.00           O
TER
ATOM      1  N   THR C   1      10.325   8.000  14.368  1.00 20.00           N
ATOM      2  CA  THR C   1      10.111   8.702  13.108  1.00 20.00           C
ATOM      3  C   THR C   1      11.313   9.570  12.750  1.00 20.00           C
ATOM      4  O   THR C   1      11.885  10.241  13.609  1.00 20.00           O
TER
ATOM      1  N   THR D   1      -0.430  15.458  18.495  1.00 20.00           N
ATOM      2  CA  THR D   1       0.086  15.020  17.204  1.00 20.00           C
ATOM      3  C   THR D   1       1.440  14.334  17.355  1.00 20.00           C
ATOM      4  O   THR D   1       2.297  14.791  18.111  1.00 20.00           O
TER
ATOM      1  N   THR E   1       1.895   8.366  19.752  1.00 20.00           N
ATOM      2  CA  THR E   1       1.682   9.068  18.491  1.00 20.00           C
ATOM      3  C   THR E   1       2.884   9.935  18.133  1.00 20.00           C
ATOM      4  O   THR E   1       3.455  10.606  18.993  1.00 20.00           O
TER
ATOM      1  N   THR F   1       8.346   7.308  15.936  1.00 20.00           N
ATOM      2  CA  THR F   1       7.054   7.796  15.467  1.00 20.00           C
ATOM      3  C   THR F   1       6.884   9.281  15.767  1.00 20.00           C
ATOM      4  O   THR F   1       7.237   9.751  16.849  1.00 20.00           O
TER
ATOM      1  N   THR G   1       0.609  -0.560  24.094  1.00 20.00           N
ATOM      2  CA  THR G   1       0.395   0.142  22.834  1.00 20.00           C
ATOM      3  C   THR G   1       1.598   1.009  22.476  1.00 20.00           C
ATOM      4  O   THR G   1       2.169   1.680  23.335  1.00 20.00           O
TER
ATOM      1  N   THR H   1       7.061  -1.617  20.279  1.00 20.00           N
ATOM      2  CA  THR H   1       5.768  -1.130  19.810  1.00 20.00           C
ATOM      3  C   THR H   1       5.599   0.356  20.109  1.00 20.00           C
ATOM      4  O   THR H   1       5.950   0.825  21.191  1.00 20.00           O
TER
ATOM      1  N   THR I   1       8.722   4.822  16.665  1.00 20.00           N
ATOM      2  CA  THR I   1       7.494   4.036  16.653  1.00 20.00           C
ATOM      3  C   THR I   1       6.628   4.350  17.868  1.00 20.00           C
ATOM      4  O   THR I   1       7.130   4.482  18.984  1.00 20.00           O
TER
ATOM      1  N   THR B   1       8.000  15.093  13.112  1.00 20.00           N
ATOM      2  CA  THR B   1       8.516  14.654  11.820  1.00 20.00           C
ATOM      3  C   THR B   1       9.870  13.968  11.972  1.00 20.00           C
ATOM      4  O   THR B   1      10.727  14.426  12.727  1.00 20.00           O
TER
"""


def exercise_00(prefix="iotbx_ncs_exercise_00"):
  pdb_file_name = "%s.pdb"%prefix
  ncs_params_str = """
ncs_group {
  master_selection = chain A
  copy_selection = chain B
  copy_selection = chain C
}
  """
  def check_result(ncs_inp, test_i):
    if test_i == 0:
      l1, l2, l3 = [0,1,2,3], [4,5,6,7], [8,9,10,11]
    elif test_i == 1:
      l1, l2, l3 = [1,2,3,4], [5,6,7,8], [9,10,11,12]
    else: assert 0
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    assert len(ncs_groups) == 1
    ncs_group = ncs_groups[0]
    assert approx_equal(ncs_group.master_iselection, l1)
    assert len(ncs_group.copies) == 2
    assert approx_equal(ncs_group.copies[0].iselection, l2)
    assert approx_equal(ncs_group.copies[1].iselection, l3)
  for test_i, pdb_str in enumerate([pdb_str_1, pdb_str_2]):
    of = open(pdb_file_name, "w")
    print >> of, pdb_str
    of.close()
    pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
    if test_i == 0: # XXX Not implemented. Fix later.
      # using pdb_inp
      ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp)
      check_result(ncs_inp,test_i)
      # using file_name
      ncs_inp = iotbx.ncs.input(file_name = pdb_file_name)
      check_result(ncs_inp,test_i)
      # using pdb string
      ncs_inp = iotbx.ncs.input(pdb_string = pdb_str)
      check_result(ncs_inp,test_i)
    # using combination of pdb_inp and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)

def exercise_01(prefix="iotbx_ncs_exercise_01"):
  """
  Make sure provided selections take precedence and are correctly respected.
  """
  pdb_file_name = "%s.pdb"%prefix
  ncs_params_str = """
ncs_group {
  master_selection = chain C
  copy_selection = chain A
}
  """
  def check_result(ncs_inp, test_i):
    if test_i == 0:
      l1, l2 = [8,9,10,11], [0,1,2,3]
    elif test_i == 1:
      l1, l2 = [9,10,11,12], [1,2,3,4]
    else: assert 0
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    assert len(ncs_groups) == 1
    ncs_group = ncs_groups[0]
    assert approx_equal(ncs_group.master_iselection, l1)
    assert len(ncs_group.copies) == 1
    assert approx_equal(ncs_group.copies[0].iselection, l2)
  for test_i, pdb_str in enumerate([pdb_str_1, pdb_str_2]):
    of = open(pdb_file_name, "w")
    print >> of, pdb_str
    of.close()
    pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
    # using combination of pdb_inp and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_phil_string = ncs_params_str)
    check_result(ncs_inp,test_i)

def exercise_02(prefix="iotbx_ncs_exercise_02"):
  """
  This is expected to fail as requested chains cannot be matched.
  """
  pdb_file_name = "%s.pdb"%prefix
  ncs_params_str = """
ncs_group {
  master_selection = chain C
  copy_selection = chain A
}
  """
  for test_i, pdb_str in enumerate([pdb_str_3, pdb_str_4]):
    of = open(pdb_file_name, "w")
    print >> of, pdb_str
    of.close()
    pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
    # using combination of pdb_inp and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp,
      ncs_phil_string = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_phil_string = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_phil_string = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()

def exercise_03(prefix="iotbx_ncs_exercise_03"):
  """
  Expect one master and 3 copies.
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_5)
  ncs_inp = iotbx.ncs.input(pdb_inp = pdb_inp)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  asc = ncs_inp.hierarchy.atom_selection_cache()
  sel_master = asc.selection(string = "chain 1")
  sel_copy_1 = asc.selection(string = "chain A")
  sel_copy_2 = asc.selection(string = "chain B")
  sel_copy_3 = asc.selection(string = "chain C")
  assert len(ncs_groups)==1
  ng = ncs_groups[0]
  # chains are sorted by name (numbers first)
  assert approx_equal(sel_master.iselection(), ng.master_iselection)
  assert approx_equal(sel_copy_1.iselection(), ng.copies[0].iselection)
  assert approx_equal(sel_copy_2.iselection(), ng.copies[1].iselection)
  assert approx_equal(sel_copy_3.iselection(), ng.copies[2].iselection)

def exercise_04(prefix="iotbx_ncs_exercise_04"):
  """
  Testing one ncs group and master ncs with two chains of different length
  and chains are not in order
  """
  ncs_inp = iotbx.ncs.input(pdb_string=pdb_str_6)
  t = ncs_inp.ncs_to_asu_selection
  assert t.keys() == ['chain A or chain B']
  assert t.values() == [['chain C or chain D', 'chain E or chain F']]

def exercise_05(prefix="iotbx_ncs_exercise_05"):
  """
  Make sure that phil selection overrides ncs grouping
  """
  phil_str = """\
ncs_group {
  master_selection = chain A
  copy_selection = chain C
}
ncs_group {
  master_selection = chain B
  copy_selection = chain D
}
"""
  ncs_inp = iotbx.ncs.input(
    pdb_string=pdb_str_6,
    ncs_phil_string=phil_str)
  expected = {'chain A': ['chain C'], 'chain B': ['chain D']}
  assert ncs_inp.ncs_to_asu_selection.keys(), expected.keys()
  assert ncs_inp.ncs_to_asu_selection.values(), expected.values()

def exercise_06(prefix="iotbx_ncs_exercise_06"):
  """
  Two groups, different number of chain in each group.
  Two chains that are NOT ncs related
  """
  ncs_inp = iotbx.ncs.input(pdb_string=pdb_str_7)
  t = ncs_inp.ncs_to_asu_selection
  assert t.keys() == ['chain D', 'chain A']
  assert t.values() == [['chain E'], ['chain B', 'chain C']]
  assert ncs_inp.ncs_group_map[1][0] == {'chain A'}
  assert ncs_inp.ncs_group_map[2][0] == {'chain D'}

def exercise_07(prefix="iotbx_ncs_exercise_07"):
  """
  Test that minimal number of chains in master ncs are selected (not the
  minimal number of transformations)
  """
  ncs_inp = iotbx.ncs.input(pdb_string=pdb_str_8)
  t = ncs_inp.ncs_to_asu_selection
  assert t.keys() == ['chain A']
  assert t.values() == \
    [['chain B', 'chain C', 'chain D', 'chain E',
      'chain F', 'chain G', 'chain H', 'chain I']]

def exercise_08(prefix="iotbx_ncs_exercise_08"):
  """
  Test that minimal number of transformations in master ncs are selected
  (not the minimal number of chains)
  """
  ncs_inp = iotbx.ncs.input(
    pdb_string=pdb_str_8,
    use_minimal_master_ncs=False)
  t = ncs_inp.ncs_to_asu_selection
  assert t.keys()==['chain A or chain B or chain C']
  assert t.values()==\
    [['chain D or chain E or chain F',
      'chain G or chain H or chain I']]

def exercise_09(prefix="iotbx_ncs_exercise_09"):
  """ ??? """
  from mmtbx.ncs import ncs_search
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_8)
  ph = pdb_inp.construct_hierarchy()
  #
  chain_match_list = ncs_search.search_ncs_relations(
    ph=ph,min_contig_length=0,min_percent=0)
  # make sure that all possible chains compared
  model  = ph.models()[0]
  chain_ids = list({x.id for x in model.chains()})
  chain_ids = sorted(chain_ids)
  n_chains = len(chain_ids)
  assert n_chains==9
  assert len(chain_match_list)==8
  #
  match_dict = ncs_search.clean_chain_matching(chain_match_list,ph)
  chains_info = ncs_search.get_chains_info(ph)
  # Test minimal master NCS
  transform_to_group,match_dict = ncs_search.minimal_master_ncs_grouping(
  match_dict=match_dict)
  group_dict = ncs_search.build_group_dict(
    transform_to_group,match_dict,chains_info)
  assert len(group_dict)==1
  gr_obj = group_dict[('A',)]
  assert len(gr_obj.transforms)==len(gr_obj.copies)
  assert len(gr_obj.iselections)==len(gr_obj.copies)
  expected = [['A'], ['B'], ['C'], ['D'], ['E'], ['F'], ['G'], ['H'], ['I']]
  assert gr_obj.copies == expected
  tr = gr_obj.transforms[0]
  assert tr.r.is_r3_identity_matrix()
  assert tr.t.is_col_zero()
  tr = gr_obj.transforms[1]
  assert not tr.r.is_r3_identity_matrix()
  assert not tr.t.is_col_zero()

def exercise_10(prefix="iotbx_ncs_exercise_10"):
  """ Test minimal NCS operators """
  from mmtbx.ncs import ncs_search
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_8)
  ph = pdb_inp.construct_hierarchy()
  #
  chain_match_list = ncs_search.search_ncs_relations(
    ph=ph,min_contig_length=0,
    min_percent=0,
    use_minimal_master_ncs=False)
  #
  match_dict = ncs_search.clean_chain_matching(chain_match_list,ph)
  chains_info = ncs_search.get_chains_info(ph)
  #
  transform_to_group,match_dict = ncs_search.minimal_ncs_operators_grouping(
  match_dict=match_dict)
  group_dict = ncs_search.build_group_dict(
    transform_to_group,match_dict,chains_info)
  assert len(group_dict)==1
  gr_obj = group_dict[('A', 'B', 'C')]
  assert len(gr_obj.transforms)==len(gr_obj.copies)
  assert len(gr_obj.iselections)==len(gr_obj.copies)
  expected = [['A', 'B', 'C'], ['D', 'E', 'F'], ['G', 'H', 'I']]
  assert gr_obj.copies==expected
  tr = gr_obj.transforms[0]
  assert tr.r.is_r3_identity_matrix()
  assert tr.t.is_col_zero()
  tr = gr_obj.transforms[1]
  assert not tr.r.is_r3_identity_matrix()
  assert not tr.t.is_col_zero()

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_02()
  exercise_03()
  exercise_04()
  exercise_05()
  exercise_06()
  exercise_07()
  exercise_08()
  exercise_09()
  exercise_10()
