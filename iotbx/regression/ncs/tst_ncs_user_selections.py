import time
from iotbx.ncs import ncs_group_master_phil
import iotbx.phil
import iotbx.ncs

def get_ncs_groups(phil_str, pdb_str, **kwargs):
  phil_groups = ncs_group_master_phil.fetch(
      iotbx.phil.parse(phil_str)).extract()
  ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_phil_groups=phil_groups.ncs_group, **kwargs)
  return ncs_inp.get_ncs_restraints_group_list()

def exercise_1():
  """
  Simple selections, no tricks. Make sure reference and selections are
  what user want. All three chains are NCS-related.
  """
  phil_str1="""
ncs_group {
  reference = chain E
  selection = chain F
}
"""
  phil_str2="""
ncs_group {
  reference = chain E
  selection = chain F
  selection = chain G
}
"""
  phil_str3="""
ncs_group {
  reference = chain G
  selection = chain E
}
"""
  phil_str4="""
ncs_group {
  reference = chain G
  selection = chain F
  selection = chain E
}
"""

  pdb_str = """\
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY E  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY E  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY E  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY E  34     127.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N   GLY F  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY F  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY F  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY F  34     249.986 144.872 173.510  1.00  0.00           O
TER
ATOM      9  N   GLY G  34     189.583 273.076 175.423  1.00  0.00           N
ATOM     10  CA  GLY G  34     188.804 273.006 174.173  1.00  0.00           C
ATOM     11  C   GLY G  34     188.920 271.622 173.510  1.00  0.00           C
ATOM     12  O   GLY G  34     189.986 271.004 173.508  1.00  0.00           O
TER
"""

  ncs_groups = get_ncs_groups(phil_str1, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3]
  assert list(ncs_groups[0].copies[0].iselection) == [4, 5, 6, 7]
  # print "master isel:", list(ncs_groups[0].master_iselection)
  # print "copy isel:", list(ncs_groups[0].copies[0].iselection)
  # print "="*80

  ncs_groups = get_ncs_groups(phil_str2, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3]
  assert list(ncs_groups[0].copies[0].iselection) == [4, 5, 6, 7]
  assert list(ncs_groups[0].copies[1].iselection) == [8, 9, 10, 11]

  # print "master isel:", list(ncs_groups[0].master_iselection)
  # print "copy isel:", list(ncs_groups[0].copies[0].iselection)
  # print "copy isel:", list(ncs_groups[0].copies[1].iselection)
  # print "="*80


  ncs_groups = get_ncs_groups(phil_str3, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [8, 9, 10, 11]
  assert list(ncs_groups[0].copies[0].iselection) == [0, 1, 2, 3]

  # print "master isel:", list(ncs_groups[0].master_iselection)
  # print "copy isel:", list(ncs_groups[0].copies[0].iselection)
  # print "="*80
  ncs_groups = get_ncs_groups(phil_str4, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [8, 9, 10, 11]
  assert list(ncs_groups[0].copies[0].iselection) == [4, 5, 6, 7]
  assert list(ncs_groups[0].copies[1].iselection) == [0, 1, 2, 3]

  # print "master isel:", list(ncs_groups[0].master_iselection)
  # print "copy isel:", list(ncs_groups[0].copies[0].iselection)
  # print "copy isel:", list(ncs_groups[0].copies[1].iselection)
  # print "="*80


def exercise_2():
  """
  User tries to define non-matching groups.
  A matches B, but C doesn't match anything
  """

  phil_str1="""
ncs_group {
  reference = chain A
  selection = chain C
}
"""

  phil_str2="""
ncs_group {
  reference = chain A
  selection = chain B
  selection = chain C
}
"""

  pdb_str = """\
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
  ncs_groups = get_ncs_groups(phil_str1, pdb_str)
  assert len(ncs_groups) == 0

  ncs_groups = get_ncs_groups(phil_str2, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3]
  assert list(ncs_groups[0].copies[0].iselection) == [4, 5, 6, 7]


def exercise_3():
  """
  User can split one chain into ncs groups
  """

  phil_str1="""
ncs_group {
  reference = chain A and resid 34
  selection = chain A and resid 35
}
"""

  pdb_str = """\
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CA  GLY A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   GLY A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  O   GLY A  34     127.304 211.628 173.503  1.00  0.00           O
ATOM      5  N   GLY A  35     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CA  GLY A  35     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   GLY A  35     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  O   GLY A  35     249.986 144.872 173.510  1.00  0.00           O
TER
"""
  ncs_groups = get_ncs_groups(phil_str1, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3]
  assert list(ncs_groups[0].copies[0].iselection) == [4, 5, 6, 7]


def exercise_4():
  """
  User can specify multiple chains in one selection
  """

  pdb_str = """\
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
ATOM      1  CA  LYS C 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS C 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS C 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS C 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS C 155       5.801   9.685  13.736  1.00 41.68           C
ATOM      6  CA  LYS C 156       4.731   9.594  12.667  1.00 55.55           C
TER
ATOM      7  CA  LYS D 157       4.334  10.965  12.119  1.00 72.27           C
ATOM      8  CA  LYS D 158       4.057  11.980  13.238  1.00 75.78           C
ATOM      9  CA  LYS D 159       3.177  11.427  14.310  1.00 75.88           C
TER
ATOM      1  CA  LYS E 151       6.987   4.106  17.432  1.00 44.22           C
ATOM      2  CA  LYS E 152       6.017   3.539  16.502  1.00 49.42           C
ATOM      3  CA  LYS E 153       6.497   3.492  15.036  1.00 46.15           C
ATOM      4  CA  LYS E 154       6.348   2.458  14.400  1.00 52.97           C
ATOM      5  CA  LYS E 155       4.647   4.221  16.634  1.00 41.68           C
ATOM      6  CA  LYS E 156       3.552   3.605  15.788  1.00 55.55           C
TER
ATOM      7  CA  LYS F 157       2.154   3.953  16.298  1.00 72.27           C
ATOM      8  CA  LYS F 158       2.014   3.732  17.811  1.00 75.78           C
ATOM      9  CA  LYS F 159       2.558   2.413  18.250  1.00 75.88           C
TER
"""

  phil_str1="""
ncs_group {
  reference = chain A or chain B
  selection = chain C or chain D
  selection = chain E or chain F
}
"""
  phil_str2="""
ncs_group {
  reference = chain E or chain F
  selection = chain C or chain D
  selection = chain A or chain B
}
"""
  phil_str3="""
ncs_group {
  reference = chain A or chain B
  selection = chain C or chain D
}
"""
  phil_str4="""
ncs_group {
  reference = chain A
  selection = chain C
  selection = chain E
}
"""

  ncs_groups = get_ncs_groups(phil_str1, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3, 4, 5, 6, 7, 8]
  assert list(ncs_groups[0].copies[0].iselection) == [9, 10, 11, 12, 13, 14, 15, 16, 17]
  assert list(ncs_groups[0].copies[1].iselection) == [18, 19, 20, 21, 22, 23, 24, 25, 26]

  ncs_groups = get_ncs_groups(phil_str2, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [18, 19, 20, 21, 22, 23, 24, 25, 26]
  assert list(ncs_groups[0].copies[0].iselection) == [9, 10, 11, 12, 13, 14, 15, 16, 17]
  assert list(ncs_groups[0].copies[1].iselection) == [0, 1, 2, 3, 4, 5, 6, 7, 8]

  ncs_groups = get_ncs_groups(phil_str3, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3, 4, 5, 6, 7, 8]
  assert list(ncs_groups[0].copies[0].iselection) == [9, 10, 11, 12, 13, 14, 15, 16, 17]

  ncs_groups = get_ncs_groups(phil_str4, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3, 4, 5]
  assert list(ncs_groups[0].copies[0].iselection) == [9, 10, 11, 12, 13, 14]
  assert list(ncs_groups[0].copies[1].iselection) == [18, 19, 20, 21, 22, 23]


def exercise_5():
  """
  User can specify selections for pdb file with blank chain ids, but with
  segids.
  Inputs are same as in exercise_1(), except chains renamed to segID.
  XXX. Currenly will not be outputted correctly, but should work fine.
  """

  phil_str1="""
ncs_group {
  reference = segid seg1
  selection = segid seg2
}
"""
  phil_str2="""
ncs_group {
  reference = segid seg1
  selection = segid seg2
  selection = segid seg3
}
"""
  phil_str3="""
ncs_group {
  reference = segid seg3
  selection = segid seg1
}
"""
  phil_str4="""
ncs_group {
  reference = segid seg3
  selection = segid seg2
  selection = segid seg1
}
"""

  pdb_str = """\
CRYST1  399.000  399.000  399.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY    34     125.208 211.886 175.417  1.00  0.00      seg1 N
ATOM      2  CA  GLY    34     125.035 211.123 174.168  1.00  0.00      seg1 C
ATOM      3  C   GLY    34     126.386 210.806 173.507  1.00  0.00      seg1 C
ATOM      4  O   GLY    34     127.304 211.628 173.503  1.00  0.00      seg1 O
TER
ATOM      5  N   GLY    34     251.532 143.432 175.422  1.00  0.00      seg2 N
ATOM      6  CA  GLY    34     252.120 143.948 174.173  1.00  0.00      seg2 C
ATOM      7  C   GLY    34     251.212 144.998 173.512  1.00  0.00      seg2 C
ATOM      8  O   GLY    34     249.986 144.872 173.510  1.00  0.00      seg2 O
TER
ATOM      9  N   GLY    34     189.583 273.076 175.423  1.00  0.00      seg3 N
ATOM     10  CA  GLY    34     188.804 273.006 174.173  1.00  0.00      seg3 C
ATOM     11  C   GLY    34     188.920 271.622 173.510  1.00  0.00      seg3 C
ATOM     12  O   GLY    34     189.986 271.004 173.508  1.00  0.00      seg3 O
TER
"""
  ncs_groups = get_ncs_groups(phil_str1, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3]
  assert list(ncs_groups[0].copies[0].iselection) == [4, 5, 6, 7]

  ncs_groups = get_ncs_groups(phil_str2, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3]
  assert list(ncs_groups[0].copies[0].iselection) == [4, 5, 6, 7]
  assert list(ncs_groups[0].copies[1].iselection) == [8, 9, 10, 11]

  ncs_groups = get_ncs_groups(phil_str3, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [8, 9, 10, 11]
  assert list(ncs_groups[0].copies[0].iselection) == [0, 1, 2, 3]

  ncs_groups = get_ncs_groups(phil_str4, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [8, 9, 10, 11]
  assert list(ncs_groups[0].copies[0].iselection) == [4, 5, 6, 7]
  assert list(ncs_groups[0].copies[1].iselection) == [0, 1, 2, 3]


def exercise_6():
  """
  Automatic water filtering. Chains A, B, C, D are NCS-related even with
  water molecule
  """

  phil_str1="""
ncs_group {
  reference = chain A
  selection = chain B
}
"""

  pdb_str = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1      11.782  12.419   4.645  1.00 10.00           N
ATOM      2  CA  THR A   1      11.671  11.061   4.125  1.00 10.00           C
ATOM      3  C   THR A   1      11.746  10.033   5.249  1.00 10.00           C
ATOM      4  O   THR A   1      12.561  10.157   6.163  1.00 10.00           O
ATOM      5  CB  THR A   1      12.772  10.760   3.092  1.00 10.00           C
ATOM      6  OG1 THR A   1      12.672  11.682   2.000  1.00 10.00           O
ATOM      7  O   HOH A   2      12.635   9.340   2.565  1.00 10.00           O
TER
ATOM      1  N   THR D   1      13.010   5.595  10.010  1.00 10.00           N
ATOM      2  CA  THR D   1      14.035   5.945   9.034  1.00 10.00           C
ATOM      3  C   THR D   1      15.310   6.423   9.720  1.00 10.00           C
ATOM      4  O   THR D   1      16.415   6.054   9.323  1.00 10.00           O
ATOM      5  CB  THR D   1      13.543   7.038   8.066  1.00 10.00           C
ATOM      6  OG1 THR D   1      13.237   8.229   8.802  1.00 10.00           O
ATOM      7  O   HOH D   2      12.300   6.572   7.323  1.00 10.00           O
TER
ATOM      1  N   THR B   1       6.768   9.093   9.237  1.00 10.00           N
ATOM      2  CA  THR B   1       7.284   8.654   7.945  1.00 10.00           C
ATOM      3  C   THR B   1       8.638   7.968   8.097  1.00 10.00           C
ATOM      4  O   THR B   1       9.495   8.426   8.852  1.00 10.00           O
ATOM      5  CB  THR B   1       7.423   9.832   6.963  1.00 10.00           C
ATOM      6  OG1 THR B   1       6.144  10.446   6.765  1.00 10.00           O
ATOM      7  O   HOH B   2       7.962   9.350   5.625  1.00 10.00           O
TER
ATOM      1  N   THR C   1       9.093   2.000  10.493  1.00 10.00           N
ATOM      2  CA  THR C   1       8.879   2.702   9.233  1.00 10.00           C
ATOM      3  C   THR C   1      10.081   3.570   8.875  1.00 10.00           C
ATOM      4  O   THR C   1      10.652   4.241   9.734  1.00 10.00           O
ATOM      5  CB  THR C   1       7.618   3.584   9.284  1.00 10.00           C
ATOM      6  OG1 THR C   1       6.472   2.770   9.559  1.00 10.00           O
ATOM      7  O   HOH C   2       7.417   4.305   7.960  1.00 10.00           O
TER
"""

  ncs_groups = get_ncs_groups(phil_str1, pdb_str)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [0,1,2,3,4,5]
  assert list(ncs_groups[0].copies[0].iselection) == [14,15,16,17,18,19]

  # Here we modifying default exclude_selection parameter to keep waters
  ncs_groups = get_ncs_groups(phil_str1, pdb_str, exclude_selection=None)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [0,1,2,3,4,5,6]
  assert list(ncs_groups[0].copies[0].iselection) == [14,15,16,17,18,19,20]

def exercise_7():
  """
  extraction from exercise_00 form tst_ncs_input.py
  """
  pdb_str1 = """
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
  pdb_str2 = """
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
  phil_str1 = """
ncs_group {
  reference = chain A
  selection = chain B
  selection = chain C
}
"""
  ncs_groups = get_ncs_groups(phil_str1, pdb_str1)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [0,1,2,3]
  assert list(ncs_groups[0].copies[0].iselection) == [4,5,6,7]
  assert list(ncs_groups[0].copies[1].iselection) == [8,9,10,11]

  ncs_groups = get_ncs_groups(phil_str1, pdb_str2)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert list(ncs_groups[0].master_iselection) == [1,2,3,4]
  assert list(ncs_groups[0].copies[0].iselection) == [5,6,7,8]
  assert list(ncs_groups[0].copies[1].iselection) == [9,10,11,12]

def exercise_8():
  """
  Here user tries to supply non-valid group along with valid group.
  """
  pass

def exercise_9():
  """
  Here user's groups should be modified to exclude a residue or atom
  """
  pass

def exercise_10():
  """
  Include water if requested by user.
  Oleg. Not anymore. Just fixing
  Previously was exercise_18 in tst_ncs_input
  """
  pdb_str_15 = """
ATOM      1  N1  XXX A  34     125.208 211.886 175.417  1.00  0.00           N
ATOM      2  CT  XXX A  34     125.035 211.123 174.168  1.00  0.00           C
ATOM      3  C   XXX A  34     126.386 210.806 173.507  1.00  0.00           C
ATOM      4  K   XXX A  34     127.304 211.628 173.503  1.00  0.00           K
ATOM      1  O   HOH A  35     135.208 211.886 175.417  1.00  0.00           O
ATOM      2  O   HOH A  36     135.035 211.123 174.168  1.00  0.00           O
ATOM      3  O   HOH A  37     136.386 210.806 173.507  1.00  0.00           O
ATOM      4  O   HOH A  38     137.304 211.628 173.503  1.00  0.00           O
TER
ATOM      5  N1  XXX B  34     251.532 143.432 175.422  1.00  0.00           N
ATOM      6  CT  XXX B  34     252.120 143.948 174.173  1.00  0.00           C
ATOM      7  C   XXX B  34     251.212 144.998 173.512  1.00  0.00           C
ATOM      8  K   XXX B  34     249.986 144.872 173.510  1.00  0.00           K
ATOM      5  O   HOH B  35     271.532 143.432 175.422  1.00  0.00           O
ATOM      6  O   HOH B  36     272.120 143.948 174.173  1.00  0.00           O
ATOM      7  O   HOH B  37     271.212 144.998 173.512  1.00  0.00           O
ATOM      8  O   HOH B  38     279.986 144.872 173.510  1.00  0.00           O
TER
"""
  phil_str="""
ncs_group {
  reference = chain A
  selection = chain B
}
"""
  asc = iotbx.pdb.input(source_info=None,
    lines=pdb_str_15).construct_hierarchy().atom_selection_cache()
  ### user-supplied
  phil_groups = ncs_group_master_phil.fetch(
      iotbx.phil.parse(phil_str)).extract()
  ncs_inp = iotbx.ncs.input(pdb_string = pdb_str_15,
      ncs_phil_groups=phil_groups.ncs_group,
      exclude_selection=None)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups)==1
  # group 1
  assert ncs_groups[0].master_iselection.all_eq(
    asc.selection(string = "chain A").iselection())
  g1_c = ncs_groups[0].copies
  assert len(g1_c)==1
  assert g1_c[0].iselection.all_eq(
    asc.selection(string = "chain B").iselection())

  ncs_inp = iotbx.ncs.input(pdb_string = pdb_str_15,
      ncs_phil_groups=phil_groups.ncs_group,
      exclude_selection=None)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups)==1
  # group 1
  assert ncs_groups[0].master_iselection.all_eq(
      asc.selection(
          string = "chain A"
      ).iselection())
  g1_c = ncs_groups[0].copies
  assert len(g1_c)==1
  assert g1_c[0].iselection.all_eq(
      asc.selection(
          string = "chain B"
      ).iselection())

if (__name__ == "__main__"):
  t0=time.time()

  exercise_1()
  exercise_2()
  exercise_3()
  exercise_4()
  # exercise_5() #temp disabled. Need to tweak selection_string_from_selection
  # in atom_selection.py
  exercise_6()
  exercise_7()
  exercise_8()
  exercise_9()
  exercise_10()


  print "Time: %6.4f"%(time.time()-t0)
  print "OK"
