from __future__ import absolute_import, division, print_function
import time
from iotbx.ncs import ncs_group_master_phil
import iotbx.phil
import iotbx.ncs
import iotbx.pdb
from libtbx.test_utils import approx_equal

pdb_str_a_lot_of_ncs = """\
ATOM      0  N   VALAa 175       4.326 -59.034 -15.595  1.00105.80           N
ATOM      1  CA  VALAa 175       5.525 -58.206 -15.617  1.00 86.22           C
ATOM      2  C   VALAa 175       6.348 -58.476 -14.358  1.00 61.80           C
ATOM      3  O   VALAa 175       5.879 -59.158 -13.445  1.00 58.55           O
TER
ATOM      4  N   THRAb 174     -15.180 -65.753 -17.549  1.00113.31           N
ATOM      5  CA  THRAb 174     -14.857 -64.334 -17.581  1.00120.94           C
ATOM      6  C   THRAb 174     -15.877 -63.524 -16.808  1.00109.70           C
ATOM      7  O   THRAb 174     -15.828 -62.302 -16.798  1.00110.07           O
TER
ATOM      8  N   LYSAc  30     -75.643   8.901 124.736  1.00205.41           N
ATOM      9  CA  LYSAc  30     -74.928  10.142 125.003  1.00200.32           C
ATOM     10  C   LYSAc  30     -75.558  11.283 124.229  1.00198.69           C
ATOM     11  O   LYSAc  30     -74.919  11.832 123.330  1.00194.89           O
ATOM     12  N   THRAc  31     -76.807  11.640 124.555  1.00163.49           N
ATOM     13  CA  THRAc  31     -77.445  12.728 123.826  1.00162.32           C
ATOM     14  C   THRAc  31     -76.491  13.902 123.640  1.00160.35           C
ATOM     15  O   THRAc  31     -76.874  15.060 123.806  1.00158.76           O
TER
ATOM     16  N   VALAd 175     148.055  93.347 -15.595  1.00105.80           N
ATOM     17  CA  VALAd 175     147.227  94.546 -15.617  1.00 86.22           C
ATOM     18  C   VALAd 175     147.496  95.369 -14.358  1.00 61.80           C
ATOM     19  O   VALAd 175     148.179  94.900 -13.445  1.00 58.55           O
TER
ATOM     20  N   THRAe 174     154.774  73.841 -17.549  1.00113.31           N
ATOM     21  CA  THRAe 174     153.355  74.164 -17.581  1.00120.94           C
ATOM     22  C   THRAe 174     152.544  73.144 -16.808  1.00109.70           C
ATOM     23  O   THRAe 174     151.323  73.193 -16.798  1.00110.07           O
TER
ATOM     24  N   LYSAf  30      80.120  13.378 124.736  1.00205.41           N
ATOM     25  CA  LYSAf  30      78.879  14.093 125.003  1.00200.32           C
ATOM     26  C   LYSAf  30      77.738  13.463 124.229  1.00198.69           C
ATOM     27  O   LYSAf  30      77.189  14.102 123.330  1.00194.89           O
ATOM     28  N   THRAf  31      77.381  12.214 124.555  1.00163.49           N
ATOM     29  CA  THRAf  31      76.293  11.576 123.826  1.00162.32           C
ATOM     30  C   THRAf  31      75.119  12.530 123.640  1.00160.35           C
ATOM     31  O   THRAf  31      73.961  12.147 123.806  1.00158.76           O
TER
ATOM     32  N   VALAg 175      29.987  84.695 -15.595  1.00105.80           N
ATOM     33  CA  VALAg 175      30.815  83.496 -15.617  1.00 86.22           C
ATOM     34  C   VALAg 175      30.545  82.672 -14.358  1.00 61.80           C
ATOM     35  O   VALAg 175      29.863  83.142 -13.445  1.00 58.55           O
TER
ATOM     36  N   THRAh 174      23.268 104.201 -17.549  1.00113.31           N
ATOM     37  CA  THRAh 174      24.687 103.878 -17.581  1.00120.94           C
ATOM     38  C   THRAh 174      25.497 104.898 -16.808  1.00109.70           C
ATOM     39  O   THRAh 174      26.719 104.849 -16.798  1.00110.07           O
TER
ATOM     40  N   LYSAi  30      97.922 164.664 124.736  1.00205.41           N
ATOM     41  CA  LYSAi  30      99.163 163.949 125.003  1.00200.32           C
ATOM     42  C   LYSAi  30     100.304 164.579 124.229  1.00198.69           C
ATOM     43  O   LYSAi  30     100.853 163.940 123.330  1.00194.89           O
ATOM     44  N   THRAi  31     100.661 165.828 124.555  1.00163.49           N
ATOM     45  CA  THRAi  31     101.749 166.466 123.826  1.00162.32           C
ATOM     46  C   THRAi  31     102.923 165.512 123.640  1.00160.35           C
ATOM     47  O   THRAi  31     104.081 165.895 123.806  1.00158.76           O
TER
ATOM     48  N   VALAj 175      93.347 148.055  15.595  1.00105.80           N
ATOM     49  CA  VALAj 175      94.546 147.227  15.617  1.00 86.22           C
ATOM     50  C   VALAj 175      95.369 147.496  14.358  1.00 61.80           C
ATOM     51  O   VALAj 175      94.900 148.179  13.445  1.00 58.55           O
TER
ATOM     52  N   THRAk 174      73.841 154.774  17.549  1.00113.31           N
ATOM     53  CA  THRAk 174      74.164 153.355  17.581  1.00120.94           C
ATOM     54  C   THRAk 174      73.144 152.544  16.808  1.00109.70           C
ATOM     55  O   THRAk 174      73.193 151.323  16.798  1.00110.07           O
TER
ATOM     56  N   LYSAl  30      13.378  80.120-124.736  1.00205.41           N
ATOM     57  CA  LYSAl  30      14.093  78.879-125.003  1.00200.32           C
ATOM     58  C   LYSAl  30      13.462  77.738-124.229  1.00198.69           C
ATOM     59  O   LYSAl  30      14.102  77.189-123.330  1.00194.89           O
ATOM     60  N   THRAl  31      12.213  77.381-124.555  1.00163.49           N
ATOM     61  CA  THRAl  31      11.576  76.293-123.826  1.00162.32           C
ATOM     62  C   THRAl  31      12.530  75.119-123.640  1.00160.35           C
ATOM     63  O   THRAl  31      12.147  73.961-123.806  1.00158.76           O
TER
ATOM     64  N   VALAm 175      84.695  29.987  15.595  1.00105.80           N
ATOM     65  CA  VALAm 175      83.496  30.815  15.617  1.00 86.22           C
ATOM     66  C   VALAm 175      82.672  30.545  14.358  1.00 61.80           C
ATOM     67  O   VALAm 175      83.142  29.863  13.445  1.00 58.55           O
TER
ATOM     68  N   THRAn 174     104.201  23.268  17.549  1.00113.31           N
ATOM     69  CA  THRAn 174     103.878  24.687  17.581  1.00120.94           C
ATOM     70  C   THRAn 174     104.898  25.497  16.808  1.00109.70           C
ATOM     71  O   THRAn 174     104.849  26.719  16.798  1.00110.07           O
TER
ATOM     72  N   LYSAo  30     164.664  97.921-124.736  1.00205.41           N
ATOM     73  CA  LYSAo  30     163.948  99.162-125.003  1.00200.32           C
ATOM     74  C   LYSAo  30     164.579 100.303-124.229  1.00198.69           C
ATOM     75  O   LYSAo  30     163.940 100.852-123.330  1.00194.89           O
ATOM     76  N   THRAo  31     165.828 100.660-124.555  1.00163.49           N
ATOM     77  CA  THRAo  31     166.466 101.748-123.826  1.00162.32           C
ATOM     78  C   THRAo  31     165.512 102.922-123.640  1.00160.35           C
ATOM     79  O   THRAo  31     165.895 104.080-123.806  1.00158.76           O
TER
ATOM     80  N   VALAp 175      -4.326  59.034 -15.595  1.00105.80           N
ATOM     81  CA  VALAp 175      -5.525  58.206 -15.617  1.00 86.22           C
ATOM     82  C   VALAp 175      -6.348  58.476 -14.358  1.00 61.80           C
ATOM     83  O   VALAp 175      -5.879  59.158 -13.445  1.00 58.55           O
TER
ATOM     84  N   THRAq 174      15.180  65.753 -17.549  1.00113.31           N
ATOM     85  CA  THRAq 174      14.857  64.334 -17.581  1.00120.94           C
ATOM     86  C   THRAq 174      15.877  63.524 -16.808  1.00109.70           C
ATOM     87  O   THRAq 174      15.828  62.302 -16.798  1.00110.07           O
TER
ATOM     88  N   LYSAr  30      75.643  -8.901 124.736  1.00205.41           N
ATOM     89  CA  LYSAr  30      74.928 -10.142 125.003  1.00200.32           C
ATOM     90  C   LYSAr  30      75.558 -11.283 124.229  1.00198.69           C
ATOM     91  O   LYSAr  30      74.919 -11.832 123.330  1.00194.89           O
ATOM     92  N   THRAr  31      76.807 -11.640 124.555  1.00163.49           N
ATOM     93  CA  THRAr  31      77.445 -12.728 123.826  1.00162.32           C
ATOM     94  C   THRAr  31      76.491 -13.902 123.640  1.00160.35           C
ATOM     95  O   THRAr  31      76.874 -15.060 123.806  1.00158.76           O
TER
ATOM     96  N   VALAs 175     -59.034   4.326  15.595  1.00105.80           N
ATOM     97  CA  VALAs 175     -58.206   5.525  15.617  1.00 86.22           C
ATOM     98  C   VALAs 175     -58.476   6.348  14.358  1.00 61.80           C
ATOM     99  O   VALAs 175     -59.158   5.879  13.445  1.00 58.55           O
TER
ATOM    100  N   THRAt 174     -65.753 -15.180  17.549  1.00113.31           N
ATOM    101  CA  THRAt 174     -64.334 -14.857  17.581  1.00120.94           C
ATOM    102  C   THRAt 174     -63.524 -15.877  16.808  1.00109.70           C
ATOM    103  O   THRAt 174     -62.302 -15.828  16.798  1.00110.07           O
TER
ATOM    104  N   LYSAu  30       8.901 -75.643-124.736  1.00205.41           N
ATOM    105  CA  LYSAu  30      10.142 -74.928-125.003  1.00200.32           C
ATOM    106  C   LYSAu  30      11.283 -75.558-124.229  1.00198.69           C
ATOM    107  O   LYSAu  30      11.832 -74.919-123.330  1.00194.89           O
ATOM    108  N   THRAu  31      11.640 -76.807-124.555  1.00163.49           N
ATOM    109  CA  THRAu  31      12.728 -77.445-123.826  1.00162.32           C
ATOM    110  C   THRAu  31      13.902 -76.491-123.640  1.00160.35           C
ATOM    111  O   THRAu  31      15.060 -76.874-123.806  1.00158.76           O
TER
ATOM    112  N   VALAv 175      59.034  -4.326  15.595  1.00105.80           N
ATOM    113  CA  VALAv 175      58.206  -5.525  15.617  1.00 86.22           C
ATOM    114  C   VALAv 175      58.476  -6.348  14.358  1.00 61.80           C
ATOM    115  O   VALAv 175      59.158  -5.879  13.445  1.00 58.55           O
TER
ATOM    116  N   THRAw 174      65.753  15.180  17.549  1.00113.31           N
ATOM    117  CA  THRAw 174      64.334  14.857  17.581  1.00120.94           C
ATOM    118  C   THRAw 174      63.524  15.877  16.808  1.00109.70           C
ATOM    119  O   THRAw 174      62.302  15.828  16.798  1.00110.07           O
TER
ATOM    120  N   LYSAx  30      -8.901  75.643-124.736  1.00205.41           N
ATOM    121  CA  LYSAx  30     -10.142  74.928-125.003  1.00200.32           C
ATOM    122  C   LYSAx  30     -11.283  75.558-124.229  1.00198.69           C
ATOM    123  O   LYSAx  30     -11.832  74.919-123.330  1.00194.89           O
ATOM    124  N   THRAx  31     -11.640  76.807-124.555  1.00163.49           N
ATOM    125  CA  THRAx  31     -12.728  77.445-123.826  1.00162.32           C
ATOM    126  C   THRAx  31     -13.902  76.491-123.640  1.00160.35           C
ATOM    127  O   THRAx  31     -15.060  76.874-123.806  1.00158.76           O
TER
"""

def get_ncs_groups(phil_str, pdb_str, **kwargs):
  phil_groups = ncs_group_master_phil.fetch(
      iotbx.phil.parse(phil_str)).extract()
  pdb_inp = iotbx.pdb.input(lines=pdb_str,source_info=None)
  ncs_inp = iotbx.ncs.input(hierarchy = pdb_inp.construct_hierarchy(),
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
  p = iotbx.ncs.input.get_default_params()
  p.ncs_search.exclude_selection=None
  ncs_groups = get_ncs_groups(phil_str1, pdb_str, params=p.ncs_search)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1
  assert list(ncs_groups[0].master_iselection) == [0,1,2,3,4,5,6], list(ncs_groups[0].master_iselection)
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
  Full list of valid NCS groups is in exercise_11, phil_str1.
  """
  phil_str1 = """\
ncs_group {
  reference = chain Aa
  selection = chain Ae
}
"""
  phil_str2 = """\
ncs_group {
  reference = chain Aa
  selection = chain Ae
}
ncs_group {
  reference = chain As
  selection = chain Av
}
"""
  pdb_h = iotbx.pdb.input(
      source_info=None,lines=pdb_str_a_lot_of_ncs.split('\n')).\
      construct_hierarchy()
  # show file:
  # pdb_h.write_pdb_file("user_ncs_phil_ex11.pdb")
  asc = pdb_h.atom_selection_cache()

  ncs_groups = get_ncs_groups(phil_str1, pdb_str_a_lot_of_ncs)
  assert len(ncs_groups)==0

  # One group is ok
  ncs_groups = get_ncs_groups(phil_str2, pdb_str_a_lot_of_ncs)
  assert len(ncs_groups)==1
  assert len(ncs_groups[0].copies) == 1
  assert approx_equal(ncs_groups[0].master_iselection,
      asc.iselection("chain As"))
  assert approx_equal(ncs_groups[0].copies[0].iselection,
      asc.iselection("chain Av"))


def exercise_9():
  """
  Here user's groups should be modified to exclude a residue or atom
  pdb_str is has good NCS, remarked atoms shows introduced discrepancies
  """
  pdb_str = """\
ATOM      1  CB  MET--   1      52.886   1.976   9.011  1.00 41.44           C
ATOM      2  CG  MET--   1      53.271   0.996  10.102  1.00 47.36           C
ATOM      3  SD  MET--   1      55.066   1.105  10.320  1.00 54.87           S
ATOM      4  CE  MET--   1      55.198   1.758  11.963  1.00 52.95           C
ATOM      5  C   MET--   1      51.257   3.002   7.473  1.00 32.90           C
ATOM      6  O   MET--   1      51.643   3.028   6.297  1.00 33.11           O
ATOM      7  N   MET--   1      51.320   0.535   7.805  1.00 38.36           N
ATOM      8  CA  MET--   1      51.486   1.834   8.434  1.00 36.79           C
ATOM      9  N   ILE--   2      50.621   4.010   8.052  1.00 27.46           N
ATOM     10  CA  ILE--   2      50.330   5.240   7.341  1.00 22.07           C
ATOM     11  CB  ILE--   2      49.096   5.952   7.960  1.00 23.08           C
ATOM     12  CG2 ILE--   2      48.824   7.256   7.244  1.00 21.31           C
remark ATOM     13  CG1 ILE--   2      47.841   5.109   7.800  1.00 23.53           C
ATOM     14  CD1 ILE--   2      46.598   5.779   8.421  1.00 26.20           C
ATOM     15  C   ILE--   2      51.562   6.118   7.509  1.00 19.91           C
ATOM     16  O   ILE--   2      51.955   6.341   8.645  1.00 17.98           O
ATOM     17  N   LYS--   3      52.103   6.642   6.398  1.00 18.18           N
ATOM     18  CA  LYS--   3      53.229   7.583   6.418  1.00 17.08           C
ATOM     19  CB  LYS--   3      54.076   7.396   5.144  1.00 19.78           C
ATOM     20  CG  LYS--   3      55.276   8.374   5.088  1.00 23.29           C
ATOM     21  CD  LYS--   3      55.948   8.330   3.717  1.00 29.04           C
ATOM     22  CE  LYS--   3      57.288   7.627   3.791  1.00 31.30           C
ATOM     23  NZ  LYS--   3      57.822   7.350   2.460  1.00 34.23           N
ATOM     24  C   LYS--   3      52.675   9.012   6.481  1.00 14.59           C
ATOM     25  O   LYS--   3      51.877   9.404   5.627  1.00 15.75           O
ATOM     26  N   VAL--   4      53.044   9.811   7.492  1.00 12.13           N
ATOM     27  CA  VAL--   4      52.594  11.201   7.600  1.00 10.72           C
ATOM     28  CB  VAL--   4      51.792  11.417   8.905  1.00 12.27           C
ATOM     29  CG1 VAL--   4      51.313  12.867   9.013  1.00 11.73           C
ATOM     30  CG2 VAL--   4      50.586  10.433   8.913  1.00 12.18           C
ATOM     31  C   VAL--   4      53.878  12.034   7.643  1.00 11.26           C
ATOM     32  O   VAL--   4      54.859  11.624   8.278  1.00 12.24           O
ATOM     33  N   GLU--   5      53.857  13.160   6.937  1.00  9.99           N
ATOM     34  CA  GLU--   5      55.027  14.050   6.987  1.00 11.22           C
ATOM     35  CB  GLU--   5      55.705  14.042   5.622  1.00 12.79           C
ATOM     36  CG  GLU--   5      56.999  14.908   5.598  1.00 16.55           C
ATOM     37  CD  GLU--   5      57.945  14.693   4.429  1.00 18.73           C
ATOM     38  OE1 GLU--   5      57.591  14.045   3.458  1.00 18.89           O
ATOM     39  OE2 GLU--   5      59.066  15.182   4.506  1.00 21.35           O
ATOM     40  C   GLU--   5      54.642  15.472   7.371  1.00 11.49           C
ATOM     41  O   GLU--   5      53.646  15.977   6.854  1.00 12.69           O
TER
ATOM      1  CB  MET G   1       0.115  20.271   3.107  1.00 41.44
ATOM      2  CG  MET G   1       0.717  21.085   1.978  1.00 47.36
ATOM      3  SD  MET G   1       1.424  22.584   2.709  1.00 54.87
ATOM      4  CE  MET G   1       3.156  22.361   2.406  1.00 52.95
ATOM      5  C   MET G   1      -1.177  18.359   3.971  1.00 32.90
ATOM      6  O   MET G   1      -2.072  18.689   4.760  1.00 33.11
ATOM      7  N   MET G   1      -1.968  19.642   1.989  1.00 38.36
ATOM      8  CA  MET G   1      -0.810  19.133   2.704  1.00 36.79
ATOM      9  N   ILE G   2      -0.401  17.300   4.154  1.00 27.46
ATOM     10  CA  ILE G   2      -0.542  16.440   5.313  1.00 22.07
ATOM     11  CB  ILE G   2      -0.020  15.011   5.001  1.00 23.08
ATOM     12  CG2 ILE G   2      -0.128  14.130   6.226  1.00 21.31
ATOM     13  CG1 ILE G   2      -0.860  14.345   3.923  1.00 23.53
ATOM     14  CD1 ILE G   2      -0.357  12.928   3.575  1.00 26.20
ATOM     15  C   ILE G   2       0.314  17.068   6.404  1.00 19.91
ATOM     16  O   ILE G   2       1.490  17.289   6.152  1.00 17.98
ATOM     17  N   LYS G   3      -0.264  17.284   7.595  1.00 18.18
ATOM     18  CA  LYS G   3       0.465  17.790   8.764  1.00 17.08
ATOM     19  CB  LYS G   3      -0.488  18.627   9.640  1.00 19.78
ATOM     20  CG  LYS G   3       0.212  19.180  10.906  1.00 23.29
ATOM     21  CD  LYS G   3      -0.806  19.795  11.864  1.00 29.04
ATOM     22  CE  LYS G   3      -0.702  21.306  11.871  1.00 31.30
ATOM     23  NZ  LYS G   3      -1.826  21.917  12.577  1.00 34.23
ATOM     24  C   LYS G   3       1.006  16.597   9.561  1.00 14.59
ATOM     25  O   LYS G   3       0.239  15.716   9.954  1.00 15.75
ATOM     26  N   VAL G   4       2.324  16.510   9.792  1.00 12.13
ATOM     27  CA  VAL G   4       2.914  15.426  10.581  1.00 10.72
ATOM     28  CB  VAL G   4       3.918  14.613   9.730  1.00 12.27
remark ATOM     29  CG1 VAL G   4       4.527  13.474  10.552  1.00 11.73
ATOM     30  CG2 VAL G   4       3.172  14.059   8.481  1.00 12.18
ATOM     31  C   VAL G   4       3.657  16.123  11.725  1.00 11.26
ATOM     32  O   VAL G   4       4.277  17.173  11.509  1.00 12.24
ATOM     33  N   GLU G   5       3.545  15.548  12.918  1.00  9.99
ATOM     34  CA  GLU G   5       4.289  16.118  14.052  1.00 11.22
ATOM     35  CB  GLU G   5       3.293  16.720  15.036  1.00 12.79
ATOM     36  CG  GLU G   5       3.995  17.409  16.243  1.00 16.55
ATOM     37  CD  GLU G   5       3.146  18.345  17.086  1.00 18.73
ATOM     38  OE1 GLU G   5       1.934  18.369  16.952  1.00 18.89
ATOM     39  OE2 GLU G   5       3.723  19.072  17.886  1.00 21.35
ATOM     40  C   GLU G   5       5.146  15.072  14.750  1.00 11.49
ATOM     41  O   GLU G   5       4.666  13.961  14.970  1.00 12.69
"""
  phil_str1 = """\
ncs_group {
  reference = chain --
  selection = chain G
}
"""
  phil_str2 = """\
ncs_group {
  reference = chain G
  selection = chain '--'
}
"""

  pdb_h = iotbx.pdb.input(
      source_info=None,lines=pdb_str.split('\n')).\
      construct_hierarchy()
  # show file:
  # pdb_h.write_pdb_file("user_ncs_phil_ex9.pdb")
  asc = pdb_h.atom_selection_cache()
  p = iotbx.ncs.input.get_default_params()
  p.ncs_search.exclude_selection = None
  ncs_groups = get_ncs_groups(phil_str1,pdb_str, params=p.ncs_search)
  assert len(ncs_groups)==1
  assert len(ncs_groups[0].copies) == 1
  assert approx_equal(ncs_groups[0].master_iselection,
      asc.iselection("chain '--' and not ((resid 2 or resid 4) and name CG1)"))
  assert approx_equal(ncs_groups[0].copies[0].iselection,
      asc.iselection("chain G and not ((resid 2 or resid 4) and name CG1)"))

  ncs_groups = get_ncs_groups(phil_str2,pdb_str, params=p.ncs_search)
  assert len(ncs_groups)==1
  assert len(ncs_groups[0].copies) == 1
  assert approx_equal(ncs_groups[0].master_iselection,
      asc.iselection("chain G and not ((resid 2 or resid 4) and name CG1)"))
  assert approx_equal(ncs_groups[0].copies[0].iselection,
      asc.iselection("chain '--' and not ((resid 2 or resid 4) and name CG1)"))


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
  ### user-supplied
  p = iotbx.ncs.input.get_default_params()
  p.ncs_search.exclude_selection = None

  ncs_groups = get_ncs_groups(phil_str,pdb_str_15, params=p.ncs_search)
  # group 1
  assert len(ncs_groups)==1
  assert len(ncs_groups[0].copies) == 1
  # chain A
  assert list(ncs_groups[0].master_iselection) == [0, 1, 2, 3, 4, 5, 6, 7]
  # chain B
  assert list(ncs_groups[0].copies[0].iselection) == [8,9,10,11,12,13,14,15]

def exercise_11():
  """
  Two-letter chain ids. Also testing user selections for multiple chains.
  """

  phil_str1 = """\
ncs_group {
  reference        = chain Aa or chain Ab or chain Ac
  selection        = chain Ad or chain Ae or chain Af
  selection        = chain Ag or chain Ah or chain Ai
  selection        = chain Aj or chain Ak or chain Al
  selection        = chain Am or chain An or chain Ao
  selection        = chain Ap or chain Aq or chain Ar
  selection        = chain As or chain At or chain Au
  selection        = chain Av or chain Aw or chain Ax
}"""

  phil_str2 = """\
ncs_group {
  reference        = chain Aa or chain Ab
  selection        = chain Ad or chain Ae
  selection        = chain Ag or chain Ah
  selection        = chain Aj or chain Ak
  selection        = chain Am or chain An
  selection        = chain Ap or chain Aq
  selection        = chain As or chain At
  selection        = chain Av or chain Aw
}"""

  phil_str3 = """\
ncs_group {
  reference        = chain Aj
  selection        = chain Am
  selection        = chain Ap
  }
"""
  phil_str4 = """\
ncs_group {
  reference        = chain Aj
  selection        = chain Am
  selection        = chain Ap
  }
ncs_group {
  reference        = chain Aa
  selection        = chain Ad
  selection        = chain Av
  }
"""

  pdb_h = iotbx.pdb.input(
      source_info=None,lines=pdb_str_a_lot_of_ncs.split('\n')).\
      construct_hierarchy()
  # show file:
  # pdb_h.write_pdb_file("user_ncs_phil_ex11.pdb")
  asc = pdb_h.atom_selection_cache()

  ncs_groups = get_ncs_groups(phil_str1, pdb_str_a_lot_of_ncs)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 7
  assert ncs_groups[0].master_iselection == asc.iselection(
      "chain Aa or chain Ab or chain Ac")
  for i, s_str in enumerate([
      "chain Ad or chain Ae or chain Af",
      "chain Ag or chain Ah or chain Ai",
      "chain Aj or chain Ak or chain Al",
      "chain Am or chain An or chain Ao",
      "chain Ap or chain Aq or chain Ar",
      "chain As or chain At or chain Au",
      "chain Av or chain Aw or chain Ax"]):
    assert approx_equal(ncs_groups[0].copies[i].iselection,
        asc.iselection(s_str))

  # Now we want only 2 chains from each selection
  ncs_groups = get_ncs_groups(phil_str2, pdb_str_a_lot_of_ncs)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 7
  assert ncs_groups[0].master_iselection == asc.iselection(
      "chain Aa or chain Ab")
  for i, s_str in enumerate([
      "chain Ad or chain Ae",
      "chain Ag or chain Ah",
      "chain Aj or chain Ak",
      "chain Am or chain An",
      "chain Ap or chain Aq",
      "chain As or chain At",
      "chain Av or chain Aw"]):
    assert approx_equal(ncs_groups[0].copies[i].iselection,
        asc.iselection(s_str))

  # User is greedy
  ncs_groups = get_ncs_groups(phil_str3, pdb_str_a_lot_of_ncs)
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 2
  assert approx_equal(ncs_groups[0].master_iselection,
      asc.iselection("chain Aj"))
  assert approx_equal(ncs_groups[0].copies[0].iselection,
      asc.iselection("chain Am"))
  assert approx_equal(ncs_groups[0].copies[1].iselection,
      asc.iselection("chain Ap"))
  # Two groups
  ncs_groups = get_ncs_groups(phil_str4, pdb_str_a_lot_of_ncs)
  assert len(ncs_groups) == 2
  assert len(ncs_groups[0].copies) == 2
  assert len(ncs_groups[1].copies) == 2
  assert approx_equal(ncs_groups[0].master_iselection,
      asc.iselection("chain Aj"))
  assert approx_equal(ncs_groups[0].copies[0].iselection,
      asc.iselection("chain Am"))
  assert approx_equal(ncs_groups[0].copies[1].iselection,
      asc.iselection("chain Ap"))
  assert approx_equal(ncs_groups[1].master_iselection,
      asc.iselection("chain Aa"))
  assert approx_equal(ncs_groups[1].copies[0].iselection,
      asc.iselection("chain Ad"))
  assert approx_equal(ncs_groups[1].copies[1].iselection,
      asc.iselection("chain Av"))

def exercise_12():
  """
  User's ncs_group is separated into several groups by validation procedure.
  """
  pdb_str = """\
CRYST1  595.620  595.620  595.620  90.00  90.00  90.00 I 2 3
SCALE1      0.001679  0.000000  0.000000        0.00000
SCALE2      0.000000  0.001679  0.000000        0.00000
SCALE3      0.000000  0.000000  0.001679        0.00000
ATOM    266  N   LEU 1  20     -63.533-187.068-206.553  1.00 47.37           N
ATOM    267  CA  LEU 1  20     -64.594-187.999-206.894  1.00 38.42           C
ATOM    268  C   LEU 1  20     -64.792-189.093-205.863  1.00 40.45           C
ATOM    269  O   LEU 1  20     -65.844-189.739-205.864  1.00 37.26           O
ATOM    270  CB  LEU 1  20     -64.312-188.680-208.230  1.00 38.09           C
ATOM    271  CG  LEU 1  20     -63.920-187.774-209.381  1.00 40.19           C
ATOM    272  CD1 LEU 1  20     -63.950-188.678-210.596  1.00 39.43           C
ATOM    273  CD2 LEU 1  20     -64.787-186.534-209.581  1.00 33.28           C
ATOM    285  N   THR 1  21     -63.827-189.307-204.987  1.00 45.58           N
ATOM    286  CA  THR 1  21     -63.813-190.421-204.050  1.00 39.34           C
ATOM    287  C   THR 1  21     -63.827-189.887-202.615  1.00 39.90           C
ATOM    288  O   THR 1  21     -63.742-188.677-202.358  1.00 36.18           O
ATOM    289  CB  THR 1  21     -62.584-191.313-204.320  1.00 49.62           C
ATOM    290  OG1 THR 1  21     -61.390-190.515-204.234  1.00 52.31           O
ATOM    291  CG2 THR 1  21     -62.664-191.988-205.729  1.00 38.87           C
ATOM    299  N   ARG 1  22     -63.976-190.808-201.669  1.00 47.15           N
ATOM    300  CA  ARG 1  22     -64.297-190.439-200.299  1.00 45.35           C
ATOM    301  C   ARG 1  22     -63.888-191.595-199.403  1.00 44.05           C
ATOM    302  O   ARG 1  22     -64.251-192.746-199.674  1.00 47.23           O
ATOM    303  CB  ARG 1  22     -65.799-190.145-200.267  1.00 50.32           C
ATOM    304  CG  ARG 1  22     -66.470-189.595-199.030  1.00 48.26           C
ATOM    305  CD  ARG 1  22     -67.916-189.242-199.482  1.00 71.42           C
ATOM    306  NE  ARG 1  22     -68.807-188.659-198.474  1.00 89.96           N
ATOM    307  CZ  ARG 1  22     -70.012-188.139-198.752  1.00 89.76           C
ATOM    308  NH1 ARG 1  22     -70.490-188.125-200.004  1.00 72.83           N
ATOM    309  NH2 ARG 1  22     -70.755-187.629-197.765  1.00 72.15           N
ATOM    323  N   ALA 1  23     -63.090-191.310-198.380  1.00 31.79           N
ATOM    324  CA  ALA 1  23     -62.661-192.357-197.458  1.00 37.78           C
ATOM    325  C   ALA 1  23     -63.777-192.758-196.478  1.00 41.29           C
ATOM    326  O   ALA 1  23     -64.556-191.929-195.997  1.00 41.52           O
ATOM    327  CB  ALA 1  23     -61.432-191.902-196.679  1.00 34.84           C
ATOM    333  N   LEU 1  24     -63.865-194.045-196.194  1.00 37.34           N
ATOM    334  CA  LEU 1  24     -64.861-194.553-195.277  1.00 27.60           C
ATOM    335  C   LEU 1  24     -64.182-195.423-194.245  1.00 35.82           C
ATOM    336  O   LEU 1  24     -63.148-196.034-194.520  1.00 42.04           O
ATOM    337  CB  LEU 1  24     -65.902-195.408-195.995  1.00 41.37           C
ATOM    338  CG  LEU 1  24     -66.637-194.703-197.121  1.00 45.57           C
ATOM    339  CD1 LEU 1  24     -67.524-195.709-197.954  1.00 28.75           C
ATOM    340  CD2 LEU 1  24     -67.430-193.557-196.482  1.00 26.24           C
ATOM    352  N   PRO 1  25     -64.769-195.542-193.072  1.00 43.07           N
ATOM    353  CA  PRO 1  25     -64.203-196.455-192.080  1.00 45.02           C
ATOM    354  C   PRO 1  25     -64.368-197.875-192.575  1.00 28.49           C
ATOM    355  O   PRO 1  25     -65.337-198.184-193.276  1.00 23.83           O
ATOM    356  CB  PRO 1  25     -65.031-196.189-190.827  1.00 42.63           C
ATOM    357  CG  PRO 1  25     -66.337-195.696-191.405  1.00 34.61           C
ATOM    358  CD  PRO 1  25     -65.960-194.851-192.567  1.00 25.43           C
TER
ATOM  A0BF2  N   THRGA  20     -52.224 -82.892-121.524  1.00 32.92           N
ATOM  A0BF3  CA  THRGA  20     -53.561 -83.049-120.966  1.00 35.79           C
ATOM  A0BF4  C   THRGA  20     -54.347 -84.185-121.629  1.00 43.91           C
ATOM  A0BF5  O   THRGA  20     -55.534 -84.024-121.903  1.00 55.90           O
ATOM  A0BF6  CB  THRGA  20     -54.363 -81.742-121.088  1.00 32.17           C
ATOM  A0BF7  OG1 THRGA  20     -54.607 -81.451-122.480  1.00 27.96           O
ATOM  A0BF8  CG2 THRGA  20     -53.683 -80.551-120.374  1.00 31.79           C
ATOM  A0BFG  N   GLUGA  21     -53.739 -85.326-121.965  1.00 56.34           N
ATOM  A0BFH  CA  GLUGA  21     -54.586 -86.501-122.243  1.00 51.10           C
ATOM  A0BFI  C   GLUGA  21     -54.919 -87.215-120.943  1.00 54.73           C
ATOM  A0BFJ  O   GLUGA  21     -54.108 -87.260-120.010  1.00 42.85           O
ATOM  A0BFK  CB  GLUGA  21     -54.002 -87.549-123.220  1.00 30.90           C
ATOM  A0BFL  CG  GLUGA  21     -54.033 -87.190-124.724  1.00 50.79           C
ATOM  A0BFM  CD  GLUGA  21     -52.748 -86.608-125.203  1.00 85.38           C
ATOM  A0BFN  OE1 GLUGA  21     -52.055 -85.975-124.369  1.00 95.33           O
ATOM  A0BFO  OE2 GLUGA  21     -52.502 -86.706-126.435  1.00 76.32           O
ATOM  A0BFV  N   GLYGA  22     -56.107 -87.813-120.913  1.00 62.38           N
ATOM  A0BFW  CA  GLYGA  22     -56.653 -88.413-119.710  1.00 46.07           C
ATOM  A0BFX  C   GLYGA  22     -57.080 -87.396-118.680  1.00 65.93           C
ATOM  A0BFY  O   GLYGA  22     -56.984 -87.663-117.475  1.00 63.56           O
ATOM  A0BG2  N   SERGA  23     -57.542 -86.227-119.121  1.00 75.27           N
ATOM  A0BG3  CA  SERGA  23     -57.933 -85.143-118.234  1.00 63.26           C
ATOM  A0BG4  C   SERGA  23     -59.451 -85.070-118.142  1.00 85.09           C
ATOM  A0BG5  O   SERGA  23     -60.193 -85.512-119.030  1.00 74.26           O
ATOM  A0BG6  CB  SERGA  23     -57.360 -83.796-118.689  1.00 62.39           C
ATOM  A0BG7  OG  SERGA  23     -55.986 -83.921-119.034  1.00 78.01           O
ATOM  A0BGD  N   THRGA  24     -59.887 -84.512-117.033  1.00103.11           N
ATOM  A0BGE  CA  THRGA  24     -61.286 -84.324-116.713  1.00 96.25           C
ATOM  A0BGF  C   THRGA  24     -61.901 -83.158-117.487  1.00 95.46           C
ATOM  A0BGG  O   THRGA  24     -63.067 -82.845-117.236  1.00 94.00           O
ATOM  A0BGH  CB  THRGA  24     -61.420 -84.124-115.185  1.00 96.30           C
ATOM  A0BGI  OG1 THRGA  24     -62.797 -83.996-114.836  1.00102.74           O
ATOM  A0BGJ  CG2 THRGA  24     -60.614 -82.887-114.622  1.00 80.83           C
ATOM  A0BGR  N   ILEGA  25     -61.158 -82.541-118.421  1.00 92.73           N
ATOM  A0BGS  CA  ILEGA  25     -61.639 -81.451-119.271  1.00 85.74           C
ATOM  A0BGT  C   ILEGA  25     -61.603 -81.911-120.740  1.00 81.39           C
ATOM  A0BGU  O   ILEGA  25     -60.740 -82.708-121.136  1.00 75.57           O
ATOM  A0BGV  CB  ILEGA  25     -60.844 -80.121-119.029  1.00 78.83           C
ATOM  A0BGW  CG1 ILEGA  25     -60.850 -79.715-117.527  1.00 81.73           C
ATOM  A0BGX  CG2 ILEGA  25     -61.406 -78.934-119.840  1.00 76.17           C
ATOM  A0BGY  CD1 ILEGA  25     -59.934 -78.530-117.162  1.00 70.88           C
TER
ATOM  A1OYJ  N   THRAB  20     -34.601-167.170-165.019  1.00 47.42           N
ATOM  A1OYK  CA  THRAB  20     -35.664-166.423-165.675  1.00 38.19           C
ATOM  A1OYL  C   THRAB  20     -36.920-166.289-164.808  1.00 43.55           C
ATOM  A1OYM  O   THRAB  20     -38.028-166.342-165.336  1.00 58.60           O
ATOM  A1OYN  CB  THRAB  20     -36.026-167.082-167.013  1.00 37.29           C
ATOM  A1OYO  OG1 THRAB  20     -36.609-168.396-166.811  1.00 34.65           O
ATOM  A1OYP  CG2 THRAB  20     -34.819-167.200-167.905  1.00 47.27           C
ATOM  A1OYX  N   GLUAB  21     -36.809-166.114-163.487  1.00 40.30           N
ATOM  A1OYY  CA  GLUAB  21     -37.997-165.653-162.772  1.00 40.72           C
ATOM  A1OYZ  C   GLUAB  21     -38.098-164.151-162.911  1.00 56.45           C
ATOM  A1OZ0  O   GLUAB  21     -37.095-163.439-163.057  1.00 48.83           O
ATOM  A1OZ1  CB  GLUAB  21     -38.090-166.019-161.270  1.00 35.24           C
ATOM  A1OZ2  CG  GLUAB  21     -38.502-167.478-160.944  1.00 54.14           C
ATOM  A1OZ3  CD  GLUAB  21     -37.345-168.335-160.677  1.00 80.61           C
ATOM  A1OZ4  OE1 GLUAB  21     -36.278-168.016-161.240  1.00 92.19           O
ATOM  A1OZ5  OE2 GLUAB  21     -37.542-169.372-160.000  1.00 74.27           O
ATOM  A1OZC  N   GLYAB  22     -39.335-163.681-162.903  1.00 58.92           N
ATOM  A1OZD  CA  GLYAB  22     -39.594-162.287-163.140  1.00 48.41           C
ATOM  A1OZE  C   GLYAB  22     -39.333-161.886-164.568  1.00 61.25           C
ATOM  A1OZF  O   GLYAB  22     -38.868-160.772-164.828  1.00 52.03           O
ATOM  A1OZJ  N   SERAB  23     -39.573-162.785-165.504  1.00 81.16           N
ATOM  A1OZK  CA  SERAB  23     -39.313-162.516-166.905  1.00 71.05           C
ATOM  A1OZL  C   SERAB  23     -40.624-162.236-167.624  1.00 79.99           C
ATOM  A1OZM  O   SERAB  23     -41.718-162.618-167.190  1.00 75.06           O
ATOM  A1OZN  CB  SERAB  23     -38.567-163.682-167.576  1.00 68.18           C
ATOM  A1OZO  OG  SERAB  23     -37.509-164.144-166.746  1.00 77.85           O
ATOM  A1OZU  N   THRAB  24     -40.481-161.521-168.708  1.00 89.58           N
ATOM  A1OZV  CA  THRAB  24     -41.586-161.179-169.573  1.00 91.58           C
ATOM  A1OZW  C   THRAB  24     -42.024-162.351-170.470  1.00 98.61           C
ATOM  A1OZX  O   THRAB  24     -42.893-162.146-171.321  1.00100.89           O
ATOM  A1OZY  CB  THRAB  24     -41.175-159.953-170.402  1.00 91.75           C
ATOM  A1OZZ  OG1 THRAB  24     -42.261-159.554-171.238  1.00103.20           O
ATOM  A1P00  CG2 THRAB  24     -39.878-160.193-171.260  1.00 77.23           C
ATOM  A1P08  N   ILEAB  25     -41.467-163.557-170.293  1.00 96.60           N
ATOM  A1P09  CA  ILEAB  25     -41.834-164.747-171.060  1.00 81.06           C
ATOM  A1P0A  C   ILEAB  25     -42.395-165.793-170.099  1.00 77.96           C
ATOM  A1P0B  O   ILEAB  25     -41.993-165.860-168.929  1.00 71.48           O
ATOM  A1P0C  CB  ILEAB  25     -40.641-165.307-171.895  1.00 82.14           C
ATOM  A1P0D  CG1 ILEAB  25     -40.019-164.206-172.805  1.00 85.32           C
ATOM  A1P0E  CG2 ILEAB  25     -41.071-166.508-172.747  1.00 75.75           C
ATOM  A1P0F  CD1 ILEAB  25     -38.682-164.562-173.521  1.00 64.61           C
TER
ATOM  A1ZO1  N   LEUFB  20    -117.661-248.006-174.897  1.00 48.46           N
ATOM  A1ZO2  CA  LEUFB  20    -116.344-247.393-174.779  1.00 31.10           C
ATOM  A1ZO3  C   LEUFB  20    -115.899-246.648-176.028  1.00 42.66           C
ATOM  A1ZO4  O   LEUFB  20    -114.898-245.929-175.971  1.00 40.76           O
ATOM  A1ZO5  CB  LEUFB  20    -115.276-248.441-174.508  1.00 31.04           C
ATOM  A1ZO6  CG  LEUFB  20    -115.569-249.438-173.423  1.00 43.61           C
ATOM  A1ZO7  CD1 LEUFB  20    -114.283-250.230-173.234  1.00 48.24           C
ATOM  A1ZO8  CD2 LEUFB  20    -116.061-248.825-172.115  1.00 38.75           C
ATOM  A1ZOK  N   THRFB  21    -116.547-246.866-177.163  1.00 55.85           N
ATOM  A1ZOL  CA  THRFB  21    -116.140-246.297-178.441  1.00 43.47           C
ATOM  A1ZOM  C   THRFB  21    -117.270-245.418-178.984  1.00 51.67           C
ATOM  A1ZON  O   THRFB  21    -118.375-245.353-178.434  1.00 50.72           O
ATOM  A1ZOO  CB  THRFB  21    -115.774-247.401-179.425  1.00 33.43           C
ATOM  A1ZOP  OG1 THRFB  21    -116.920-248.240-179.614  1.00 50.93           O
ATOM  A1ZOQ  CG2 THRFB  21    -114.581-248.247-178.851  1.00 26.17           C
ATOM  A1ZOY  N   ARGFB  22    -116.980-244.715-180.072  1.00 49.13           N
ATOM  A1ZOZ  CA  ARGFB  22    -117.818-243.600-180.478  1.00 35.40           C
ATOM  A1ZP0  C   ARGFB  22    -117.553-243.328-181.949  1.00 44.20           C
ATOM  A1ZP1  O   ARGFB  22    -116.395-243.146-182.343  1.00 41.59           O
ATOM  A1ZP2  CB  ARGFB  22    -117.441-242.430-179.578  1.00 49.27           C
ATOM  A1ZP3  CG  ARGFB  22    -118.165-241.112-179.677  1.00 59.06           C
ATOM  A1ZP4  CD  ARGFB  22    -117.649-240.290-178.480  1.00 81.73           C
ATOM  A1ZP5  NE  ARGFB  22    -118.235-238.966-178.278  1.00 93.03           N
ATOM  A1ZP6  CZ  ARGFB  22    -118.016-238.230-177.181  1.00 94.32           C
ATOM  A1ZP7  NH1 ARGFB  22    -117.217-238.667-176.191  1.00 76.59           N
ATOM  A1ZP8  NH2 ARGFB  22    -118.595-237.036-177.072  1.00 78.59           N
ATOM  A1ZPM  N   ALAFB  23    -118.604-243.290-182.762  1.00 39.74           N
ATOM  A1ZPN  CA  ALAFB  23    -118.401-243.076-184.189  1.00 40.06           C
ATOM  A1ZPO  C   ALAFB  23    -118.078-241.613-184.495  1.00 42.35           C
ATOM  A1ZPP  O   ALAFB  23    -118.587-240.687-183.861  1.00 44.55           O
ATOM  A1ZPQ  CB  ALAFB  23    -119.627-243.535-184.974  1.00 39.71           C
ATOM  A1ZPW  N   LEUFB  24    -117.174-241.404-185.439  1.00 44.81           N
ATOM  A1ZPX  CA  LEUFB  24    -116.769-240.066-185.835  1.00 35.11           C
ATOM  A1ZPY  C   LEUFB  24    -116.826-239.893-187.335  1.00 40.09           C
ATOM  A1ZPZ  O   LEUFB  24    -116.660-240.858-188.085  1.00 51.80           O
ATOM  A1ZQ0  CB  LEUFB  24    -115.345-239.756-185.414  1.00 42.69           C
ATOM  A1ZQ1  CG  LEUFB  24    -115.154-239.849-183.909  1.00 55.89           C
ATOM  A1ZQ2  CD1 LEUFB  24    -113.661-239.714-183.568  1.00 33.93           C
ATOM  A1ZQ3  CD2 LEUFB  24    -116.049-238.781-183.220  1.00 35.88           C
ATOM  A1ZQF  N   PROFB  25    -117.004-238.671-187.801  1.00 40.39           N
ATOM  A1ZQG  CA  PROFB  25    -117.007-238.457-189.243  1.00 49.95           C
ATOM  A1ZQH  C   PROFB  25    -115.622-238.648-189.848  1.00 35.64           C
ATOM  A1ZQI  O   PROFB  25    -114.604-238.253-189.273  1.00 21.90           O
ATOM  A1ZQJ  CB  PROFB  25    -117.529-237.028-189.377  1.00 36.12           C
ATOM  A1ZQK  CG  PROFB  25    -117.210-236.413-188.096  1.00 26.49           C
ATOM  A1ZQL  CD  PROFB  25    -117.424-237.464-187.089  1.00 32.95           C
TER
END
"""
  phil_str1 = """\
ncs_group {
  reference        = chain 1
  selection        = chain FB
  selection        = chain AB
  selection        = chain GA
  }
"""
  pdb_h = iotbx.pdb.input(
      source_info=None,lines=pdb_str.split('\n')).\
      construct_hierarchy()
  # show file:
  # pdb_h.write_pdb_file("user_ncs_phil_ex12.pdb")
  asc = pdb_h.atom_selection_cache()


  ncs_groups = get_ncs_groups(phil_str1, pdb_str)
  # ncs_groups._show()
  assert len(ncs_groups) == 2
  assert len(ncs_groups[0].copies) == 1
  assert len(ncs_groups[1].copies) == 1
  # from iotbx.pdb.atom_selection import selection_string_from_selection
  # print list(ncs_groups[0].master_iselection)
  # print list(asc.iselection("chain 1"))
  # print selection_string_from_selection(pdb_h, ncs_groups[0].master_iselection)
  # print selection_string_from_selection(pdb_h, ncs_groups[0].copies[0].iselection)
  # STOP()
  assert approx_equal(ncs_groups[1].master_iselection,
      asc.iselection("chain AB"))
  assert approx_equal(ncs_groups[1].copies[0].iselection,
      asc.iselection("chain GA"))
  assert approx_equal(ncs_groups[0].master_iselection,
      asc.iselection("chain 1"))
  assert approx_equal(ncs_groups[0].copies[0].iselection,
      asc.iselection("chain FB"))



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
  exercise_11()
  exercise_12()

  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
