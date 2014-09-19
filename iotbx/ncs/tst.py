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
      ncs_selection_params = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_selection_params = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_selection_params = ncs_params_str)
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
      ncs_selection_params = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_selection_params = ncs_params_str)
    check_result(ncs_inp,test_i)
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_selection_params = ncs_params_str)
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
      ncs_selection_params = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    # using combination of pdb file name and Phil parameter string
    ncs_inp = iotbx.ncs.input(file_name = pdb_file_name,
      ncs_selection_params = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    # using combination of pdb string and Phil parameter string
    ncs_inp = iotbx.ncs.input(pdb_string = pdb_str,
      ncs_selection_params = ncs_params_str)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_02()
