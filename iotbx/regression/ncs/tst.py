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

if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_02()
  exercise_03()
