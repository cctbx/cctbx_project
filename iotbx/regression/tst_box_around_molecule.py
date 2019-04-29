from __future__ import division, print_function
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file

pdb_str = """\
ATOM      1  N   SER 1  71      10.595  29.158 114.534  1.00300.00           N
ATOM      2  CA  SER 1  71      10.156  30.582 114.477  1.00300.00           C
ATOM      3  C   SER 1  71      11.003  31.474 115.392  1.00300.00           C
ATOM      4  O   SER 1  71      10.892  31.400 116.622  1.00300.00           O
"""

def exercise(prefix="tst_pdb_box_around_mol"):
  with open("%s.pdb" % prefix, 'w') as f:
    f.write(pdb_str)
  assert not easy_run.call("iotbx.pdb.box_around_molecule %s.pdb > %s_boxed.pdb" % (
      prefix, prefix))
  assert_lines_in_file(file_name="%s_boxed.pdb" % prefix,
    lines="""\
CRYST1   10.847   12.316   12.145  90.00  90.00  90.00 P 1
SCALE1      0.092191  0.000000  0.000000        0.00000
SCALE2      0.000000  0.081195  0.000000        0.00000
SCALE3      0.000000  0.000000  0.082338        0.00000
ATOM      1  N   SER 1  71       5.439   5.000   5.057  1.00300.00           N
ATOM      2  CA  SER 1  71       5.000   6.424   5.000  1.00300.00           C
ATOM      3  C   SER 1  71       5.847   7.316   5.915  1.00300.00           C
ATOM      4  O   SER 1  71       5.736   7.242   7.145  1.00300.00           O
TER""")
  print("OK")

if __name__=="__main__":
  exercise()
