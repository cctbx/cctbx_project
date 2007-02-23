from iotbx.pdb.atom_name_interpretation import interpreters
from cStringIO import StringIO
from libtbx.test_utils import show_diff
import sys

def exercise():
  atom_names = ["CA", "C", "C", "X", "HB1", "HB2", "HB3"]
  matched = interpreters["LEU"].match_atom_names(atom_names)
  s = StringIO()
  matched.show_problems(out=s, prefix=">")
  assert not show_diff(s.getvalue(), """\
>unexpected atom names: "X"
>multiple matches: expected pattern=C  names="C", "C"
>mutually exclusive: 1hB 3hB
""")
  print "OK"

if (__name__ == "__main__"):
  exercise()
