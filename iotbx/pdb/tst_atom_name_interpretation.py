from iotbx.pdb.atom_name_interpretation import interpreters
from cStringIO import StringIO
from libtbx.test_utils import show_diff

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

  ala_names = ["N", "CA", "C", "O", "CB"]
  ala_missing_names = ["CA", "CB", "O"]
  ala_with_h = ["N","CA","C","O","CB","HA","HB1","HB2","HB3","H"]
  ala_matched = interpreters["ALA"].match_atom_names(ala_names)
  ala_missing = interpreters["ALA"].match_atom_names(ala_missing_names)
  ala_h = interpreters["ALA"].match_atom_names(ala_with_h)
  assert ala_matched.missing_atom_names(ignore_hydrogen=True) == set(())
  assert ala_missing.missing_atom_names(ignore_hydrogen=True) == set(("C","N"))
  assert ala_h.missing_atom_names(ignore_hydrogen=False) == set(())


  print "OK"

if (__name__ == "__main__"):
  exercise()
