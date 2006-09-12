from cctbx.eltbx import chemical_elements
from scitbx import stl
import scitbx.stl.set

def exercise():
  lc = chemical_elements.proper_caps_list()
  assert len(lc) == 111
  assert lc[:3] == ["H", "He", "Li"]
  lu = chemical_elements.proper_upper_list()
  assert len(lu) == len(lc)
  assert lu[:3] == ["H", "HE", "LI"]
  for c,u in zip(lc,lu): assert c.upper() == u
  sc = chemical_elements.proper_caps_set()
  assert len(sc) == len(lc)
  assert list(sc) == list(stl.set.stl_string(lc))
  su = chemical_elements.proper_upper_set()
  assert len(su) == len(lc)
  assert list(su) == list(stl.set.stl_string(lu))
  su = chemical_elements.proper_and_isotopes_upper_set()
  assert len(su) == len(lc) + 2
  assert list(su) == list(stl.set.stl_string(lu+["D", "T"]))
  print "OK"

if (__name__ == "__main__"):
  exercise()
