"""mmtbx.model.manager.add_terminal_op3: OP3 on the phosphate of a 5' end
that lacks it (P, OP1, OP2 and O5' present, P not bonded to another
residue); nothing on a P bonded to the previous O3', nor where OP3 exists.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from scitbx import matrix
from mmtbx.regression.tst_reduce2_capping import pdb_dna3

def get_model(text):
  return mmtbx.model.manager(model_input=iotbx.pdb.input(
    source_info=None, lines=text.splitlines()))

def atoms_of(model, resseq):
  return dict((a.name.strip(), a) for a in model.get_hierarchy().atoms()
              if a.parent().parent().resseq.strip() == str(resseq))

def test_bonded_p_untouched():
  # DC 1 is a 5'-OH end (no P), DG 2 and DC 3 are bonded to the previous O3'
  m = get_model(pdb_dna3)
  assert m.add_terminal_op3() == 0
  assert "OP3" not in atoms_of(m, 2) and "OP3" not in atoms_of(m, 3)

def test_5p_phosphate():
  # DC 1 removed: DG 2 is a 5' end with P, OP1, OP2 and no OP3
  text = "\n".join(l for l in pdb_dna3.splitlines()
    if not (l.startswith("ATOM") and l[22:26].strip() == "1")) + "\n"
  m = get_model(text)
  n_atoms = m.get_number_of_atoms()
  assert m.add_terminal_op3() == 1
  assert m.get_number_of_atoms() == n_atoms + 1
  a = atoms_of(m, 2)
  assert "OP3" in a and "OP3" not in atoms_of(m, 3), sorted(a)
  op3, p = a["OP3"], a["P"]
  assert abs(op3.distance(p) - 1.48) < 1e-3, op3.distance(p)
  # tetrahedral: ~109 degrees to the three other O
  for name in ("O5'", "OP1", "OP2"):
    u = matrix.col(a[name].xyz) - matrix.col(p.xyz)
    v = matrix.col(op3.xyz) - matrix.col(p.xyz)
    assert 95 < u.angle(v, deg=True) < 125, (name, u.angle(v, deg=True))
  assert op3.occ == a["OP2"].occ and op3.b == a["OP2"].b
  assert op3.element.strip() == "O"
  # a second call adds nothing
  assert m.add_terminal_op3() == 0

def test_missing_atoms_skipped():
  text = "\n".join(l for l in pdb_dna3.splitlines()
    if not (l.startswith("ATOM") and l[22:26].strip() == "1")
    and not (l.startswith("ATOM") and l[22:26].strip() == "2"
             and l[12:16].strip() == "OP1")) + "\n"
  m = get_model(text)
  assert m.add_terminal_op3() == 0

if __name__ == "__main__":
  test_bonded_p_untouched()
  test_5p_phosphate()
  test_missing_atoms_skipped()
  print("OK")
