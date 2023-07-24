from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb

pdb_str = """
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
"""

def get_h_and_sel():
  inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
  h = inp.construct_hierarchy()
  asc=h.atom_selection_cache()
  sel = asc.selection(string = "all")
  return h, sel


# In these tests we will use id_str() to show if the information is available in the hierarchy
# These tests show how .select() can corrupt resulting hierarchy if not done carefully.

def test1():
  """Normal operations  """
  h, sel = get_h_and_sel()
  assert h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % h.atoms()[0].id_str()
  sel_h = h.select(sel)
  assert sel_h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % sel_h.atoms()[0].id_str()

def test2():
  """ Still works"""
  h, sel = get_h_and_sel()
  assert h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % h.atoms()[0].id_str()
  copy_h = h.deep_copy()
  sel_h = copy_h.select(sel)
  assert sel_h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % sel_h.atoms()[0].id_str()

def test3():
  """ Does not work"""
  h, sel = get_h_and_sel()
  assert h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % h.atoms()[0].id_str()
  sel_h = h.deep_copy().select(sel)
  print("There is no info:", sel_h.atoms()[0].id_str())
  print("And no parent:", sel_h.atoms()[0].parent())
  # assert sel_h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % sel_h.atoms()[0].id_str()

def test4():
  """ What happened in test3(): """
  h, sel = get_h_and_sel()
  assert h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % h.atoms()[0].id_str()
  copy_h = h.deep_copy()
  sel_h = copy_h.select(sel)
  # Everything is good for now
  assert sel_h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % sel_h.atoms()[0].id_str()
  del copy_h
  # but when the copy_h is deleted (or garbage collector removes it) the info is gone:
  # That's exactly what happens when one does sel_h = h.deep_copy().select(), the result of
  # h.deep_copy() is dropped immediately.
  print("Info is gone:", sel_h.atoms()[0].id_str())
  # assert sel_h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % sel_h.atoms()[0].id_str()

def test5():
  """ That's how test3 should be written, note copy_atoms=True """
  h, sel = get_h_and_sel()
  assert h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % h.atoms()[0].id_str()
  sel_h = h.deep_copy().select(sel, copy_atoms=True)
  print("Info is still here:", sel_h.atoms()[0].id_str())
  assert sel_h.atoms()[0].id_str() == 'pdb=" N   GLY A   1 "', "'%s'" % sel_h.atoms()[0].id_str()

if (__name__ == "__main__"):
  t0 = time.time()
  test1()
  test2()
  test3()
  test4()
  test5()
  print("OK. Time: %8.3f"%(time.time()-t0))
