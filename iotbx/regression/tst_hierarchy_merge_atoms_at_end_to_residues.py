from __future__ import absolute_import, division, print_function
import iotbx.pdb

def exercise_1():
  pdb_str_1 = """\
ATOM   1492  N   TYR 1   7     -13.645   6.683   6.147  1.00 14.70           N
ATOM   1493  CA  TYR 1   7     -12.778   7.010   7.299  1.00 15.18           C
ATOM   1494  C   TYR 1   7     -11.334   7.197   6.885  1.00 15.91           C
ATOM   1495  O   TYR 1   7     -10.896   6.677   5.855  1.00 15.76           O
ATOM   1496  CB  TYR 1   7     -12.876   5.931   8.369  1.00 15.35           C
ATOM   1497  CG  TYR 1   7     -14.272   5.795   8.902  1.00 14.45           C
ATOM   1498  CD1 TYR 1   7     -14.727   6.622   9.920  1.00 14.80           C
ATOM   1499  CD2 TYR 1   7     -15.166   4.887   8.327  1.00 15.68           C
ATOM   1500  CE1 TYR 1   7     -16.033   6.515  10.416  1.00 14.33           C
ATOM   1501  CE2 TYR 1   7     -16.457   4.772   8.796  1.00 13.46           C
ATOM   1502  CZ  TYR 1   7     -16.890   5.595   9.831  1.00 15.09           C
ATOM   1503  OH  TYR 1   7     -18.171   5.455  10.291  1.00 14.39           O
ATOM   1504  OXT TYR 1   7     -10.579   7.865   7.612  1.00 17.49           O
ATOM   1505  H   TYR 1   7     -13.863   5.691   6.055  0.00 14.70           H
ATOM   1506  HA  TYR 1   7     -13.122   7.951   7.728  0.00 15.18           H
ATOM   1507  HB2 TYR 1   7     -12.217   6.188   9.199  0.00 15.35           H
ATOM   1508  HB3 TYR 1   7     -12.579   4.974   7.941  0.00 15.35           H
ATOM   1509  HD1 TYR 1   7     -14.061   7.362  10.338  0.00 14.80           H
ATOM   1510  HD2 TYR 1   7     -14.841   4.267   7.504  0.00 15.68           H
ATOM   1511  HE1 TYR 1   7     -16.365   7.135  11.236  0.00 14.33           H
ATOM   1512  HE2 TYR 1   7     -17.129   4.047   8.362  0.00 13.46           H
ATOM   1513  HH  TYR 1   7     -18.330   6.073  11.035  0.00 14.39           H
TER
ATOM   1560  N   GLN 3   4     -20.932   4.661  -3.598  1.00 10.29           N
ATOM   1561  CA  GLN 3   4     -22.321   4.321  -3.199  1.00 10.53           C
ATOM   1562  C   GLN 3   4     -23.372   5.039  -4.088  1.00 10.24           C
ATOM   1563  O   GLN 3   4     -23.484   6.276  -4.115  1.00  8.86           O
ATOM   1564  CB  GLN 3   4     -22.593   4.581  -1.711  1.00  9.80           C
ATOM   1565  CG  GLN 3   4     -23.881   3.891  -1.213  1.00 10.25           C
ATOM   1566  CD  GLN 3   4     -24.441   4.477   0.089  1.00 12.43           C
ATOM   1567  OE1 GLN 3   4     -24.681   5.701   0.190  1.00 14.62           O
ATOM   1568  NE2 GLN 3   4     -24.687   3.594   1.091  1.00  9.05           N
ATOM   1569  H   GLN 3   4     -20.775   5.669  -3.622  0.00 10.29           H
ATOM   1570  HA  GLN 3   4     -22.442   3.248  -3.346  0.00 10.53           H
ATOM   1571  HB2 GLN 3   4     -22.701   5.654  -1.552  0.00  9.80           H
ATOM   1572  HB3 GLN 3   4     -21.757   4.199  -1.125  0.00  9.80           H
ATOM   1573  HG2 GLN 3   4     -24.650   3.992  -1.979  0.00 10.25           H
ATOM   1574  HG3 GLN 3   4     -23.667   2.837  -1.036  0.00 10.25           H
ATOM   1575 HE21 GLN 3   4     -25.060   3.924   1.981  0.00  9.05           H
ATOM   1576 HE22 GLN 3   4     -24.499   2.601   0.953  0.00  9.05           H
TER
ATOM         H2  GLN 1   4     -23.584   6.695   2.944  1.00 10.29           H
TER
ATOM         HC  GLN 3   4     -23.988   4.503  -4.665  1.00 10.24           H
TER
"""
  h = iotbx.pdb.input(lines=pdb_str_1, source_info=None).construct_hierarchy()
  atoms_in_chains = [(c.id, c.atoms_size()) for c in h.only_model().chains()]
  print (atoms_in_chains)
  assert atoms_in_chains == [('1', 22), ('3', 17), ('1', 1), ('3', 1)], atoms_in_chains

  h.merge_atoms_at_end_to_residues()
  atoms_in_chains2 = [(c.id, c.atoms_size()) for c in h.only_model().chains()]
  # Note that the second chain '3' is removed and there's no chain with zero atoms.
  print (atoms_in_chains2)
  assert atoms_in_chains2 == [('1', 22), ('3', 18), ('1', 1)], atoms_in_chains

if (__name__ == "__main__"):
  exercise_1()
  print("OK")
