from __future__ import absolute_import, division, print_function


def exercise():
  from mmtbx.refinement.print_statistics import coordinate_shifts
  import iotbx.pdb
  hierarchy_start = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   GLY A   1      -9.009  14.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052  14.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015  13.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523  12.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656  12.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522  12.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241  12.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978  13.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346  11.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584  11.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025  10.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204  12.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438  11.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193  11.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955  11.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872  10.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259  11.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006  11.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702  12.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271  10.715   7.306  1.00 13.48           N
END""").construct_hierarchy()
  hierarchy_end = iotbx.pdb.input(source_info=None, lines="""\
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
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
TER
HETATM   60  O   HOH     8      -6.471   5.227   7.124  1.00 20.00           O
HETATM   61  O   HOH     9      10.431   1.858   3.216  1.00 20.00           O
HETATM   62  O   HOH    10     -11.286   1.756  -1.468  1.00 20.00           O
END""").construct_hierarchy()
  cs = coordinate_shifts(hierarchy_start, hierarchy_end)
  #hierarchy_shifted = cs.get_shifts()
  mmm = cs.min_max_mean()
  assert (mmm.min==10) and (mmm.max==10) and (mmm.mean==10)
  cs.save_pdb_file("shifted.pdb")
  with open("shifted.pdb") as f:
    lines = f.readlines()
  h3 = iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()
  a3 = h3.atoms()
  shifts = a3.extract_b()
  assert (shifts[-1] == -1.0)

if (__name__ == "__main__"):
  exercise()
  print("OK")
