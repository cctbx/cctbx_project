from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import iotbx.ncs
from libtbx.utils import null_out

# ---------------------------------------------------------------------------
# Make sure iotbx.ncs.input (NCS search procedure) preserves iseqs on the passed
# hierarchy.
# ---------------------------------------------------------------------------

pdb_str = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1       5.043   9.706  12.193  1.00 10.00           N
ATOM      2  CA  GLY A   1       5.000   9.301  10.742  1.00 10.00           C
ATOM      3  C   GLY A   1       6.037   8.234  10.510  1.00 10.00           C
ATOM      4  O   GLY A   1       6.529   7.615  11.472  1.00 10.00           O
ATOM      5  H1  GLY A   1       5.852  10.022  12.389  1.00 10.00           H
ATOM      6  H2  GLY A   1       4.434  10.337  12.343  1.00 10.00           H
ATOM      7  H3  GLY A   1       4.869   8.998  12.704  1.00 10.00           H
ATOM      8  HA2 GLY A   1       5.191  10.064  10.175  1.00 10.00           H
ATOM      9  HA3 GLY A   1       4.123   8.952  10.517  1.00 10.00           H
ATOM     10  N   ASN A   2       6.396   8.017   9.246  1.00 10.00           N
ATOM     11  CA  ASN A   2       7.530   7.132   8.922  1.00 10.00           C
ATOM     12  C   ASN A   2       8.811   7.631   9.518  1.00 10.00           C
ATOM     13  O   ASN A   2       9.074   8.836   9.517  1.00 10.00           O
ATOM     14  CB  ASN A   2       7.706   6.975   7.432  1.00 10.00           C
ATOM     15  CG  ASN A   2       6.468   6.436   6.783  1.00 10.00           C
ATOM     16  OD1 ASN A   2       6.027   5.321   7.107  1.00 10.00           O
ATOM     17  ND2 ASN A   2       5.848   7.249   5.922  1.00 10.00           N
ATOM     18  H   ASN A   2       6.007   8.363   8.561  1.00 10.00           H
ATOM     19  HA  ASN A   2       7.333   6.259   9.297  1.00 10.00           H
ATOM     20  HB2 ASN A   2       7.904   7.840   7.040  1.00 10.00           H
ATOM     21  HB3 ASN A   2       8.434   6.358   7.259  1.00 10.00           H
ATOM     22 HD21 ASN A   2       5.133   6.987   5.522  1.00 10.00           H
ATOM     23 HD22 ASN A   2       6.164   8.034   5.768  1.00 10.00           H
ATOM     24  N   TYR A   3       9.614   6.684   9.996  1.00 10.00           N
ATOM     25  CA  TYR A   3      10.859   6.998  10.680  1.00 10.00           C
ATOM     26  C   TYR A   3      12.097   6.426   9.986  1.00 10.00           C
ATOM     27  O   TYR A   3      12.180   5.213   9.739  1.00 10.00           O
ATOM     28  CB  TYR A   3      10.793   6.472  12.133  1.00 10.00           C
ATOM     29  CG  TYR A   3      12.046   6.833  12.952  1.00 10.00           C
ATOM     30  CD1 TYR A   3      12.176   8.072  13.565  1.00 10.00           C
ATOM     31  CD2 TYR A   3      13.056   5.900  13.140  1.00 10.00           C
ATOM     32  CE1 TYR A   3      13.297   8.409  14.301  1.00 10.00           C
ATOM     33  CE2 TYR A   3      14.180   6.229  13.874  1.00 10.00           C
ATOM     34  CZ  TYR A   3      14.297   7.480  14.457  1.00 10.00           C
ATOM     35  OH  TYR A   3      15.417   7.809  15.189  1.00 10.00           O
ATOM     36  H   TYR A   3       9.456   5.841   9.935  1.00 10.00           H
ATOM     37  HA  TYR A   3      10.965   7.962  10.677  1.00 10.00           H
ATOM     38  HB2 TYR A   3      10.022   6.861  12.575  1.00 10.00           H
ATOM     39  HB3 TYR A   3      10.713   5.505  12.116  1.00 10.00           H
ATOM     40  HD1 TYR A   3      11.489   8.693  13.478  1.00 10.00           H
ATOM     41  HD2 TYR A   3      12.976   5.050  12.771  1.00 10.00           H
ATOM     42  HE1 TYR A   3      13.373   9.253  14.685  1.00 10.00           H
ATOM     43  HE2 TYR A   3      14.863   5.606  13.977  1.00 10.00           H
ATOM     44  HH  TYR A   3      15.352   8.596  15.475  1.00 10.00           H
END
"""

def exercise():
  h = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n")).construct_hierarchy()
  h.atoms().reset_i_seq()
  n = h.atoms_size()
  before = list(h.atoms().extract_i_seq())
  assert len(set(before)) == n, (len(set(before)), n)  # sanity: unique to start

  iotbx.ncs.input(hierarchy=h, log=null_out())

  after = list(h.atoms().extract_i_seq())
  # iotbx.ncs.input must not touch the caller's hierarchy.
  assert len(set(after)) == n, \
    "iotbx.ncs.input corrupted the caller's atom i_seqs: %d unique of %d atoms" % (
      len(set(after)), n)
  assert after == before, "iotbx.ncs.input changed the caller's atom i_seqs"
  print("OK")

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK. Time: %8.3f" % (time.time() - t0))
