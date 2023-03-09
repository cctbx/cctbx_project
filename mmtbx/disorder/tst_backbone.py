
from __future__ import absolute_import, division, print_function
from mmtbx.disorder import backbone
import iotbx.pdb
from six.moves import cStringIO as StringIO

pdb_raw = """
CRYST1   21.937    6.000   23.477  90.00 107.08  90.00 P 1 21 1      2
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
ATOM     13  N  AASN A   3      -4.438   1.590   3.905  0.50 12.26           N
ATOM     14  CA AASN A   3      -3.193   1.904   4.589  0.50 11.74           C
ATOM     15  C  AASN A   3      -1.955   1.332   3.895  0.50 11.10           C
ATOM     16  O  AASN A   3      -1.872   0.119   3.648  0.50 10.42           O
ATOM     17  CB AASN A   3      -3.259   1.378   6.042  0.50 12.15           C
ATOM     18  CG AASN A   3      -2.006   1.739   6.861  0.50 12.82           C
ATOM     19  OD1AASN A   3      -1.702   2.925   7.072  0.50 15.05           O
ATOM     20  ND2AASN A   3      -1.271   0.715   7.306  0.50 13.48           N
ATOM     13  N  BASN A   3      -4.438   1.590   3.905  0.50 12.26           N
ATOM     14  CA BASN A   3      -3.193   1.904   4.589  0.50 11.74           C
ATOM     15  C  BASN A   3      -1.993   1.011   4.266  0.50 11.10           C
ATOM     16  O  BASN A   3      -2.000  -0.192   4.568  0.50 10.42           O
ATOM     17  CB BASN A   3      -3.259   1.378   6.042  0.50 12.15           C
ATOM     18  CG BASN A   3      -2.006   1.739   6.861  0.50 12.82           C
ATOM     19  OD1BASN A   3      -1.702   2.925   7.072  0.50 15.05           O
ATOM     20  ND2BASN A   3      -1.271   0.715   7.306  0.50 13.48           N
ATOM     21  N  AGLN A   4      -1.005   2.228   3.598  0.50 10.29           N
ATOM     22  CA AGLN A   4       0.384   1.888   3.199  0.50 10.53           C
ATOM     23  C  AGLN A   4       1.435   2.606   4.088  0.50 10.24           C
ATOM     24  O  AGLN A   4       1.547   3.843   4.115  0.50  8.86           O
ATOM     25  CB AGLN A   4       0.656   2.148   1.711  0.50  9.80           C
ATOM     26  CG AGLN A   4       1.944   1.458   1.213  0.50 10.25           C
ATOM     27  CD AGLN A   4       2.504   2.044  -0.089  0.50 12.43           C
ATOM     28  OE1AGLN A   4       2.744   3.268  -0.190  0.50 14.62           O
ATOM     29  NE2AGLN A   4       2.750   1.161  -1.091  0.50  9.05           N
ATOM     21  N  BGLN A   4      -0.972   1.629   3.658  0.50 10.29           N
ATOM     22  CA BGLN A   4       0.392   1.063   3.508  0.50 10.53           C
ATOM     23  C  BGLN A   4       1.482   2.027   4.050  0.50 10.24           C
ATOM     24  O  BGLN A   4       1.689   3.145   3.548  0.50  8.86           O
ATOM     25  CB BGLN A   4       0.705   0.638   2.067  0.50  9.80           C
ATOM     26  CG BGLN A   4       1.943  -0.280   1.974  0.50 10.25           C
ATOM     27  CD BGLN A   4       2.565  -0.350   0.574  0.50 12.43           C
ATOM     28  OE1BGLN A   4       2.900   0.692  -0.032  0.50 14.62           O
ATOM     29  NE2BGLN A   4       2.757  -1.593   0.061  0.50  9.05           N
ATOM     30  N  AGLN A   5       2.154   1.821   4.871  0.50 10.38           N
ATOM     31  CA AGLN A   5       3.270   2.361   5.640  0.50 11.39           C
ATOM     32  C  AGLN A   5       4.594   1.768   5.172  0.50 11.52           C
ATOM     33  O  AGLN A   5       4.768   0.546   5.054  0.50 12.05           O
ATOM     34  CB AGLN A   5       3.056   2.183   7.147  0.50 11.96           C
ATOM     35  CG AGLN A   5       1.829   2.950   7.647  0.50 10.81           C
ATOM     36  CD AGLN A   5       1.344   2.414   8.954  0.50 13.10           C
ATOM     37  OE1AGLN A   5       0.774   1.325   9.002  0.50 10.65           O
ATOM     38  NE2AGLN A   5       1.549   3.187  10.039  0.50 12.30           N
ATOM     30  N  BGLN A   5       2.127   1.614   5.127  0.50 10.38           N
ATOM     31  CA BGLN A   5       3.270   2.361   5.640  0.50 11.39           C
ATOM     32  C  BGLN A   5       4.594   1.768   5.172  0.50 11.52           C
ATOM     33  O  BGLN A   5       4.768   0.546   5.054  0.50 12.05           O
ATOM     34  CB BGLN A   5       3.056   2.183   7.147  0.50 11.96           C
ATOM     35  CG BGLN A   5       1.829   2.950   7.647  0.50 10.81           C
ATOM     36  CD BGLN A   5       1.344   2.414   8.954  0.50 13.10           C
ATOM     37  OE1BGLN A   5       0.774   1.325   9.002  0.50 10.65           O
ATOM     38  NE2BGLN A   5       1.549   3.187  10.039  0.50 12.30           N
ATOM     39  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
TER      60      TYR A   7
HETATM   61  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   62  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   63  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   64  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM   65  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   66  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   67  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
END
"""

def exercise():
  hierarchy = iotbx.pdb.input(source_info=None, lines=pdb_raw).construct_hierarchy()
  with open("tmp.pdb", "w") as f:
    f.write(hierarchy.as_pdb_string())
  backrubs = backbone.find_backrubs(
    pdb_hierarchy=hierarchy)
  assert (len(backrubs) == 1)
  out = StringIO()
  backrubs[0].show(out=out)
  assert (out.getvalue() == """backrub  A   4  (A,B): angle=-26.0\n"""), \
    out.getvalue()

if (__name__ == "__main__"):
  exercise()
  print("OK")
