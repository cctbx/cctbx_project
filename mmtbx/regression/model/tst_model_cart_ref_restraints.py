from __future__ import absolute_import, division, print_function
import mmtbx.model
import libtbx.load_env
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
import iotbx.pdb
from six.moves import range

# 1yjp
pdb_str = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
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
ATOM     21  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     22  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     23  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     24  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.944   1.458   1.213  1.00 10.25           C
ATOM     27  CD  GLN A   4       2.504   2.044  -0.089  1.00 12.43           C
ATOM     28  OE1 GLN A   4       2.744   3.268  -0.190  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.750   1.161  -1.091  1.00  9.05           N
ATOM     30  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     31  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     32  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     33  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     34  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     35  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     36  CD  GLN A   5       1.344   2.414   8.954  1.00 13.10           C
ATOM     37  OE1 GLN A   5       0.774   1.325   9.002  1.00 10.65           O
ATOM     38  NE2 GLN A   5       1.549   3.187  10.039  1.00 12.30           N
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
"""

pdb_str_h = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      0  H1  GLY A   1      -9.510   5.337   6.223  1.00 16.77           H   new
ATOM      0  H2  GLY A   1      -9.322   3.948   6.605  1.00 16.77           H   new
ATOM      0  H3  GLY A   1      -8.169   4.790   6.336  1.00 16.77           H   new
ATOM      0  HA2 GLY A   1      -9.934   3.875   4.420  1.00 16.57           H   new
ATOM      0  HA3 GLY A   1      -8.879   4.974   4.082  1.00 16.57           H   new
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM      0  H   ASN A   2      -8.047   3.272   2.473  1.00 15.02           H   new
ATOM      0  HA  ASN A   2      -6.733   1.174   3.217  1.00 14.10           H   new
ATOM      0  HB2 ASN A   2      -6.122   2.739   0.949  1.00 15.38           H   new
ATOM      0  HB3 ASN A   2      -5.601   1.285   1.164  1.00 15.38           H   new
ATOM      0 HD21 ASN A   2      -8.947   1.914  -0.529  1.00 11.72           H   new
ATOM      0 HD22 ASN A   2      -7.860   2.919  -0.363  1.00 11.72           H   new
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
ATOM      0  H   ASN A   3      -4.602   0.748   3.841  1.00 12.26           H   new
ATOM      0  HA  ASN A   3      -3.100   2.869   4.573  1.00 11.74           H   new
ATOM      0  HB2 ASN A   3      -4.044   1.744   6.479  1.00 12.15           H   new
ATOM      0  HB3 ASN A   3      -3.366   0.414   6.029  1.00 12.15           H   new
ATOM      0 HD21 ASN A   3      -0.555   0.864   7.759  1.00 13.48           H   new
ATOM      0 HD22 ASN A   3      -1.514  -0.093   7.140  1.00 13.48           H   new
ATOM     21  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     22  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     23  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     24  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.944   1.458   1.213  1.00 10.25           C
ATOM     27  CD  GLN A   4       2.504   2.044  -0.089  1.00 12.43           C
ATOM     28  OE1 GLN A   4       2.744   3.268  -0.190  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.750   1.161  -1.091  1.00  9.05           N
ATOM      0  H   GLN A   4      -1.149   3.076   3.622  1.00 10.29           H   new
ATOM      0  HA  GLN A   4       0.474   0.933   3.341  1.00 10.53           H   new
ATOM      0  HB2 GLN A   4      -0.098   1.834   1.188  1.00  9.80           H   new
ATOM      0  HB3 GLN A   4       0.726   3.104   1.561  1.00  9.80           H   new
ATOM      0  HG2 GLN A   4       2.622   1.523   1.904  1.00 10.25           H   new
ATOM      0  HG3 GLN A   4       1.763   0.514   1.080  1.00 10.25           H   new
ATOM      0 HE21 GLN A   4       2.571   0.327  -0.983  1.00  9.05           H   new
ATOM      0 HE22 GLN A   4       3.085   1.436  -1.834  1.00  9.05           H   new
ATOM     30  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     31  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     32  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     33  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     34  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     35  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     36  CD  GLN A   5       1.344   2.414   8.954  1.00 13.10           C
ATOM     37  OE1 GLN A   5       0.774   1.325   9.002  1.00 10.65           O
ATOM     38  NE2 GLN A   5       1.549   3.187  10.039  1.00 12.30           N
ATOM      0  H   GLN A   5       2.017   0.978   4.974  1.00 10.38           H   new
ATOM      0  HA  GLN A   5       3.309   3.316   5.477  1.00 11.39           H   new
ATOM      0  HB2 GLN A   5       2.951   1.240   7.349  1.00 11.96           H   new
ATOM      0  HB3 GLN A   5       3.844   2.489   7.623  1.00 11.96           H   new
ATOM      0  HG2 GLN A   5       2.050   3.889   7.743  1.00 10.81           H   new
ATOM      0  HG3 GLN A   5       1.119   2.892   6.989  1.00 10.81           H   new
ATOM      0 HE21 GLN A   5       1.953   3.942   9.959  1.00 12.30           H   new
ATOM      0 HE22 GLN A   5       1.276   2.925  10.811  1.00 12.30           H   new
ATOM     39  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM      0  H   ASN A   6       5.395   3.511   4.947  1.00 11.99           H   new
ATOM      0  HA  ASN A   6       6.892   1.355   4.162  1.00 12.30           H   new
ATOM      0  HB2 ASN A   6       7.129   3.972   3.143  1.00 12.13           H   new
ATOM      0  HB3 ASN A   6       7.914   2.730   2.620  1.00 12.13           H   new
ATOM      0 HD21 ASN A   6       4.545   3.635   1.128  1.00 10.07           H   new
ATOM      0 HD22 ASN A   6       5.342   4.524   2.018  1.00 10.07           H   new
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
ATOM      0  H   TYR A   7       8.106   0.981   6.066  1.00 14.70           H   new
ATOM      0  HA  TYR A   7       8.843   2.985   7.664  1.00 15.18           H   new
ATOM      0  HB2 TYR A   7       9.349   0.216   7.999  1.00 15.35           H   new
ATOM      0  HB3 TYR A   7       9.666   1.277   9.097  1.00 15.35           H   new
ATOM      0  HD1 TYR A   7       7.051  -0.513   7.619  1.00 15.68           H   new
ATOM      0  HD2 TYR A   7       7.783   2.394  10.280  1.00 14.80           H   new
ATOM      0  HE1 TYR A   7       4.901  -0.719   8.422  1.00 13.46           H   new
ATOM      0  HE2 TYR A   7       5.618   2.183  11.122  1.00 14.33           H   new
ATOM      0  HH  TYR A   7       3.629   1.138  10.911  1.00 14.39           H   new
"""

pdb_str_base = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
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
"""

pdb_str_water = """\
HETATM   61  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   62  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   63  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   64  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM   65  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
"""

def exercise_adopting_coord_restraints():
  inp_1 = iotbx.pdb.input(lines=pdb_str_h, source_info=None)
  h_model = mmtbx.model.manager(model_input = inp_1)
  h_model.process(make_restraints=True)

  inp_2 = iotbx.pdb.input(lines=pdb_str, source_info=None)
  model = mmtbx.model.manager(model_input = inp_2)
  model.process(make_restraints=True)

  # case 1: big model adopting small model.
  h_model.set_reference_coordinate_restraints(ref_model = model)
  # print h_model.restraints_as_geo()
  rcp = h_model.get_restraints_manager().geometry.reference_coordinate_proxies
  assert len(rcp) == 59
  # Note iseqs, hydrogens are skipped
  answer = [ [(0,) , (-9.009, 4.612, 6.102), 25.0],
             [(1,) , (-9.052, 4.207, 4.651), 25.0],
             [(2,) , (-8.015, 3.14,  4.419), 25.0],
             [(3,) , (-7.523, 2.521, 5.381), 25.0],
             [(9,) , (-7.656, 2.923, 3.155), 25.0],
             [(10,), (-6.522, 2.038, 2.831), 25.0]]
  # for p in rcp:
  #   print p.i_seqs, p.ref_sites, p.weight

  for i in range(6):
    assert rcp[i].i_seqs == answer[i][0], "%s %s" % (rcp[i].i_seqs, answer[i][0])
    assert approx_equal(rcp[i].ref_sites, answer[i][1])
    assert approx_equal(rcp[i].weight, answer[i][2]), "%s, %s" % (rcp[i].weight, answer[i][2])

  # case 2: small model adopting big model.
  model.set_reference_coordinate_restraints(ref_model=h_model)
  rcp = model.get_restraints_manager().geometry.reference_coordinate_proxies
  assert len(rcp) == 59
  # Note iseqs, hydrogens are skipped
  answer = [ [(0,), (-9.009, 4.612, 6.102), 25.0],
             [(1,), (-9.052, 4.207, 4.651), 25.0],
             [(2,), (-8.015, 3.14,  4.419), 25.0],
             [(3,), (-7.523, 2.521, 5.381), 25.0],
             [(4,), (-7.656, 2.923, 3.155), 25.0],
             [(5,),  (-6.522, 2.038, 2.831), 25.0]]
  for i in range(6):
    assert rcp[i].i_seqs == answer[i][0]
    assert approx_equal(rcp[i].ref_sites, answer[i][1])
    assert approx_equal(rcp[i].weight, answer[i][2]), "%s, %s" % (rcp[i].weight, answer[i][2])

def exercise_adopting_coord_restraints_water():
  inp_1 = iotbx.pdb.input(lines=pdb_str_base, source_info=None)
  nowat_model = mmtbx.model.manager(model_input = inp_1)
  nowat_model.process(make_restraints=True)
  # Changing coordinates a little to make sure right ones are used for
  # reference restraints
  for a in nowat_model.get_hierarchy().atoms():
    a.xyz = (a.xyz[0]+1, a.xyz[1], a.xyz[2])
  nowat_model.set_sites_cart_from_hierarchy()

  inp_2 = iotbx.pdb.input(lines=pdb_str_base+pdb_str_water, source_info=None)
  wat_model = mmtbx.model.manager(model_input = inp_2)
  wat_model.process(make_restraints=True)
  nowat_model.set_reference_coordinate_restraints(
      ref_model = wat_model)
  rcp = nowat_model.get_restraints_manager().geometry.reference_coordinate_proxies
  assert len(rcp) == 20, len(rcp)
  answer = [[(0,) , (-9.009, 4.612, 6.102), 25.0],
            [(1,) , (-9.052, 4.207, 4.651), 25.0],
            [(2,) , (-8.015, 3.14,  4.419), 25.0],
            [(3,) , (-7.523, 2.521, 5.381), 25.0],
            [(4,) , (-7.656, 2.923, 3.155), 25.0],
            [(5,) , (-6.522, 2.038, 2.831), 25.0],
            [(6,) , (-5.241, 2.537, 3.427), 25.0],
            [(7,) , (-4.978, 3.742, 3.426), 25.0],
            [(8,) , (-6.346, 1.881, 1.341), 25.0],
            [(9,) , (-7.584, 1.342, 0.692), 25.0],
            [(10,), (-8.025, 0.227, 1.016), 25.0],
            [(11,), (-8.204, 2.155, -0.169), 25.0],
            [(12,), (-4.438, 1.59,  3.905), 25.0],
            [(13,), (-3.193, 1.904, 4.589), 25.0],
            [(14,), (-1.955, 1.332, 3.895), 25.0],
            [(15,), (-1.872, 0.119, 3.648), 25.0],
            [(16,), (-3.259, 1.378, 6.042), 25.0],
            [(17,), (-2.006, 1.739, 6.861), 25.0],
            [(18,), (-1.702, 2.925, 7.072), 25.0],
            [(19,), (-1.271, 0.715, 7.306), 25.0]]
  for i in range(20):
    assert rcp[i].i_seqs == answer[i][0]
    assert approx_equal(rcp[i].ref_sites, answer[i][1])
    assert approx_equal(rcp[i].weight, answer[i][2]), "%s, %s" % (rcp[i].weight, answer[i][2])

  # print nowat_model.model_as_pdb()
  wat_model.set_reference_coordinate_restraints(ref_model=nowat_model)
  rcp = wat_model.get_restraints_manager().geometry.reference_coordinate_proxies
  assert len(rcp) == 20, len(rcp)
  # Note different x coord here
  answer = [[(0,) , (-8.009, 4.612, 6.102), 25.0],
            [(1,) , (-8.052, 4.207, 4.651), 25.0],
            [(2,) , (-7.015, 3.14,  4.419), 25.0],
            [(3,) , (-6.523, 2.521, 5.381), 25.0],
            [(4,) , (-6.656, 2.923, 3.155), 25.0],
            [(5,) , (-5.522, 2.038, 2.831), 25.0],
            [(6,) , (-4.241, 2.537, 3.427), 25.0],
            [(7,) , (-3.978, 3.742, 3.426), 25.0],
            [(8,) , (-5.346, 1.881, 1.341), 25.0],
            [(9,) , (-6.584, 1.342, 0.692), 25.0],
            [(10,), (-7.025, 0.227, 1.016), 25.0],
            [(11,), (-7.204, 2.155, -0.169), 25.0],
            [(12,), (-3.438, 1.59,  3.905), 25.0],
            [(13,), (-2.193, 1.904, 4.589), 25.0],
            [(14,), (-0.955, 1.332, 3.895), 25.0],
            [(15,), (-0.872, 0.119, 3.648), 25.0],
            [(16,), (-2.259, 1.378, 6.042), 25.0],
            [(17,), (-1.006, 1.739, 6.861), 25.0],
            [(18,), (-0.702, 2.925, 7.072), 25.0],
            [(19,), (-0.271, 0.715, 7.306), 25.0]]
  for i in range(20):
    assert rcp[i].i_seqs == answer[i][0]
    assert approx_equal(rcp[i].ref_sites, answer[i][1])
    assert approx_equal(rcp[i].weight, answer[i][2]), "%s, %s" % (rcp[i].weight, answer[i][2])
  # for p in rcp:
  #   print p.i_seqs, p.ref_sites, p.weight

def run():
  if (not libtbx.env.has_module("reduce")):
    print("Reduce not installed.")
    return
  exercise_adopting_coord_restraints()
  exercise_adopting_coord_restraints_water()

  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
