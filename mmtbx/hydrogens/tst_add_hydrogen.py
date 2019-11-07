from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
import mmtbx.hydrogens
from libtbx.utils import null_out
#from libtbx.test_utils import approx_equal


def run():
  incomplete_peptide_unit()
  CD1_TYR_missing()
  no_H_on_disulfides()
  add_H_to_lone_CYS()
  peptide_unit_N_missing()
  no_H_on_coordinating_CYS()
  sym_related_disulfide()


def get_model(pdb_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp,
                              log         = null_out())
  model_with_h = mmtbx.hydrogens.add(model = model)
  return model_with_h


def compare_models(pdb_str,
                   contains     = None,
                   not_contains = None,
                   number_h     = None):
  #
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  xrs = model_initial.get_xray_structure()
  hd_sel_initial = model_initial.get_hd_selection()
  number_h_expected = hd_sel_initial.count(True)

  model_without_h = model_initial.select(~hd_sel_initial)
  hd_sel_without_h = model_without_h.get_hd_selection()
  assert (hd_sel_without_h is not None)
  assert (hd_sel_without_h.count(True) == 0)

  model_h_added = mmtbx.hydrogens.add(model = model_without_h)
  hd_sel_h_added = model_h_added.get_hd_selection()
  number_h_added = hd_sel_h_added.count(True)
  ph_initial = model_initial.get_hierarchy()
  ph_h_added = model_h_added.get_hierarchy()
  assert ph_initial.is_similar_hierarchy(other=ph_h_added)

  if number_h:
    assert(number_h == number_h_added)

  if not_contains:
    h_atoms_added = model_h_added.get_hierarchy().select(hd_sel_h_added).atoms()
    h_names_added = list(h_atoms_added.extract_name())
    assert (not_contains not in h_names_added)

  if contains:
    h_atoms_added = model_h_added.get_hierarchy().select(hd_sel_h_added).atoms()
    h_names_added = list(h_atoms_added.extract_name())
    assert (contains in h_names_added)


def incomplete_peptide_unit():
  """Incomplete peptide unit (N-terminal) --> H should not be placed."""
  compare_models(pdb_str = pdb_str_1)


def CD1_TYR_missing():
  """Atom CD1 is missing --> HE1 should not be placed

     H, HD1, HE1 cannot be parameterized and are not added.
  """
  compare_models(pdb_str = pdb_str_2,
                 number_h = 6,
                 not_contains = ' HE1')


def no_H_on_disulfides():
  """Don't put H on disulfides."""
  model_h_added = compare_models(pdb_str = pdb_str_3,
                                 number_h = 6,
                                 not_contains = ' HG ')


def add_H_to_lone_CYS():
  """Make sure to add HG on lone Cys."""
  compare_models(pdb_str = pdb_str_4,
                 number_h = 4,
                 contains = ' HG ')


def peptide_unit_N_missing():
  """N atom is missing, make sure peptide H is not placed"""
  compare_models(pdb_str = pdb_str_5,
                 number_h = 8,
                 not_contains = ' H  ')


def no_H_on_coordinating_CYS():
  """Don't put H on Cys coordinating a metal"""
  compare_models(pdb_str = pdb_str_6,
                 number_h = 3,
                 not_contains = ' HG ')


def sym_related_disulfide():
  """Don't put H on Cys forming a disulfide bond with symmetry related molecule"""
  compare_models(pdb_str = pdb_str_7,
                 number_h = 3,
                 not_contains = ' HG ')




pdb_str_1 = """
CRYST1   14.074   14.385   15.410  90.00  90.00  90.00 P 1
SCALE1      0.071053  0.000000  0.000000        0.00000
SCALE2      0.000000  0.069517  0.000000        0.00000
SCALE3      0.000000  0.000000  0.064893        0.00000
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  HA1 GLY A   1      -8.861   4.970   4.084  1.00 16.57           H
ATOM      6  HA2 GLY A   1      -9.929   3.858   4.426  1.00 16.57           H
ATOM      8  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      9  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM     10  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM     11  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM     12  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     13  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     14  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     15  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     16 HD22 ASN A   2      -7.888   2.940  -0.323  1.00 11.72           H
ATOM     17 HD21 ASN A   2      -8.919   1.893  -0.569  1.00 11.72           H
ATOM     18  H   ASN A   2      -8.045   3.269   2.470  1.00 15.02           H
ATOM     19  HA  ASN A   2      -6.719   1.165   3.206  1.00 14.10           H
ATOM     20  HB1 ASN A   2      -6.148   2.746   0.949  1.00 15.38           H
ATOM     21  HB2 ASN A   2      -5.618   1.264   1.168  1.00 15.38           H
TER
"""

pdb_str_2 = """
CRYST1   17.955   13.272   13.095  90.00  90.00  90.00 P 1
SCALE1      0.055695  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075347  0.000000        0.00000
SCALE3      0.000000  0.000000  0.076365        0.00000
ATOM      1  N   TYR A 139      10.241   7.920   5.000  1.00 10.00           N
ATOM      2  CA  TYR A 139      10.853   7.555   6.271  1.00 10.00           C
ATOM      3  C   TYR A 139      12.362   7.771   6.227  1.00 10.00           C
ATOM      4  O   TYR A 139      12.955   8.272   7.181  1.00 10.00           O
ATOM      5  CB  TYR A 139      10.540   6.098   6.617  1.00 10.00           C
ATOM      6  CG  TYR A 139       9.063   5.805   6.749  1.00 10.00           C
ATOM      7  CD2 TYR A 139       8.414   5.943   7.969  1.00 10.00           C
ATOM      8  CE1 TYR A 139       6.966   5.122   5.770  1.00 10.00           C
ATOM      9  CE2 TYR A 139       7.064   5.676   8.095  1.00 10.00           C
ATOM     10  CZ  TYR A 139       6.345   5.266   6.993  1.00 10.00           C
ATOM     11  OH  TYR A 139       5.000   5.000   7.113  1.00 10.00           O
ATOM     12  HE2 TYR A 139       6.643   5.772   8.919  1.00 10.00           H
ATOM     13  HD2 TYR A 139       8.896   6.220   8.714  1.00 10.00           H
ATOM     17  HH  TYR A 139       4.565   5.718   7.120  1.00 10.00           H
ATOM     18  HA  TYR A 139      10.487   8.115   6.973  1.00 10.00           H
ATOM     19  HB1 TYR A 139      10.961   5.881   7.464  1.00 10.00           H
ATOM     20  HB2 TYR A 139      10.893   5.529   5.916  1.00 10.00           H
TER
"""

pdb_str_3 = """
CRYST1   14.197   15.507   16.075  90.00  90.00  90.00 P 1
SCALE1      0.070437  0.000000  0.000000        0.00000
SCALE2      0.000000  0.064487  0.000000        0.00000
SCALE3      0.000000  0.000000  0.062208        0.00000
ATOM      1  N   CYS A 128      14.558 -30.378 -19.478  1.00 11.73           N
ATOM      2  CA  CYS A 128      13.499 -29.511 -19.901  1.00 11.51           C
ATOM      3  C   CYS A 128      12.454 -30.316 -20.640  1.00 11.99           C
ATOM      4  O   CYS A 128      12.765 -31.329 -21.265  1.00 16.86           O
ATOM      5  CB  CYS A 128      14.059 -28.357 -20.691  1.00 12.94           C
ATOM      6  SG  CYS A 128      15.213 -27.350 -19.760  1.00 12.36           S
ATOM      8  HA  CYS A 128      13.044 -29.102 -19.148  1.00 11.51           H
ATOM      9  HB1 CYS A 128      14.526 -28.706 -21.466  1.00 12.94           H
ATOM     11  HB2 CYS A 128      13.327 -27.787 -20.974  1.00 12.94           H
ATOM     12  N   CYS A 232      13.105 -25.822 -15.371  1.00  4.75           N
ATOM     13  CA  CYS A 232      14.298 -26.646 -15.573  1.00  5.50           C
ATOM     14  C   CYS A 232      15.596 -25.838 -15.551  1.00  5.70           C
ATOM     15  O   CYS A 232      16.651 -26.346 -15.190  1.00  8.25           O
ATOM     16  CB  CYS A 232      14.178 -27.497 -16.824  1.00  6.69           C
ATOM     17  SG  CYS A 232      14.060 -26.477 -18.360  1.00  7.70           S
ATOM     19  HA  CYS A 232      14.363 -27.258 -14.823  1.00  5.50           H
ATOM     20  HB1 CYS A 232      14.961 -28.065 -16.898  1.00  6.69           H
ATOM     22  HB2 CYS A 232      13.378 -28.042 -16.761  1.00  6.69           H
TER
"""

pdb_str_4 = """
CRYST1   14.197   15.507   16.075  90.00  90.00  90.00 P 1
SCALE1      0.070437  0.000000  0.000000        0.00000
SCALE2      0.000000  0.064487  0.000000        0.00000
SCALE3      0.000000  0.000000  0.062208        0.00000
ATOM      1  N   CYS A 128      14.558 -30.378 -19.478  1.00 11.73           N
ATOM      2  CA  CYS A 128      13.499 -29.511 -19.901  1.00 11.51           C
ATOM      3  C   CYS A 128      12.454 -30.316 -20.640  1.00 11.99           C
ATOM      4  O   CYS A 128      12.765 -31.329 -21.265  1.00 16.86           O
ATOM      5  CB  CYS A 128      14.059 -28.357 -20.691  1.00 12.94           C
ATOM      6  SG  CYS A 128      15.213 -27.350 -19.760  1.00 12.36           S
ATOM      8  HA  CYS A 128      13.044 -29.102 -19.148  1.00 11.51           H
ATOM      9  HB1 CYS A 128      14.526 -28.706 -21.466  1.00 12.94           H
ATOM     10  HG  CYS A 128      15.585 -27.968 -18.801  1.00 12.36           H
ATOM     11  HB2 CYS A 128      13.327 -27.787 -20.974  1.00 12.94           H
TER
"""

pdb_str_5 = """
CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 1
SCALE1      0.066667  0.000000  0.000000        0.00000
SCALE2      0.000000  0.066667  0.000000        0.00000
SCALE3      0.000000  0.000000  0.066667        0.00000
ATOM      1  CA  PHE H   1       7.000   7.000   7.000  1.00 15.00           C
ATOM      2  C   PHE H   1       7.956   7.811   6.133  1.00 15.00           C
ATOM      3  O   PHE H   1       8.506   7.237   5.169  1.00 15.00           O
ATOM      4  CB  PHE H   1       7.767   5.853   7.671  1.00 15.00           C
ATOM      5  CG  PHE H   1       6.935   5.032   8.622  1.00 15.00           C
ATOM      6  CD1 PHE H   1       5.918   4.176   8.140  1.00 15.00           C
ATOM      7  CD2 PHE H   1       7.161   5.107  10.012  1.00 15.00           C
ATOM      8  CE1 PHE H   1       5.126   3.395   9.038  1.00 15.00           C
ATOM      9  CE2 PHE H   1       6.382   4.336  10.930  1.00 15.00           C
ATOM     10  CZ  PHE H   1       5.360   3.476  10.439  1.00 15.00           C
ATOM     11  HZ  PHE H   1       4.849   2.971  11.029  1.00 15.00           H
ATOM     12  HE2 PHE H   1       6.542   4.396  11.844  1.00 15.00           H
ATOM     13  HD2 PHE H   1       7.828   5.668  10.337  1.00 15.00           H
ATOM     14  HD1 PHE H   1       5.762   4.120   7.225  1.00 15.00           H
ATOM     16  HE1 PHE H   1       4.460   2.836   8.708  1.00 15.00           H
ATOM     17  HA  PHE H   1       6.622   7.575   7.683  1.00 15.00           H
ATOM     18  HB1 PHE H   1       8.101   5.258   6.982  1.00 15.00           H
ATOM     19  HB2 PHE H   1       8.507   6.227   8.174  1.00 15.00           H
TER
"""

pdb_str_6 = """
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
SCALE1      0.009182  0.000000  0.000000        0.00000
SCALE2      0.000000  0.009182  0.000000        0.00000
SCALE3      0.000000  0.000000  0.009182        0.00000
ATOM      1  N   CYS A  48     -40.561 -16.057   6.830  1.00 99.72           N
ATOM      2  CA  CYS A  48     -41.351 -14.915   6.380  1.00106.82           C
ATOM      3  C   CYS A  48     -41.146 -13.719   7.305  1.00104.36           C
ATOM      4  O   CYS A  48     -40.407 -13.794   8.283  1.00104.56           O
ATOM      5  CB  CYS A  48     -42.835 -15.275   6.303  1.00105.79           C
ATOM      6  SG  CYS A  48     -43.707 -15.186   7.879  1.00101.00           S
ATOM      8  HA  CYS A  48     -41.063 -14.664   5.488  1.00106.82           H
ATOM      9  HB1 CYS A  48     -43.271 -14.661   5.692  1.00105.79           H
ATOM     11  HB2 CYS A  48     -42.916 -16.184   5.974  1.00105.79           H
TER
HETATM   12 ZN    ZN A 102     -44.322 -17.446   8.351  1.00 98.67          Zn
END
"""

pdb_str_7 = """
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
SCALE1      0.009182  0.000000  0.000000        0.00000
SCALE2      0.000000  0.009182  0.000000        0.00000
SCALE3      0.000000  0.000000  0.009182        0.00000
ATOM      1  N   CYS A  12     -50.777 -13.205 -12.246  1.00 95.50           N
ATOM      2  CA  CYS A  12     -51.774 -13.797 -13.124  1.00103.68           C
ATOM      3  C   CYS A  12     -51.135 -14.327 -14.399  1.00101.77           C
ATOM      4  O   CYS A  12     -50.987 -13.597 -15.382  1.00106.09           O
ATOM      5  CB  CYS A  12     -52.861 -12.775 -13.466  1.00103.91           C
ATOM      6  SG  CYS A  12     -54.064 -13.357 -14.680  1.00106.39           S
ATOM      8  HA  CYS A  12     -52.194 -14.542 -12.667  1.00103.68           H
ATOM      9  HB1 CYS A  12     -52.438 -11.981 -13.828  1.00103.91           H
ATOM     11  HB2 CYS A  12     -53.344 -12.553 -12.655  1.00103.91           H
TER
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
