from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

# ------------------------------------------------------------------------------

def run():
  test_000()
  test_001()
  test_002()
  test_003()
  test_004()
  test_005()
  test_006()
  test_007()
  test_008()
  test_009()

# ------------------------------------------------------------------------------

def compare_models(pdb_str,
                   contains     = None,
                   not_contains = None):
  '''
    Function to compare model with new H to the known answer (pdb_str)
  '''
  #
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  # initial model
  model_initial = mmtbx.model.manager(model_input = pdb_inp,
                                      log         = null_out())
  ph_initial = model_initial.get_hierarchy()
  hd_sel_initial = model_initial.get_hd_selection()
  h_atoms_initial = ph_initial.select(hd_sel_initial).atoms()
  h_names_initial = list(h_atoms_initial.extract_name())
  # number of H in pdb string (right answer)
  number_h_expected = hd_sel_initial.count(True)

  # get model obj without H atoms
  model_without_h = model_initial.select(~hd_sel_initial)
  # make sure model without H indeed has no H atoms
  hd_sel_without_h = model_without_h.get_hd_selection()
  assert (hd_sel_without_h is not None)
  assert (hd_sel_without_h.count(True) == 0)

  #model_h_added = reduce.add(model = model_without_h)
  # place H atoms again
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model_without_h)
  try:
    reduce_add_h_obj.run()
  except Exception as e:
    assert 0
  #
  model_h_added = reduce_add_h_obj.get_model()
  hd_sel_h_added = model_h_added.get_hd_selection()

  # For debugging
  if 0:
    f = open("m_initial.pdb","w")
    f.write(model_initial.model_as_pdb())
    f.close()
    f = open("m_initial.cif","w")
    f.write(model_initial.model_as_mmcif())
    f.close()
    f = open("m_added.pdb","w")
    f.write(model_h_added.model_as_pdb())
    f.close()
    f = open("m_added.cif","w")
    f.write(model_h_added.model_as_mmcif())
    f.close()

  ph_h_added = model_h_added.get_hierarchy()
  assert ph_initial.is_similar_hierarchy(other=ph_h_added)

  number_h_added = hd_sel_h_added.count(True)
  assert(number_h_expected == number_h_added)

  h_atoms_added = ph_h_added.select(hd_sel_h_added).atoms()
  h_names_added = list(h_atoms_added.extract_name())

  if not_contains:
    assert (not_contains not in h_names_added)
  if contains:
    assert (contains in h_names_added)

  sc_h_initial = model_initial.select(hd_sel_initial).get_sites_cart()
  sc_h_added = model_h_added.select(hd_sel_h_added).get_sites_cart()

  d1 = {h_names_initial[i]: sc_h_initial[i] for i in range(len(h_names_initial))}
  d2 = {h_names_added[i]: sc_h_added[i] for i in range(len(h_names_added))}

  # check if coordinates are correct
  for name, sc in d2.items():
    assert(name in d1)
    assert approx_equal(sc, d1[name], 0.01)

# ------------------------------------------------------------------------------

def test_000():
  '''
    Incomplete peptide unit at the N-terminal --> H should not be placed.
  '''
  compare_models(pdb_str = pdb_str_000)

# ------------------------------------------------------------------------------

def test_001():
  '''
    CD1 is missing --> H, HD1, HE1 can't be parameterized --> should not be added.
  '''
  compare_models(pdb_str      = pdb_str_001,
                 not_contains = ' HE1')

# ------------------------------------------------------------------------------

def test_002():
  '''
    If there is a disulfide bond --> don't add H on CYS
  '''
  compare_models(pdb_str      = pdb_str_002,
                 not_contains = ' HG ')

# ------------------------------------------------------------------------------

def test_003():
  '''
    If CYS is not in disulfide or does not coordinate --> place HG
  '''
  compare_models(pdb_str  = pdb_str_003,
                 contains = ' HG ')

# ------------------------------------------------------------------------------

def test_004():
  '''
    N atom is missing --> make sure peptide H is not placed
  '''
  compare_models(pdb_str      = pdb_str_004,
                 not_contains = ' H  ')

# ------------------------------------------------------------------------------

def test_005():
  '''
    CYS coordinates a metal --> HG should not be placed
  '''
  compare_models(pdb_str      = pdb_str_005,
                 not_contains = ' HG ')

# ------------------------------------------------------------------------------

def test_006():
  '''
    CYS forms a disulfid with symmetry mate --> don't place HG
  '''
  compare_models(pdb_str      = pdb_str_006,
                 not_contains = ' HG ')

# ------------------------------------------------------------------------------

def test_007():
  '''
    CYS forms a disulfid with symmetry mate --> don't place HG
  '''
  compare_models(pdb_str      = pdb_str_007,
                 not_contains = ' HG ')

# ------------------------------------------------------------------------------

def test_008():
  '''
    disulfide bridge between neighboring residues --> don't place HG atoms
  '''
  compare_models(pdb_str      = pdb_str_008,
                 not_contains = ' HG ')

# ------------------------------------------------------------------------------

def test_009():
  '''
    Incomplete peptide unit at the N-terminal --> H should not be placed.
  '''
  compare_models(pdb_str = pdb_str_009)

# ------------------------------------------------------------------------------

pdb_str_000 = """
REMARK two AA fragments, GLY peptide H atom is not expected to be placed
CRYST1   14.074   14.385   15.410  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  H1  GLY A   1      -8.818   5.479   6.162  1.00 16.77           H
ATOM      6  H2  GLY A   1      -9.801   4.456   6.477  1.00 16.77           H
ATOM      7  H3  GLY A   1      -8.383   4.140   6.523  1.00 16.77           H
ATOM      8  HA3 GLY A   1      -9.929   3.858   4.426  1.00 16.57           H
ATOM      9  HA2 GLY A   1      -8.861   4.970   4.084  1.00 16.57           H
ATOM     10  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM     11  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM     12  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM     13  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM     14  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     15  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     16  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     17  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     18  H   ASN A   2      -8.045   3.269   2.470  1.00 15.02           H
ATOM     19  HA  ASN A   2      -6.719   1.165   3.206  1.00 14.10           H
ATOM     20  HB2 ASN A   2      -6.148   2.746   0.949  1.00 15.38           H
ATOM     21  HB3 ASN A   2      -5.618   1.264   1.168  1.00 15.38           H
ATOM     22 HD21 ASN A   2      -8.919   1.893  -0.569  1.00 11.72           H
ATOM     23 HD22 ASN A   2      -7.888   2.940  -0.323  1.00 11.72           H
TER
"""

pdb_str_001 = """
REMARK CD1 is missing --> HD1, HE1 cannot not be placed; H is not placed
CRYST1   17.955   13.272   13.095  90.00  90.00  90.00 P 1
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
ATOM     12  H   TYR A 139       9.382   7.879   5.001  1.00 10.00           H
ATOM     13  HA  TYR A 139      10.487   8.115   6.973  1.00 10.00           H
ATOM     14  HB2 TYR A 139      10.961   5.881   7.464  1.00 10.00           H
ATOM     15  HB3 TYR A 139      10.893   5.529   5.916  1.00 10.00           H
ATOM     16  HD2 TYR A 139       8.896   6.220   8.714  1.00 10.00           H
ATOM     17  HE2 TYR A 139       6.643   5.772   8.919  1.00 10.00           H
ATOM     18  HH  TYR A 139       4.752   5.127   7.905  1.00 10.00           H
TER
"""

pdb_str_002 = """
REMARK fragment of disulfide bridge --> don't place H on Cys SG atom
CRYST1   14.197   15.507   16.075  90.00  90.00  90.00 P 1
ATOM      1  N   CYS A 128      14.558 -30.378 -19.478  1.00 11.73           N
ATOM      2  CA  CYS A 128      13.499 -29.511 -19.901  1.00 11.51           C
ATOM      3  C   CYS A 128      12.454 -30.316 -20.640  1.00 11.99           C
ATOM      4  O   CYS A 128      12.765 -31.329 -21.265  1.00 16.86           O
ATOM      5  CB  CYS A 128      14.059 -28.357 -20.691  1.00 12.94           C
ATOM      6  SG  CYS A 128      15.213 -27.350 -19.760  1.00 12.36           S
ATOM      7  H   CYS A 128      14.993 -30.091 -18.794  1.00 11.73           H
ATOM      8  HA  CYS A 128      13.044 -29.102 -19.148  1.00 11.51           H
ATOM      9  HB2 CYS A 128      14.526 -28.706 -21.466  1.00 12.94           H
ATOM     10  HB3 CYS A 128      13.327 -27.787 -20.974  1.00 12.94           H
ATOM     11  N   CYS A 232      13.105 -25.822 -15.371  1.00  4.75           N
ATOM     12  CA  CYS A 232      14.298 -26.646 -15.573  1.00  5.50           C
ATOM     14  C   CYS A 232      15.596 -25.838 -15.551  1.00  5.70           C
ATOM     15  O   CYS A 232      16.651 -26.346 -15.190  1.00  8.25           O
ATOM     16  CB  CYS A 232      14.178 -27.497 -16.824  1.00  6.69           C
ATOM     17  SG  CYS A 232      14.060 -26.477 -18.360  1.00  7.70           S
ATOM     18  H   CYS A 232      13.275 -24.997 -15.197  1.00  4.75           H
ATOM     19  HA  CYS A 232      14.363 -27.258 -14.823  1.00  5.50           H
ATOM     20  HB2 CYS A 232      14.961 -28.065 -16.898  1.00  6.69           H
ATOM     21  HB3 CYS A 232      13.378 -28.042 -16.761  1.00  6.69           H
TER
"""

pdb_str_003 = """
REMARK fragment of lone CYS --> place the HG atom
CRYST1   14.197   15.507   16.075  90.00  90.00  90.00 P 1
ATOM      1  N   CYS A 128      14.558 -30.378 -19.478  1.00 11.73           N
ATOM      2  CA  CYS A 128      13.499 -29.511 -19.901  1.00 11.51           C
ATOM      3  C   CYS A 128      12.454 -30.316 -20.640  1.00 11.99           C
ATOM      4  O   CYS A 128      12.765 -31.329 -21.265  1.00 16.86           O
ATOM      5  CB  CYS A 128      14.059 -28.357 -20.691  1.00 12.94           C
ATOM      6  SG  CYS A 128      15.213 -27.350 -19.760  1.00 12.36           S
ATOM      7  H   CYS A 128      14.993 -30.091 -18.794  1.00 11.73           H
ATOM      8  HA  CYS A 128      13.044 -29.102 -19.148  1.00 11.51           H
ATOM      9  HB2 CYS A 128      14.526 -28.706 -21.466  1.00 12.94           H
ATOM     10  HB3 CYS A 128      13.327 -27.787 -20.974  1.00 12.94           H
ATOM     11  HG  CYS A 128      16.134 -28.034 -19.410  1.00 12.36           H
TER
"""

pdb_str_004 = """
REMARK peptide N is missing --> peptide H atom should not be placed
CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 1
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
ATOM     18  HB3 PHE H   1       8.101   5.258   6.982  1.00 15.00           H
ATOM     19  HB2 PHE H   1       8.507   6.227   8.174  1.00 15.00           H
TER
"""

pdb_str_005 = """
REMARK CYS coordinates a metal --> HG should not be placed
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
ATOM      1  N   CYS A  48     -40.561 -16.057   6.830  1.00 99.72           N
ATOM      2  CA  CYS A  48     -41.351 -14.915   6.380  1.00106.82           C
ATOM      3  C   CYS A  48     -41.146 -13.719   7.305  1.00104.36           C
ATOM      4  O   CYS A  48     -40.407 -13.794   8.283  1.00104.56           O
ATOM      5  CB  CYS A  48     -42.835 -15.275   6.303  1.00105.79           C
ATOM      6  SG  CYS A  48     -43.707 -15.186   7.879  1.00101.00           S
ATOM      7  H   CYS A  48     -39.952 -15.857   7.403  1.00 99.72           H
ATOM      8  HA  CYS A  48     -41.063 -14.664   5.488  1.00106.82           H
ATOM      9  HB2 CYS A  48     -43.271 -14.661   5.692  1.00105.79           H
ATOM     10  HB3 CYS A  48     -42.916 -16.184   5.974  1.00105.79           H
TER
HETATM   12 ZN    ZN A 102     -44.322 -17.446   8.351  1.00 98.67          Zn
END
"""

pdb_str_006 = """
REMARK CYS forms a disulfide bridge with symmetry mate --> don't place HG
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
ATOM      1  N   CYS A  12     -50.777 -13.205 -12.246  1.00 95.50           N
ATOM      2  CA  CYS A  12     -51.774 -13.797 -13.124  1.00103.68           C
ATOM      3  C   CYS A  12     -51.135 -14.327 -14.399  1.00101.77           C
ATOM      4  O   CYS A  12     -50.987 -13.597 -15.382  1.00106.09           O
ATOM      5  CB  CYS A  12     -52.861 -12.775 -13.466  1.00103.91           C
ATOM      6  SG  CYS A  12     -54.064 -13.357 -14.680  1.00106.39           S
ATOM      7  H   CYS A  12     -50.242 -13.776 -11.889  1.00 95.50           H
ATOM      8  HA  CYS A  12     -52.194 -14.542 -12.667  1.00103.68           H
ATOM      9  HB2 CYS A  12     -52.438 -11.981 -13.828  1.00103.91           H
ATOM     10  HB3 CYS A  12     -53.344 -12.553 -12.655  1.00103.91           H
TER
"""

pdb_str_007 = """
REMARK CYS forms a disulfide with symmetry mate --> don't place HG
REMARK from PDB 5c11
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
ATOM      1  N   PRO A   1     -50.110 -12.340  -8.990  1.00 87.52           N
ATOM      2  CA  PRO A   1     -49.518 -12.918 -10.200  1.00 88.41           C
ATOM      3  C   PRO A   1     -50.549 -13.676 -11.025  1.00 97.43           C
ATOM      4  O   PRO A   1     -51.129 -14.656 -10.556  1.00 93.89           O
ATOM      5  CB  PRO A   1     -48.446 -13.865  -9.649  1.00 78.46           C
ATOM      6  CG  PRO A   1     -48.110 -13.314  -8.317  1.00 80.77           C
ATOM      7  CD  PRO A   1     -49.402 -12.783  -7.777  1.00 92.02           C
ATOM      8  H2  PRO A   1     -51.029 -12.612  -8.932  1.00 87.52           H
ATOM      9  H3  PRO A   1     -50.063 -11.383  -9.046  1.00 87.52           H
ATOM     10  HA  PRO A   1     -49.092 -12.241 -10.748  1.00 88.41           H
ATOM     11  HB2 PRO A   1     -48.804 -14.764  -9.574  1.00 78.46           H
ATOM     12  HB3 PRO A   1     -47.671 -13.864 -10.232  1.00 78.46           H
ATOM     13  HG2 PRO A   1     -47.455 -12.604  -8.409  1.00 80.77           H
ATOM     14  HG3 PRO A   1     -47.760 -14.017  -7.748  1.00 80.77           H
ATOM     15  HD2 PRO A   1     -49.899 -13.480  -7.321  1.00 92.02           H
ATOM     16  HD3 PRO A   1     -49.246 -12.038  -7.175  1.00 92.02           H
ATOM     17  N   CYS A   2     -50.777 -13.205 -12.246  1.00 95.50           N
ATOM     18  CA  CYS A   2     -51.774 -13.797 -13.124  1.00103.68           C
ATOM     19  C   CYS A   2     -51.135 -14.327 -14.399  1.00101.77           C
ATOM     20  O   CYS A   2     -50.987 -13.597 -15.382  1.00106.09           O
ATOM     21  CB  CYS A   2     -52.861 -12.775 -13.466  1.00103.91           C
ATOM     22  SG  CYS A   2     -54.064 -13.357 -14.680  1.00106.39           S
ATOM     23  H   CYS A   2     -50.361 -12.536 -12.590  1.00 95.50           H
ATOM     24  HA  CYS A   2     -52.147 -14.558 -12.653  1.00103.68           H
ATOM     25  HB2 CYS A   2     -52.438 -11.981 -13.828  1.00103.91           H
ATOM     26  HB3 CYS A   2     -53.344 -12.553 -12.655  1.00103.91           H
ATOM     28  N   LYS A   3     -50.756 -15.600 -14.375  1.00105.98           N
ATOM     29  CA  LYS A   3     -50.135 -16.230 -15.528  1.00115.21           C
ATOM     30  C   LYS A   3     -50.110 -17.740 -15.357  1.00105.26           C
ATOM     31  O   LYS A   3     -50.331 -18.255 -14.264  1.00118.65           O
ATOM     32  CB  LYS A   3     -48.714 -15.704 -15.737  1.00110.03           C
ATOM     33  CG  LYS A   3     -47.784 -15.954 -14.563  1.00100.14           C
ATOM     34  CD  LYS A   3     -46.379 -15.452 -14.853  1.00101.77           C
ATOM     35  CE  LYS A   3     -45.456 -15.678 -13.663  1.00112.21           C
ATOM     36  NZ  LYS A   3     -44.063 -15.229 -13.941  1.00110.68           N
ATOM     37  H   LYS A   3     -50.849 -16.122 -13.698  1.00105.98           H
ATOM     38  HA  LYS A   3     -50.648 -16.015 -16.323  1.00115.21           H
ATOM     39  HB2 LYS A   3     -48.755 -14.746 -15.882  1.00110.03           H
ATOM     40  HB3 LYS A   3     -48.332 -16.141 -16.514  1.00110.03           H
ATOM     41  HG2 LYS A   3     -48.120 -15.488 -13.782  1.00100.14           H
ATOM     42  HG3 LYS A   3     -47.738 -16.907 -14.386  1.00100.14           H
ATOM     43  HD2 LYS A   3     -46.409 -14.501 -15.041  1.00101.77           H
ATOM     44  HD3 LYS A   3     -46.017 -15.929 -15.616  1.00101.77           H
ATOM     45  HE2 LYS A   3     -45.789 -15.178 -12.902  1.00112.21           H
ATOM     46  HE3 LYS A   3     -45.432 -16.625 -13.453  1.00112.21           H
ATOM     47  HZ1 LYS A   3     -43.731 -15.679 -14.633  1.00110.68           H
ATOM     48  HZ2 LYS A   3     -43.549 -15.374 -13.229  1.00110.68           H
ATOM     49  HZ3 LYS A   3     -44.056 -14.360 -14.131  1.00110.68           H
TER
"""

pdb_str_007_old = """
REMARK CYS forms a disulfide with symmetry mate --> don't place HG
REMARK from PDB 5c11
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
ATOM      1  N   PRO A   1     -50.110 -12.340  -8.990  1.00 87.52           N
ATOM      2  CA  PRO A   1     -49.518 -12.918 -10.200  1.00 88.41           C
ATOM      3  C   PRO A   1     -50.549 -13.676 -11.025  1.00 97.43           C
ATOM      4  O   PRO A   1     -51.129 -14.656 -10.556  1.00 93.89           O
ATOM      5  CB  PRO A   1     -48.446 -13.865  -9.649  1.00 78.46           C
ATOM      6  CG  PRO A   1     -48.110 -13.314  -8.317  1.00 80.77           C
ATOM      7  CD  PRO A   1     -49.402 -12.783  -7.777  1.00 92.02           C
ATOM      8  H2  PRO A   1     -50.010 -11.446  -9.004  1.00 87.52           H
ATOM      9  H3  PRO A   1     -50.986 -12.543  -8.958  1.00 87.52           H
ATOM     10  HA  PRO A   1     -49.092 -12.241 -10.748  1.00 88.41           H
ATOM     11  HB2 PRO A   1     -48.804 -14.764  -9.574  1.00 78.46           H
ATOM     12  HB3 PRO A   1     -47.671 -13.864 -10.232  1.00 78.46           H
ATOM     13  HG2 PRO A   1     -47.455 -12.604  -8.409  1.00 80.77           H
ATOM     14  HG3 PRO A   1     -47.760 -14.017  -7.748  1.00 80.77           H
ATOM     15  HD2 PRO A   1     -49.899 -13.480  -7.321  1.00 92.02           H
ATOM     16  HD3 PRO A   1     -49.246 -12.038  -7.175  1.00 92.02           H
ATOM     17  N   CYS A   2     -50.777 -13.205 -12.246  1.00 95.50           N
ATOM     18  CA  CYS A   2     -51.774 -13.797 -13.124  1.00103.68           C
ATOM     19  C   CYS A   2     -51.135 -14.327 -14.399  1.00101.77           C
ATOM     20  O   CYS A   2     -50.987 -13.597 -15.382  1.00106.09           O
ATOM     21  CB  CYS A   2     -52.861 -12.775 -13.466  1.00103.91           C
ATOM     22  SG  CYS A   2     -54.064 -13.357 -14.680  1.00106.39           S
ATOM     23  H   CYS A   2     -50.361 -12.536 -12.590  1.00 95.50           H
ATOM     24  HA  CYS A   2     -52.147 -14.558 -12.653  1.00103.68           H
ATOM     25  HB2 CYS A   2     -52.438 -11.981 -13.828  1.00103.91           H
ATOM     26  HB3 CYS A   2     -53.344 -12.553 -12.655  1.00103.91           H
ATOM     28  N   LYS A   3     -50.756 -15.600 -14.375  1.00105.98           N
ATOM     29  CA  LYS A   3     -50.135 -16.230 -15.528  1.00115.21           C
ATOM     30  C   LYS A   3     -50.110 -17.740 -15.357  1.00105.26           C
ATOM     31  O   LYS A   3     -50.331 -18.255 -14.264  1.00118.65           O
ATOM     32  CB  LYS A   3     -48.714 -15.704 -15.737  1.00110.03           C
ATOM     33  CG  LYS A   3     -47.784 -15.954 -14.563  1.00100.14           C
ATOM     34  CD  LYS A   3     -46.379 -15.452 -14.853  1.00101.77           C
ATOM     35  CE  LYS A   3     -45.456 -15.678 -13.663  1.00112.21           C
ATOM     36  NZ  LYS A   3     -44.063 -15.229 -13.941  1.00110.68           N
ATOM     37  H   LYS A   3     -50.849 -16.122 -13.698  1.00105.98           H
ATOM     38  HA  LYS A   3     -50.648 -16.015 -16.323  1.00115.21           H
ATOM     39  HB2 LYS A   3     -48.755 -14.746 -15.882  1.00110.03           H
ATOM     40  HB3 LYS A   3     -48.332 -16.141 -16.514  1.00110.03           H
ATOM     41  HG2 LYS A   3     -48.120 -15.488 -13.782  1.00100.14           H
ATOM     42  HG3 LYS A   3     -47.738 -16.907 -14.386  1.00100.14           H
ATOM     43  HD2 LYS A   3     -46.409 -14.501 -15.041  1.00101.77           H
ATOM     44  HD3 LYS A   3     -46.017 -15.929 -15.616  1.00101.77           H
ATOM     45  HE2 LYS A   3     -45.789 -15.178 -12.902  1.00112.21           H
ATOM     46  HE3 LYS A   3     -45.432 -16.625 -13.453  1.00112.21           H
ATOM     47  HZ1 LYS A   3     -43.731 -15.679 -14.633  1.00110.68           H
ATOM     48  HZ2 LYS A   3     -43.549 -15.374 -13.229  1.00110.68           H
ATOM     49  HZ3 LYS A   3     -44.056 -14.360 -14.131  1.00110.68           H
TER
"""

pdb_str_008 = """
REMARK disulfide bridge between neighboring residues --> don't place HG atoms
REMARK from PDB 3cu9
CRYST1   86.141   89.579   76.041  90.00  90.00  90.00 C 2 2 21
ATOM      1  N   CYS A 221     -25.409   6.530   5.604  1.00  8.50           N
ATOM      2  CA  CYS A 221     -25.611   7.909   5.995  1.00  8.62           C
ATOM      3  C   CYS A 221     -25.469   8.835   4.795  1.00  8.84           C
ATOM      4  O   CYS A 221     -25.131   8.426   3.672  1.00  9.51           O
ATOM      5  CB  CYS A 221     -24.546   8.272   7.047  1.00  9.34           C
ATOM      6  SG  CYS A 221     -22.908   8.520   6.302  1.00 10.59           S
ATOM      7  H   CYS A 221     -26.072   6.187   5.176  1.00  8.50           H
ATOM      8  HA  CYS A 221     -26.500   8.030   6.365  1.00  8.62           H
ATOM      9  HB2 CYS A 221     -24.805   9.094   7.491  1.00  9.34           H
ATOM     10  HB3 CYS A 221     -24.479   7.552   7.694  1.00  9.34           H
ATOM     11  N   CYS A 222     -25.649  10.119   5.071  1.00  8.82           N
ATOM     12  CA  CYS A 222     -25.005  11.177   4.321  1.00  8.99           C
ATOM     14  C   CYS A 222     -25.532  11.323   2.911  1.00  9.12           C
ATOM     15  O   CYS A 222     -24.803  11.847   2.036  1.00 11.09           O
ATOM     16  CB  CYS A 222     -23.463  10.894   4.313  1.00 11.23           C
ATOM     17  SG  CYS A 222     -22.792  10.598   5.958  1.00 11.91           S
ATOM     18  H   CYS A 222     -26.153  10.408   5.705  1.00  8.82           H
ATOM     19  HA  CYS A 222     -25.188  12.029   4.746  1.00  8.99           H
ATOM     20  HB2 CYS A 222     -23.004  11.661   3.936  1.00 11.23           H
ATOM     21  HB3 CYS A 222     -23.291  10.107   3.773  1.00 11.23           H
TER
"""

pdb_str_009 = """
CRYST1   22.251   49.068   23.478  90.00  90.00  90.00 P 1
ATOM      1  N   HIS A  44       5.727  38.590   5.000  1.00  7.72           N
ATOM      2  CA  HIS A  44       6.991  38.641   5.686  1.00  8.00           C
ATOM      3  C   HIS A  44       7.274  37.329   6.472  1.00  8.78           C
ATOM      4  O   HIS A  44       6.403  36.603   6.844  1.00  9.64           O
ATOM      5  CB  HIS A  44       6.917  39.719   6.794  1.00 11.03           C
ATOM      6  CG  HIS A  44       6.477  41.058   6.353  1.00  9.11           C
ATOM      7  ND1 HIS A  44       5.698  41.854   7.162  1.00  9.47           N
ATOM      8  CD2 HIS A  44       6.675  41.756   5.230  1.00 10.66           C
ATOM      9  CE1 HIS A  44       5.476  43.005   6.544  1.00 11.24           C
ATOM     10  NE2 HIS A  44       6.065  42.976   5.358  1.00  9.90           N
ATOM     11  H   HIS A  44       5.757  38.876   4.190  1.00  7.72           H
ATOM     12  HA  HIS A  44       7.678  38.811   5.022  1.00  8.00           H
ATOM     13  HB2 HIS A  44       7.801  39.819   7.181  1.00 11.03           H
ATOM     14  HB3 HIS A  44       6.290  39.418   7.471  1.00 11.03           H
ATOM     15  HD1 HIS A  44       5.404  41.639   7.941  1.00  9.47           H
ATOM     16  HD2 HIS A  44       7.148  41.462   4.485  1.00 10.66           H
ATOM     17  HE1 HIS A  44       4.986  43.717   6.887  1.00 11.24           H
ATOM     18  HE2 HIS A  44       6.065  43.608   4.775  1.00  9.90           H
ATOM     19  N   SER A  60       8.719  39.746  15.915  1.00  9.28           N
ATOM     20  CA  SER A  60       8.995  40.756  14.888  1.00  8.34           C
ATOM     21  C   SER A  60       9.663  40.135  13.685  1.00  9.31           C
ATOM     22  O   SER A  60       9.754  38.932  13.513  1.00  9.78           O
ATOM     23  CB ASER A  60       9.760  41.944  15.459  0.50 11.89           C
ATOM     24  OG ASER A  60       9.077  42.510  16.590  0.33 10.99           O
ATOM     25  H  ASER A  60       7.932  39.801  16.258  1.00  9.28           H
ATOM     26  HA ASER A  60       8.159  41.130  14.569  1.00  8.34           H
ATOM     27  HB2ASER A  60      10.639  41.645  15.740  0.50 11.89           H
ATOM     28  HB3ASER A  60       9.844  42.623  14.771  0.50 11.89           H
ATOM     29  HG ASER A  60       9.512  43.162  16.891  0.33 10.99           H
ATOM     30  CB BSER A  60       9.680  41.954  15.550  0.50 12.25           C
ATOM     31  OG BSER A  60      10.939  41.568  16.145  0.66 17.46           O
ATOM     32  H  BSER A  60       7.882  39.617  16.063  1.00  9.28           H
ATOM     33  HA BSER A  60       8.196  41.138  14.493  1.00  8.34           H
ATOM     34  HB2BSER A  60       9.099  42.306  16.243  0.50 12.25           H
ATOM     35  HB3BSER A  60       9.845  42.634  14.878  0.50 12.25           H
ATOM     36  HG BSER A  60      11.303  42.234  16.504  0.66 17.46           H
ATOM     37  N   HIS A  61      10.095  40.969  12.694  1.00  7.19           N
ATOM     38  CA  HIS A  61      10.687  40.443  11.504  1.00  7.24           C
ATOM     39  C   HIS A  61      12.059  39.821  11.742  1.00  7.38           C
ATOM     40  O   HIS A  61      12.874  40.275  12.498  1.00  9.97           O
ATOM     41  CB  HIS A  61      10.837  41.605  10.486  1.00  9.38           C
ATOM     42  CG  HIS A  61       9.512  42.014   9.909  1.00  8.78           C
ATOM     43  ND1 HIS A  61       9.298  43.039   9.081  1.00  9.06           N
ATOM     44  CD2 HIS A  61       8.304  41.455  10.105  1.00 11.08           C
ATOM     45  CE1 HIS A  61       7.997  43.118   8.789  1.00 12.84           C
ATOM     46  NE2 HIS A  61       7.381  42.061   9.398  1.00 11.66           N
ATOM     47  H   HIS A  61      10.043  41.827  12.717  1.00  7.19           H
ATOM     48  HA  HIS A  61      10.119  39.738  11.156  1.00  7.24           H
ATOM     49  HB2 HIS A  61      11.225  42.374  10.933  1.00  9.38           H
ATOM     50  HB3 HIS A  61      11.411  41.320   9.758  1.00  9.38           H
ATOM     51  HD1 HIS A  61       9.906  43.568   8.781  1.00  9.06           H
ATOM     52  HD2 HIS A  61       8.144  40.733  10.669  1.00 11.08           H
ATOM     53  HE1 HIS A  61       7.587  43.770   8.268  1.00 12.84           H
ATOM     54  HE2 HIS A  61       6.553  41.838   9.331  1.00 11.66           H
ATOM     55  N   PHE A  62      12.300  38.781  10.945  1.00  7.37           N
ATOM     56  CA  PHE A  62      13.573  38.090  10.989  1.00  7.69           C
ATOM     57  C   PHE A  62      14.659  39.009  10.437  1.00  8.66           C
ATOM     58  O   PHE A  62      14.567  39.469   9.302  1.00  8.67           O
ATOM     59  CB  PHE A  62      13.468  36.801  10.201  1.00  9.12           C
ATOM     60  CG  PHE A  62      14.736  35.956  10.204  1.00  8.46           C
ATOM     61  CD1 PHE A  62      15.452  35.723  11.369  1.00 10.92           C
ATOM     62  CD2 PHE A  62      15.186  35.393   9.017  1.00 12.32           C
ATOM     63  CE1 PHE A  62      16.602  34.930  11.343  1.00 13.49           C
ATOM     64  CE2 PHE A  62      16.297  34.631   8.988  1.00 14.14           C
ATOM     65  CZ  PHE A  62      16.979  34.320  10.172  1.00 13.67           C
ATOM     66  H   PHE A  62      11.741  38.462  10.375  1.00  7.37           H
ATOM     67  HA  PHE A  62      13.822  37.850  11.895  1.00  7.69           H
ATOM     68  HB2 PHE A  62      13.263  37.018   9.278  1.00  9.12           H
ATOM     69  HB3 PHE A  62      12.756  36.263  10.581  1.00  9.12           H
ATOM     70  HD1 PHE A  62      15.165  36.096  12.171  1.00 10.92           H
ATOM     71  HD2 PHE A  62      14.714  35.544   8.230  1.00 12.32           H
ATOM     72  HE1 PHE A  62      17.109  34.816  12.114  1.00 13.49           H
ATOM     73  HE2 PHE A  62      16.614  34.309   8.175  1.00 14.14           H
ATOM     74  HZ  PHE A  62      17.679  33.708  10.165  1.00 13.67           H
TER
HETATM   75 CU    CU A   1       5.421  41.446   9.097  1.00 13.22          Cu
HETATM   76 ZN    ZN A 152      10.669  44.068   7.982  1.00  9.34          Zn
HETATM   77  O   HOH A 167       5.000  43.153  10.522  1.00 24.06           O
HETATM   78  O   HOH A 204      13.077  42.283  14.236  1.00 35.24           O
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
