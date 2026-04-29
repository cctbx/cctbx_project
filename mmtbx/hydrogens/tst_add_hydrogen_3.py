from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models
from libtbx.utils import null_out

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

# ------------------------------------------------------------------------------

def test_000():
  '''
    Make sure reduce does not crash for single_atom_residue models
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_000.split("\n"), source_info=None)
  # initial model (has no H atoms)
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  number_h_expected = model_initial.get_hd_selection().count(True)
  assert(number_h_expected == 0)
  # place H atoms
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model_initial)
  reduce_add_h_obj.run()
  # We don't expect H atoms to be placed
  # (not enough restraints for single atom residues)
  model_h_added = reduce_add_h_obj.get_model()
  number_h_placed = model_h_added.get_hd_selection().count(True)
  assert(number_h_placed == 0)

# ------------------------------------------------------------------------------

def test_001():
  '''
    Check keyword n_terminal_charge:
    NH3 on resseq 1, first residue in chain, or no NH3 at all
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_001.split("\n"), source_info=None)
  # initial model
  model_initial = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  #
  # place H atoms: NH3 at resseq 1 only
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model_initial)
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  #
  hd_sel_h_added = model_h_added.get_hd_selection()
  ph_h_added     = model_h_added.get_hierarchy()
  h_atoms_added = ph_h_added.select(hd_sel_h_added).atoms()
  h_names_added = list(h_atoms_added.extract_name())
  assert(h_names_added.count(' H1 ')==1)
  assert(h_names_added.count(' H2 ')==1)
  assert(h_names_added.count(' H3 ')==1)
  #
  # place H atoms: NH3 at first residue in chain
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    n_terminal_charge = 'first_in_chain')
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  #
  hd_sel_h_added = model_h_added.get_hd_selection()
  ph_h_added     = model_h_added.get_hierarchy()
  h_atoms_added = ph_h_added.select(hd_sel_h_added).atoms()
  h_names_added = list(h_atoms_added.extract_name())
  assert(h_names_added.count(' H1 ')==3)
  assert(h_names_added.count(' H2 ')==3)
  assert(h_names_added.count(' H3 ')==3)
  #
  # place H atoms: no NH3
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model = model_initial,
    n_terminal_charge = 'no_charge')
  reduce_add_h_obj.run()
  model_h_added = reduce_add_h_obj.get_model()
  #
  hd_sel_h_added = model_h_added.get_hd_selection()
  ph_h_added     = model_h_added.get_hierarchy()
  h_atoms_added = ph_h_added.select(hd_sel_h_added).atoms()
  h_names_added = list(h_atoms_added.extract_name())
  assert(h_names_added.count(' H1 ')==0)
  assert(h_names_added.count(' H2 ')==0)
  assert(h_names_added.count(' H3 ')==0)

# ------------------------------------------------------------------------------

def test_002():
  '''
    SIN forms covalent link to GLU 1. Make sure the default NH3  at the
    N-terminal becomes a single peptide H in this particular scenario.
  '''
  compare_models(pdb_str = pdb_str_002)

# ------------------------------------------------------------------------------

def test_003():
  '''
    Carbohydrate forms covlalent link to ASN. Make sure valence is correct for
    the NAG and protein links. Both parts are in double conformation.
  '''
  compare_models(pdb_str = pdb_str_003)

# ------------------------------------------------------------------------------

def test_004():
  '''
    Carbohydrate forms covlalent link to ASN. Make sure valence is correct for
    the NAG and protein links.
  '''
  compare_models(pdb_str = pdb_str_004)

# ------------------------------------------------------------------------------

def test_005():
  '''
    Carbohydrate (BMA) forms covlalent link to ASN. Make sure valence is correct for
    the NAG and protein links. This example is from Russel (6vw1).
  '''
  compare_models(pdb_str = pdb_str_005)

# ------------------------------------------------------------------------------

def test_006():
  '''
    Ligand in double conformation. If atom v3 naming is not restricted to aa's,
    then two atoms in the ligand will be renamed, leading to seemingly missing
    atoms; then pH dependent restraints kick in and they will create
    nonsensical restraints involving A and B conformers, which throws off riding H.
  '''
  compare_models(pdb_str = pdb_str_006)

# ------------------------------------------------------------------------------

def test_007():
  '''
    MET in triple conformation, with split at CAlpha (i.e. N-C-O have altloc
    blank). Make sure that H atoms are added at each alt conf.
  '''
  compare_models(pdb_str = pdb_str_007)

# ------------------------------------------------------------------------------

pdb_str_000 = """
REMARK Make sure reduce does not crash for single_atom_residue models
CRYST1   22.029   33.502   24.035  90.00  90.00  90.00 P 1
ATOM      6  P     A A  10     -62.272  56.445  13.820  1.00 15.00           P
ATOM      7  P     G A  11     -63.673  51.410  11.026  1.00 15.00           P
ATOM      8  P     U A  12     -62.888  45.926   9.711  1.00 15.00           P
ATOM      9  P     U A  13     -60.326  41.305  11.244  1.00 15.00           P
ATOM     10  P     U A  14     -57.909  36.481  13.207  1.00 15.00           P
ATOM     11  P     G A  15     -62.106  32.943  15.800  1.00 15.00           P
ATOM     12  P     A A  16     -65.446  37.240  15.291  1.00 15.00           P
ATOM     13  P     U A  17     -66.286  42.354  18.232  1.00 15.00           P
ATOM     14  P     C A  18     -64.629  46.517  21.258  1.00 15.00           P
ATOM     15  P     A A  19     -60.460  50.019  23.746  1.00 15.00           P
ATOM     16  P     U A  20     -54.257  51.133  23.481  1.00 15.00           P
"""

pdb_str_001 = """
REMARK Make sure NH3 is applied correctly depending on keyword
CRYST1   30.200   47.800   61.300  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   TYR A   5       8.831  48.837  54.788  1.00  9.67           N
ATOM      2  CA  TYR A   5       7.706  49.436  54.084  1.00  9.24           C
ATOM      3  C   TYR A   5       6.456  48.981  54.822  1.00 10.02           C
ATOM      4  O   TYR A   5       6.310  47.784  55.139  1.00 10.62           O
ATOM      5  CB  TYR A   5       7.599  48.942  52.633  1.00  9.83           C
ATOM      6  CG  TYR A   5       8.692  49.472  51.736  1.00 12.79           C
ATOM      7  CD1 TYR A   5       9.959  48.898  51.781  1.00 13.60           C
ATOM      8  CD2 TYR A   5       8.415  50.537  50.880  1.00 12.13           C
ATOM      9  CE1 TYR A   5      10.961  49.400  50.960  1.00 14.77           C
ATOM     10  CE2 TYR A   5       9.426  51.032  50.065  1.00 14.21           C
ATOM     11  CZ  TYR A   5      10.685  50.467  50.116  1.00 14.05           C
ATOM     12  OH  TYR A   5      11.708  50.978  49.318  1.00 17.48           O
ATOM     13  H1  TYR A   5       9.284  48.311  54.231  1.00  9.67           H
ATOM     14  H2  TYR A   5       9.367  49.479  55.092  1.00  9.67           H
ATOM     15  H3  TYR A   5       8.530  48.354  55.472  1.00  9.67           H
ATOM     16  HA  TYR A   5       7.805  50.400  54.049  1.00  9.24           H
ATOM     17  HB2 TYR A   5       7.653  47.974  52.627  1.00  9.83           H
ATOM     18  HB3 TYR A   5       6.748  49.229  52.266  1.00  9.83           H
ATOM     19  HD1 TYR A   5      10.133  48.187  52.354  1.00 13.60           H
ATOM     20  HD2 TYR A   5       7.564  50.912  50.855  1.00 12.13           H
ATOM     21  HE1 TYR A   5      11.811  49.024  50.976  1.00 14.77           H
ATOM     22  HE2 TYR A   5       9.255  51.741  49.488  1.00 14.21           H
ATOM     23  HH  TYR A   5      12.418  50.550  49.452  1.00 17.48           H
TER
ATOM     24  N   GLU B   1      14.684  52.510  58.829  1.00 14.47           N
ATOM     25  CA  GLU B   1      13.950  53.593  58.225  1.00 15.28           C
ATOM     26  C   GLU B   1      12.566  53.009  57.937  1.00 14.74           C
ATOM     27  O   GLU B   1      12.459  51.916  57.371  1.00 14.27           O
ATOM     28  CB  GLU B   1      14.667  53.949  56.967  1.00 18.89           C
ATOM     29  CG  GLU B   1      13.973  54.945  56.083  1.00 27.57           C
ATOM     30  CD  GLU B   1      14.729  55.326  54.802  1.00 32.66           C
ATOM     31  OE1 GLU B   1      15.802  54.776  54.508  1.00 34.50           O
ATOM     32  OE2 GLU B   1      14.224  56.201  54.090  1.00 36.48           O
ATOM     33  H1  GLU B   1      15.398  52.322  58.332  1.00 14.47           H
ATOM     34  H2  GLU B   1      14.945  52.748  59.646  1.00 14.47           H
ATOM     35  H3  GLU B   1      14.162  51.791  58.882  1.00 14.47           H
ATOM     36  HA  GLU B   1      13.870  54.399  58.759  1.00 15.28           H
ATOM     37  HB2 GLU B   1      15.529  54.326  57.204  1.00 18.89           H
ATOM     38  HB3 GLU B   1      14.790  53.139  56.447  1.00 18.89           H
ATOM     39  HG2 GLU B   1      13.835  55.760  56.590  1.00 27.57           H
ATOM     40  HG3 GLU B   1      13.119  54.573  55.814  1.00 27.57           H
TER
ATOM     41  N   PHE C  -3       6.984  40.342  58.778  1.00  8.80           N
ATOM     42  CA  PHE C  -3       7.384  40.247  60.166  1.00  8.32           C
ATOM     43  C   PHE C  -3       7.719  38.788  60.513  1.00  9.15           C
ATOM     44  O   PHE C  -3       8.710  38.555  61.201  1.00  9.31           O
ATOM     45  CB  PHE C  -3       6.284  40.754  61.091  1.00  8.94           C
ATOM     46  CG  PHE C  -3       6.705  40.642  62.560  1.00  9.84           C
ATOM     47  CD1 PHE C  -3       7.825  41.288  63.022  1.00 10.62           C
ATOM     48  CD2 PHE C  -3       5.989  39.828  63.426  1.00 12.63           C
ATOM     49  CE1 PHE C  -3       8.229  41.132  64.328  1.00 11.20           C
ATOM     50  CE2 PHE C  -3       6.398  39.679  64.737  1.00 13.74           C
ATOM     51  CZ  PHE C  -3       7.527  40.326  65.197  1.00 12.55           C
ATOM     52  H1  PHE C  -3       6.155  40.662  58.729  1.00  8.80           H
ATOM     53  H2  PHE C  -3       7.538  40.888  58.346  1.00  8.80           H
ATOM     54  H3  PHE C  -3       7.013  39.534  58.405  1.00  8.80           H
ATOM     55  HA  PHE C  -3       8.167  40.801  60.311  1.00  8.32           H
ATOM     56  HB2 PHE C  -3       6.101  41.686  60.894  1.00  8.94           H
ATOM     57  HB3 PHE C  -3       5.483  40.224  60.960  1.00  8.94           H
ATOM     58  HD1 PHE C  -3       8.313  41.834  62.449  1.00 10.62           H
ATOM     59  HD2 PHE C  -3       5.232  39.381  63.123  1.00 12.63           H
ATOM     60  HE1 PHE C  -3       8.988  41.578  64.629  1.00 11.20           H
ATOM     61  HE2 PHE C  -3       5.909  39.138  65.314  1.00 13.74           H
ATOM     62  HZ  PHE C  -3       7.808  40.220  66.077  1.00 12.55           H
TER
END
"""

pdb_str_002 = """
REMARK Make sure linking and valence are correct (keep one H of NH3 terminal)
CRYST1   33.386   33.386   65.691  90.00  90.00 120.00 H 3 2
SCALE1      0.029953  0.017293  0.000000        0.00000
SCALE2      0.000000  0.034586  0.000000        0.00000
SCALE3      0.000000  0.000000  0.015223        0.00000
ATOM      1  N   GLU A   1      10.117  25.200   2.571  1.00 -1.00           N
ANISOU    1  N   GLU A   1     3000   2835   1127    208   -142   -228       N
ATOM      2  CA  GLU A   1       9.338  24.612   3.642  1.00 -1.00           C
ANISOU    2  CA  GLU A   1     2595   2691   1161    236   -238   -217       C
ATOM      3  C   GLU A   1      10.179  23.909   4.698  1.00 -1.00           C
ANISOU    3  C   GLU A   1     2462   2375    957    271   -124   -310       C
ATOM      4  O   GLU A   1       9.927  24.071   5.895  1.00 -1.00           O
ANISOU    4  O   GLU A   1     2439   2116    850    143   -268   -383       O
ATOM      5  CB  GLU A   1       8.308  23.601   3.105  1.00 -1.00           C
ANISOU    5  CB  GLU A   1     2772   2664   1227    243   -283   -303       C
ATOM      6  CG  GLU A   1       7.502  22.878   4.179  1.00 -1.00           C
ANISOU    6  CG  GLU A   1     2867   3180   1472     61   -157   -268       C
ATOM      7  CD  GLU A   1       6.505  23.776   4.881  1.00 -1.00           C
ANISOU    7  CD  GLU A   1     3055   3241   1610    122   -165   -363       C
ATOM      8  OE1 GLU A   1       6.017  24.742   4.259  1.00 -1.00           O
ANISOU    8  OE1 GLU A   1     3535   3430   1924    346   -221   -295       O
ATOM      9  OE2 GLU A   1       6.197  23.520   6.058  1.00 -1.00           O
ANISOU    9  OE2 GLU A   1     3239   3344   1625     40   -119   -371       O
ATOM     10  H1  GLU A   1       9.897  25.002   1.732  1.00 18.32           H
ATOM     13  HA  GLU A   1       8.880  25.353   4.068  1.00 16.97           H
ATOM     14  HB2 GLU A   1       7.679  24.074   2.537  1.00 17.54           H
ATOM     15  HB3 GLU A   1       8.778  22.927   2.590  1.00 17.54           H
ATOM     16  HG2 GLU A   1       8.111  22.529   4.848  1.00 19.79           H
ATOM     17  HG3 GLU A   1       7.010  22.151   3.767  1.00 19.79           H
TER
HETATM   18  C1  SIN A   0      11.143  26.020   2.776  1.00 -1.00           C
ANISOU   18  C1  SIN A   0     3119   2972   1277    124   -163   -240       C
HETATM   19  C2  SIN A   0      11.845  26.557   1.542  1.00 -1.00           C
ANISOU   19  C2  SIN A   0     3459   3368   1388    102    -39   -151       C
HETATM   20  C3  SIN A   0      12.465  27.926   1.781  1.00 -1.00           C
ANISOU   20  C3  SIN A   0     3640   3427   1555     45     -9   -143       C
HETATM   21  C4  SIN A   0      11.485  29.081   1.820  1.00 -1.00           C
ANISOU   21  C4  SIN A   0     3756   3638   1693    179      8   -113       C
HETATM   22  O1  SIN A   0      11.536  26.296   3.914  1.00 -1.00           O
ANISOU   22  O1  SIN A   0     3338   3021   1204     43   -125   -241       O
HETATM   23  O3  SIN A   0      11.966  30.230   1.941  1.00 -1.00           O
ANISOU   23  O3  SIN A   0     3949   3621   1603    176    -78    -76       O
HETATM   24  O4  SIN A   0      10.281  28.844   2.062  1.00 -1.00           O
ANISOU   24  O4  SIN A   0     3699   3930   2011    285      2   -324       O
HETATM   25  H21 SIN A   0      11.207  26.621   0.815  1.00 21.62           H
HETATM   26  H22 SIN A   0      12.541  25.936   1.276  1.00 21.62           H
HETATM   27  H31 SIN A   0      12.954  27.930   2.619  1.00 22.69           H
HETATM   28  H32 SIN A   0      13.119  28.122   1.092  1.00 22.69           H
END
"""


pdb_str_003 = """
REMARK scenario for linking, valence and a double conformation
CRYST1   39.670   48.940   71.610  88.74  97.15 108.44 P 1
SCALE1      0.025208  0.008405  0.003318        0.00000
SCALE2      0.000000  0.021539  0.000398        0.00000
SCALE3      0.000000  0.000000  0.014076        0.00000
ATOM      1  N  AASN A  74      12.741   6.714   3.033  0.88  5.87           N
ATOM      2  CA AASN A  74      13.110   7.769   2.091  0.88  6.12           C
ATOM      3  C  AASN A  74      11.884   8.130   1.219  0.88  5.86           C
ATOM      4  O  AASN A  74      10.858   7.490   1.271  0.88  6.44           O
ATOM      5  CB AASN A  74      14.293   7.285   1.218  0.88  7.27           C
ATOM      6  CG AASN A  74      15.066   8.432   0.603  0.88  8.82           C
ATOM      7  OD1AASN A  74      15.487   9.375   1.229  0.88 11.36           O
ATOM      8  ND2AASN A  74      15.252   8.317  -0.703  0.88  8.84           N
ATOM      9  H  AASN A  74      13.382   6.467   3.550  0.88  5.87           H
ATOM     10  HA AASN A  74      13.396   8.571   2.556  0.88  6.12           H
ATOM     11  HB2AASN A  74      14.904   6.771   1.769  0.88  7.27           H
ATOM     12  HB3AASN A  74      13.951   6.732   0.498  0.88  7.27           H
ATOM     13 HD21AASN A  74      14.860   7.684  -1.133  0.88  8.84           H
ATOM     15  N  BASN A  74      12.696   6.684   3.030  0.12  5.92           N
ATOM     16  CA BASN A  74      13.040   7.695   2.058  0.12  6.02           C
ATOM     17  C  BASN A  74      11.859   8.100   1.181  0.12  5.93           C
ATOM     18  O  BASN A  74      10.802   7.435   1.187  0.12  7.07           O
ATOM     19  CB BASN A  74      14.200   7.209   1.214  0.12  5.18           C
ATOM     20  CG BASN A  74      15.493   7.906   1.579  0.12  3.11           C
ATOM     21  OD1BASN A  74      15.680   8.391   2.668  0.12  3.08           O
ATOM     22  ND2BASN A  74      16.411   7.915   0.608  0.12  3.22           N
ATOM     23  H  BASN A  74      13.365   6.396   3.487  0.12  5.92           H
ATOM     24  HA BASN A  74      13.309   8.501   2.527  0.12  6.02           H
ATOM     25  HB2BASN A  74      14.012   7.388   0.279  0.12  5.18           H
ATOM     26  HB3BASN A  74      14.319   6.256   1.352  0.12  5.18           H
ATOM     27 HD21BASN A  74      16.235   7.488  -0.118  0.12  3.22           H
TER
HETATM   29  C1 ANAG C   1      16.113   9.249  -1.424  0.74  8.24           C
HETATM   30  C2 ANAG C   1      16.970   8.487  -2.451  0.74  7.99           C
HETATM   31  C3 ANAG C   1      17.836   9.477  -3.177  0.74  8.76           C
HETATM   32  C4 ANAG C   1      17.107  10.640  -3.663  0.74 10.13           C
HETATM   33  C5 ANAG C   1      16.092  11.231  -2.643  0.74  8.62           C
HETATM   34  C6 ANAG C   1      15.151  12.287  -3.168  0.74  9.76           C
HETATM   35  C7 ANAG C   1      17.696   6.180  -1.895  0.74  8.89           C
HETATM   36  C8 ANAG C   1      18.571   5.397  -0.988  0.74  9.89           C
HETATM   37  N2 ANAG C   1      17.767   7.513  -1.768  0.74  8.57           N
HETATM   38  O3 ANAG C   1      18.515   8.763  -4.246  0.74 10.65           O
HETATM   39  O4 ANAG C   1      18.042  11.683  -3.981  0.74 13.35           O
HETATM   40  O5 ANAG C   1      15.310  10.151  -2.140  0.74  7.25           O
HETATM   41  O6 ANAG C   1      14.360  11.822  -4.218  0.74  8.68           O
HETATM   42  O7 ANAG C   1      16.861   5.728  -2.714  0.74 10.29           O
HETATM   43  H1 ANAG C   1      16.696   9.715  -0.805  0.74  8.24           H
HETATM   44  H2 ANAG C   1      16.401   8.030  -3.090  0.74  7.99           H
HETATM   45  H3 ANAG C   1      18.456   9.829  -2.519  0.74  8.76           H
HETATM   46  H4 ANAG C   1      16.631  10.272  -4.424  0.74 10.13           H
HETATM   47  H5 ANAG C   1      16.602  11.677  -1.949  0.74  8.62           H
HETATM   48  H61ANAG C   1      15.688  13.047  -3.442  0.74  9.76           H
HETATM   49  H62ANAG C   1      14.606  12.589  -2.424  0.74  9.76           H
HETATM   50  H81ANAG C   1      19.036   6.002  -0.389  0.74  9.89           H
HETATM   51  H82ANAG C   1      19.216   4.900  -1.515  0.74  9.89           H
HETATM   52  H83ANAG C   1      18.028   4.781  -0.471  0.74  9.89           H
HETATM   53  HN2ANAG C   1      18.352   7.824  -1.219  0.74  8.57           H
HETATM   55  HO3ANAG C   1      18.959   8.129  -3.895  0.74 10.65           H
HETATM   57  HO6ANAG C   1      14.872  11.499  -4.815  0.74  8.68           H
HETATM   58  C1 BNAG C   1      17.701   8.609   0.704  0.26  6.07           C
HETATM   59  C2 BNAG C   1      18.523   8.170  -0.427  0.26  8.31           C
HETATM   60  C3 BNAG C   1      19.761   8.820  -0.183  0.26  6.55           C
HETATM   61  C4 BNAG C   1      19.463  10.425  -0.498  0.26 10.04           C
HETATM   62  C5 BNAG C   1      18.496  10.802   0.555  0.26  8.70           C
HETATM   63  C6 BNAG C   1      17.911  12.176   0.535  0.26 13.23           C
HETATM   64  C7 BNAG C   1      18.720   5.989  -1.367  0.26 10.97           C
HETATM   65  C8 BNAG C   1      19.249   4.641  -1.192  0.26 12.42           C
HETATM   66  N2 BNAG C   1      18.962   6.846  -0.336  0.26  9.44           N
HETATM   67  O3 BNAG C   1      20.699   8.721  -1.435  0.26 14.71           O
HETATM   68  O4 BNAG C   1      20.552  11.276  -0.389  0.26 12.04           O
HETATM   69  O5 BNAG C   1      17.320   9.980   0.458  0.26  7.16           O
HETATM   70  O6 BNAG C   1      17.693  12.418  -0.804  0.26 19.92           O
HETATM   71  O7 BNAG C   1      18.181   6.220  -2.430  0.26 11.69           O
HETATM   72  H1 BNAG C   1      18.137   8.480   1.561  0.26  6.07           H
HETATM   73  H2 BNAG C   1      18.012   8.327  -1.237  0.26  8.31           H
HETATM   74  H3 BNAG C   1      20.062   8.479   0.674  0.26  6.55           H
HETATM   75  H4 BNAG C   1      19.176  10.456  -1.424  0.26 10.04           H
HETATM   76  H5 BNAG C   1      19.021  10.704   1.364  0.26  8.70           H
HETATM   77  H61BNAG C   1      18.525  12.807   0.942  0.26 13.23           H
HETATM   78  H62BNAG C   1      17.099  12.201   1.065  0.26 13.23           H
HETATM   79  H81BNAG C   1      19.656   4.567  -0.315  0.26 12.42           H
HETATM   80  H82BNAG C   1      19.915   4.466  -1.875  0.26 12.42           H
HETATM   81  H83BNAG C   1      18.526   4.000  -1.271  0.26 12.42           H
HETATM   82  HN2BNAG C   1      19.380   6.582   0.368  0.26  9.44           H
HETATM   84  HO3BNAG C   1      20.316   9.121  -2.080  0.26 14.71           H
HETATM   86  HO6BNAG C   1      18.423  12.288  -1.220  0.26 19.92           H
HETATM   87  C1 ANAG C   2      17.995  12.199  -5.260  0.74 19.86           C
HETATM   88  C2 ANAG C   2      18.677  13.574  -5.253  0.74 23.29           C
HETATM   89  C3 ANAG C   2      18.830  14.151  -6.591  0.74 29.68           C
HETATM   90  C4 ANAG C   2      19.318  13.037  -7.482  0.74 30.27           C
HETATM   91  C5 ANAG C   2      18.591  11.736  -7.460  0.74 29.19           C
HETATM   92  C6 ANAG C   2      19.149  10.584  -8.227  0.74 33.79           C
HETATM   93  C7 ANAG C   2      18.471  14.805  -3.050  0.74 27.91           C
HETATM   94  C8 ANAG C   2      17.772  15.886  -2.243  0.74 37.24           C
HETATM   95  N2 ANAG C   2      18.165  14.566  -4.332  0.74 22.31           N
HETATM   96  O3 ANAG C   2      19.864  15.164  -6.531  0.74 40.55           O
HETATM   97  O4 ANAG C   2      19.327  13.494  -8.818  0.74 46.76           O
HETATM   98  O5 ANAG C   2      18.644  11.227  -6.086  0.74 21.90           O
HETATM   99  O6 ANAG C   2      20.567  10.602  -8.018  0.74 41.31           O
HETATM  100  O7 ANAG C   2      19.396  14.260  -2.622  0.74 24.10           O
HETATM  101  H1 ANAG C   2      17.076  12.314  -5.548  0.74 19.86           H
HETATM  102  H2 ANAG C   2      19.547  13.331  -4.899  0.74 23.29           H
HETATM  103  H3 ANAG C   2      17.981  14.517  -6.883  0.74 29.68           H
HETATM  104  H4 ANAG C   2      20.179  12.882  -7.064  0.74 30.27           H
HETATM  105  H5 ANAG C   2      17.721  11.946  -7.834  0.74 29.19           H
HETATM  106  H61ANAG C   2      18.738   9.762  -7.916  0.74 33.79           H
HETATM  107  H62ANAG C   2      18.910  10.675  -9.163  0.74 33.79           H
HETATM  108  H81ANAG C   2      17.079  16.294  -2.785  0.74 37.24           H
HETATM  109  H82ANAG C   2      18.419  16.560  -1.981  0.74 37.24           H
HETATM  110  H83ANAG C   2      17.375  15.490  -1.451  0.74 37.24           H
HETATM  111  HN2ANAG C   2      17.567  15.083  -4.671  0.74 22.31           H
HETATM  113  HO3ANAG C   2      20.590  14.783  -6.307  0.74 40.55           H
HETATM  114  HO4ANAG C   2      19.521  14.322  -8.826  0.74 46.76           H
HETATM  115  HO6ANAG C   2      20.857  11.373  -8.226  0.74 41.31           H
HETATM  116  C1 BNAG C   2      20.710  12.370  -1.239  0.26 14.76           C
HETATM  117  C2 BNAG C   2      21.943  13.194  -0.806  0.26 19.56           C
HETATM  118  C3 BNAG C   2      22.215  14.268  -1.785  0.26 24.38           C
HETATM  119  C4 BNAG C   2      22.458  13.627  -3.140  0.26 21.94           C
HETATM  120  C5 BNAG C   2      21.254  12.760  -3.461  0.26 20.16           C
HETATM  121  C6 BNAG C   2      21.355  12.027  -4.829  0.26 23.66           C
HETATM  122  C7 BNAG C   2      22.490  13.708   1.511  0.26 32.27           C
HETATM  123  C8 BNAG C   2      21.906  14.269   2.792  0.26 36.29           C
HETATM  124  N2 BNAG C   2      21.644  13.753   0.496  0.26 26.94           N
HETATM  125  O3 BNAG C   2      23.271  15.027  -1.295  0.26 27.97           O
HETATM  126  O4 BNAG C   2      22.468  14.619  -4.109  0.26 26.97           O
HETATM  127  O5 BNAG C   2      20.976  11.868  -2.415  0.26 11.40           O
HETATM  128  O6 BNAG C   2      22.466  11.207  -5.128  0.26 29.08           O
HETATM  129  O7 BNAG C   2      23.600  13.383   1.407  0.26 21.94           O
HETATM  130  H1 BNAG C   2      19.925  12.940  -1.237  0.26 14.76           H
HETATM  131  H2 BNAG C   2      22.730  12.628  -0.768  0.26 19.56           H
HETATM  132  H3 BNAG C   2      21.443  14.847  -1.886  0.26 24.38           H
HETATM  133  H4 BNAG C   2      23.296  13.144  -3.071  0.26 21.94           H
HETATM  134  H5 BNAG C   2      20.480  13.340  -3.539  0.26 20.16           H
HETATM  135  H61BNAG C   2      21.284  12.729  -5.495  0.26 23.66           H
HETATM  136  H62BNAG C   2      20.541  11.503  -4.893  0.26 23.66           H
HETATM  137  H81BNAG C   2      20.980  14.517   2.641  0.26 36.29           H
HETATM  138  H82BNAG C   2      22.413  15.051   3.060  0.26 36.29           H
HETATM  139  H83BNAG C   2      21.953  13.595   3.488  0.26 36.29           H
HETATM  140  HN2BNAG C   2      20.879  14.130   0.606  0.26 26.94           H
HETATM  142  HO3BNAG C   2      23.951  14.520  -1.236  0.26 27.97           H
HETATM  143  HO4BNAG C   2      22.798  15.328  -3.776  0.26 26.97           H
HETATM  144  HO6BNAG C   2      23.173  11.663  -5.010  0.26 29.08           H
END
"""

pdb_str_004 = """
REMARK Linking + valence scenario for a carbohydrate
CRYST1   39.670   48.940   71.610  88.74  97.15 108.44 P 1
SCALE1      0.025208  0.008405  0.003318        0.00000
SCALE2      0.000000  0.021539  0.000398        0.00000
SCALE3      0.000000  0.000000  0.014076        0.00000
ATOM      1  N   ASN B  74      -3.026  34.995  30.688  1.00  9.05           N
ATOM      2  CA  ASN B  74      -2.493  35.978  31.595  1.00  9.88           C
ATOM      3  C   ASN B  74      -1.547  35.327  32.616  1.00  9.26           C
ATOM      4  O   ASN B  74      -1.420  34.120  32.689  1.00  9.22           O
ATOM      5  CB  ASN B  74      -3.624  36.727  32.324  1.00 12.24           C
ATOM      6  CG  ASN B  74      -3.216  38.094  32.845  1.00 14.21           C
ATOM      7  OD1 ASN B  74      -2.472  38.834  32.203  1.00 18.21           O
ATOM      8  ND2 ASN B  74      -3.692  38.409  34.028  1.00 16.96           N
ATOM      9  H   ASN B  74      -2.549  34.876  29.982  1.00  9.05           H
ATOM     10  HA  ASN B  74      -1.988  36.631  31.085  1.00  9.88           H
ATOM     11  HB2 ASN B  74      -3.912  36.195  33.082  1.00 12.24           H
ATOM     12  HB3 ASN B  74      -4.362  36.854  31.708  1.00 12.24           H
ATOM     13 HD21 ASN B  74      -4.198  37.842  34.432  1.00 16.96           H
TER
HETATM   15  C1  NAG D   1      -3.389  39.731  34.716  1.00 15.06           C
HETATM   16  C2  NAG D   1      -4.468  40.105  35.528  1.00 16.43           C
HETATM   17  C3  NAG D   1      -4.201  41.525  36.083  1.00 15.29           C
HETATM   18  C4  NAG D   1      -2.754  41.564  36.747  1.00 15.66           C
HETATM   19  C5  NAG D   1      -1.708  40.987  35.892  1.00 17.17           C
HETATM   20  C6  NAG D   1      -0.490  40.596  36.940  1.00 25.12           C
HETATM   21  C7  NAG D   1      -6.706  39.271  34.932  1.00 17.46           C
HETATM   22  C8  NAG D   1      -7.943  39.400  34.022  1.00 20.06           C
HETATM   23  N2  NAG D   1      -5.748  40.154  34.750  1.00 16.04           N
HETATM   24  O3  NAG D   1      -5.243  41.846  36.992  1.00 17.10           O
HETATM   25  O4  NAG D   1      -2.597  43.022  36.916  1.00 18.40           O
HETATM   26  O5  NAG D   1      -2.177  39.700  35.427  1.00 15.81           O
HETATM   27  O6  NAG D   1       0.621  40.445  35.919  1.00 26.15           O
HETATM   28  O7  NAG D   1      -6.670  38.382  35.716  1.00 17.36           O
HETATM   29  H1  NAG D   1      -3.303  40.393  34.012  1.00 15.06           H
HETATM   30  H2  NAG D   1      -4.545  39.445  36.235  1.00 16.43           H
HETATM   31  H3  NAG D   1      -4.201  42.160  35.350  1.00 15.29           H
HETATM   32  H4  NAG D   1      -2.726  41.025  37.553  1.00 15.66           H
HETATM   33  H5  NAG D   1      -1.456  41.561  35.151  1.00 17.17           H
HETATM   34  H61 NAG D   1      -0.689  39.788  37.439  1.00 25.12           H
HETATM   35  H62 NAG D   1      -0.340  41.295  37.596  1.00 25.12           H
HETATM   36  H81 NAG D   1      -7.800  40.119  33.386  1.00 20.06           H
HETATM   37  H82 NAG D   1      -8.722  39.596  34.566  1.00 20.06           H
HETATM   38  H83 NAG D   1      -8.079  38.566  33.546  1.00 20.06           H
HETATM   39  HN2 NAG D   1      -5.852  40.780  34.170  1.00 16.04           H
HETATM   41  HO3 NAG D   1      -5.213  41.286  37.631  1.00 17.10           H
HETATM   43  HO6 NAG D   1       0.364  39.894  35.325  1.00 26.15           H
HETATM   44  C1  NAG D   2      -2.391  43.486  38.203  1.00 21.11           C
HETATM   45  C2  NAG D   2      -1.775  44.823  37.997  1.00 26.59           C
HETATM   46  C3  NAG D   2      -1.585  45.513  39.294  1.00 26.26           C
HETATM   47  C4  NAG D   2      -2.959  45.457  39.969  1.00 29.92           C
HETATM   48  C5  NAG D   2      -3.462  44.035  40.104  1.00 30.15           C
HETATM   49  C6  NAG D   2      -4.756  43.848  40.745  1.00 35.86           C
HETATM   50  C7  NAG D   2      -0.234  44.953  35.978  1.00 30.77           C
HETATM   51  C8  NAG D   2       1.206  44.728  35.397  1.00 33.11           C
HETATM   52  N2  NAG D   2      -0.459  44.736  37.311  1.00 26.48           N
HETATM   53  O3  NAG D   2      -1.072  46.826  39.192  1.00 28.25           O
HETATM   54  O4  NAG D   2      -2.849  45.859  41.286  1.00 32.37           O
HETATM   55  O5  NAG D   2      -3.648  43.578  38.774  1.00 23.11           O
HETATM   56  O6  NAG D   2      -5.481  44.858  40.253  1.00 30.27           O
HETATM   57  O7  NAG D   2      -1.176  45.443  35.293  1.00 29.59           O
HETATM   58  H1  NAG D   2      -1.798  42.910  38.711  1.00 21.11           H
HETATM   59  H2  NAG D   2      -2.392  45.317  37.434  1.00 26.59           H
HETATM   60  H3  NAG D   2      -0.902  45.033  39.787  1.00 26.26           H
HETATM   61  H4  NAG D   2      -3.496  46.020  39.389  1.00 29.92           H
HETATM   62  H5  NAG D   2      -2.821  43.547  40.644  1.00 30.15           H
HETATM   63  H61 NAG D   2      -4.662  43.876  41.710  1.00 35.86           H
HETATM   64  H62 NAG D   2      -5.117  42.974  40.529  1.00 35.86           H
HETATM   65  H81 NAG D   2       1.786  44.392  36.098  1.00 33.11           H
HETATM   66  H82 NAG D   2       1.552  45.570  35.063  1.00 33.11           H
HETATM   67  H83 NAG D   2       1.163  44.083  34.673  1.00 33.11           H
HETATM   68  HN2 NAG D   2       0.216  44.531  37.802  1.00 26.48           H
HETATM   70  HO3 NAG D   2      -1.634  47.292  38.756  1.00 28.25           H
HETATM   71  HO4 NAG D   2      -2.248  46.457  41.343  1.00 32.37           H
HETATM   72  HO6 NAG D   2      -5.060  45.582  40.396  1.00 30.27           H
END
"""

pdb_str_005 = """
CRYST1   80.435  118.034  112.075  90.00  93.12  90.00 P 1 21 1
SCALE1      0.012432  0.000000  0.000678        0.00000
SCALE2      0.000000  0.008472  0.000000        0.00000
SCALE3      0.000000  0.000000  0.008936        0.00000
ATOM      1  N   ASN A  90      76.017 -20.785 155.731  1.00 -1.00           N
ANISOU    1  N   ASN A  90    17156  21970  15197  -6703    124    870       N
ATOM      2  CA  ASN A  90      75.942 -19.415 156.227  1.00 -1.00           C
ANISOU    2  CA  ASN A  90    17577  22797  16129  -6189     29   1252       C
ATOM      3  C   ASN A  90      77.350 -18.831 156.235  1.00 -1.00           C
ANISOU    3  C   ASN A  90    17170  21823  15645  -5768    103   1134       C
ATOM      4  O   ASN A  90      78.239 -19.359 156.912  1.00 -1.00           O
ANISOU    4  O   ASN A  90    17802  21588  16095  -5569    283    881       O
ATOM      5  CB  ASN A  90      75.320 -19.390 157.626  1.00 -1.00           C
ANISOU    5  CB  ASN A  90    18758  23893  17672  -5918    105   1419       C
ATOM      6  CG  ASN A  90      75.068 -17.982 158.153  1.00 -1.00           C
ANISOU    6  CG  ASN A  90    19640  25239  19120  -5436     46   1805       C
ATOM      7  OD1 ASN A  90      75.922 -17.100 158.061  1.00 -1.00           O
ANISOU    7  OD1 ASN A  90    18832  24282  18399  -5077     50   1839       O
ATOM      8  ND2 ASN A  90      73.886 -17.776 158.730  1.00 -1.00           N
ANISOU    8  ND2 ASN A  90    21686  27829  21572  -5430     16   2088       N
ATOM      9  H   ASN A  90      76.270 -20.847 154.911  1.00142.97           H
ATOM     10  HA  ASN A  90      75.375 -18.865 155.664  1.00148.71           H
ATOM     11  HB2 ASN A  90      75.921 -19.835 158.244  1.00158.76           H
ATOM     12  HB3 ASN A  90      74.468 -19.854 157.599  1.00158.76           H
ATOM     13 HD21 ASN A  90      73.319 -18.422 158.759  1.00187.09           H
TER
HETATM   15  C1  NAG C   1      73.523 -16.504 159.309  1.00103.52           C
HETATM   16  C2  NAG C   1      72.427 -15.749 158.564  1.00114.27           C
HETATM   17  C3  NAG C   1      72.112 -14.434 159.280  1.00117.99           C
HETATM   18  C4  NAG C   1      71.829 -14.658 160.763  1.00124.42           C
HETATM   19  C5  NAG C   1      72.930 -15.509 161.403  1.00121.57           C
HETATM   20  C6  NAG C   1      72.619 -15.937 162.820  1.00122.81           C
HETATM   21  C7  NAG C   1      71.992 -15.651 156.147  1.00125.03           C
HETATM   22  C8  NAG C   1      72.570 -15.344 154.798  1.00127.15           C
HETATM   23  N2  NAG C   1      72.819 -15.497 157.186  1.00121.51           N
HETATM   24  O3  NAG C   1      70.990 -13.814 158.662  1.00116.89           O
HETATM   25  O4  NAG C   1      71.791 -13.386 161.404  1.00133.35           O
HETATM   26  O5  NAG C   1      73.123 -16.717 160.654  1.00114.93           O
HETATM   27  O6  NAG C   1      73.745 -15.787 163.673  1.00121.65           O
HETATM   28  O7  NAG C   1      70.831 -16.020 156.289  1.00125.47           O
HETATM   29  H1  NAG C   1      74.310 -15.937 159.309  1.00103.52           H
HETATM   30  H2  NAG C   1      71.632 -16.305 158.561  1.00114.27           H
HETATM   31  H3  NAG C   1      72.914 -13.895 159.195  1.00117.99           H
HETATM   32  H4  NAG C   1      70.979 -15.124 160.802  1.00124.42           H
HETATM   33  H5  NAG C   1      73.736 -14.970 161.406  1.00121.57           H
HETATM   34  H61 NAG C   1      71.867 -15.409 163.130  1.00122.81           H
HETATM   35  H62 NAG C   1      72.316 -16.858 162.792  1.00122.81           H
HETATM   36  H81 NAG C   1      73.509 -15.120 154.894  1.00127.15           H
HETATM   37  H82 NAG C   1      72.095 -14.593 154.409  1.00127.15           H
HETATM   38  H83 NAG C   1      72.477 -16.121 154.225  1.00127.15           H
HETATM   39  HN2 NAG C   1      73.625 -15.236 157.040  1.00121.51           H
HETATM   41  HO3 NAG C   1      70.851 -13.074 159.056  1.00116.89           H
HETATM   43  HO6 NAG C   1      74.400 -16.202 163.325  1.00121.65           H
HETATM   44  C1  NAG C   2      70.691 -13.209 162.330  1.00143.22           C
HETATM   45  C2  NAG C   2      71.046 -12.013 163.217  1.00149.90           C
HETATM   46  C3  NAG C   2      69.929 -11.751 164.223  1.00153.90           C
HETATM   47  C4  NAG C   2      68.594 -11.600 163.506  1.00157.09           C
HETATM   48  C5  NAG C   2      68.339 -12.806 162.605  1.00152.27           C
HETATM   49  C6  NAG C   2      67.082 -12.672 161.777  1.00151.75           C
HETATM   50  C7  NAG C   2      73.401 -11.485 163.670  1.00150.52           C
HETATM   51  C8  NAG C   2      74.621 -11.843 164.465  1.00149.56           C
HETATM   52  N2  NAG C   2      72.312 -12.226 163.901  1.00151.39           N
HETATM   53  O3  NAG C   2      70.226 -10.571 164.961  1.00153.21           O
HETATM   54  O4  NAG C   2      67.540 -11.490 164.456  1.00163.16           O
HETATM   55  O5  NAG C   2      69.429 -12.972 161.684  1.00147.83           O
HETATM   56  O6  NAG C   2      67.313 -11.906 160.602  1.00151.36           O
HETATM   57  O7  NAG C   2      73.402 -10.566 162.856  1.00149.17           O
HETATM   58  H1  NAG C   2      70.620 -13.994 162.895  1.00143.22           H
HETATM   59  H2  NAG C   2      71.140 -11.238 162.641  1.00149.90           H
HETATM   60  H3  NAG C   2      69.893 -12.528 164.802  1.00153.90           H
HETATM   61  H4  NAG C   2      68.684 -10.788 162.984  1.00157.09           H
HETATM   62  H5  NAG C   2      68.256 -13.586 163.175  1.00152.27           H
HETATM   63  H61 NAG C   2      66.399 -12.268 162.336  1.00151.75           H
HETATM   64  H62 NAG C   2      66.770 -13.565 161.563  1.00151.75           H
HETATM   65  H81 NAG C   2      74.432 -12.623 165.009  1.00149.56           H
HETATM   66  H82 NAG C   2      74.862 -11.097 165.037  1.00149.56           H
HETATM   67  H83 NAG C   2      75.353 -12.038 163.859  1.00149.56           H
HETATM   68  HN2 NAG C   2      72.354 -12.860 164.480  1.00151.39           H
HETATM   70  HO3 NAG C   2      70.256  -9.922 164.413  1.00153.21           H
HETATM   72  HO6 NAG C   2      67.657 -11.163 160.832  1.00151.36           H
HETATM   73  C1  BMA C   3      67.036 -10.135 164.465  1.00164.84           C
HETATM   74  C2  BMA C   3      65.513 -10.187 164.205  1.00164.51           C
HETATM   75  C3  BMA C   3      64.912  -8.784 164.326  1.00164.84           C
HETATM   76  C4  BMA C   3      65.353  -8.085 165.626  1.00165.75           C
HETATM   77  C5  BMA C   3      66.886  -8.117 165.760  1.00164.34           C
HETATM   78  C6  BMA C   3      67.377  -7.517 167.067  1.00159.95           C
HETATM   79  O2  BMA C   3      64.864 -10.998 165.174  1.00162.15           O
HETATM   80  O3  BMA C   3      63.491  -8.819 164.251  1.00162.71           O
HETATM   81  O4  BMA C   3      64.902  -6.737 165.632  1.00164.98           O
HETATM   82  O5  BMA C   3      67.321  -9.489 165.704  1.00164.62           O
HETATM   83  O6  BMA C   3      66.718  -6.271 167.265  1.00154.88           O
HETATM   84  H1  BMA C   3      67.445  -9.638 163.740  1.00164.84           H
HETATM   85  H2  BMA C   3      65.448 -10.536 163.302  1.00164.51           H
HETATM   86  H3  BMA C   3      65.260  -8.317 163.550  1.00164.84           H
HETATM   87  H4  BMA C   3      64.914  -8.536 166.365  1.00165.75           H
HETATM   88  H5  BMA C   3      67.297  -7.591 165.056  1.00164.34           H
HETATM   89  H61 BMA C   3      67.191  -8.125 167.800  1.00159.95           H
HETATM   90  H62 BMA C   3      68.339  -7.396 167.035  1.00159.95           H
HETATM   92  HO2 BMA C   3      65.200 -10.810 165.932  1.00162.15           H
HETATM   93  HO3 BMA C   3      63.182  -9.145 164.973  1.00162.71           H
HETATM   94  HO4 BMA C   3      65.175  -6.362 166.344  1.00164.98           H
HETATM   95  HO6 BMA C   3      65.881  -6.399 167.190  1.00154.88           H
END
"""

pdb_str_006 = '''
CRYST1   56.620   56.620  181.790  90.00  90.00  90.00 P 41 21 2
SCALE1      0.017662  0.000000  0.000000        0.00000
SCALE2      0.000000  0.017662  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005501        0.00000
HETATM    1  CB1ACNZ A 301      21.699  31.831  83.042  0.42 31.90           C
HETATM    2  CB2ACNZ A 301      20.840  35.774  78.508  0.42 36.67           C
HETATM    3  CG1ACNZ A 301      20.748  32.695  82.223  0.42 32.30           C
HETATM    4  SG2ACNZ A 301      21.831  34.484  77.712  0.42 38.01           S
HETATM    5  CD1ACNZ A 301      21.433  33.671  81.296  0.42 32.98           C
HETATM    6  OE1ACNZ A 301      22.657  33.747  81.222  0.42 32.67           O
HETATM    7  C1 ACNZ A 301      20.264  29.771  83.305  0.42 30.73           C
HETATM    8  C13ACNZ A 301      23.480  35.240  77.602  0.42 38.76           C
HETATM    9  C14ACNZ A 301      24.306  34.544  76.542  0.42 39.44           C
HETATM   10  C2 ACNZ A 301      20.164  36.824  80.647  0.42 35.25           C
HETATM   11  C3 ACNZ A 301      20.811  40.295  81.232  0.42 36.72           C
HETATM   12  CA1ACNZ A 301      20.988  30.898  84.024  0.42 31.32           C
HETATM   13  CA2ACNZ A 301      21.019  35.721  80.028  0.42 35.21           C
HETATM   14  CA3ACNZ A 301      20.125  39.024  81.660  0.42 36.22           C
HETATM   15  CB4ACNZ A 301      30.740  35.079  78.050  0.42 38.90           C
HETATM   16  CB5ACNZ A 301      27.945  33.837  75.863  0.42 38.95           C
HETATM   17  CD4ACNZ A 301      32.849  34.182  77.341  0.42 38.95           C
HETATM   18  CD5ACNZ A 301      25.785  34.478  76.724  0.42 39.23           C
HETATM   19  CE4ACNZ A 301      32.199  33.407  76.406  0.42 38.70           C
HETATM   20  CE5ACNZ A 301      26.419  34.990  77.854  0.42 39.19           C
HETATM   21  CG4ACNZ A 301      32.120  35.018  78.161  0.42 39.01           C
HETATM   22  CG5ACNZ A 301      26.569  33.894  75.734  0.42 39.12           C
HETATM   23  CH4ACNZ A 301      30.067  34.300  77.110  0.42 38.71           C
HETATM   24  CH5ACNZ A 301      28.589  34.355  76.985  0.42 38.89           C
HETATM   25  CZ4ACNZ A 301      30.821  33.462  76.288  0.42 38.64           C
HETATM   26  CZ5ACNZ A 301      27.797  34.927  77.981  0.42 39.11           C
HETATM   27  N1 ACNZ A 301      21.971  30.290  84.963  0.42 31.34           N
HETATM   28  N2 ACNZ A 301      20.617  34.430  80.560  0.42 33.98           N
HETATM   29  N3 ACNZ A 301      20.820  37.838  81.211  0.42 35.82           N
HETATM   30  O11ACNZ A 301      18.995  30.026  83.078  0.42 30.27           O
HETATM   31  O12ACNZ A 301      20.816  28.737  83.002  0.42 30.57           O
HETATM   32  O2 ACNZ A 301      18.939  36.766  80.595  0.42 34.92           O
HETATM   33  O31ACNZ A 301      21.252  40.487  80.111  0.42 36.85           O
HETATM   34  O32ACNZ A 301      20.885  41.163  82.211  0.42 36.84           O
HETATM   35  O5 ACNZ A 301      23.775  34.065  75.554  0.42 39.70           O
HETATM   36  H1 ACNZ A 301      22.793  30.250  84.624  0.42 31.34           H
HETATM   37  H2 ACNZ A 301      19.822  34.147  80.392  0.42 33.98           H
HETATM   38  H3 ACNZ A 301      21.672  37.827  81.328  0.42 35.82           H
HETATM   39  HA1ACNZ A 301      20.347  31.361  84.585  0.42 31.32           H
HETATM   40  HA2ACNZ A 301      21.958  35.873  80.218  0.42 35.21           H
HETATM   41 HA31ACNZ A 301      20.051  39.024  82.627  0.42 36.22           H
HETATM   42 HA32ACNZ A 301      19.218  39.033  81.316  0.42 36.22           H
HETATM   43 HB11ACNZ A 301      22.306  32.400  83.541  0.42 31.90           H
HETATM   44 HB12ACNZ A 301      22.245  31.294  82.446  0.42 31.90           H
HETATM   45 HB21ACNZ A 301      21.113  36.646  78.181  0.42 36.67           H
HETATM   46 HB22ACNZ A 301      19.903  35.650  78.289  0.42 36.67           H
HETATM   47 HG11ACNZ A 301      20.164  33.218  82.794  0.42 32.30           H
HETATM   48 HG12ACNZ A 301      20.163  32.153  81.671  0.42 32.30           H
HETATM   49  HB4ACNZ A 301      30.261  35.647  78.610  0.42 38.90           H
HETATM   50  HB5ACNZ A 301      28.443  33.444  75.183  0.42 38.95           H
HETATM   51  HD4ACNZ A 301      33.775  34.142  77.419  0.42 38.95           H
HETATM   52  HE4ACNZ A 301      32.685  32.843  75.849  0.42 38.70           H
HETATM   53  HE5ACNZ A 301      25.931  35.382  78.542  0.42 39.19           H
HETATM   54  HG4ACNZ A 301      32.553  35.545  78.793  0.42 39.01           H
HETATM   55  HG5ACNZ A 301      26.167  33.537  74.975  0.42 39.12           H
HETATM   56  HZ4ACNZ A 301      30.397  32.932  75.652  0.42 38.64           H
HETATM   57  HZ5ACNZ A 301      28.195  35.274  78.747  0.42 39.11           H
HETATM   58 H11NACNZ A 301      22.044  30.742  85.726  0.42 31.34           H
HETATM   59 H12NACNZ A 301      21.760  29.456  85.192  0.42 31.34           H
HETATM   60 H131ACNZ A 301      23.934  35.181  78.457  0.42 38.76           H
HETATM   61 H132ACNZ A 301      23.403  36.183  77.386  0.42 38.76           H
HETATM   62  CB1BCNZ A 301      21.628  31.747  83.072  0.58 24.44           C
HETATM   63  CB2BCNZ A 301      20.928  35.615  78.523  0.58 32.96           C
HETATM   64  CG1BCNZ A 301      20.751  32.485  82.072  0.58 25.13           C
HETATM   65  SG2BCNZ A 301      21.870  34.251  77.786  0.58 35.21           S
HETATM   66  CD1BCNZ A 301      21.399  33.719  81.489  0.58 26.54           C
HETATM   67  OE1BCNZ A 301      22.558  34.030  81.757  0.58 26.24           O
HETATM   68  C1 BCNZ A 301      20.084  29.772  83.397  0.58 21.85           C
HETATM   69  C13BCNZ A 301      23.611  34.725  78.008  0.58 36.08           C
HETATM   70  C14BCNZ A 301      24.428  34.553  76.745  0.58 36.66           C
HETATM   71  C2 BCNZ A 301      20.268  36.859  80.558  0.58 31.00           C
HETATM   72  C3 BCNZ A 301      20.953  40.442  81.245  0.58 35.20           C
HETATM   73  CA1BCNZ A 301      20.838  30.909  84.075  0.58 23.32           C
HETATM   74  CA2BCNZ A 301      21.083  35.672  80.046  0.58 30.56           C
HETATM   75  CA3BCNZ A 301      20.362  39.080  81.531  0.58 33.89           C
HETATM   76  CB4BCNZ A 301      29.376  30.173  77.597  0.58 27.67           C
HETATM   77  CB5BCNZ A 301      26.920  32.007  75.501  0.58 33.10           C
HETATM   78  CD4BCNZ A 301      30.221  28.218  76.489  0.58 26.53           C
HETATM   79  CD5BCNZ A 301      25.519  33.534  76.736  0.58 35.16           C
HETATM   80  CE4BCNZ A 301      29.225  28.329  75.545  0.58 26.68           C
HETATM   81  CE5BCNZ A 301      26.041  33.009  77.917  0.58 34.37           C
HETATM   82  CG4BCNZ A 301      30.291  29.135  77.520  0.58 26.89           C
HETATM   83  CG5BCNZ A 301      25.993  33.034  75.528  0.58 34.22           C
HETATM   84  CH4BCNZ A 301      28.371  30.309  76.639  0.58 28.91           C
HETATM   85  CH5BCNZ A 301      27.424  31.454  76.677  0.58 31.80           C
HETATM   86  CZ4BCNZ A 301      28.306  29.363  75.614  0.58 27.81           C
HETATM   87  CZ5BCNZ A 301      26.976  31.988  77.885  0.58 33.11           C
HETATM   88  N1 BCNZ A 301      21.748  30.353  85.116  0.58 23.71           N
HETATM   89  N2 BCNZ A 301      20.625  34.442  80.675  0.58 28.24           N
HETATM   90  N3 BCNZ A 301      20.934  38.001  80.752  0.58 32.37           N
HETATM   91  O11BCNZ A 301      20.774  28.662  83.286  0.58 21.59           O
HETATM   92  O12BCNZ A 301      18.936  29.904  82.989  0.58 20.54           O
HETATM   93  O2 BCNZ A 301      19.059  36.753  80.739  0.58 30.28           O
HETATM   94  O31BCNZ A 301      20.309  41.472  81.325  0.58 35.67           O
HETATM   95  O32BCNZ A 301      22.229  40.420  80.923  0.58 35.87           O
HETATM   96  O5 BCNZ A 301      24.172  35.209  75.746  0.58 37.96           O
HETATM   97  H1 BCNZ A 301      22.583  30.252  84.826  0.58 23.71           H
HETATM   98  H2 BCNZ A 301      19.824  34.173  80.514  0.58 28.24           H
HETATM   99  H3 BCNZ A 301      21.718  38.138  80.426  0.58 32.37           H
HETATM  100  HA1BCNZ A 301      20.186  31.424  84.576  0.58 23.32           H
HETATM  101  HA2BCNZ A 301      22.025  35.797  80.242  0.58 30.56           H
HETATM  102 HA31BCNZ A 301      20.474  38.898  82.477  0.58 33.89           H
HETATM  103 HA32BCNZ A 301      19.406  39.129  81.373  0.58 33.89           H
HETATM  104 HB11BCNZ A 301      22.242  31.164  82.598  0.58 24.44           H
HETATM  105 HB12BCNZ A 301      22.170  32.385  83.562  0.58 24.44           H
HETATM  106 HB21BCNZ A 301      21.231  36.455  78.144  0.58 32.96           H
HETATM  107 HB22BCNZ A 301      19.988  35.509  78.306  0.58 32.96           H
HETATM  108 HG11BCNZ A 301      20.506  31.917  81.325  0.58 25.13           H
HETATM  109 HG12BCNZ A 301      19.915  32.770  82.473  0.58 25.13           H
HETATM  110  HB4BCNZ A 301      29.436  30.783  78.296  0.58 27.67           H
HETATM  111  HB5BCNZ A 301      27.207  31.685  74.677  0.58 33.10           H
HETATM  112  HD4BCNZ A 301      30.844  27.529  76.433  0.58 26.53           H
HETATM  113  HE4BCNZ A 301      29.166  27.708  74.855  0.58 26.68           H
HETATM  114  HE5BCNZ A 301      25.772  33.333  78.746  0.58 34.37           H
HETATM  115  HG4BCNZ A 301      30.954  29.059  78.168  0.58 26.89           H
HETATM  116  HG5BCNZ A 301      25.686  33.391  74.726  0.58 34.22           H
HETATM  117  HZ4BCNZ A 301      27.640  29.424  74.968  0.58 27.81           H
HETATM  118  HZ5BCNZ A 301      27.307  31.657  78.688  0.58 33.11           H
HETATM  119 H11NBCNZ A 301      21.801  30.869  85.839  0.58 23.71           H
HETATM  120 H12NBCNZ A 301      21.491  29.552  85.406  0.58 23.71           H
HETATM  121 H131BCNZ A 301      24.007  34.190  78.713  0.58 36.08           H
HETATM  122 H132BCNZ A 301      23.663  35.651  78.292  0.58 36.08           H
END'''

pdb_str_007 = '''
REMARK scenario for 3 alt confs with split at Calpha (1ssx)
CRYST1   14.055   14.336   16.437  90.00  90.00  90.00 P 1
ATOM      1  N   MET A 190       8.147   8.228   5.841  1.00  5.11           N
ATOM      2  C   MET A 190       6.079   7.001   6.145  1.00  4.90           C
ATOM      3  O   MET A 190       5.489   8.085   6.030  1.00  5.59           O
ATOM      4  CA AMET A 190       7.607   6.937   6.231  0.64  5.00           C
ATOM      5  CB AMET A 190       8.086   6.534   7.644  0.64  5.98           C
ATOM      6  CG AMET A 190       7.600   7.386   8.781  0.64  8.55           C
ATOM      7  SD AMET A 190       5.950   7.100   9.364  0.64  9.70           S
ATOM      8  CE AMET A 190       6.057   5.778  10.487  0.64 19.81           C
ATOM      9  H  AMET A 190       8.372   8.277   5.013  1.00  5.11           H
ATOM     10  HA AMET A 190       7.935   6.240   5.642  0.64  5.00           H
ATOM     11  HB2AMET A 190       7.785   5.628   7.818  0.64  5.98           H
ATOM     12  HB3AMET A 190       9.055   6.569   7.655  0.64  5.98           H
ATOM     13  HG2AMET A 190       8.193   7.243   9.535  0.64  8.55           H
ATOM     14  HG3AMET A 190       7.641   8.313   8.498  0.64  8.55           H
ATOM     15  HE1AMET A 190       6.408   5.000  10.026  0.64 19.81           H
ATOM     16  HE2AMET A 190       5.172   5.586  10.834  0.64 19.81           H
ATOM     17  HE3AMET A 190       6.649   6.028  11.213  0.64 19.81           H
ATOM     18  CA BMET A 190       7.599   6.926   6.172  0.21  4.35           C
ATOM     19  CB BMET A 190       8.059   6.429   7.549  0.21  3.88           C
ATOM     20  CG BMET A 190       7.945   7.470   8.646  0.21  3.81           C
ATOM     21  SD BMET A 190       6.333   8.029   9.240  0.21  5.37           S
ATOM     22  CE BMET A 190       5.873   6.543  10.201  0.21  7.00           C
ATOM     23  H  BMET A 190       8.283   8.347   5.000  1.00  5.11           H
ATOM     24  HA BMET A 190       7.929   6.282   5.526  0.21  4.35           H
ATOM     25  HB2BMET A 190       7.513   5.669   7.804  0.21  3.88           H
ATOM     26  HB3BMET A 190       8.990   6.163   7.489  0.21  3.88           H
ATOM     27  HG2BMET A 190       8.403   7.114   9.423  0.21  3.81           H
ATOM     28  HG3BMET A 190       8.400   8.266   8.329  0.21  3.81           H
ATOM     29  HE1BMET A 190       5.845   5.778   9.605  0.21  7.00           H
ATOM     30  HE2BMET A 190       5.000   6.683  10.601  0.21  7.00           H
ATOM     31  HE3BMET A 190       6.535   6.398  10.895  0.21  7.00           H
ATOM     32  CA CMET A 190       7.592   6.970   6.290  0.15  4.69           C
ATOM     33  CB CMET A 190       7.939   6.723   7.773  0.15  5.99           C
ATOM     34  CG CMET A 190       7.637   7.940   8.598  0.15  7.61           C
ATOM     35  SD CMET A 190       7.421   7.776  10.371  0.15  6.55           S
ATOM     36  CE CMET A 190       5.808   8.422  10.506  0.15  6.32           C
ATOM     37  H  CMET A 190       8.470   8.203   5.044  1.00  5.11           H
ATOM     38  HA CMET A 190       7.972   6.249   5.764  0.15  4.69           H
ATOM     39  HB2CMET A 190       8.884   6.521   7.855  0.15  5.99           H
ATOM     40  HB3CMET A 190       7.411   5.982   8.110  0.15  5.99           H
ATOM     41  HG2CMET A 190       8.368   8.564   8.471  0.15  7.61           H
ATOM     42  HG3CMET A 190       6.812   8.322   8.259  0.15  7.61           H
ATOM     43  HE1CMET A 190       5.536   8.401  11.437  0.15  6.32           H
ATOM     44  HE2CMET A 190       5.205   7.880   9.973  0.15  6.32           H
ATOM     45  HE3CMET A 190       5.806   9.336  10.181  0.15  6.32           H
END
'''

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
