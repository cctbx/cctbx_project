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
  #test_009()

# ------------------------------------------------------------------------------

def compare_models(pdb_str,
                   contains     = None,
                   not_contains = None):
  #
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  # initial model
  model_initial = mmtbx.model.manager(model_input = pdb_inp,
                                      log         = null_out())
  hd_sel_initial = model_initial.get_hd_selection()
  number_h_expected = hd_sel_initial.count(True)
  ph_initial = model_initial.get_hierarchy()
  h_atoms_initial = ph_initial.select(hd_sel_initial).atoms()
  h_names_initial = list(h_atoms_initial.extract_name())
  # remove H atoms
  model_without_h = model_initial.select(~hd_sel_initial)
  hd_sel_without_h = model_without_h.get_hd_selection()
  assert (hd_sel_without_h is not None)
  assert (hd_sel_without_h.count(True) == 0)
  # place H atoms again
  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(model = model_without_h)
  reduce_add_h_obj.run()
  #
  model_h_added = reduce_add_h_obj.get_model()
  hd_sel_h_added = model_h_added.get_hd_selection()
  ph_h_added = model_h_added.get_hierarchy()
  h_atoms_added = ph_h_added.select(hd_sel_h_added).atoms()
  h_names_added = list(h_atoms_added.extract_name())
  number_h_added = hd_sel_h_added.count(True)

  #f = open("bla_intermediate.pdb","w")
  #f.write(model_h_added.model_as_pdb())
  #f.close()

  assert ph_initial.is_similar_hierarchy(other=ph_h_added)

  assert(number_h_expected == number_h_added)

  if not_contains:
    assert (not_contains not in h_names_added)
  if contains:
    assert (contains in h_names_added)

  sc_h_initial = model_initial.select(hd_sel_initial).get_sites_cart()
  sc_h_added   = model_h_added.select(hd_sel_h_added).get_sites_cart()

  d1 = {h_names_initial[i]: sc_h_initial[i] for i in range(len(h_names_initial))}
  d2 = {h_names_added[i]: sc_h_added[i] for i in range(len(h_names_added))}

  for name, sc in d2.items():
    assert(name in d1)
    assert approx_equal(sc, d1[name], 0.01), name

# ------------------------------------------------------------------------------

def test_000():
  '''
    CDL is true by default. Will crash for this example if is False.
  '''
  compare_models(pdb_str = pdb_str_000)

# ------------------------------------------------------------------------------

def test_001():
  '''
  clash between LYS 246 NZ & Tam 2 C5 --> fails in REDUCE

  Removing either Lys 246 NZ, or TAM 2 C5, or both of Lys 246 HZ1 and HZ3 will
  allow Reduce to run to completion and produce a usable result.
  However, Lys 246 HZ1 and HZ3 will not be added back by Reduce.
  '''
  compare_models(pdb_str = pdb_str_001)

# ------------------------------------------------------------------------------

def test_002():
  '''
    Test if multi-model file is supported
  '''
  compare_models(pdb_str = pdb_str_002)

# ------------------------------------------------------------------------------

def test_003():
  '''
    Check if normal and modified nucleic acid work.
  '''
  compare_models(pdb_str = pdb_str_003)

# ------------------------------------------------------------------------------

def test_004():
  '''
    Check if dc is processed correctly: TYR dc with every atom in either A or B
  '''
  compare_models(pdb_str = pdb_str_004)

# ------------------------------------------------------------------------------

def test_005():
  '''
    Check if dc is processed correctly: Glu dc with atoms in A or B or '' (blank)
  '''
  compare_models(pdb_str = pdb_str_005)

# ------------------------------------------------------------------------------

def test_006():
  '''
    Check if model without H to be placed is correctly processed
  '''
  compare_models(pdb_str = pdb_str_006)

# ------------------------------------------------------------------------------

def test_007():
  '''
    Check if mmCIF file format works
  '''
  compare_models(pdb_str = pdb_str_007)

# ------------------------------------------------------------------------------

def test_008():
  '''
    Test CYS-CYS disulfide with double conformation
  '''
  compare_models(pdb_str = pdb_str_008)

# ------------------------------------------------------------------------------

def test_009():
  '''
    1KYC: SIN linked to GLU1, but later in file. Goal: no NH3 to GLU.
  '''
  compare_models(pdb_str = pdb_str_009)
  assert 0

# ------------------------------------------------------------------------------

pdb_str_000 = """
REMARK This will crash if CDL is set to FALSE
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   PRO H   1      52.628 -74.147  33.427  1.00 20.43           N
ATOM      2  CA  PRO H   1      53.440 -73.630  34.533  1.00 20.01           C
ATOM      3  C   PRO H   1      54.482 -72.584  34.124  1.00 20.76           C
ATOM      4  O   PRO H   1      55.025 -72.627  33.021  1.00 16.34           O
ATOM      5  CB  PRO H   1      54.055 -74.895  35.134  1.00 22.06           C
ATOM      6  CG  PRO H   1      54.084 -75.862  33.972  1.00 25.16           C
ATOM      7  CD  PRO H   1      52.770 -75.608  33.294  1.00 17.36           C
ATOM      8  H2  PRO H   1      51.760 -73.969  33.583  1.00 20.43           H
ATOM      9  H3  PRO H   1      52.881 -73.753  32.659  1.00 20.43           H
ATOM     10  HA  PRO H   1      52.877 -73.212  35.203  1.00 20.01           H
ATOM     11  HB2 PRO H   1      54.949 -74.711  35.461  1.00 22.06           H
ATOM     12  HB3 PRO H   1      53.499 -75.228  35.856  1.00 22.06           H
ATOM     13  HG2 PRO H   1      54.830 -75.664  33.385  1.00 25.16           H
ATOM     14  HG3 PRO H   1      54.147 -76.775  34.292  1.00 25.16           H
ATOM     15  HD3 PRO H   1      52.801 -75.872  32.361  1.00 17.36           H
ATOM     16  HD2 PRO H   1      52.048 -76.072  33.746  1.00 17.36           H
ATOM     17  N   SER H   2      54.727 -71.646  35.038  1.00 21.70           N
ATOM     18  CA  SER H   2      55.670 -70.537  34.874  1.00 25.33           C
ATOM     19  C   SER H   2      55.049 -69.401  34.057  1.00 24.78           C
ATOM     20  O   SER H   2      55.581 -68.291  34.023  1.00 27.51           O
ATOM     21  CB  SER H   2      56.982 -71.005  34.219  1.00 25.20           C
ATOM     22  OG  SER H   2      56.914 -70.938  32.802  1.00 28.91           O
ATOM     23  H   SER H   2      54.335 -71.634  35.803  1.00 21.70           H
ATOM     24  HA  SER H   2      55.899 -70.163  35.739  1.00 25.33           H
ATOM     25  HB2 SER H   2      57.705 -70.434  34.524  1.00 25.20           H
ATOM     26  HB3 SER H   2      57.151 -71.924  34.481  1.00 25.20           H
ATOM     27  HG  SER H   2      56.769 -70.147  32.558  1.00 28.91           H
ATOM     28  N   GLN H   3      53.918 -69.678  33.412  1.00 24.55           N
ATOM     29  CA  GLN H   3      53.224 -68.673  32.611  1.00 29.39           C
ATOM     30  C   GLN H   3      52.340 -67.778  33.475  1.00 28.13           C
ATOM     31  O   GLN H   3      52.234 -67.987  34.681  1.00 26.35           O
ATOM     32  CB  GLN H   3      52.371 -69.346  31.533  1.00 31.67           C
ATOM     33  CG  GLN H   3      53.196 -70.112  30.524  1.00 44.80           C
ATOM     34  CD  GLN H   3      54.379 -69.303  30.030  1.00 48.55           C
ATOM     35  OE1 GLN H   3      54.213 -68.269  29.386  1.00 52.45           O
ATOM     36  NE2 GLN H   3      55.584 -69.766  30.342  1.00 55.07           N
ATOM     37  H   GLN H   3      53.530 -70.445  33.423  1.00 24.55           H
ATOM     38  HA  GLN H   3      53.888 -68.112  32.179  1.00 29.39           H
ATOM     39  HB2 GLN H   3      51.871 -68.665  31.056  1.00 31.67           H
ATOM     40  HB3 GLN H   3      51.761 -69.970  31.957  1.00 31.67           H
ATOM     41  HG2 GLN H   3      53.533 -70.922  30.937  1.00 44.80           H
ATOM     42  HG3 GLN H   3      52.641 -70.335  29.761  1.00 44.80           H
ATOM     43 HE21 GLN H   3      56.287 -69.343  30.085  1.00 55.07           H
ATOM     44 HE22 GLN H   3      55.661 -70.489  30.801  1.00 55.07           H
"""

pdb_str_001 = """
REMARK this fails in REDUCE because of clash between LYS 246 NZ & Tam 2 C5
CRYST1   24.984   25.729   23.590  90.00  90.00  90.00 P 1
ATOM      1  N   PRO A 245      13.194  10.192  16.658  1.00 41.32           N
ATOM      2  CA  PRO A 245      12.939  11.276  15.705  1.00 43.09           C
ATOM      3  C   PRO A 245      13.983  11.305  14.601  1.00 46.16           C
ATOM      4  O   PRO A 245      15.086  10.768  14.728  1.00 44.69           O
ATOM      5  CB  PRO A 245      13.007  12.541  16.569  1.00 42.50           C
ATOM      6  CG  PRO A 245      13.795  12.147  17.772  1.00 47.69           C
ATOM      7  CD  PRO A 245      13.504  10.698  18.006  1.00 45.37           C
REMARK ATOM      8  H2  PRO A 245      12.414   9.635  16.709  1.00 41.32           H
REMARK ATOM      9  H3  PRO A 245      13.938   9.673  16.344  1.00 41.32           H
ATOM     10  HA  PRO A 245      12.051  11.211  15.320  1.00 43.09           H
ATOM     11  HB2 PRO A 245      13.452  13.251  16.080  1.00 42.50           H
ATOM     12  HB3 PRO A 245      12.112  12.820  16.817  1.00 42.50           H
ATOM     13  HG2 PRO A 245      14.740  12.283  17.601  1.00 47.69           H
ATOM     14  HG3 PRO A 245      13.515  12.679  18.533  1.00 47.69           H
ATOM     15  HD2 PRO A 245      14.279  10.246  18.375  1.00 45.37           H
ATOM     16  HD3 PRO A 245      12.744  10.590  18.598  1.00 45.37           H
ATOM     17  N   LYS A 246      13.611  11.942  13.495  1.00 42.99           N
ATOM     18  CA  LYS A 246      14.551  12.132  12.404  1.00 43.44           C
ATOM     19  C   LYS A 246      15.633  13.122  12.832  1.00 46.25           C
ATOM     20  O   LYS A 246      15.332  14.110  13.510  1.00 45.93           O
ATOM     21  CB  LYS A 246      13.837  12.652  11.156  1.00 49.81           C
ATOM     22  CG  LYS A 246      12.652  11.809  10.713  1.00 54.70           C
ATOM     23  CD  LYS A 246      13.071  10.649   9.829  1.00 62.40           C
ATOM     24  CE  LYS A 246      11.928   9.661   9.640  1.00 71.25           C
ATOM     25  NZ  LYS A 246      10.594  10.329   9.556  1.00 76.52           N
ATOM     26  H   LYS A 246      12.828  12.269  13.355  1.00 42.99           H
ATOM     27  HA  LYS A 246      14.966  11.285  12.179  1.00 43.44           H
ATOM     28  HB2 LYS A 246      14.472  12.674  10.423  1.00 49.81           H
ATOM     29  HB3 LYS A 246      13.509  13.547  11.338  1.00 49.81           H
ATOM     30  HG2 LYS A 246      12.208  11.447  11.496  1.00 54.70           H
ATOM     31  HG3 LYS A 246      12.036  12.365  10.210  1.00 54.70           H
ATOM     32  HD2 LYS A 246      13.815  10.182  10.241  1.00 62.40           H
ATOM     33  HD3 LYS A 246      13.332  10.985   8.958  1.00 62.40           H
ATOM     34  HE2 LYS A 246      11.910   9.050  10.393  1.00 71.25           H
ATOM     35  HE3 LYS A 246      12.069   9.168   8.816  1.00 71.25           H
REMARK ATOM     36  HZ1 LYS A 246       9.955   9.720   9.446  1.00 76.52           H
REMARK ATOM     37  HZ2 LYS A 246      10.578  10.892   8.867  1.00 76.52           H
REMARK ATOM     38  HZ3 LYS A 246      10.434  10.784  10.304  1.00 76.52           H
ATOM     39  N   PRO A 247      16.897  12.890  12.459  1.00 44.07           N
ATOM     40  CA  PRO A 247      17.961  13.810  12.905  1.00 39.69           C
ATOM     41  C   PRO A 247      17.659  15.272  12.622  1.00 41.96           C
ATOM     42  O   PRO A 247      17.848  16.127  13.497  1.00 43.51           O
ATOM     43  CB  PRO A 247      19.188  13.316  12.126  1.00 43.94           C
ATOM     44  CG  PRO A 247      18.912  11.884  11.852  1.00 47.67           C
ATOM     45  CD  PRO A 247      17.430  11.781  11.650  1.00 44.76           C
ATOM     46  HA  PRO A 247      18.124  13.716  13.857  1.00 39.69           H
ATOM     47  HB3 PRO A 247      19.279  13.817  11.300  1.00 43.94           H
ATOM     48  HB2 PRO A 247      19.987  13.419  12.667  1.00 43.94           H
ATOM     49  HG2 PRO A 247      19.193  11.346  12.609  1.00 47.67           H
ATOM     50  HG3 PRO A 247      19.387  11.607  11.053  1.00 47.67           H
ATOM     51  HD2 PRO A 247      17.201  11.896  10.714  1.00 44.76           H
ATOM     52  HD3 PRO A 247      17.097  10.929  11.971  1.00 44.76           H
TER
HETATM   48  N   TAM H   2       9.323  12.496   7.335  1.00 20.00           N
HETATM   49  C   TAM H   2       8.060  12.492   8.002  1.00 20.00           C
HETATM   50  C1  TAM H   2       7.540  13.901   8.071  1.00 20.00           C
HETATM   51  C2  TAM H   2       8.386  11.881   9.335  1.00 20.00           C
HETATM   52  C3  TAM H   2       7.035  11.686   7.294  1.00 20.00           C
HETATM   53  C4  TAM H   2       7.128  14.539   6.744  1.00 20.00           C
HETATM   54  C5  TAM H   2       8.930  10.458   9.271  1.00 20.00           C
HETATM   55  C6  TAM H   2       5.660  11.992   7.821  1.00 20.00           C
HETATM   56  O4  TAM H   2       5.710  14.391   6.585  1.00 20.00           O
HETATM   57  O5  TAM H   2       7.872   9.487   9.299  1.00 20.00           O
HETATM   58  O6  TAM H   2       5.714  12.262   9.200  1.00 20.00           O
HETATM   59  H11 TAM H   2       6.779  13.910   8.672  1.00 20.00           H
HETATM   60  H12 TAM H   2       8.223  14.452   8.485  1.00 20.00           H
HETATM   61  H21 TAM H   2       7.589  11.879   9.888  1.00 20.00           H
HETATM   62  H22 TAM H   2       9.038  12.437   9.790  1.00 20.00           H
HETATM   63  H31 TAM H   2       7.225  10.740   7.392  1.00 20.00           H
HETATM   64  H32 TAM H   2       7.071  11.859   6.340  1.00 20.00           H
HETATM   65  H41 TAM H   2       7.600  14.117   6.009  1.00 20.00           H
HETATM   66  H42 TAM H   2       7.377  15.476   6.732  1.00 20.00           H
HETATM   67  H61 TAM H   2       5.264  12.763   7.385  1.00 20.00           H
HETATM   68  H62 TAM H   2       5.044  11.254   7.689  1.00 20.00           H
HETATM   69  HN1 TAM H   2       9.235  12.859   6.527  1.00 20.00           H
HETATM   70  HN2 TAM H   2       9.942  12.134   7.862  1.00 20.00           H
"""

pdb_str_002 = """
REMARK This is a multi model file --> check if this works
CRYST1   16.760   20.171   17.648  90.00  90.00  90.00 P 1
MODEL        1
ATOM      1  N   GLY A  -3      14.573   7.304   5.082  1.00 23.20           N
ATOM      2  CA  GLY A  -3      15.503   6.752   6.050  1.00 43.12           C
ATOM      3  C   GLY A  -3      16.822   7.516   6.092  1.00 10.33           C
ATOM      4  O   GLY A  -3      17.833   7.000   5.608  1.00 34.23           O
ATOM      5  HA2 GLY A  -3      15.692   5.828   5.822  1.00 43.12           H
ATOM      6  HA3 GLY A  -3      15.103   6.784   6.933  1.00 43.12           H
ATOM      7  N   PRO A  -2      16.855   8.759   6.667  1.00 43.42           N
ATOM      8  CA  PRO A  -2      18.084   9.573   6.756  1.00 72.12           C
ATOM      9  C   PRO A  -2      19.050   9.085   7.843  1.00 24.33           C
ATOM     10  O   PRO A  -2      20.269   9.138   7.662  1.00  1.42           O
ATOM     11  CB  PRO A  -2      17.569  10.986   7.098  1.00 65.32           C
ATOM     12  CG  PRO A  -2      16.078  10.915   7.005  1.00 23.11           C
ATOM     13  CD  PRO A  -2      15.713   9.483   7.257  1.00  3.14           C
ATOM     14  HA  PRO A  -2      18.527   9.578   5.893  1.00 72.12           H
ATOM     15  HB2 PRO A  -2      17.846  11.226   7.996  1.00 65.32           H
ATOM     16  HB3 PRO A  -2      17.923  11.625   6.460  1.00 65.32           H
ATOM     17  HG2 PRO A  -2      15.793  11.190   6.120  1.00 23.11           H
ATOM     18  HG3 PRO A  -2      15.682  11.493   7.676  1.00 23.11           H
ATOM     19  HD2 PRO A  -2      15.639   9.304   8.207  1.00  3.14           H
ATOM     20  HD3 PRO A  -2      14.884   9.250   6.810  1.00  3.14           H
ATOM     21  N   SER A  -1      18.485   8.603   8.974  1.00 54.30           N
ATOM     22  CA  SER A  -1      19.258   8.090  10.133  1.00 74.21           C
ATOM     23  C   SER A  -1      20.126   9.182  10.770  1.00 44.31           C
ATOM     24  O   SER A  -1      20.919   9.834  10.087  1.00 24.24           O
ATOM     25  CB  SER A  -1      20.129   6.881   9.741  1.00 35.31           C
ATOM     26  OG  SER A  -1      20.678   6.247  10.885  1.00 34.24           O
ATOM     27  H   SER A  -1      17.635   8.561   9.098  1.00 54.30           H
ATOM     28  HA  SER A  -1      18.623   7.790  10.802  1.00 74.21           H
ATOM     29  HB2 SER A  -1      20.854   7.187   9.174  1.00 35.31           H
ATOM     30  HB3 SER A  -1      19.580   6.241   9.261  1.00 35.31           H
ATOM     31  HG  SER A  -1      21.157   6.789  11.312  1.00 34.24           H
TER
ENDMDL
MODEL        2
ATOM     35  N   GLY A  -3      15.598  12.155  12.730  1.00 23.20           N
ATOM     36  CA  GLY A  -3      15.801  13.217  11.761  1.00 43.12           C
ATOM     37  C   GLY A  -3      15.603  12.746  10.322  1.00 10.33           C
ATOM     38  O   GLY A  -3      14.940  11.727  10.105  1.00 34.23           O
ATOM     39  HA3 GLY A  -3      15.172  13.934  11.936  1.00 43.12           H
ATOM     40  HA2 GLY A  -3      16.704  13.561  11.848  1.00 43.12           H
ATOM     41  N   PRO A  -2      16.164  13.467   9.301  1.00 43.42           N
ATOM     42  CA  PRO A  -2      16.028  13.090   7.880  1.00 72.12           C
ATOM     43  C   PRO A  -2      16.904  11.891   7.491  1.00 24.33           C
ATOM     44  O   PRO A  -2      16.479  11.039   6.705  1.00  1.42           O
ATOM     45  CB  PRO A  -2      16.479  14.351   7.114  1.00 65.32           C
ATOM     46  CG  PRO A  -2      16.651  15.418   8.147  1.00 23.11           C
ATOM     47  CD  PRO A  -2      16.953  14.708   9.432  1.00  3.14           C
ATOM     48  HA  PRO A  -2      15.094  12.906   7.692  1.00 72.12           H
ATOM     49  HB2 PRO A  -2      17.317  14.174   6.659  1.00 65.32           H
ATOM     50  HB3 PRO A  -2      15.799  14.602   6.470  1.00 65.32           H
ATOM     51  HG3 PRO A  -2      15.832  15.932   8.224  1.00 23.11           H
ATOM     52  HG2 PRO A  -2      17.386  15.999   7.896  1.00 23.11           H
ATOM     53  HD2 PRO A  -2      17.900  14.513   9.508  1.00  3.14           H
ATOM     54  HD3 PRO A  -2      16.659  15.229  10.196  1.00  3.14           H
ATOM     55  N   SER A  -1      18.123  11.839   8.047  1.00 54.30           N
ATOM     56  CA  SER A  -1      19.064  10.752   7.768  1.00 74.21           C
ATOM     57  C   SER A  -1      19.650  10.188   9.066  1.00 44.31           C
ATOM     58  O   SER A  -1      19.683   8.969   9.256  1.00 24.24           O
ATOM     59  CB  SER A  -1      20.187  11.243   6.843  1.00 35.31           C
ATOM     60  OG  SER A  -1      20.983  10.164   6.382  1.00 34.24           O
ATOM     61  H   SER A  -1      18.429  12.428   8.594  1.00 54.30           H
ATOM     62  HA  SER A  -1      18.597  10.032   7.316  1.00 74.21           H
ATOM     63  HB2 SER A  -1      20.750  11.861   7.334  1.00 35.31           H
ATOM     64  HB3 SER A  -1      19.791  11.690   6.079  1.00 35.31           H
ATOM     65  HG  SER A  -1      21.333   9.767   7.034  1.00 34.24           H
TER
ENDMDL
"""

pdb_str_003 = """
REMARK PDB snippet with a normal and a modified nucleic acid
CRYST1   17.826   22.060   19.146  90.00  90.00  90.00 P 1
ATOM      1  P     U A   2       7.236  16.525   9.726  1.00 37.21           P
ATOM      2  OP1   U A   2       6.663  17.060  10.993  1.00 38.75           O
ATOM      3  OP2   U A   2       8.650  16.805   9.385  1.00 34.84           O
ATOM      4  O5'   U A   2       7.029  14.947   9.685  1.00 34.31           O
ATOM      5  C5'   U A   2       5.756  14.354   9.974  1.00 37.02           C
ATOM      6  C4'   U A   2       5.821  12.860   9.763  1.00 36.81           C
ATOM      7  O4'   U A   2       5.957  12.552   8.350  1.00 36.96           O
ATOM      8  C3'   U A   2       7.025  12.178  10.388  1.00 36.55           C
ATOM      9  O3'   U A   2       6.846  11.969  11.775  1.00 36.22           O
ATOM     10  C2'   U A   2       7.106  10.884   9.592  1.00 36.94           C
ATOM     11  O2'   U A   2       6.138   9.930   9.980  1.00 38.55           O
ATOM     12  C1'   U A   2       6.755  11.383   8.191  1.00 36.22           C
ATOM     13  N1    U A   2       7.938  11.735   7.391  1.00 34.47           N
ATOM     14  C2    U A   2       8.652  10.705   6.795  1.00 33.67           C
ATOM     15  O2    U A   2       8.375   9.526   6.954  1.00 33.54           O
ATOM     16  N3    U A   2       9.706  11.110   6.014  1.00 33.13           N
ATOM     17  C4    U A   2      10.119  12.408   5.784  1.00 32.43           C
ATOM     18  O4    U A   2      11.043  12.620   5.000  1.00 31.05           O
ATOM     19  C5    U A   2       9.352  13.409   6.464  1.00 33.80           C
ATOM     20  C6    U A   2       8.316  13.048   7.223  1.00 33.46           C
ATOM     21  H4*   U A   2       5.000  12.459  10.088  1.00 36.81           H
ATOM     22  H3*   U A   2       7.826  12.702  10.232  1.00 36.55           H
ATOM     23  H2*   U A   2       7.997  10.505   9.648  1.00 36.94           H
ATOM     24 HO2*   U A   2       5.578   9.841   9.360  1.00 38.55           H
ATOM     25  H1*   U A   2       6.256  10.694   7.724  1.00 36.22           H
ATOM     26  H5    U A   2       9.576  14.307   6.377  1.00 33.80           H
ATOM     27 H5*2   U A   2       5.515  14.540  10.895  1.00 37.02           H
ATOM     28 H5*1   U A   2       5.082  14.734   9.388  1.00 37.02           H
ATOM     29  H3    U A   2      10.157  10.489   5.626  1.00 33.13           H
ATOM     30  H6    U A   2       7.829  13.711   7.657  1.00 33.46           H
HETATM   31  P   UMS A   3       8.115  12.049  12.757  1.00 35.65           P
HETATM   32  OP1 UMS A   3       7.608  12.006  14.146  1.00 37.56           O
HETATM   33  OP2 UMS A   3       8.986  13.173  12.329  1.00 35.79           O
HETATM   34  O5' UMS A   3       8.946  10.722  12.449  1.00 37.12           O
HETATM   35  C5' UMS A   3       8.323   9.447  12.591  1.00 36.88           C
HETATM   36  C4' UMS A   3       9.203   8.363  12.041  1.00 39.11           C
HETATM   37  O4' UMS A   3       9.187   8.277  10.585  1.00 38.36           O
HETATM   38  C3' UMS A   3      10.640   8.278  12.458  1.00 40.33           C
HETATM   39  O3' UMS A   3      10.895   8.001  13.828  1.00 40.80           O
HETATM   40  C2' UMS A   3      11.291   7.467  11.374  1.00 40.93           C
HETATM   41  C1' UMS A   3      10.496   7.904  10.133  1.00 38.75           C
HETATM   42  N1  UMS A   3      11.106   9.042   9.429  1.00 36.80           N
HETATM   43  C2  UMS A   3      12.013   8.736   8.436  1.00 35.25           C
HETATM   44  O2  UMS A   3      12.328   7.589   8.163  1.00 36.36           O
HETATM   45  N3  UMS A   3      12.557   9.822   7.786  1.00 33.98           N
HETATM   46  C4  UMS A   3      12.275  11.159   7.996  1.00 32.51           C
HETATM   47  O4  UMS A   3      12.826  12.001   7.302  1.00 33.01           O
HETATM   48  C5  UMS A   3      11.323  11.402   9.046  1.00 33.90           C
HETATM   49  C6  UMS A   3      10.817  10.354   9.743  1.00 34.62           C
HETATM   50  CA' UMS A   3      11.776   5.000  13.480  1.00 46.26           C
HETATM   51 SE2' UMS A   3      10.815   5.400  11.905  1.00 48.17          Se
HETATM   53  H3* UMS A   3      11.022   9.152  12.279  1.00 40.33           H
HETATM   54  H3  UMS A   3      13.140   9.648   7.179  1.00 33.98           H
HETATM   55  H1* UMS A   3      10.387   7.156   9.525  1.00 38.75           H
HETATM   56  H6  UMS A   3      10.252  10.527  10.462  1.00 34.62           H
HETATM   57  H4* UMS A   3       8.782   7.553  12.369  1.00 39.11           H
HETATM   58 H5*2 UMS A   3       7.481   9.448  12.109  1.00 36.88           H
HETATM   59  H5  UMS A   3      11.056  12.270   9.246  1.00 33.90           H
HETATM   60  H5* UMS A   3       8.158   9.276  13.531  1.00 36.88           H
HETATM   64  H2* UMS A   3      12.022   6.833  11.440  1.00 40.93           H
"""

pdb_str_004 = '''
REMARK TYR double conformation where *every* atom is in either A or B
CRYST1   15.639   15.148   16.657  90.00  90.00  90.00 P 1
ATOM      1  N  ATYR A  59       5.624   5.492   5.997  0.63  5.05           N
ATOM      2  CA ATYR A  59       6.283   5.821   7.250  0.63  5.48           C
ATOM      3  C  ATYR A  59       5.451   6.841   8.030  0.63  6.01           C
ATOM      4  O  ATYR A  59       5.000   7.863   7.506  0.63  6.38           O
ATOM      5  CB ATYR A  59       7.724   6.421   6.963  0.63  5.57           C
ATOM      6  CG ATYR A  59       8.212   7.215   8.170  0.63  6.71           C
ATOM      7  CD1ATYR A  59       8.690   6.541   9.297  0.63  7.05           C
ATOM      8  CD2ATYR A  59       8.071   8.583   8.242  0.63  8.31           C
ATOM      9  CE1ATYR A  59       9.100   7.172  10.481  0.63  7.99           C
ATOM     10  CE2ATYR A  59       8.408   9.230   9.447  0.63  9.07           C
ATOM     11  CZ ATYR A  59       8.919   8.547  10.507  0.63  9.01           C
ATOM     12  OH ATYR A  59       9.211   9.255  11.657  0.63 12.31           O
ATOM     13  HA ATYR A  59       6.384   5.020   7.788  0.63  5.48           H
ATOM     14  HB2ATYR A  59       7.683   7.014   6.196  0.63  5.57           H
ATOM     15  HB3ATYR A  59       8.349   5.699   6.792  0.63  5.57           H
ATOM     16  HD1ATYR A  59       8.740   5.613   9.260  0.63  7.05           H
ATOM     17  HD2ATYR A  59       7.760   9.068   7.512  0.63  8.31           H
ATOM     18  HE1ATYR A  59       9.466   6.703  11.196  0.63  7.99           H
ATOM     19  HE2ATYR A  59       8.278  10.148   9.520  0.63  9.07           H
ATOM     20  HH ATYR A  59       9.058  10.072  11.538  0.63 12.31           H
ATOM     21  N  BTYR A  59       5.613   5.513   5.963  0.37  5.75           N
ATOM     22  CA BTYR A  59       6.322   5.809   7.211  0.37  5.49           C
ATOM     23  C  BTYR A  59       5.795   6.953   8.094  0.37  5.14           C
ATOM     24  O  BTYR A  59       5.668   8.090   7.641  0.37  6.42           O
ATOM     25  CB BTYR A  59       7.798   6.076   6.900  0.37  7.77           C
ATOM     26  CG BTYR A  59       8.556   6.722   8.038  0.37  5.20           C
ATOM     27  CD1BTYR A  59       9.162   5.951   9.021  0.37  8.94           C
ATOM     28  CD2BTYR A  59       8.665   8.103   8.129  0.37  6.25           C
ATOM     29  CE1BTYR A  59       9.856   6.537  10.063  0.37 11.97           C
ATOM     30  CE2BTYR A  59       9.357   8.699   9.167  0.37  9.52           C
ATOM     31  CZ BTYR A  59       9.950   7.911  10.131  0.37 12.68           C
ATOM     32  OH BTYR A  59      10.639   8.500  11.166  0.37 26.50           O
ATOM     33  HA BTYR A  59       6.185   5.019   7.758  0.37  5.49           H
ATOM     34  HB2BTYR A  59       7.853   6.669   6.134  0.37  7.77           H
ATOM     35  HB3BTYR A  59       8.231   5.232   6.698  0.37  7.77           H
ATOM     36  HD1BTYR A  59       9.100   5.024   8.978  0.37  8.94           H
ATOM     37  HD2BTYR A  59       8.266   8.637   7.480  0.37  6.25           H
ATOM     38  HE1BTYR A  59      10.257   6.008  10.714  0.37 11.97           H
ATOM     39  HE2BTYR A  59       9.422   9.625   9.215  0.37  9.52           H
ATOM     40  HH BTYR A  59      10.618   9.336  11.086  0.37 26.50           H
'''

pdb_str_005 = '''
REMARK Glu double conformation where atoms are either A, B or '' (blank)
CRYST1   13.702   13.985   14.985  90.00  90.00  90.00 P 1
ATOM      1  N   GLU A  78       8.702   8.360   5.570  1.00 35.65           N
ATOM      2  C   GLU A  78       6.379   7.842   5.202  1.00 35.59           C
ATOM      3  O   GLU A  78       5.975   8.985   5.000  1.00 35.38           O
ATOM      4  CA AGLU A  78       7.598   7.571   6.076  0.70 35.57           C
ATOM      5  CB AGLU A  78       7.301   7.887   7.536  0.70 35.75           C
ATOM      6  CG AGLU A  78       6.481   6.798   8.188  0.70 36.10           C
ATOM      7  CD AGLU A  78       5.833   7.232   9.476  0.70 37.70           C
ATOM      8  OE1AGLU A  78       6.155   8.333   9.982  0.70 38.74           O
ATOM      9  OE2AGLU A  78       5.000   6.456   9.985  0.70 37.65           O
ATOM     10  HA AGLU A  78       7.819   6.627   6.051  0.70 35.57           H
ATOM     11  HB2AGLU A  78       6.802   8.717   7.588  0.70 35.75           H
ATOM     12  HB3AGLU A  78       8.137   7.971   8.021  0.70 35.75           H
ATOM     13  HG2AGLU A  78       5.778   6.526   7.578  0.70 36.10           H
ATOM     14  HG3AGLU A  78       7.059   6.044   8.385  0.70 36.10           H
ATOM     15  CA BGLU A  78       7.581   7.608   6.115  0.30 35.61           C
ATOM     16  CB BGLU A  78       7.269   8.093   7.534  0.30 35.70           C
ATOM     17  CG BGLU A  78       6.166   7.322   8.245  0.30 36.05           C
ATOM     18  CD BGLU A  78       6.683   6.115   9.003  0.30 36.79           C
ATOM     19  OE1BGLU A  78       7.585   6.285   9.856  0.30 37.19           O
ATOM     20  OE2BGLU A  78       6.173   5.000   8.760  0.30 37.31           O
ATOM     21  HA BGLU A  78       7.779   6.660   6.170  0.30 35.61           H
ATOM     22  HB2BGLU A  78       8.073   8.013   8.070  0.30 35.70           H
ATOM     23  HB3BGLU A  78       6.992   9.022   7.488  0.30 35.70           H
ATOM     24  HG3BGLU A  78       5.525   7.010   7.587  0.30 36.05           H
ATOM     25  HG2BGLU A  78       5.730   7.910   8.881  0.30 36.05           H
'''

pdb_str_006 = '''
REMARK Hg, Sr and HOH: no H atoms are expected to be placed
CRYST1   10.286   24.260   13.089  90.00  90.00  90.00 P 1
HETATM    1 HG    HG A 101       5.000   7.056   5.951  0.60 17.71          HG
HETATM    2 SR    SR A 102       5.182  10.793   8.089  0.85 18.78          SR
HETATM    3  O   HOH A 201       5.093   5.000   5.000  1.00 25.34           O
HETATM    4  O   HOH A 202       5.286  19.260   7.818  1.00 28.43           O
'''

# check if mmCIF file works
pdb_str_007 = '''
data_1UBQ
#
_entry.id   1UBQ
#
_cell.entry_id           1UBQ
_cell.length_a           50.840
_cell.length_b           42.770
_cell.length_c           28.950
_cell.angle_alpha        90.00
_cell.angle_beta         90.00
_cell.angle_gamma        90.00
_cell.Z_PDB              4
#
_symmetry.entry_id                         1UBQ
_symmetry.space_group_name_H-M             'P 21 21 21'
_symmetry.Int_Tables_number                19
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   594 N N   . GLY A 1 75 ? 41.165 35.531 31.898 0.25 36.31 ? 75  GLY A N   1
ATOM   595 C CA  . GLY A 1 75 ? 41.845 36.550 32.686 0.25 36.07 ? 75  GLY A CA  1
ATOM   596 C C   . GLY A 1 75 ? 41.251 37.941 32.588 0.25 36.16 ? 75  GLY A C   1
ATOM   597 O O   . GLY A 1 75 ? 41.102 38.523 31.500 0.25 36.26 ? 75  GLY A O   1
ATOM   598 H HA2 . GLY A 1 75 ? 42.768 36.603 32.393 0.25 36.07 ? 75  GLY A HA3 1
ATOM   599 H HA3 . GLY A 1 75 ? 41.823 36.286 33.619 0.25 36.07 ? 75  GLY A HA2 1
ATOM   600 N N   . GLY A 1 76 ? 40.946 38.472 33.757 0.25 36.05 ? 76  GLY A N   1
ATOM   601 C CA  . GLY A 1 76 ? 40.373 39.813 33.944 0.25 36.19 ? 76  GLY A CA  1
ATOM   602 C C   . GLY A 1 76 ? 40.031 39.992 35.432 0.25 36.20 ? 76  GLY A C   1
ATOM   603 O O   . GLY A 1 76 ? 38.933 40.525 35.687 0.25 36.13 ? 76  GLY A O   1
ATOM   604 O OXT . GLY A 1 76 ? 40.862 39.575 36.251 0.25 36.27 ? 76  GLY A OXT 1
ATOM   605 H H   . GLY A 1 76 ? 41.063 38.062 34.504 0.25 36.05 ? 76  GLY A H   1
ATOM   606 H HA3 . GLY A 1 76 ? 39.566 39.910 33.413 0.25 36.19 ? 76  GLY A HA3 1
ATOM   607 H HA2 . GLY A 1 76 ? 41.011 40.491 33.675 0.25 36.19 ? 76  GLY A HA2 1
HETATM 611 O O   . HOH B 2 .  ? 45.747 30.081 19.708 1.00 12.43 ? 77  HOH A O   1
'''

pdb_str_008 = '''
REMARK disulfide bond with altloc
CRYST1   13.626   15.799   17.617  90.00  90.00  90.00 P 1
ATOM      1  N   CYS A  27      17.541   4.439  12.897  1.00 13.99           N
ATOM      2  CA  CYS A  27      16.566   5.527  12.862  1.00 14.57           C
ATOM      3  C   CYS A  27      16.236   6.026  11.467  1.00 14.53           C
ATOM      4  O   CYS A  27      15.254   6.760  11.351  1.00 16.95           O
ATOM      5  CB  CYS A  27      17.114   6.698  13.662  1.00 15.77           C
ATOM      6  SG  CYS A  27      17.230   6.332  15.443  1.00 17.57           S
ATOM      7  HA  CYS A  27      15.739   5.186  13.237  1.00 14.57           H
ATOM      8  HB2 CYS A  27      16.526   7.461  13.549  1.00 15.77           H
ATOM      9  HB3 CYS A  27      18.003   6.913  13.340  1.00 15.77           H
ATOM     10  CB  CYS A 123      14.607   7.591  16.260  1.00 24.16           C
ATOM     11  SG  CYS A 123      15.316   5.939  15.946  1.00 20.05           S
ATOM     12  N  ACYS A 123      15.023   7.279  18.624  0.58 26.40           N
ATOM     14  CA ACYS A 123      15.266   8.190  17.491  0.58 25.69           C
ATOM     15  C  ACYS A 123      14.764   9.599  17.776  0.58 26.33           C
ATOM     16  O  ACYS A 123      14.197  10.238  16.886  0.58 28.70           O
ATOM     17  OXTACYS A 123      14.975  10.081  18.878  0.58 28.31           O
ATOM     18  HA ACYS A 123      16.217   8.287  17.324  0.58 25.69           H
ATOM     19  HB2ACYS A 123      13.652   7.502  16.408  1.00 24.16           H
ATOM     20  HB3ACYS A 123      14.772   8.157  15.490  1.00 24.16           H
ATOM     22  N  BCYS A 123      15.023   7.288  18.685  0.42 25.68           N
ATOM     23  CA BCYS A 123      15.108   8.205  17.548  0.42 25.86           C
ATOM     24  C  BCYS A 123      14.270   9.460  17.813  0.42 26.42           C
ATOM     26  O  BCYS A 123      13.915  10.125  16.837  0.42 27.75           O
ATOM     27  OXTBCYS A 123      13.981   9.728  18.968  0.42 28.04           O
ATOM     28  HA BCYS A 123      16.045   8.426  17.432  0.42 25.86           H
ATOM     29  HB2BCYS A 123      13.642   7.500  16.307  1.00 24.16           H
ATOM     30  HB3BCYS A 123      14.850   8.168  15.519  1.00 24.16           H
'''

pdb_str_009 = '''
REMARK 1KYC: SIN is linked to GLU 1, but at end of model file
CRYST1   19.769   30.000   32.270  90.00  90.00  90.00 P 1
ATOM      1  N   GLU A   1       9.100  10.903   6.846  1.00 18.32           N
ATOM      2  CA  GLU A   1       8.321  10.315   7.917  1.00 16.97           C
ATOM      3  C   GLU A   1       9.162   9.612   8.973  1.00 15.25           C
ATOM      4  O   GLU A   1       8.910   9.774  10.170  1.00 14.23           O
ATOM      5  CB  GLU A   1       7.291   9.304   7.380  1.00 17.54           C
ATOM      6  CG  GLU A   1       6.485   8.581   8.454  1.00 19.79           C
ATOM      7  CD  GLU A   1       5.488   9.479   9.156  1.00 20.81           C
ATOM      8  OE1 GLU A   1       5.000  10.445   8.534  1.00 23.39           O
ATOM      9  OE2 GLU A   1       5.180   9.223  10.333  1.00 21.60           O
ATOM     10  HA  GLU A   1       7.863  11.056   8.343  1.00 16.97           H
ATOM     11  HB2 GLU A   1       6.662   9.777   6.812  1.00 17.54           H
ATOM     12  HB3 GLU A   1       7.761   8.630   6.865  1.00 17.54           H
ATOM     13  HG2 GLU A   1       7.094   8.232   9.123  1.00 19.79           H
ATOM     14  HG3 GLU A   1       5.993   7.854   8.042  1.00 19.79           H
ATOM     15  N   GLU A   2      10.111   8.769   8.574  1.00 15.28           N
ATOM     16  CA  GLU A   2      10.916   8.018   9.519  1.00 14.99           C
ATOM     17  C   GLU A   2      11.541   8.882  10.615  1.00 13.85           C
ATOM     18  O   GLU A   2      11.362   8.667  11.822  1.00 13.88           O
ATOM     19  CB  GLU A   2      12.024   7.270   8.757  1.00 16.20           C
ATOM     20  CG  GLU A   2      12.869   6.426   9.713  1.00 17.22           C
ATOM     21  CD  GLU A   2      14.152   5.903   9.124  1.00 17.43           C
ATOM     22  OE1 GLU A   2      14.556   6.396   8.052  1.00 18.62           O
ATOM     23  OE2 GLU A   2      14.769   5.000   9.734  1.00 17.26           O
ATOM     24  H   GLU A   2      10.307   8.616   7.751  1.00 15.28           H
ATOM     25  HA  GLU A   2      10.334   7.388   9.972  1.00 14.99           H
ATOM     26  HB2 GLU A   2      12.604   7.912   8.318  1.00 16.20           H
ATOM     27  HB3 GLU A   2      11.622   6.681   8.099  1.00 16.20           H
ATOM     28  HG2 GLU A   2      13.101   6.970  10.482  1.00 17.22           H
ATOM     29  HG3 GLU A   2      12.344   5.661   9.995  1.00 17.22           H
HETATM   30  O1  SO4 A 101      11.629  25.000  26.048  0.50 20.11           O
HETATM   31  O2  SO4 A 101      10.541  22.822  26.297  0.50 21.85           O
HETATM   32  O3  SO4 A 101       9.524  24.537  24.875  0.50 20.45           O
HETATM   33  O4  SO4 A 101       9.520  24.852  27.270  0.50 19.22           O
HETATM   34  S   SO4 A 101      10.313  24.291  26.126  0.50 19.91           S
HETATM   35  C1  SIN A   0      10.126  11.723   7.051  1.00 19.39           C
HETATM   36  C2  SIN A   0      10.828  12.260   5.817  1.00 21.62           C
HETATM   37  C3  SIN A   0      11.448  13.629   6.056  1.00 22.69           C
HETATM   38  C4  SIN A   0      10.468  14.784   6.095  1.00 23.92           C
HETATM   39  O1  SIN A   0      10.519  11.999   8.189  1.00 19.91           O
HETATM   40  O3  SIN A   0      10.949  15.933   6.216  1.00 24.14           O
HETATM   41  O4  SIN A   0       9.264  14.547   6.337  1.00 25.37           O
HETATM   42  H21 SIN A   0      10.110  12.330   5.000  1.00 21.62           H
HETATM   43  H22 SIN A   0      11.609  11.561   5.518  1.00 21.62           H
HETATM   44  H31 SIN A   0      11.985  13.602   7.004  1.00 22.69           H
HETATM   45  H32 SIN A   0      12.173  13.820   5.265  1.00 22.69           H
'''

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
