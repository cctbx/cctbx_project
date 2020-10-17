from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce
from libtbx.utils import null_out


def run():
  correct_H_position_with_cdl()
  snippet_with_clashing_H()


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

  model_h_added = reduce.add(model = model_without_h)
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


def correct_H_position_with_cdl():
  compare_models(pdb_str = pdb_str_1)

def snippet_with_clashing_H():
  '''
  This is an example for a file that fails with REDUCE.

  The NZ atom of Lys 246 and the C5 atom of Tam 2 are extremely close together,
  which seems to give Reduce a problem.

  Removing either Lys 246 NZ, or TAM 2 C5, or both of Lys 246 HZ1 and HZ3 will
  allow Reduce to run to completion and produce a usable result.
  However, Lys 246 HZ1 and HZ3 will not be added back by Reduce.
  '''
  compare_models(pdb_str = pdb_str_2)


pdb_str_1 = """
CRYST1   72.240   72.010   86.990  90.00  90.00  90.00 P 21 21 21
SCALE1      0.013843  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013887  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011496        0.00000
ATOM      1  N   PRO H  14      52.628 -74.147  33.427  1.00 20.43           N
ATOM      2  CA  PRO H  14      53.440 -73.630  34.533  1.00 20.01           C
ATOM      3  C   PRO H  14      54.482 -72.584  34.124  1.00 20.76           C
ATOM      4  O   PRO H  14      55.025 -72.627  33.021  1.00 16.34           O
ATOM      5  CB  PRO H  14      54.055 -74.895  35.134  1.00 22.06           C
ATOM      6  CG  PRO H  14      54.084 -75.862  33.972  1.00 25.16           C
ATOM      7  CD  PRO H  14      52.770 -75.608  33.294  1.00 17.36           C
ATOM      8  HD2 PRO H  14      52.801 -75.872  32.361  1.00 17.36           H
ATOM      9  HD1 PRO H  14      52.048 -76.072  33.746  1.00 17.36           H
ATOM     10  HG2 PRO H  14      54.830 -75.664  33.385  1.00 25.16           H
ATOM     11  HG1 PRO H  14      54.147 -76.775  34.292  1.00 25.16           H
ATOM     12  HA  PRO H  14      52.877 -73.212  35.203  1.00 20.01           H
ATOM     13  HB1 PRO H  14      54.949 -74.711  35.461  1.00 22.06           H
ATOM     14  HB2 PRO H  14      53.499 -75.228  35.856  1.00 22.06           H
ATOM     15  N   SER H  15      54.727 -71.646  35.038  1.00 21.70           N
ATOM     16  CA  SER H  15      55.670 -70.537  34.874  1.00 25.33           C
ATOM     17  C   SER H  15      55.049 -69.401  34.057  1.00 24.78           C
ATOM     18  O   SER H  15      55.581 -68.291  34.023  1.00 27.51           O
ATOM     19  CB  SER H  15      56.982 -71.005  34.219  1.00 25.20           C
ATOM     20  OG  SER H  15      56.914 -70.938  32.802  1.00 28.91           O
ATOM     21  H   SER H  15      54.335 -71.634  35.803  1.00 21.70           H
ATOM     22  HA  SER H  15      55.899 -70.163  35.739  1.00 25.33           H
ATOM     23  HB1 SER H  15      57.705 -70.434  34.524  1.00 25.20           H
ATOM     24  HG  SER H  15      56.225 -70.518  32.567  1.00 28.91           H
ATOM     25  HB2 SER H  15      57.151 -71.924  34.481  1.00 25.20           H
ATOM     26  N   GLN H  16      53.918 -69.678  33.412  1.00 24.55           N
ATOM     27  CA  GLN H  16      53.224 -68.673  32.611  1.00 29.39           C
ATOM     28  C   GLN H  16      52.340 -67.778  33.475  1.00 28.13           C
ATOM     29  O   GLN H  16      52.234 -67.987  34.681  1.00 26.35           O
ATOM     30  CB  GLN H  16      52.371 -69.346  31.533  1.00 31.67           C
ATOM     31  CG  GLN H  16      53.196 -70.112  30.524  1.00 44.80           C
ATOM     32  CD  GLN H  16      54.379 -69.303  30.030  1.00 48.55           C
ATOM     33  OE1 GLN H  16      54.213 -68.269  29.386  1.00 52.45           O
ATOM     34  NE2 GLN H  16      55.584 -69.766  30.342  1.00 55.07           N
ATOM     35  H   GLN H  16      53.530 -70.445  33.423  1.00 24.55           H
ATOM     36  HG2 GLN H  16      53.533 -70.922  30.937  1.00 44.80           H
ATOM     37  HG1 GLN H  16      52.641 -70.335  29.761  1.00 44.80           H
ATOM     38 HE22 GLN H  16      55.661 -70.489  30.801  1.00 55.07           H
ATOM     39 HE21 GLN H  16      56.287 -69.343  30.085  1.00 55.07           H
ATOM     40  HA  GLN H  16      53.888 -68.112  32.179  1.00 29.39           H
ATOM     41  HB1 GLN H  16      51.871 -68.665  31.056  1.00 31.67           H
ATOM     42  HB2 GLN H  16      51.761 -69.970  31.957  1.00 31.67           H
TER
"""

pdb_str_2 = """
CRYST1   24.984   25.729   23.590  90.00  90.00  90.00 P 1
ATOM      1  N   PRO A 245      13.194  10.192  16.658  1.00 41.32           N
ANISOU    1  N   PRO A 245     5445   5107   5149   -761   -337   -188       N
ATOM      2  CA  PRO A 245      12.939  11.276  15.705  1.00 43.09           C
ANISOU    2  CA  PRO A 245     5616   5296   5462   -703   -345   -192       C
ATOM      3  C   PRO A 245      13.983  11.305  14.601  1.00 46.16           C
ANISOU    3  C   PRO A 245     6029   5612   5898   -683   -379   -133       C
ATOM      4  O   PRO A 245      15.086  10.768  14.728  1.00 44.69           O
ANISOU    4  O   PRO A 245     5890   5400   5691   -698   -398   -108       O
ATOM      5  CB  PRO A 245      13.007  12.541  16.569  1.00 42.50           C
ANISOU    5  CB  PRO A 245     5520   5200   5429   -665   -331   -270       C
ATOM      6  CG  PRO A 245      13.795  12.147  17.772  1.00 47.69           C
ANISOU    6  CG  PRO A 245     6219   5884   6018   -710   -338   -300       C
ATOM      7  CD  PRO A 245      13.504  10.698  18.006  1.00 45.37           C
ANISOU    7  CD  PRO A 245     5958   5649   5630   -763   -332   -248       C
ATOM      8  HD2 PRO A 245      14.279  10.246  18.375  1.00 45.37           H
ATOM      9  HD1 PRO A 245      12.744  10.590  18.598  1.00 45.37           H
ATOM     10  HG2 PRO A 245      14.740  12.283  17.601  1.00 47.69           H
ATOM     11  HG1 PRO A 245      13.515  12.679  18.533  1.00 47.69           H
ATOM     12  HA  PRO A 245      12.051  11.211  15.320  1.00 43.09           H
ATOM     13  HB1 PRO A 245      13.452  13.251  16.080  1.00 42.50           H
ATOM     14  HB2 PRO A 245      12.112  12.820  16.817  1.00 42.50           H
ATOM     15  N   LYS A 246      13.611  11.942  13.495  1.00 42.99           N
ANISOU   15  N   LYS A 246     5586   5194   5553   -642   -384   -103       N
ATOM     16  CA  LYS A 246      14.551  12.132  12.404  1.00 43.44           C
ANISOU   16  CA  LYS A 246     5657   5196   5651   -626   -410    -46       C
ATOM     17  C   LYS A 246      15.633  13.122  12.832  1.00 46.25           C
ANISOU   17  C   LYS A 246     6029   5476   6066   -612   -418    -67       C
ATOM     18  O   LYS A 246      15.332  14.110  13.510  1.00 45.93           O
ANISOU   18  O   LYS A 246     5978   5405   6070   -590   -399   -123       O
ATOM     19  CB  LYS A 246      13.837  12.652  11.156  1.00 49.81           C
ANISOU   19  CB  LYS A 246     6412   6022   6492   -584   -414      4       C
ATOM     20  CG  LYS A 246      12.652  11.809  10.713  1.00 54.70           C
ANISOU   20  CG  LYS A 246     6990   6749   7043   -615   -407      9       C
ATOM     21  CD  LYS A 246      13.071  10.649   9.829  1.00 62.40           C
ANISOU   21  CD  LYS A 246     8002   7745   7964   -674   -416     37       C
ATOM     22  CE  LYS A 246      11.928   9.661   9.640  1.00 71.25           C
ANISOU   22  CE  LYS A 246     9098   8967   9007   -746   -399     13       C
ATOM     23  NZ  LYS A 246      10.594  10.329   9.556  1.00 76.52           N
ANISOU   23  NZ  LYS A 246     9661   9743   9672   -714   -398     10       N
ATOM     24  HE1 LYS A 246      11.910   9.050  10.393  1.00 71.25           H
ATOM     25  HE2 LYS A 246      12.069   9.168   8.816  1.00 71.25           H
ATOM     26  H   LYS A 246      12.828  12.269  13.355  1.00 42.99           H
ATOM     27  HG2 LYS A 246      12.208  11.447  11.496  1.00 54.70           H
ATOM     28  HG1 LYS A 246      12.036  12.365  10.210  1.00 54.70           H
ATOM     29  HD1 LYS A 246      13.815  10.182  10.241  1.00 62.40           H
ATOM     30  HD2 LYS A 246      13.332  10.985   8.958  1.00 62.40           H
ATOM     31  HZ1 LYS A 246      10.608  11.094  10.010  1.00 76.52           H
ATOM     32  HZ3 LYS A 246       9.966   9.800   9.900  1.00 76.52           H
ATOM     33  HZ2 LYS A 246      10.393  10.501   8.706  1.00 76.52           H
ATOM     34  HA  LYS A 246      14.966  11.285  12.179  1.00 43.44           H
ATOM     35  HB1 LYS A 246      14.472  12.674  10.423  1.00 49.81           H
ATOM     36  HB2 LYS A 246      13.509  13.547  11.338  1.00 49.81           H
ATOM     37  N   PRO A 247      16.897  12.890  12.459  1.00 44.07           N
ANISOU   37  N   PRO A 247     5777   5176   5792   -629   -439    -35       N
ATOM     38  CA  PRO A 247      17.961  13.810  12.905  1.00 39.69           C
ANISOU   38  CA  PRO A 247     5226   4573   5283   -643   -446    -65       C
ATOM     39  C   PRO A 247      17.659  15.272  12.622  1.00 41.96           C
ANISOU   39  C   PRO A 247     5505   4774   5665   -619   -425    -77       C
ATOM     40  O   PRO A 247      17.848  16.127  13.497  1.00 43.51           O
ANISOU   40  O   PRO A 247     5712   4925   5896   -637   -407   -153       O
ATOM     41  CB  PRO A 247      19.188  13.316  12.126  1.00 43.94           C
ANISOU   41  CB  PRO A 247     5766   5120   5808   -656   -468     -8       C
ATOM     42  CG  PRO A 247      18.912  11.884  11.852  1.00 47.67           C
ANISOU   42  CG  PRO A 247     6260   5641   6212   -648   -469     24       C
ATOM     43  CD  PRO A 247      17.430  11.781  11.650  1.00 44.76           C
ANISOU   43  CD  PRO A 247     5884   5287   5837   -641   -451     19       C
ATOM     44  HD2 PRO A 247      17.201  11.896  10.714  1.00 44.76           H
ATOM     45  HD1 PRO A 247      17.097  10.929  11.971  1.00 44.76           H
ATOM     46  HG2 PRO A 247      19.193  11.346  12.609  1.00 47.67           H
ATOM     47  HG1 PRO A 247      19.387  11.607  11.053  1.00 47.67           H
ATOM     48  HA  PRO A 247      18.124  13.716  13.857  1.00 39.69           H
ATOM     49  HB1 PRO A 247      19.279  13.817  11.300  1.00 43.94           H
ATOM     50  HB2 PRO A 247      19.987  13.419  12.667  1.00 43.94           H
TER
HETATM   51  N   TAM H   2       9.323  12.496   7.335  1.00 20.00           N
HETATM   52  C   TAM H   2       8.060  12.492   8.002  1.00 20.00           C
HETATM   53  C1  TAM H   2       7.540  13.901   8.071  1.00 20.00           C
HETATM   54  C2  TAM H   2       8.386  11.881   9.335  1.00 20.00           C
HETATM   55  C3  TAM H   2       7.035  11.686   7.294  1.00 20.00           C
HETATM   56  C4  TAM H   2       7.128  14.539   6.744  1.00 20.00           C
HETATM   57  C5  TAM H   2       8.930  10.458   9.271  1.00 20.00           C
HETATM   58  C6  TAM H   2       5.660  11.992   7.821  1.00 20.00           C
HETATM   59  O4  TAM H   2       5.710  14.391   6.585  1.00 20.00           O
HETATM   60  O5  TAM H   2       7.872   9.487   9.299  1.00 20.00           O
HETATM   61  O6  TAM H   2       5.714  12.262   9.200  1.00 20.00           O
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
