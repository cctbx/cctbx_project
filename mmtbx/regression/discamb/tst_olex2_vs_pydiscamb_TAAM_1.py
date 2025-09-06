from __future__ import absolute_import, division, print_function
import os, time
from libtbx.utils import null_out
import mmtbx.model
import iotbx.pdb
import libtbx.load_env
from cctbx.xray import structure
from iotbx import reflection_file_reader
import pydiscamb
from mmtbx.regression.discamb import tst_IAM_discamb_vs_cctbx_1 as regression_discamb

verbose = False

# ------------------------------------------------------------------------------

def has_unassigned(atom_type_assignment):
  has_000 = any(value[0].endswith('000')
    for value in atom_type_assignment.values())
  has_IAM = any(value[0] == '' for value in atom_type_assignment.values())
  return has_000 or has_IAM

# ------------------------------------------------------------------------------

def run(pdb_str, mtz_fn, target_scores):
  regression_mtz = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/discamb/mtz/"+mtz_fn,
    test=os.path.isfile)
  if (regression_mtz is None):
    print("Skipping tst_olex2_vs_pydiscamb_TAAM_1: input mtz not available: ", mtz_fn)
    return


  ma = reflection_file_reader.any_reflection_file(file_name =
    regression_mtz).as_miller_arrays()
  fcalc_olex_full = ma[0]
  # mtz has data up to 0.7 A resolution, so cut at 2.0 to speed up computations
  tmp_sel = fcalc_olex_full.d_spacings().data() > 2.0
  fcalc_olex_selected = fcalc_olex_full.select(tmp_sel)
  #print(fcalc_olex_full.d_min())

  # save test model as file
  model_fn = "tst_olex2_vs_pydiscamb_TAAM_1.res"
  with open(model_fn, "w") as f:
    f.write(pdb_str)
  # get xray_structure obj
  xrs = structure.from_shelx(filename=model_fn)
#  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
#  model = mmtbx.model.manager(
#    model_input = pdb_inp,
#    log         = null_out())
#  xrs = model.get_xray_structure()
  # this crashes the wrapper with table=n_gaussian
  xrs.scattering_type_registry(table='wk1995')
  atom_labels = xrs.scatterers().extract_labels()

  w = pydiscamb.DiscambWrapper(
    xrs,
    method=pydiscamb.FCalcMethod.TAAM,
    assignment_info="atom_type_assignment.log")
  # use hkl from olex file
  w.set_indices(fcalc_olex_selected.indices())
  #fcalc = w.f_calc(2.0)
  fcalc_wrapper = w.f_calc()

  fcalc_olex = fcalc_olex_selected.data()
  assert(fcalc_olex.size() == fcalc_wrapper.size())

  score, mean_diff, max_diff = regression_discamb.compare_structure_factors(
    x=fcalc_wrapper, y=fcalc_olex)

  if verbose:
    print('Score     mean diff max diff')
    print(round(score, 7), round(mean_diff, 7), round(max_diff,7))
    print()
  assert (score < target_scores[0])
  assert (mean_diff < target_scores[1])
  assert (max_diff < target_scores[2])
#  print('Score from Viljars ipynb')
#  print(
#    sum(map(lambda a, b: abs(abs(a) - abs(b)), fcalc_olex, fcalc_wrapper))
#    / sum(map(lambda a, b: abs(a) * abs(b), fcalc_wrapper, fcalc_wrapper))
#    )

  ata = w.atom_type_assignment
  # check if assignment dictionary and xrs content is the same
  assert set(atom_labels)==set(ata.keys())
  # make sure no type is unassigned
  assert not has_unassigned(ata)

  os.remove(model_fn)
  os.remove('atom_type_assignment.log')

# ------------------------------------------------------------------------------

pdb_str_3 = '''
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
SCALE1      0.009182  0.000000  0.000000        0.00000
SCALE2      0.000000  0.009182  0.000000        0.00000
SCALE3      0.000000  0.000000  0.009182        0.00000
ATOM      1  N   PRO A  11     -50.099 -12.405  -8.949  1.00 93.57           N
ANISOU    1  N   PRO A  11    15200  10950   9403   2061  -1279    939       N
ATOM      2  CA  PRO A  11     -49.528 -12.967 -10.176  1.00 94.56           C
ANISOU    2  CA  PRO A  11    14433  11580   9914   3081   2193   3834       C
ATOM      3  C   PRO A  11     -50.591 -13.707 -10.969  1.00100.49           C
ANISOU    3  C   PRO A  11    15775  10884  11521   2605   -315   3797       C
ATOM      4  O   PRO A  11     -51.235 -14.637 -10.473  1.00100.63           O
ANISOU    4  O   PRO A  11    14729   9501  14004   3762     76   3742       O
ATOM      5  CB  PRO A  11     -48.440 -13.919  -9.662  1.00 84.07           C
ANISOU    5  CB  PRO A  11    11697  12697   7547   3799   3975   2401       C
ATOM      6  CG  PRO A  11     -48.136 -13.471  -8.321  1.00 90.09           C
ANISOU    6  CG  PRO A  11    12973  10953  10303   2960   1344   1360       C
ATOM      7  CD  PRO A  11     -49.420 -12.928  -7.755  1.00 98.86           C
ANISOU    7  CD  PRO A  11    14002   9616  13946   2347   -699   1417       C
ATOM      8  H2  PRO A  11     -51.030 -12.634  -8.904  1.00 93.57           H
ATOM      9  H3  PRO A  11     -50.005 -11.450  -8.969  1.00 93.57           H
ATOM     10  HA  PRO A  11     -49.104 -12.293 -10.729  1.00 94.56           H
ATOM     11  HB2 PRO A  11     -48.777 -14.829  -9.654  1.00 84.07           H
ATOM     12  HB3 PRO A  11     -47.657 -13.860 -10.231  1.00 84.07           H
ATOM     13  HG2 PRO A  11     -47.458 -12.778  -8.355  1.00 90.09           H
ATOM     14  HG3 PRO A  11     -47.819 -14.220  -7.792  1.00 90.09           H
ATOM     15  HD2 PRO A  11     -49.942 -13.632  -7.338  1.00 98.86           H
ATOM     16  HD3 PRO A  11     -49.246 -12.220  -7.115  1.00 98.86           H
ATOM     17  N   CYS A  12     -50.789 -13.266 -12.204  1.00103.83           N
ANISOU   17  N   CYS A  12    17280   9800  12371    708    451   3961       N
ATOM     18  CA  CYS A  12     -51.831 -13.789 -13.088  1.00116.34           C
ANISOU   18  CA  CYS A  12    17459  12703  14041   1832    700   3834       C
ATOM     19  C   CYS A  12     -51.157 -14.344 -14.339  1.00107.18           C
ANISOU   19  C   CYS A  12    16138  11838  12749    946    545   4232       C
ATOM     20  O   CYS A  12     -51.033 -13.657 -15.353  1.00109.11           O
ANISOU   20  O   CYS A  12    16094  11046  14315  -1900    527   1716       O
ATOM     21  CB  CYS A  12     -52.851 -12.712 -13.426  1.00125.22           C
ANISOU   21  CB  CYS A  12    18286  15938  13352   3474    240   3605       C
ATOM     22  SG  CYS A  12     -54.098 -13.254 -14.597  1.00120.89           S
ANISOU   22  SG  CYS A  12    20261  14589  11084   4721   -996   4480       S
ATOM     23  H   CYS A  12     -50.320 -12.644 -12.569  1.00103.83           H
ATOM     24  HA  CYS A  12     -52.334 -14.498 -12.657  1.00116.34           H
ATOM     25  HB2 CYS A  12     -52.388 -11.953 -13.814  1.00125.22           H
ATOM     26  HB3 CYS A  12     -53.304 -12.443 -12.612  1.00125.22           H
ATOM     27  N   LYS A  13     -50.734 -15.601 -14.263  1.00109.38           N
ANISOU   27  N   LYS A  13    15794  13544  12222   3570    477   4467       N
ATOM     28  CA  LYS A  13     -50.115 -16.290 -15.386  1.00117.72           C
ANISOU   28  CA  LYS A  13    17031  15495  12201   3524  -1585   6315       C
ATOM     29  C   LYS A  13     -50.232 -17.792 -15.161  1.00116.82           C
ANISOU   29  C   LYS A  13    17271  13627  13488   6093  -1719   4334       C
ATOM     30  O   LYS A  13     -50.435 -18.250 -14.034  1.00131.07           O
ANISOU   30  O   LYS A  13    19695  14881  15223   7808  -1127   2418       O
ATOM     31  CB  LYS A  13     -48.645 -15.876 -15.563  1.00120.45           C
ANISOU   31  CB  LYS A  13    16090  15490  14185   1188   -465   7692       C
ATOM     32  CG  LYS A  13     -47.762 -15.990 -14.322  1.00116.42           C
ANISOU   32  CG  LYS A  13    14396  14920  14919   1953    -91   8714       C
ATOM     33  CD  LYS A  13     -46.345 -15.545 -14.642  1.00124.35           C
ANISOU   33  CD  LYS A  13    16425  16041  14783   2857   1187   8695       C
ATOM     34  CE  LYS A  13     -45.378 -15.850 -13.516  1.00130.78           C
ANISOU   34  CE  LYS A  13    18779  16955  13957   2389   3276   8665       C
ATOM     35  NZ  LYS A  13     -43.978 -15.507 -13.891  1.00137.62           N
ANISOU   35  NZ  LYS A  13    20535  17795  13958   5370   2906   4712       N
ATOM     36  H   LYS A  13     -50.796 -16.086 -13.556  1.00109.38           H
ATOM     37  HA  LYS A  13     -50.562 -16.050 -16.213  1.00117.72           H
ATOM     38  HB2 LYS A  13     -48.624 -14.948 -15.845  1.00120.45           H
ATOM     39  HB3 LYS A  13     -48.252 -16.440 -16.247  1.00120.45           H
ATOM     40  HG2 LYS A  13     -48.114 -15.423 -13.618  1.00116.42           H
ATOM     41  HG3 LYS A  13     -47.737 -16.912 -14.023  1.00116.42           H
ATOM     42  HD2 LYS A  13     -46.339 -14.587 -14.793  1.00124.35           H
ATOM     43  HD3 LYS A  13     -46.038 -16.008 -15.437  1.00124.35           H
ATOM     44  HE2 LYS A  13     -45.621 -15.329 -12.734  1.00130.78           H
ATOM     45  HE3 LYS A  13     -45.415 -16.797 -13.309  1.00130.78           H
ATOM     46  HZ1 LYS A  13     -43.918 -14.640 -14.081  1.00137.62           H
ATOM     47  HZ2 LYS A  13     -43.730 -15.977 -14.605  1.00137.62           H
ATOM     48  HZ3 LYS A  13     -43.428 -15.694 -13.217  1.00137.62           H
ATOM     49  OXT LYS A  13     -50.126 -18.580 -16.101  1.00 30.00           O
TER      50      LYS A  13
END
'''

pdb_str_2 = '''
CRYST1   80.435  118.034  112.075  90.00  93.12  90.00 P 1 21 1
SCALE1      0.012432  0.000000  0.000678        0.00000
SCALE2      0.000000  0.008472  0.000000        0.00000
SCALE3      0.000000  0.000000  0.008936        0.00000
ATOM      1  N   ASN A  90      76.017 -20.785 155.731  1.00142.97           N
ANISOU    1  N   ASN A  90    17156  21970  15197  -6703    124    870       N
ATOM      2  CA  ASN A  90      75.942 -19.415 156.227  1.00148.71           C
ANISOU    2  CA  ASN A  90    17577  22797  16129  -6189     29   1252       C
ATOM      3  C   ASN A  90      77.350 -18.831 156.235  1.00143.80           C
ANISOU    3  C   ASN A  90    17170  21823  15645  -5768    103   1134       C
ATOM      4  O   ASN A  90      78.239 -19.359 156.912  1.00146.03           O
ANISOU    4  O   ASN A  90    17802  21588  16095  -5569    283    881       O
ATOM      5  CB  ASN A  90      75.320 -19.390 157.626  1.00158.76           C
ANISOU    5  CB  ASN A  90    18758  23893  17672  -5918    105   1419       C
ATOM      6  CG  ASN A  90      75.068 -17.982 158.153  1.00168.44           C
ANISOU    6  CG  ASN A  90    19640  25239  19120  -5436     46   1805       C
ATOM      7  OD1 ASN A  90      75.922 -17.100 158.061  1.00161.90           O
ANISOU    7  OD1 ASN A  90    18832  24282  18399  -5077     50   1839       O
ATOM      8  ND2 ASN A  90      73.886 -17.776 158.730  1.00187.09           N
ANISOU    8  ND2 ASN A  90    21686  27829  21572  -5430     16   2088       N
ATOM      9  H1  ASN A  90      76.441 -20.787 154.781  1.00142.97           H
ATOM     10  H2  ASN A  90      76.604 -21.355 156.374  1.00142.97           H
ATOM     11  H3  ASN A  90      75.060 -21.190 155.685  1.00142.97           H
ATOM     12  HA  ASN A  90      75.305 -18.797 155.594  1.00148.71           H
ATOM     13  HB2 ASN A  90      75.995 -19.890 158.320  1.00158.76           H
ATOM     14  HB3 ASN A  90      74.363 -19.911 157.595  1.00158.76           H
ATOM     15 HD21 ASN A  90      73.214 -18.543 158.764  1.00187.09           H
ATOM     16  OXT ASN A  90      77.621 -17.830 155.571  1.00 30.00           O
TER      17      ASN A  90
HETATM   16  C1  NAG C   1      73.523 -16.504 159.309  1.00103.52           C
HETATM   17  C2  NAG C   1      72.427 -15.749 158.564  1.00114.27           C
HETATM   18  C3  NAG C   1      72.112 -14.434 159.280  1.00117.99           C
HETATM   19  C4  NAG C   1      71.829 -14.658 160.763  1.00124.42           C
HETATM   20  C5  NAG C   1      72.930 -15.509 161.403  1.00121.57           C
HETATM   21  C6  NAG C   1      72.619 -15.937 162.820  1.00122.81           C
HETATM   22  C7  NAG C   1      71.992 -15.651 156.147  1.00125.03           C
HETATM   23  C8  NAG C   1      72.570 -15.344 154.798  1.00127.15           C
HETATM   24  N2  NAG C   1      72.819 -15.497 157.186  1.00121.51           N
HETATM   25  O3  NAG C   1      70.990 -13.814 158.662  1.00116.89           O
HETATM   26  O4  NAG C   1      71.791 -13.386 161.404  1.00133.35           O
HETATM   27  O5  NAG C   1      73.123 -16.717 160.654  1.00114.93           O
HETATM   28  O6  NAG C   1      73.745 -15.787 163.673  1.00121.65           O
HETATM   29  O7  NAG C   1      70.831 -16.020 156.289  1.00125.47           O
HETATM   30  H1  NAG C   1      74.408 -15.867 159.309  1.00103.52           H
HETATM   31  H2  NAG C   1      71.534 -16.374 158.560  1.00114.27           H
HETATM   32  H3  NAG C   1      73.013 -13.828 159.184  1.00117.99           H
HETATM   33  H4  NAG C   1      70.874 -15.182 160.807  1.00124.42           H
HETATM   34  H5  NAG C   1      73.836 -14.903 161.407  1.00121.57           H
HETATM   35  H61 NAG C   1      71.774 -15.344 163.169  1.00122.81           H
HETATM   36  H62 NAG C   1      72.279 -16.972 162.788  1.00122.81           H
HETATM   37  H81 NAG C   1      73.010 -16.249 154.379  1.00127.15           H
HETATM   38  H82 NAG C   1      73.338 -14.578 154.901  1.00127.15           H
HETATM   39  H83 NAG C   1      71.779 -14.983 154.140  1.00127.15           H
HETATM   40  HN2 NAG C   1      73.775 -15.187 157.013  1.00121.51           H
HETATM   41  HO3 NAG C   1      70.461 -14.515 158.228  1.00116.89           H
HETATM   42  HO6 NAG C   1      73.838 -14.837 163.894  1.00121.65           H
HETATM   43  C1  NAG C   2      70.691 -13.209 162.330  1.00143.22           C
HETATM   44  C2  NAG C   2      71.046 -12.013 163.217  1.00149.90           C
HETATM   45  C3  NAG C   2      69.929 -11.751 164.223  1.00153.90           C
HETATM   46  C4  NAG C   2      68.594 -11.600 163.506  1.00157.09           C
HETATM   47  C5  NAG C   2      68.339 -12.806 162.605  1.00152.27           C
HETATM   48  C6  NAG C   2      67.082 -12.672 161.777  1.00151.75           C
HETATM   49  C7  NAG C   2      73.401 -11.485 163.670  1.00150.52           C
HETATM   50  C8  NAG C   2      74.621 -11.843 164.465  1.00149.56           C
HETATM   51  N2  NAG C   2      72.312 -12.226 163.901  1.00151.39           N
HETATM   52  O3  NAG C   2      70.226 -10.571 164.961  1.00153.21           O
HETATM   53  O4  NAG C   2      67.540 -11.490 164.456  1.00163.16           O
HETATM   54  O5  NAG C   2      69.429 -12.972 161.684  1.00147.83           O
HETATM   55  O6  NAG C   2      67.313 -11.906 160.602  1.00151.36           O
HETATM   56  O7  NAG C   2      73.402 -10.566 162.856  1.00149.17           O
HETATM   57  H1  NAG C   2      70.612 -14.091 162.965  1.00143.22           H
HETATM   58  H2  NAG C   2      71.151 -11.142 162.570  1.00149.90           H
HETATM   59  H3  NAG C   2      69.889 -12.624 164.874  1.00153.90           H
HETATM   60  H4  NAG C   2      68.695 -10.687 162.919  1.00157.09           H
HETATM   61  H5  NAG C   2      68.245 -13.683 163.246  1.00152.27           H
HETATM   62  H61 NAG C   2      66.315 -12.218 162.405  1.00151.75           H
HETATM   63  H62 NAG C   2      66.732 -13.676 161.537  1.00151.75           H
HETATM   64  H81 NAG C   2      74.625 -11.274 165.395  1.00149.56           H
HETATM   65  H82 NAG C   2      75.512 -11.604 163.885  1.00149.56           H
HETATM   66  H83 NAG C   2      74.605 -12.910 164.689  1.00149.56           H
HETATM   67  HN2 NAG C   2      72.362 -12.979 164.588  1.00151.39           H
HETATM   68  HO3 NAG C   2      69.769  -9.821 164.526  1.00153.21           H
HETATM   69  HO6 NAG C   2      67.822 -11.107 160.852  1.00151.36           H
HETATM   70  C1  BMA C   3      67.036 -10.135 164.465  1.00164.84           C
HETATM   71  C2  BMA C   3      65.513 -10.187 164.205  1.00164.51           C
HETATM   72  C3  BMA C   3      64.912  -8.784 164.326  1.00164.84           C
HETATM   73  C4  BMA C   3      65.353  -8.085 165.626  1.00165.75           C
HETATM   74  C5  BMA C   3      66.886  -8.117 165.760  1.00164.34           C
HETATM   75  C6  BMA C   3      67.377  -7.517 167.067  1.00159.95           C
HETATM   76  O2  BMA C   3      64.864 -10.998 165.174  1.00162.15           O
HETATM   77  O3  BMA C   3      63.491  -8.819 164.251  1.00162.71           O
HETATM   78  O4  BMA C   3      64.902  -6.737 165.632  1.00164.98           O
HETATM   79  O5  BMA C   3      67.321  -9.489 165.704  1.00164.62           O
HETATM   80  O6  BMA C   3      66.718  -6.271 167.265  1.00154.88           O
HETATM   81  H1  BMA C   3      67.496  -9.576 163.650  1.00164.84           H
HETATM   82  H2  BMA C   3      65.440 -10.579 163.191  1.00164.51           H
HETATM   83  H3  BMA C   3      65.303  -8.259 163.454  1.00164.84           H
HETATM   84  H4  BMA C   3      64.860  -8.591 166.456  1.00165.75           H
HETATM   85  H5  BMA C   3      67.348  -7.526 164.969  1.00164.34           H
HETATM   86  H61 BMA C   3      67.168  -8.200 167.891  1.00159.95           H
HETATM   87  H62 BMA C   3      68.458  -7.381 167.030  1.00159.95           H
HETATM   88  HO2 BMA C   3      65.215 -10.745 166.053  1.00162.15           H
HETATM   89  HO3 BMA C   3      63.149  -9.448 164.920  1.00162.71           H
HETATM   90  HO4 BMA C   3      64.777  -6.452 166.561  1.00164.98           H
HETATM   91  HO6 BMA C   3      66.477  -5.909 166.387  1.00154.88           H
END
'''

pdb_str_1 = '''
CRYST1   39.153   32.696   29.501  90.00  90.00  90.00 P 1
SCALE1      0.025541  0.000000  0.000000        0.00000
SCALE2      0.000000  0.030585  0.000000        0.00000
SCALE3      0.000000  0.000000  0.033897        0.00000
ATOM      1  N   ASP A   8      15.273  21.826   9.207  1.00 14.56           N
ATOM      2  CA  ASP A   8      14.438  22.654  10.069  1.00 16.10           C
ATOM      3  C   ASP A   8      13.125  21.951  10.395  1.00 15.76           C
ATOM      4  O   ASP A   8      12.396  21.528   9.498  1.00 15.20           O
ATOM      5  CB  ASP A   8      14.163  24.007   9.409  1.00 18.35           C
ATOM      6  CG  ASP A   8      13.956  25.117  10.421  1.00 22.97           C
ATOM      7  OD1 ASP A   8      13.547  24.817  11.562  1.00 22.44           O
ATOM      8  OD2 ASP A   8      14.202  26.292  10.075  1.00 27.92           O
ATOM      9  H1  ASP A   8      15.450  20.912   9.670  1.00 14.56           H
ATOM     10  H2  ASP A   8      14.787  21.668   8.302  1.00 14.56           H
ATOM     11  H3  ASP A   8      16.178  22.308   9.035  1.00 14.56           H
ATOM     12  HA  ASP A   8      14.964  22.837  11.006  1.00 16.10           H
ATOM     13  HB2 ASP A   8      15.013  24.276   8.781  1.00 18.35           H
ATOM     14  HB3 ASP A   8      13.261  23.930   8.802  1.00 18.35           H
ATOM     15  N   SER A   9      12.829  21.829  11.688  1.00 16.08           N
ATOM     16  CA  SER A   9      11.611  21.183  12.154  1.00 17.80           C
ATOM     17  C   SER A   9      10.503  22.173  12.488  1.00 19.90           C
ATOM     18  O   SER A   9       9.436  21.755  12.949  1.00 19.71           O
ATOM     19  CB  SER A   9      11.911  20.314  13.380  1.00 16.95           C
ATOM     20  OG  SER A   9      12.676  19.175  13.023  1.00 17.71           O
ATOM     21  H   SER A   9      13.422  22.173  12.443  1.00 16.08           H
ATOM     22  HA  SER A   9      11.239  20.553  11.346  1.00 17.80           H
ATOM     23  HB2 SER A   9      12.473  20.905  14.103  1.00 16.95           H
ATOM     24  HB3 SER A   9      10.969  19.984  13.819  1.00 16.95           H
ATOM     25  HG  SER A   9      12.572  18.475  13.701  1.00 17.71           H
ATOM     26  N   ASN A  10      10.728  23.470  12.269  1.00 20.66           N
ATOM     27  CA  ASN A  10       9.704  24.463  12.570  1.00 22.74           C
ATOM     28  C   ASN A  10       8.677  24.586  11.451  1.00 22.30           C
ATOM     29  O   ASN A  10       7.532  24.971  11.708  1.00 24.46           O
ATOM     30  CB  ASN A  10      10.353  25.823  12.836  1.00 27.34           C
ATOM     31  CG  ASN A  10      11.152  25.845  14.124  1.00 31.14           C
ATOM     32  OD1 ASN A  10      10.662  25.443  15.179  1.00 34.86           O
ATOM     33  ND2 ASN A  10      12.391  26.315  14.043  1.00 32.87           N
ATOM     34  H   ASN A  10      11.594  23.854  11.891  1.00 20.66           H
ATOM     35  HA  ASN A  10       9.169  24.147  13.465  1.00 22.74           H
ATOM     36  HB2 ASN A  10       9.573  26.581  12.909  1.00 27.34           H
ATOM     37  HB3 ASN A  10      11.029  26.061  12.015  1.00 27.34           H
ATOM     38 HD21 ASN A  10      12.977  26.354  14.877  1.00 32.87           H
ATOM     39 HD22 ASN A  10      12.757  26.637  13.147  1.00 32.87           H
ATOM     40  N   ILE A  11       9.061  24.266  10.217  1.00 20.31           N
ATOM     41  CA  ILE A  11       8.152  24.348   9.079  1.00 18.78           C
ATOM     42  C   ILE A  11       7.553  22.974   8.815  1.00 18.68           C
ATOM     43  O   ILE A  11       6.328  22.811   8.791  1.00 18.26           O
ATOM     44  CB  ILE A  11       8.869  24.888   7.828  1.00 18.40           C
ATOM     45  CG1 ILE A  11       9.332  26.328   8.059  1.00 18.41           C
ATOM     46  CG2 ILE A  11       7.957  24.805   6.614  1.00 17.87           C
ATOM     47  CD1 ILE A  11      10.420  26.779   7.110  1.00 20.89           C
ATOM     48  H   ILE A  11       9.998  23.946   9.973  1.00 20.31           H
ATOM     49  HA  ILE A  11       7.346  25.043   9.313  1.00 18.78           H
ATOM     50  HB  ILE A  11       9.747  24.270   7.641  1.00 18.40           H
ATOM     51 HG12 ILE A  11       9.720  26.413   9.074  1.00 18.41           H
ATOM     52 HG13 ILE A  11       8.480  26.995   7.929  1.00 18.41           H
ATOM     53 HG21 ILE A  11       7.063  25.401   6.800  1.00 17.87           H
ATOM     54 HG22 ILE A  11       7.681  23.764   6.448  1.00 17.87           H
ATOM     55 HG23 ILE A  11       8.488  25.193   5.744  1.00 17.87           H
ATOM     56 HD11 ILE A  11      10.048  26.714   6.087  1.00 20.89           H
ATOM     57 HD12 ILE A  11      11.288  26.131   7.232  1.00 20.89           H
ATOM     58 HD13 ILE A  11      10.689  27.809   7.343  1.00 20.89           H
ATOM     59  N   HIS A  12       8.412  21.977   8.614  1.00 17.15           N
ATOM     60  CA  HIS A  12       7.991  20.609   8.351  1.00 15.83           C
ATOM     61  C   HIS A  12       8.585  19.690   9.407  1.00 16.05           C
ATOM     62  O   HIS A  12       9.796  19.723   9.654  1.00 15.95           O
ATOM     63  CB  HIS A  12       8.420  20.159   6.951  1.00 14.81           C
ATOM     64  CG  HIS A  12       7.831  18.849   6.530  1.00 16.05           C
ATOM     65  ND1 HIS A  12       6.721  18.758   5.718  1.00 17.40           N
ATOM     66  CD2 HIS A  12       8.198  17.575   6.808  1.00 16.62           C
ATOM     67  CE1 HIS A  12       6.429  17.486   5.514  1.00 17.35           C
ATOM     68  NE2 HIS A  12       7.310  16.747   6.164  1.00 19.34           N
ATOM     69  H   HIS A  12       9.426  22.090   8.628  1.00 17.15           H
ATOM     70  HA  HIS A  12       6.904  20.543   8.385  1.00 15.83           H
ATOM     71  HB2 HIS A  12       9.505  20.057   6.933  1.00 14.81           H
ATOM     72  HB3 HIS A  12       8.105  20.912   6.229  1.00 14.81           H
ATOM     73  HD1 HIS A  12       6.205  19.550   5.334  1.00 17.40           H
ATOM     74  HD2 HIS A  12       9.032  17.267   7.421  1.00 16.62           H
ATOM     75  HE1 HIS A  12       5.609  17.114   4.917  1.00 17.35           H
ATOM     76  HE2 HIS A  12       7.328  15.727   6.185  1.00 19.34           H
ATOM     77  N   LYS A  13       7.735  18.874  10.026  1.00 15.54           N
ATOM     78  CA  LYS A  13       8.178  17.945  11.058  1.00 16.13           C
ATOM     79  C   LYS A  13       8.782  16.708  10.405  1.00 14.55           C
ATOM     80  O   LYS A  13       8.110  16.013   9.635  1.00 14.66           O
ATOM     81  CB  LYS A  13       7.011  17.563  11.966  1.00 19.94           C
ATOM     82  CG  LYS A  13       7.415  16.752  13.187  0.50 26.76           C
ATOM     83  CD  LYS A  13       6.228  16.509  14.105  0.50 30.57           C
ATOM     84  CE  LYS A  13       6.678  16.255  15.535  0.50 34.22           C
ATOM     85  NZ  LYS A  13       5.524  16.138  16.468  0.50 36.77           N
ATOM     86  H   LYS A  13       6.734  18.834   9.834  1.00 15.54           H
ATOM     87  HA  LYS A  13       8.936  18.411  11.687  1.00 16.13           H
ATOM     88  HB2 LYS A  13       6.530  18.475  12.319  1.00 19.94           H
ATOM     89  HB3 LYS A  13       6.303  16.966  11.391  1.00 19.94           H
ATOM     90  HG2 LYS A  13       8.177  17.295  13.747  0.50 26.76           H
ATOM     91  HG3 LYS A  13       7.805  15.786  12.867  0.50 26.76           H
ATOM     92  HD2 LYS A  13       5.581  17.386  14.099  0.50 30.57           H
ATOM     93  HD3 LYS A  13       5.676  15.636  13.759  0.50 30.57           H
ATOM     94  HE2 LYS A  13       7.302  17.084  15.868  0.50 34.22           H
ATOM     95  HE3 LYS A  13       7.243  15.324  15.574  0.50 34.22           H
ATOM     96  HZ1 LYS A  13       4.912  15.348  16.179  0.50 36.77           H
ATOM     97  HZ2 LYS A  13       4.968  17.017  16.458  0.50 36.77           H
ATOM     98  HZ3 LYS A  13       5.865  15.965  17.435  0.50 36.77           H
ATOM     99  N   CYS A  14      10.047  16.434  10.711  1.00 12.28           N
ATOM    100  CA  CYS A  14      10.764  15.291  10.166  1.00 11.10           C
ATOM    101  C   CYS A  14      11.140  14.338  11.291  1.00 10.09           C
ATOM    102  O   CYS A  14      11.423  14.768  12.414  1.00 10.97           O
ATOM    103  CB  CYS A  14      12.022  15.734   9.411  1.00 10.22           C
ATOM    104  SG  CYS A  14      12.793  14.437   8.414  1.00 11.10           S
ATOM    105  H   CYS A  14      10.612  16.997  11.346  1.00 12.28           H
ATOM    106  HA  CYS A  14      10.126  14.771   9.451  1.00 11.10           H
ATOM    107  HB2 CYS A  14      12.761  16.075  10.137  1.00 10.22           H
ATOM    108  HB3 CYS A  14      11.755  16.550   8.740  1.00 10.22           H
ATOM    109  N   GLY A  15      11.141  13.043  10.984  1.00 10.65           N
ATOM    110  CA  GLY A  15      11.473  12.029  11.954  1.00 11.27           C
ATOM    111  C   GLY A  15      12.968  11.858  12.127  1.00 10.83           C
ATOM    112  O   GLY A  15      13.561  12.321  13.106  1.00 11.43           O
ATOM    113  H   GLY A  15      10.913  12.672  10.062  1.00 10.65           H
ATOM    114  HA2 GLY A  15      11.045  12.300  12.919  1.00 11.27           H
ATOM    115  HA3 GLY A  15      11.053  11.075  11.635  1.00 11.27           H
ATOM    116  N   PRO A  16      13.611  11.180  11.170  1.00 10.74           N
ATOM    117  CA  PRO A  16      15.065  10.968  11.264  1.00  9.60           C
ATOM    118  C   PRO A  16      15.886  12.240  11.125  1.00  8.87           C
ATOM    119  O   PRO A  16      17.095  12.201  11.389  1.00  8.86           O
ATOM    120  CB  PRO A  16      15.357   9.997  10.110  1.00 10.32           C
ATOM    121  CG  PRO A  16      14.025   9.413   9.734  1.00 10.84           C
ATOM    122  CD  PRO A  16      13.025  10.487  10.012  1.00 11.34           C
ATOM    123  HA  PRO A  16      15.312  10.499  12.216  1.00  9.60           H
ATOM    124  HB2 PRO A  16      15.794  10.540   9.272  1.00 10.32           H
ATOM    125  HB3 PRO A  16      16.041   9.218  10.447  1.00 10.32           H
ATOM    126  HG2 PRO A  16      13.824   8.530  10.341  1.00 10.84           H
ATOM    127  HG3 PRO A  16      14.024   9.147   8.677  1.00 10.84           H
ATOM    128  HD2 PRO A  16      12.929  11.161   9.161  1.00 11.34           H
ATOM    129  HD3 PRO A  16      12.053  10.061  10.262  1.00 11.34           H
ATOM    130  N   CYS A  17      15.278  13.358  10.721  1.00  8.07           N
ATOM    131  CA  CYS A  17      16.035  14.597  10.576  1.00  8.42           C
ATOM    132  C   CYS A  17      16.450  15.167  11.927  1.00  9.31           C
ATOM    133  O   CYS A  17      17.458  15.877  12.013  1.00 10.92           O
ATOM    134  CB  CYS A  17      15.215  15.625   9.797  1.00  8.17           C
ATOM    135  SG  CYS A  17      14.702  15.079   8.152  1.00 10.96           S
ATOM    136  H   CYS A  17      14.287  13.433  10.493  1.00  8.07           H
ATOM    137  HA  CYS A  17      16.948  14.388  10.019  1.00  8.42           H
ATOM    138  HB2 CYS A  17      15.815  16.526   9.673  1.00  8.17           H
ATOM    139  HB3 CYS A  17      14.313  15.853  10.365  1.00  8.17           H
ATOM    140  N   GLU A  18      15.695  14.869  12.986  1.00  9.20           N
ATOM    141  CA  GLU A  18      16.038  15.375  14.309  1.00 10.73           C
ATOM    142  C   GLU A  18      17.095  14.519  14.994  1.00  9.89           C
ATOM    143  O   GLU A  18      17.810  15.013  15.874  1.00 11.78           O
ATOM    144  CB  GLU A  18      14.785  15.457  15.183  1.00 30.00           C
ATOM    145  CG  GLU A  18      13.812  16.550  14.772  1.00 30.00           C
ATOM    146  CD  GLU A  18      12.571  16.582  15.642  1.00 30.00           C
ATOM    147  OE1 GLU A  18      12.356  15.620  16.409  1.00 30.00           O
ATOM    148  OE2 GLU A  18      11.811  17.570  15.560  1.00 30.00           O
ATOM    149  H   GLU A  18      14.856  14.290  12.957  1.00  9.20           H
ATOM    150  HA  GLU A  18      16.455  16.376  14.198  1.00 10.73           H
ATOM    151  HB2 GLU A  18      15.089  15.654  16.211  1.00 30.00           H
ATOM    152  HB3 GLU A  18      14.257  14.505  15.126  1.00 30.00           H
ATOM    153  HG2 GLU A  18      14.308  17.517  14.855  1.00 30.00           H
ATOM    154  HG3 GLU A  18      13.498  16.379  13.742  1.00 30.00           H
ATOM    155  N   GLN A  19      17.211  13.248  14.611  1.00  9.22           N
ATOM    156  CA  GLN A  19      18.196  12.358  15.211  1.00 10.53           C
ATOM    157  C   GLN A  19      19.540  12.399  14.498  1.00 10.38           C
ATOM    158  O   GLN A  19      20.561  12.045  15.098  1.00 11.57           O
ATOM    159  CB  GLN A  19      17.670  10.920  15.224  1.00 30.00           C
ATOM    160  CG  GLN A  19      16.457  10.711  16.116  1.00 30.00           C
ATOM    161  CD  GLN A  19      15.899   9.304  16.021  1.00 30.00           C
ATOM    162  OE1 GLN A  19      16.084   8.618  15.016  1.00 30.00           O
ATOM    163  NE2 GLN A  19      15.212   8.868  17.070  1.00 30.00           N
ATOM    164  H   GLN A  19      16.638  12.809  13.890  1.00  9.22           H
ATOM    165  HA  GLN A  19      18.365  12.691  16.235  1.00 10.53           H
ATOM    166  HB2 GLN A  19      18.462  10.262  15.582  1.00 30.00           H
ATOM    167  HB3 GLN A  19      17.386  10.643  14.209  1.00 30.00           H
ATOM    168  HG2 GLN A  19      16.742  10.892  17.152  1.00 30.00           H
ATOM    169  HG3 GLN A  19      15.672  11.406  15.818  1.00 30.00           H
ATOM    170 HE21 GLN A  19      14.813   7.929  17.064  1.00 30.00           H
ATOM    171 HE22 GLN A  19      15.084   9.472  17.882  1.00 30.00           H
ATOM    172  N   ALA A  20      19.564  12.822  13.232  1.00  9.54           N
ATOM    173  CA  ALA A  20      20.819  12.886  12.493  1.00  9.47           C
ATOM    174  C   ALA A  20      21.645  14.109  12.870  1.00 10.34           C
ATOM    175  O   ALA A  20      22.868  14.105  12.685  1.00 10.59           O
ATOM    176  CB  ALA A  20      20.546  12.881  10.989  1.00  9.78           C
ATOM    177  H   ALA A  20      18.744  13.121  12.704  1.00  9.54           H
ATOM    178  HA  ALA A  20      21.407  12.004  12.747  1.00  9.47           H
ATOM    179  HB1 ALA A  20      19.934  13.747  10.737  1.00  9.78           H
ATOM    180  HB2 ALA A  20      20.019  11.964  10.727  1.00  9.78           H
ATOM    181  HB3 ALA A  20      21.495  12.929  10.456  1.00  9.78           H
ATOM    182  N   MET A  21      21.006  15.156  13.393  1.00 11.01           N
ATOM    183  CA  MET A  21      21.736  16.358  13.777  1.00 12.37           C
ATOM    184  C   MET A  21      22.379  16.225  15.151  1.00 11.41           C
ATOM    185  O   MET A  21      23.430  16.827  15.398  1.00 12.52           O
ATOM    186  CB  MET A  21      20.804  17.571  13.750  1.00 30.00           C
ATOM    187  CG  MET A  21      20.295  17.932  12.364  1.00 30.00           C
ATOM    188  SD  MET A  21      21.626  18.312  11.210  1.00 30.00           S
ATOM    189  CE  MET A  21      22.285  19.812  11.935  2.00 30.00           C
ATOM    190  H   MET A  21      20.001  15.199  13.559  1.00 11.01           H
ATOM    191  HA  MET A  21      22.532  16.518  13.050  1.00 12.37           H
ATOM    192  HB2 MET A  21      21.343  18.434  14.140  1.00 30.00           H
ATOM    193  HB3 MET A  21      19.937  17.359  14.376  1.00 30.00           H
ATOM    194  HG2 MET A  21      19.653  18.810  12.439  1.00 30.00           H
ATOM    195  HG3 MET A  21      19.731  17.090  11.962  1.00 30.00           H
ATOM    196  HE1 MET A  21      22.629  19.596  12.946  2.00 30.00           H
ATOM    197  HE2 MET A  21      21.500  20.567  11.963  2.00 30.00           H
ATOM    198  HE3 MET A  21      23.118  20.164  11.326  2.00 30.00           H
ATOM    199  N   ARG A  22      21.774  15.451  16.048  1.00 12.13           N
ATOM    200  CA  ARG A  22      22.302  15.253  17.390  1.00 13.74           C
ATOM    201  C   ARG A  22      23.261  14.073  17.484  1.00 13.22           C
ATOM    202  O   ARG A  22      23.793  13.812  18.568  1.00 12.80           O
ATOM    203  CB  ARG A  22      21.153  15.061  18.386  1.00 17.56           C
ATOM    204  CG  ARG A  22      20.591  16.358  18.944  1.00 21.02           C
ATOM    205  CD  ARG A  22      19.078  16.408  18.812  1.00 25.74           C
ATOM    206  NE  ARG A  22      18.562  17.764  18.953  1.00 27.54           N
ATOM    207  CZ  ARG A  22      17.273  18.077  18.972  1.00 28.60           C
ATOM    208  NH1 ARG A  22      16.335  17.151  18.861  1.00 26.59           N
ATOM    209  NH2 ARG A  22      16.917  19.352  19.106  1.00 31.75           N
ATOM    210  H   ARG A  22      20.907  14.943  15.872  1.00 12.13           H
ATOM    211  HA  ARG A  22      22.869  16.144  17.658  1.00 13.74           H
ATOM    212  HB2 ARG A  22      20.340  14.537  17.884  1.00 17.56           H
ATOM    213  HB3 ARG A  22      21.515  14.469  19.226  1.00 17.56           H
ATOM    214  HG2 ARG A  22      21.011  17.200  18.395  1.00 21.02           H
ATOM    215  HG3 ARG A  22      20.846  16.438  20.001  1.00 21.02           H
ATOM    216  HD2 ARG A  22      18.791  16.035  17.829  1.00 25.74           H
ATOM    217  HD3 ARG A  22      18.629  15.790  19.590  1.00 25.74           H
ATOM    218  HE  ARG A  22      19.237  18.523  19.043  1.00 27.54           H
ATOM    219 HH11 ARG A  22      15.351  17.418  18.878  1.00 26.59           H
ATOM    220 HH12 ARG A  22      16.596  16.170  18.758  1.00 26.59           H
ATOM    221 HH21 ARG A  22      15.929  19.607  19.122  1.00 31.75           H
ATOM    222 HH22 ARG A  22      17.632  20.075  19.193  1.00 31.75           H
ATOM    223  N   LEU A  23      23.493  13.358  16.382  1.00 11.70           N
ATOM    224  CA  LEU A  23      24.403  12.218  16.415  1.00 11.05           C
ATOM    225  C   LEU A  23      25.860  12.662  16.344  1.00 10.67           C
ATOM    226  O   LEU A  23      26.717  12.095  17.030  1.00 10.75           O
ATOM    227  CB  LEU A  23      24.082  11.256  15.271  1.00 11.60           C
ATOM    228  CG  LEU A  23      24.910   9.971  15.207  1.00 14.28           C
ATOM    229  CD1 LEU A  23      24.776   9.181  16.500  1.00 15.73           C
ATOM    230  CD2 LEU A  23      24.498   9.125  14.012  1.00 13.45           C
ATOM    231  H   LEU A  23      23.073  13.542  15.471  1.00 11.70           H
ATOM    232  HA  LEU A  23      24.272  11.686  17.357  1.00 11.05           H
ATOM    233  HB2 LEU A  23      24.239  11.785  14.331  1.00 11.60           H
ATOM    234  HB3 LEU A  23      23.037  10.960  15.362  1.00 11.60           H
ATOM    235  HG  LEU A  23      25.960  10.234  15.083  1.00 14.28           H
ATOM    236 HD11 LEU A  23      23.728   8.923  16.651  1.00 15.73           H
ATOM    237 HD12 LEU A  23      25.131   9.794  17.329  1.00 15.73           H
ATOM    238 HD13 LEU A  23      25.375   8.273  16.425  1.00 15.73           H
ATOM    239 HD21 LEU A  23      24.658   9.698  13.099  1.00 13.45           H
ATOM    240 HD22 LEU A  23      23.444   8.866  14.108  1.00 13.45           H
ATOM    241 HD23 LEU A  23      25.104   8.219  13.993  1.00 13.45           H
ATOM    242  N   PHE A  24      26.156  13.671  15.525  1.00  9.90           N
ATOM    243  CA  PHE A  24      27.519  14.166  15.386  1.00  9.78           C
ATOM    244  C   PHE A  24      27.964  15.025  16.562  1.00 10.19           C
ATOM    245  O   PHE A  24      29.165  15.287  16.695  1.00 11.78           O
ATOM    246  CB  PHE A  24      27.658  14.965  14.088  1.00 30.00           C
ATOM    247  CG  PHE A  24      27.431  14.150  12.847  1.00 30.00           C
ATOM    248  CD1 PHE A  24      28.452  13.385  12.308  1.00 30.00           C
ATOM    249  CD2 PHE A  24      26.196  14.149  12.218  1.00 30.00           C
ATOM    250  CE1 PHE A  24      28.246  12.633  11.166  1.00 30.00           C
ATOM    251  CE2 PHE A  24      25.984  13.400  11.076  1.00 30.00           C
ATOM    252  CZ  PHE A  24      27.010  12.641  10.549  1.00 30.00           C
ATOM    253  H   PHE A  24      25.474  14.163  14.948  1.00  9.90           H
ATOM    254  HA  PHE A  24      28.184  13.303  15.355  1.00  9.78           H
ATOM    255  HB2 PHE A  24      28.666  15.376  14.036  1.00 30.00           H
ATOM    256  HB3 PHE A  24      26.926  15.773  14.094  1.00 30.00           H
ATOM    257  HD1 PHE A  24      29.421  13.376  12.785  1.00 30.00           H
ATOM    258  HD2 PHE A  24      25.389  14.741  12.625  1.00 30.00           H
ATOM    259  HE1 PHE A  24      29.051  12.040  10.757  1.00 30.00           H
ATOM    260  HE2 PHE A  24      25.016  13.408  10.596  1.00 30.00           H
ATOM    261  HZ  PHE A  24      26.847  12.055   9.657  1.00 30.00           H
ATOM    262  N   THR A  25      27.036  15.469  17.411  1.00 11.80           N
ATOM    263  CA  THR A  25      27.409  16.293  18.556  1.00 13.89           C
ATOM    264  C   THR A  25      28.067  15.464  19.653  1.00 14.79           C
ATOM    265  O   THR A  25      29.052  15.902  20.258  1.00 15.74           O
ATOM    266  CB  THR A  25      26.179  17.019  19.104  1.00 14.78           C
ATOM    267  OG1 THR A  25      25.522  17.720  18.041  1.00 18.27           O
ATOM    268  CG2 THR A  25      26.581  18.012  20.185  1.00 14.92           C
ATOM    269  H   THR A  25      26.037  15.277  17.332  1.00 11.80           H
ATOM    270  HA  THR A  25      28.144  17.040  18.255  1.00 13.89           H
ATOM    271  HB  THR A  25      25.494  16.292  19.541  1.00 14.78           H
ATOM    272  HG1 THR A  25      24.934  17.109  17.550  1.00 18.27           H
ATOM    273 HG21 THR A  25      27.268  18.752  19.775  1.00 14.92           H
ATOM    274 HG22 THR A  25      27.072  17.490  21.006  1.00 14.92           H
ATOM    275 HG23 THR A  25      25.698  18.523  20.568  1.00 14.92           H
ATOM    276  N   VAL A  26      27.541  14.267  19.918  1.00 15.02           N
ATOM    277  CA  VAL A  26      28.115  13.420  20.960  1.00 16.72           C
ATOM    278  C   VAL A  26      29.456  12.852  20.511  1.00 16.09           C
ATOM    279  O   VAL A  26      30.409  12.779  21.296  1.00 16.42           O
ATOM    280  CB  VAL A  26      27.125  12.305  21.343  1.00 18.10           C
ATOM    281  CG1 VAL A  26      27.687  11.460  22.477  1.00 20.04           C
ATOM    282  CG2 VAL A  26      25.780  12.900  21.730  1.00 18.95           C
ATOM    283  H   VAL A  26      26.735  13.866  19.438  1.00 15.02           H
ATOM    284  HA  VAL A  26      28.306  14.019  21.851  1.00 16.72           H
ATOM    285  HB  VAL A  26      26.973  11.655  20.481  1.00 18.10           H
ATOM    286 HG11 VAL A  26      27.860  12.099  23.343  1.00 20.04           H
ATOM    287 HG12 VAL A  26      28.626  11.010  22.153  1.00 20.04           H
ATOM    288 HG13 VAL A  26      26.968  10.680  22.728  1.00 20.04           H
ATOM    289 HG21 VAL A  26      25.382  13.458  20.883  1.00 18.95           H
ATOM    290 HG22 VAL A  26      25.918  13.566  22.582  1.00 18.95           H
ATOM    291 HG23 VAL A  26      25.098  12.092  21.996  1.00 18.95           H
ATOM    292  N   TRP A  27      29.555  12.448  19.246  1.00 15.34           N
ATOM    293  CA  TRP A  27      30.787  11.882  18.712  1.00 15.59           C
ATOM    294  C   TRP A  27      31.824  12.939  18.353  1.00 16.30           C
ATOM    295  O   TRP A  27      32.947  12.575  17.986  1.00 16.26           O
ATOM    296  CB  TRP A  27      30.480  11.025  17.480  1.00 30.00           C
ATOM    297  CG  TRP A  27      29.784   9.738  17.803  1.00 30.00           C
ATOM    298  CD1 TRP A  27      28.439   9.511  17.802  1.00 30.00           C
ATOM    299  CD2 TRP A  27      30.400   8.499  18.174  1.00 30.00           C
ATOM    300  NE1 TRP A  27      28.179   8.208  18.150  1.00 30.00           N
ATOM    301  CE2 TRP A  27      29.366   7.565  18.383  1.00 30.00           C
ATOM    302  CE3 TRP A  27      31.725   8.088  18.349  1.00 30.00           C
ATOM    303  CZ2 TRP A  27      29.615   6.246  18.759  1.00 30.00           C
ATOM    304  CZ3 TRP A  27      31.971   6.779  18.722  1.00 30.00           C
ATOM    305  CH2 TRP A  27      30.921   5.874  18.923  1.00 30.00           C
ATOM    306  H   TRP A  27      28.796  12.501  18.566  1.00 15.34           H
ATOM    307  HA  TRP A  27      31.229  11.258  19.489  1.00 15.59           H
ATOM    308  HB2 TRP A  27      31.417  10.781  16.980  1.00 30.00           H
ATOM    309  HB3 TRP A  27      29.836  11.593  16.809  1.00 30.00           H
ATOM    310  HD1 TRP A  27      27.688  10.249  17.562  1.00 30.00           H
ATOM    311  HE1 TRP A  27      27.253   7.787  18.223  1.00 30.00           H
ATOM    312  HE3 TRP A  27      32.541   8.779  18.196  1.00 30.00           H
ATOM    313  HZ2 TRP A  27      28.808   5.545  18.915  1.00 30.00           H
ATOM    314  HZ3 TRP A  27      32.990   6.449  18.861  1.00 30.00           H
ATOM    315  HH2 TRP A  27      31.148   4.859  19.214  1.00 30.00           H
ATOM    316  N   TYR A  28      31.477  14.224  18.450  1.00 15.46           N
ATOM    317  CA  TYR A  28      32.376  15.336  18.135  1.00 17.48           C
ATOM    318  C   TYR A  28      32.897  15.222  16.698  1.00 17.25           C
ATOM    319  O   TYR A  28      34.037  14.833  16.437  1.00 16.98           O
ATOM    320  CB  TYR A  28      33.531  15.409  19.141  1.00 30.00           C
ATOM    321  CG  TYR A  28      33.098  15.747  20.550  1.00 30.00           C
ATOM    322  CD1 TYR A  28      32.852  17.062  20.925  1.00 30.00           C
ATOM    323  CD2 TYR A  28      32.938  14.753  21.506  1.00 30.00           C
ATOM    324  CE1 TYR A  28      32.457  17.376  22.212  1.00 30.00           C
ATOM    325  CE2 TYR A  28      32.543  15.057  22.795  1.00 30.00           C
ATOM    326  CZ  TYR A  28      32.304  16.370  23.143  1.00 30.00           C
ATOM    327  OH  TYR A  28      31.911  16.677  24.425  1.00 30.00           O
ATOM    328  H   TYR A  28      30.554  14.537  18.752  1.00 15.46           H
ATOM    329  HA  TYR A  28      31.829  16.276  18.214  1.00 17.48           H
ATOM    330  HB2 TYR A  28      34.230  16.179  18.815  1.00 30.00           H
ATOM    331  HB3 TYR A  28      34.030  14.441  19.171  1.00 30.00           H
ATOM    332  HD1 TYR A  28      32.971  17.852  20.199  1.00 30.00           H
ATOM    333  HD2 TYR A  28      33.126  13.724  21.238  1.00 30.00           H
ATOM    334  HE1 TYR A  28      32.269  18.403  22.487  1.00 30.00           H
ATOM    335  HE2 TYR A  28      32.422  14.271  23.526  1.00 30.00           H
ATOM    336  HH  TYR A  28      31.828  17.648  24.527  1.00 30.00           H
ATOM    337  N   CYS A  29      32.016  15.579  15.769  1.00 15.90           N
ATOM    338  CA  CYS A  29      32.342  15.537  14.348  1.00 17.10           C
ATOM    339  C   CYS A  29      31.892  16.813  13.644  1.00 16.88           C
ATOM    340  O   CYS A  29      30.779  17.293  13.862  1.00 17.69           O
ATOM    341  CB  CYS A  29      31.700  14.316  13.687  1.00 30.00           C
ATOM    342  SG  CYS A  29      32.353  12.728  14.255  1.00 30.00           S
ATOM    343  OXT CYS A  29      32.632  17.387  12.845  1.00 30.00           O
ATOM    344  H   CYS A  29      31.069  15.902  15.968  1.00 15.90           H
ATOM    345  HA  CYS A  29      33.423  15.453  14.235  1.00 17.10           H
ATOM    346  HB2 CYS A  29      30.631  14.328  13.900  1.00 30.00           H
ATOM    347  HB3 CYS A  29      31.867  14.375  12.611  1.00 30.00           H
ATOM    348  HG  CYS A  29      33.116  12.929  15.288  1.00 30.00           H
TER
END
'''

olex_res_str_3 = '''
TITL ???? in I 4 3 2
CELL 0.7107 108.91 108.91 108.91 90 90 90
ZERR 48 0 0 0 0 0 0
LATT -2
SYMM +X,-Z,+Y
SYMM +X,+Z,-Y
SYMM +Z,+Y,-X
SYMM -Z,+Y,+X
SYMM -Y,+X,+Z
SYMM +Y,-X,+Z
SYMM +Z,+X,+Y
SYMM +Y,+Z,+X
SYMM -Y,-Z,+X
SYMM +Z,-X,-Y
SYMM -Y,+Z,-X
SYMM -Z,-X,+Y
SYMM -Z,+X,-Y
SYMM +Y,-Z,-X
SYMM +X,-Y,-Z
SYMM -X,+Y,-Z
SYMM -X,-Y,+Z
SYMM +Y,+X,-Z
SYMM -Y,-X,-Z
SYMM +Z,-Y,+X
SYMM -Z,-Y,-X
SYMM -X,+Z,+Y
SYMM -X,-Z,-Y
SFAC C H N O S
DISP C 0 0 0
DISP H 0 0 0
DISP N 0 0 0
DISP O 0 0 0
DISP S 0 0 0
UNIT 672 1248 192 192 48

L.S. 0
PLAN  1872
CONF
list 6
ACTA
MORE -1
BOND $H
WGHT 0.1
FVAR 0.5
REM <olex2.extras>
REM <HklSrc "%.\\olex_test.hkl">
REM </olex2.extras>

N001  3    -0.460004 -0.113901 -0.082169  11.00000  1.52000  1.09500 =
 0.94030  0.09390 -0.12790  0.20610
C001  1    -0.454761 -0.119062 -0.093435  11.00000  1.44330  1.15800 =
 0.99140  0.38340  0.21930  0.30810
C002  1    -0.464521 -0.125856 -0.100716  11.00000  1.57750  1.08840 =
 1.15210  0.37970 -0.03150  0.26050
O001  4    -0.470434 -0.134395 -0.096162  11.00000  1.47290  0.95010 =
 1.40040  0.37420  0.00760  0.37620
C003  1    -0.444771 -0.127803 -0.088715  11.00000  1.16970  1.26970 =
 0.75470  0.24010  0.39750  0.37990
C004  1    -0.441980 -0.123689 -0.076403  11.00000  1.29730  1.09530 =
 1.03030  0.13600  0.13440  0.29600
C005  1    -0.453769 -0.118704 -0.071206  11.00000  1.40020  0.96160 =
 1.39460  0.14170 -0.06990  0.23470
H001  2    -0.468552 -0.116004 -0.081756  11.00000  1.18508
H002  2    -0.459141 -0.105133 -0.082352  11.00000  1.18508
H003  2    -0.450868 -0.112873 -0.098513  11.00000  1.19762
H004  2    -0.447865 -0.136158 -0.088642  11.00000  1.06476
H005  2    -0.437581 -0.127261 -0.093940  11.00000  1.06476
H006  2    -0.435754 -0.117326 -0.076715  11.00000  1.14100
H007  2    -0.439069 -0.130567 -0.071545  11.00000  1.14100
H008  2    -0.458562 -0.125168 -0.067377  11.00000  1.25208
H009  2    -0.452172 -0.112203 -0.065329  11.00000  1.25208
N002  3    -0.466339 -0.121807 -0.112056  11.00000  1.72800  0.98000 =
 1.23710  0.39610  0.04510  0.07080
C006  1    -0.475907 -0.126609 -0.120173  11.00000  1.74590  1.27030 =
 1.40410  0.38340  0.07000  0.18320
C007  1    -0.469718 -0.131705 -0.131659  11.00000  1.61380  1.18380 =
 1.27490  0.42320  0.05450  0.09460
O002  4    -0.468580 -0.125397 -0.140970  11.00000  1.60940  1.10460 =
 1.43150  0.17160  0.05270 -0.19000
C008  1    -0.485272 -0.116720 -0.123276  11.00000  1.82860  1.59380 =
 1.33520  0.36050  0.02400  0.34740
S001  5    -0.496722 -0.121697 -0.134028  11.00000  2.02610  1.45890 =
 1.10840  0.44800 -0.09960  0.47210
H010  2    -0.462033 -0.116096 -0.115407  11.00000  1.31502
H011  2    -0.480525 -0.133119 -0.116215  11.00000  1.47346
H012  2    -0.481021 -0.109751 -0.126839  11.00000  1.58593
H013  2    -0.489432 -0.114250 -0.115802  11.00000  1.58593
N003  3    -0.465834 -0.143247 -0.130961  11.00000  1.57940  1.35440 =
 1.22220  0.44670  0.04770  0.35700
C009  1    -0.460151 -0.149573 -0.141273  11.00000  1.70310  1.54950 =
 1.22010  0.63150 -0.15850  0.35240
C010  1    -0.461225 -0.163364 -0.139207  11.00000  1.72710  1.36270 =
 1.34880  0.43340 -0.17190  0.60930
O003  4    -0.463089 -0.167570 -0.128859  11.00000  1.96950  1.48810 =
 1.52230  0.24180 -0.11270  0.78080
C011  1    -0.446653 -0.145772 -0.142898  11.00000  1.60900  1.54900 =
 1.41850  0.76920 -0.04650  0.11880
C012  1    -0.438546 -0.146818 -0.131503  11.00000  1.43960  1.49200 =
 1.49190  0.87140 -0.00910  0.19530
C013  1    -0.425535 -0.142733 -0.134441  11.00000  1.64250  1.60410 =
 1.47830  0.86950  0.11870  0.28570
C014  1    -0.416656 -0.145533 -0.124102  11.00000  1.87790  1.69550 =
 1.39570  0.86650  0.32760  0.23890
N004  3    -0.403801 -0.142384 -0.127546  11.00000  2.05350  1.77950 =
 1.39580  0.47120  0.29060  0.53700
O004  4    -0.460252 -0.170600 -0.147838  11.00000  0.37995
H014  2    -0.466403 -0.147700 -0.124470  11.00000  1.38531
H015  2    -0.464255 -0.147369 -0.148866  11.00000  1.49094
H016  2    -0.446460 -0.137251 -0.145487  11.00000  1.52552
H017  2    -0.443045 -0.150950 -0.149178  11.00000  1.52552
H018  2    -0.441778 -0.141612 -0.125039  11.00000  1.47448
H019  2    -0.438316 -0.155284 -0.128758  11.00000  1.47448
H020  2    -0.425480 -0.133936 -0.135828  11.00000  1.57491
H021  2    -0.422716 -0.146984 -0.141741  11.00000  1.57491
H022  2    -0.418887 -0.140749 -0.116922  11.00000  1.65635
H023  2    -0.416996 -0.154228 -0.122202  11.00000  1.65635
H024  2    -0.403250 -0.134423 -0.129290  11.00000  1.74298
H025  2    -0.401524 -0.146699 -0.134102  11.00000  1.74298
H026  2    -0.398751 -0.144101 -0.121357  11.00000  1.74298
HKLF 4

END
'''

olex_res_str_2 = '''
TITL ???? in P 1 21 1
CELL 0.7107 80.435 118.034 112.075 90 93.12 90
ZERR 2 0 0 0 0 0 0
LATT -1
SYMM -X,0.5+Y,-Z
SFAC C H N O
DISP C 0 0 0
DISP H 0 0 0
DISP N 0 0 0
DISP O 0 0 0
UNIT 52 88 8 36

L.S. 0
PLAN  37441
CONF
list 6
ACTA
MORE -1
BOND $H
WGHT 0.1
FVAR 0.5
REM <olex2.extras>
REM <HklSrc "%.\\olex_test.hkl">
REM </olex2.extras>

N001  3     1.050607 -0.176093  1.391588  11.00000  1.71637  2.19700 =
 1.51970  0.08700  0.09509 -0.66457
C001  1     1.050011 -0.164487  1.396020  11.00000  1.75759  2.27970 =
 1.61290  0.12520  0.09068 -0.61117
C002  1     1.067521 -0.159539  1.396091  11.00000  1.71767  2.18230 =
 1.56450  0.11340  0.09544 -0.56977
O001  4     1.079033 -0.164012  1.402141  11.00000  1.78277  2.15880 =
 1.60950  0.08810  0.11586 -0.55128
C003  1     1.043226 -0.164275  1.408521  11.00000  1.87662  2.38930 =
 1.76720  0.14190  0.10667 -0.58320
C004  1     1.040450 -0.152346  1.413230  11.00000  1.96435  2.52390 =
 1.91200  0.18050  0.10866 -0.53297
O002  4     1.051005 -0.144874  1.412408  11.00000  1.88362  2.42820 =
 1.83990  0.18390  0.10513 -0.49694
N002  3     1.026146 -0.150601  1.418386  11.00000  2.16874  2.78290 =
 2.15720  0.20880  0.11901 -0.53083
O003  4     1.070441 -0.151058  1.390158  11.00000  0.37995
H001  2     1.055235 -0.176110  1.383099  11.00000  1.81074
H002  2     1.058341 -0.180922  1.397333  11.00000  1.81074
H003  2     1.038678 -0.179525  1.391177  11.00000  1.81074
H004  2     1.041663 -0.159251  1.390363  11.00000  1.88343
H005  2     1.052088 -0.168511  1.414722  11.00000  2.01072
H006  2     1.031307 -0.168689  1.408244  11.00000  2.01072
H007  2     1.017815 -0.157099  1.418690  11.00000  2.36952
C005  1     1.022026 -0.139824  1.423560  11.00000  1.31110
C006  1     1.007895 -0.133428  1.416903  11.00000  1.44725
C007  1     1.004464 -0.122287  1.423301  11.00000  1.49436
C008  1     1.001951 -0.124185  1.436553  11.00000  1.57580
C009  1     1.016072 -0.131394  1.442272  11.00000  1.53970
C010  1     1.013166 -0.135020  1.454934  11.00000  1.55541
C011  1     1.000849 -0.132597  1.395305  11.00000  1.58352
C012  1     1.007121 -0.129996  1.383250  11.00000  1.61037
N003  3     1.011835 -0.131293  1.404589  11.00000  1.53894
O004  4     0.990096 -0.117034  1.417779  11.00000  1.48043
O005  4     1.001913 -0.113408  1.442281  11.00000  1.68890
O006  4     1.017964 -0.141629  1.435579  11.00000  1.45561
O007  4     1.027743 -0.133750  1.462556  11.00000  1.54072
O008  4     0.986511 -0.135724  1.396574  11.00000  1.58910
H008  2     1.033028 -0.134427  1.423560  11.00000  1.31110
H009  2     0.996790 -0.138723  1.416867  11.00000  1.44725
H010  2     1.015601 -0.117153  1.422443  11.00000  1.49436
H011  2     0.990107 -0.128624  1.436946  11.00000  1.57580
H012  2     1.027339 -0.126260  1.442307  11.00000  1.53970
H013  2     1.002897 -0.129996  1.458052  11.00000  1.55541
H014  2     1.008917 -0.143789  1.454648  11.00000  1.55541
H015  2     1.012307 -0.137664  1.379506  11.00000  1.61037
H016  2     1.016739 -0.123507  1.384171  11.00000  1.61037
H017  2     0.996841 -0.126938  1.377371  11.00000  1.61037
H018  2     1.023603 -0.128666  1.403043  11.00000  1.53894
H019  2     0.983225 -0.122973  1.413900  11.00000  1.48043
H020  2     1.029049 -0.125701  1.464531  11.00000  1.54072
C013  1     0.988864 -0.111908  1.450555  11.00000  1.81390
C014  1     0.993879 -0.101776  1.458481  11.00000  1.89851
C015  1     0.980674 -0.099556  1.467471  11.00000  1.94917
C016  1     0.963591 -0.098277  1.461064  11.00000  1.98957
C017  1     0.959810 -0.108494  1.453013  11.00000  1.92852
C018  1     0.943621 -0.107359  1.445614  11.00000  1.92194
C019  1     1.023464 -0.097302  1.462529  11.00000  1.90636
C020  1     1.039171 -0.100335  1.469633  11.00000  1.89420
N004  3     1.010082 -0.103580  1.464593  11.00000  1.91738
O009  4     0.984866 -0.089559  1.474065  11.00000  1.94043
O010  4     0.951131 -0.097345  1.469553  11.00000  2.06645
O011  4     0.972737 -0.109901  1.444783  11.00000  1.87229
O012  4     0.945697 -0.100869  1.435114  11.00000  1.91700
O013  4     1.022925 -0.089517  1.455255  11.00000  1.88926
H021  2     0.988313 -0.119381  1.456229  11.00000  1.81390
H022  2     0.994746 -0.094397  1.452700  11.00000  1.89851
H023  2     0.980618 -0.106952  1.473288  11.00000  1.94917
H024  2     0.964448 -0.090542  1.455818  11.00000  1.98957
H025  2     0.959076 -0.115924  1.458740  11.00000  1.92852
H026  2     0.934511 -0.103513  1.451225  11.00000  1.92194
H027  2     0.939107 -0.115865  1.443469  11.00000  1.92194
H028  2     1.039851 -0.095515  1.477944  11.00000  1.89420
H029  2     1.049855 -0.098311  1.464450  11.00000  1.89420
H030  2     1.039123 -0.109375  1.471635  11.00000  1.89420
H031  2     1.011169 -0.109960  1.470732  11.00000  1.91738
H032  2     0.978890 -0.083205  1.470178  11.00000  1.94043
H033  2     0.952194 -0.094100  1.437348  11.00000  1.91700
C021  1     0.944871 -0.085865  1.469633  11.00000  2.08772
C022  1     0.925760 -0.086306  1.467310  11.00000  2.08354
C023  1     0.918370 -0.074419  1.468391  11.00000  2.08772
C024  1     0.924734 -0.068497  1.480008  11.00000  2.09925
C025  1     0.943884 -0.068768  1.481205  11.00000  2.08139
C026  1     0.950874 -0.063685  1.492884  11.00000  2.02579
O014  4     0.918348 -0.093177  1.475969  11.00000  2.05365
O015  4     0.900653 -0.074716  1.467721  11.00000  2.06075
O016  4     0.919131 -0.057077  1.480061  11.00000  2.08950
O017  4     0.949254 -0.080392  1.480705  11.00000  2.08494
O018  4     0.942815 -0.053129  1.494654  11.00000  1.96158
H034  2     0.950037 -0.081129  1.462351  11.00000  2.08772
H035  2     0.924165 -0.089627  1.458249  11.00000  2.08354
H036  2     0.922640 -0.069971  1.460599  11.00000  2.08772
H037  2     0.919167 -0.072784  1.487424  11.00000  2.09925
H038  2     0.949091 -0.063761  1.474137  11.00000  2.08139
H039  2     0.948834 -0.069472  1.500247  11.00000  2.02579
H040  2     0.964288 -0.062533  1.492554  11.00000  2.02579
H041  2     0.923308 -0.091033  1.483823  11.00000  2.05365
H042  2     0.896854 -0.080045  1.473699  11.00000  2.06075
H043  2     0.918206 -0.054662  1.488363  11.00000  2.08950
H044  2     0.939224 -0.050062  1.486808  11.00000  1.96158
HKLF 4
'''

olex_res_str_1 = '''
TITL ???? in P 1
CELL 0.7107 39.153 32.696 29.501 90 90 90
ZERR 1 0 0 0 0 0 0
LATT -1
SFAC C H N O S
DISP C 0 0 0
DISP H 0 0 0
DISP N 0 0 0
DISP O 0 0 0
DISP S 0 0 0
UNIT 112.5 166.5 30.5 32 4

L.S. 0
PLAN  2452
CONF
list 6
ACTA
MORE -1
BOND $H
WGHT 0.1
FVAR 0.5
REM <olex2.extras>
REM <HklSrc "%.\\olex_test.hkl">
REM </olex2.extras>

N001  3     0.390085  0.667543  0.312091  11.00000  0.18440
C001  1     0.368758  0.692868  0.341310  11.00000  0.20391
C002  1     0.335223  0.671367  0.352361  11.00000  0.19960
O001  4     0.316604  0.658429  0.321955  11.00000  0.19251
C003  1     0.361735  0.734249  0.318938  11.00000  0.23241
C004  1     0.356448  0.768198  0.353242  11.00000  0.29092
O002  4     0.346002  0.759023  0.391919  11.00000  0.28421
O003  4     0.362731  0.804135  0.341514  11.00000  0.35361
H001  2     0.394606  0.639589  0.327785  11.00000  0.18440
H002  2     0.377672  0.662711  0.281414  11.00000  0.18440
H003  2     0.413199  0.682285  0.306261  11.00000  0.18440
H004  2     0.382193  0.698465  0.373072  11.00000  0.20391
H005  2     0.383444  0.742476  0.297651  11.00000  0.23241
H006  2     0.338697  0.731894  0.298363  11.00000  0.23241
N002  3     0.327663  0.667635  0.396190  11.00000  0.20366
C005  1     0.296555  0.647877  0.411986  11.00000  0.22544
C006  1     0.268255  0.678156  0.423308  11.00000  0.25204
O004  4     0.241003  0.665372  0.438934  11.00000  0.24963
C007  1     0.304217  0.621299  0.453544  11.00000  0.21467
O005  4     0.323756  0.586463  0.441443  11.00000  0.22430
H007  2     0.342809  0.678156  0.421782  11.00000  0.20366
H008  2     0.287053  0.628609  0.384597  11.00000  0.22544
H009  2     0.318571  0.639375  0.478052  11.00000  0.21467
H010  2     0.280157  0.611206  0.468425  11.00000  0.21467
H011  2     0.321099  0.565054  0.464425  11.00000  0.22430
N003  3     0.274002  0.717825  0.415884  11.00000  0.26166
C008  1     0.247848  0.748195  0.426087  11.00000  0.28801
C009  1     0.221618  0.751957  0.388156  11.00000  0.28243
O006  4     0.192374  0.763733  0.396868  11.00000  0.30979
C010  1     0.264424  0.789791  0.435104  11.00000  0.34627
C011  1     0.284831  0.790464  0.478763  11.00000  0.39439
O007  4     0.272316  0.778169  0.514525  11.00000  0.44151
N004  3     0.316476  0.804839  0.476018  11.00000  0.41630
H012  2     0.296120  0.729569  0.403071  11.00000  0.26166
H013  2     0.234184  0.738531  0.456425  11.00000  0.28801
H014  2     0.244502  0.812974  0.437578  11.00000  0.34627
H015  2     0.281690  0.797070  0.407274  11.00000  0.34627
H016  2     0.331443  0.806031  0.504288  11.00000  0.41630
H017  2     0.325824  0.814687  0.445646  11.00000  0.41630
N005  3     0.231425  0.742170  0.346327  11.00000  0.25723
C012  1     0.208209  0.744678  0.307752  11.00000  0.23785
C013  1     0.192910  0.702655  0.298803  11.00000  0.23658
O008  4     0.161622  0.697669  0.297990  11.00000  0.23127
C014  1     0.226522  0.761194  0.265347  11.00000  0.23304
C015  1     0.238347  0.805236  0.273177  11.00000  0.23317
C016  1     0.203228  0.758655  0.224196  11.00000  0.22633
C017  1     0.266135  0.819030  0.241009  11.00000  0.26457
H018  2     0.255357  0.732383  0.338056  11.00000  0.25723
H019  2     0.187623  0.765935  0.315684  11.00000  0.23785
H020  2     0.248946  0.742293  0.259008  11.00000  0.23304
H021  2     0.248257  0.807836  0.307583  11.00000  0.23317
H022  2     0.216586  0.825636  0.268771  11.00000  0.23317
H023  2     0.180395  0.776884  0.230501  11.00000  0.22633
H024  2     0.196179  0.726817  0.218569  11.00000  0.22633
H025  2     0.216791  0.770522  0.194705  11.00000  0.22633
H026  2     0.256634  0.817042  0.206332  11.00000  0.26457
H027  2     0.288305  0.799211  0.245144  11.00000  0.26457
H028  2     0.273006  0.850532  0.248907  11.00000  0.26457
N006  3     0.214849  0.672162  0.291990  11.00000  0.21721
C018  1     0.204097  0.630322  0.283075  11.00000  0.20049
C019  1     0.219268  0.602214  0.318871  11.00000  0.20328
O009  4     0.250198  0.603224  0.327243  11.00000  0.20201
C020  1     0.215054  0.616559  0.235619  11.00000  0.18757
C021  1     0.200010  0.576493  0.221348  11.00000  0.20328
N007  3     0.171660  0.573709  0.193824  11.00000  0.22037
C022  1     0.209384  0.537528  0.230772  11.00000  0.21049
C023  1     0.164202  0.534805  0.186909  11.00000  0.21974
N008  3     0.186703  0.512203  0.208942  11.00000  0.24494
H029  2     0.240748  0.675618  0.292465  11.00000  0.21721
H030  2     0.176334  0.628303  0.284228  11.00000  0.20049
H031  2     0.242766  0.613439  0.235009  11.00000  0.18757
H032  2     0.207008  0.639589  0.211145  11.00000  0.18757
H033  2     0.158481  0.597932  0.180807  11.00000  0.22037
H034  2     0.230685  0.528107  0.251551  11.00000  0.21049
H035  2     0.143258  0.523428  0.166672  11.00000  0.21974
H036  2     0.187163  0.481007  0.209654  11.00000  0.24494
N009  3     0.197558  0.577257  0.339853  11.00000  0.19682
C024  1     0.208873  0.548844  0.374835  11.00000  0.20429
C025  1     0.224300  0.511011  0.352700  11.00000  0.18428
O010  4     0.207136  0.489754  0.326599  11.00000  0.18567
C026  1     0.179067  0.537161  0.405613  11.00000  0.25254
C027  1     0.189385  0.512356  0.447002  10.50000  0.33892
C028  1     0.159068  0.504924  0.478119  10.50000  0.38717
C029  1     0.170562  0.497156  0.526592  10.50000  0.43340
N010  3     0.141088  0.493577  0.558218  10.50000  0.46570
H037  2     0.171992  0.576034  0.333345  11.00000  0.19682
H038  2     0.228233  0.563096  0.396156  11.00000  0.20429
H039  2     0.166782  0.565054  0.417579  11.00000  0.25254
H040  2     0.160984  0.518901  0.386123  11.00000  0.25254
H041  2     0.208847  0.528964  0.465984  10.50000  0.33892
H042  2     0.199346  0.482811  0.436155  10.50000  0.33892
H043  2     0.142543  0.531747  0.477916  10.50000  0.38717
H044  2     0.144970  0.478224  0.466391  10.50000  0.38717
H045  2     0.186499  0.522510  0.537880  10.50000  0.43340
H046  2     0.184992  0.468681  0.527914  10.50000  0.43340
H047  2     0.125457  0.469415  0.548422  10.50000  0.46570
H048  2     0.126887  0.520461  0.557879  10.50000  0.46570
H049  2     0.149797  0.488286  0.590997  10.50000  0.46570
N011  3     0.256609  0.502630  0.363072  11.00000  0.15553
C030  1     0.274921  0.467672  0.344598  11.00000  0.14058
C031  1     0.284525  0.438525  0.382733  11.00000  0.12779
O011  4     0.291753  0.451676  0.420799  11.00000  0.13894
C032  1     0.307052  0.481221  0.319006  11.00000  0.12944
S001  5     0.326744  0.441552  0.285211  11.00000  0.14058
H050  2     0.271039  0.519850  0.384597  11.00000  0.15553
H051  2     0.258626  0.451768  0.320362  11.00000  0.14058
H052  2     0.325926  0.491650  0.343615  11.00000  0.12944
H053  2     0.300232  0.506178  0.296261  11.00000  0.12944
N012  3     0.284550  0.398917  0.372326  11.00000  0.13488
C033  1     0.293030  0.367904  0.405207  11.00000  0.14274
C034  1     0.331213  0.362674  0.411071  11.00000  0.13716
O012  4     0.346359  0.376835  0.444256  11.00000  0.14476
H054  2     0.278727  0.387570  0.341073  11.00000  0.13488
H055  2     0.282098  0.376193  0.437917  11.00000  0.14274
H056  2     0.282303  0.338726  0.394393  11.00000  0.14274
N013  3     0.347636  0.341938  0.378631  11.00000  0.13602
C035  1     0.384773  0.335454  0.381818  11.00000  0.12159
C036  1     0.405742  0.374358  0.377106  11.00000  0.11234
O013  4     0.436620  0.373165  0.386055  11.00000  0.11221
C037  1     0.392230  0.305756  0.342700  11.00000  0.13070
C038  1     0.358210  0.287895  0.329955  11.00000  0.13729
C039  1     0.332669  0.320743  0.339378  11.00000  0.14362
H057  2     0.391081  0.321110  0.414088  11.00000  0.12159
H058  2     0.403392  0.322364  0.314294  11.00000  0.13070
H059  2     0.409700  0.281931  0.354124  11.00000  0.13070
H060  2     0.353076  0.260888  0.350530  11.00000  0.13729
H061  2     0.358185  0.279759  0.294126  11.00000  0.13729
H062  2     0.330217  0.341357  0.310532  11.00000  0.14362
H063  2     0.307844  0.307713  0.347853  11.00000  0.14362
N014  3     0.390213  0.408552  0.363411  11.00000  0.10221
C040  1     0.409547  0.446446  0.358496  11.00000  0.10664
C041  1     0.420147  0.463879  0.404291  11.00000  0.11791
O014  4     0.445892  0.485595  0.407207  11.00000  0.13830
C042  1     0.388604  0.477887  0.332090  11.00000  0.10347
S002  5     0.375501  0.461188  0.276330  11.00000  0.13881
H064  2     0.364902  0.410845  0.355683  11.00000  0.10221
H065  2     0.432866  0.440054  0.339616  11.00000  0.10664
H066  2     0.403928  0.505444  0.327887  11.00000  0.10347
H067  2     0.365566  0.484861  0.351344  11.00000  0.10347
N015  3     0.400863  0.454765  0.440188  11.00000  0.11652
C043  1     0.409624  0.470241  0.485034  11.00000  0.13590
C044  1     0.436620  0.444060  0.508254  11.00000  0.12526
O015  4     0.454882  0.459169  0.538083  11.00000  0.14920
C045  1     0.377621  0.472749  0.514661  11.00000  0.37995
C046  1     0.352770  0.506178  0.500729  11.00000  0.37995
C047  1     0.321074  0.507157  0.530219  11.00000  0.37995
O016  4     0.315582  0.477734  0.556218  11.00000  0.37995
O017  4     0.301663  0.537375  0.527440  11.00000  0.37995
H068  2     0.379435  0.437057  0.439205  11.00000  0.11652
H069  2     0.420274  0.500856  0.481272  11.00000  0.13590
H070  2     0.385386  0.478774  0.549507  11.00000  0.37995
H071  2     0.364136  0.443632  0.512728  11.00000  0.37995
H072  2     0.365438  0.535754  0.503542  11.00000  0.37995
H073  2     0.344750  0.500948  0.465815  11.00000  0.37995
N016  3     0.439583  0.405187  0.495271  11.00000  0.11677
C048  1     0.464741  0.377967  0.515610  11.00000  0.13336
C049  1     0.499068  0.379221  0.491441  11.00000  0.13146
O018  4     0.525145  0.368394  0.511779  11.00000  0.14654
C050  1     0.451306  0.333986  0.516050  11.00000  0.37995
C051  1     0.420325  0.327594  0.546287  11.00000  0.37995
C052  1     0.406074  0.284561  0.543066  11.00000  0.37995
O019  4     0.410799  0.263580  0.509000  11.00000  0.37995
N017  3     0.388527  0.271226  0.578624  11.00000  0.37995
H074  2     0.424948  0.391760  0.470831  11.00000  0.11677
H075  2     0.469057  0.388151  0.550320  11.00000  0.13336
H076  2     0.471535  0.313861  0.528185  11.00000  0.37995
H077  2     0.444053  0.325514  0.481645  11.00000  0.37995
H078  2     0.427605  0.333129  0.581404  11.00000  0.37995
H079  2     0.400276  0.348850  0.536185  11.00000  0.37995
H080  2     0.378336  0.242507  0.578421  11.00000  0.37995
H081  2     0.385258  0.289699  0.606149  11.00000  0.37995
N018  3     0.499681  0.392158  0.448527  11.00000  0.12083
C053  1     0.531734  0.394115  0.423477  11.00000  0.11994
C054  1     0.552831  0.431521  0.436256  11.00000  0.13096
O020  4     0.584068  0.431398  0.429985  11.00000  0.13412
C055  1     0.524762  0.393963  0.372496  11.00000  0.12387
H082  2     0.478737  0.401303  0.430629  11.00000  0.12083
H083  2     0.546752  0.367140  0.432087  11.00000  0.11994
H084  2     0.509131  0.420449  0.363954  11.00000  0.12387
H085  2     0.511302  0.365916  0.363615  11.00000  0.12387
H086  2     0.549000  0.395431  0.354429  11.00000  0.12387
N019  3     0.536511  0.463543  0.453985  11.00000  0.13944
C056  1     0.555155  0.500306  0.467001  11.00000  0.15667
C057  1     0.571578  0.496238  0.513576  11.00000  0.14451
O021  4     0.598422  0.514650  0.521948  11.00000  0.15857
C058  1     0.531351  0.537405  0.466086  11.00000  0.37995
C059  1     0.518351  0.548446  0.419104  11.00000  0.37995
S003  5     0.552346  0.560069  0.379987  11.00000  0.37995
C060  1     0.569177  0.605946  0.404563  12.00000  0.37995
H087  2     0.510842  0.464858  0.459612  11.00000  0.13944
H088  2     0.575486  0.505199  0.442358  11.00000  0.15667
H089  2     0.545118  0.563800  0.479306  11.00000  0.37995
H090  2     0.509207  0.530921  0.487306  11.00000  0.37995
H091  2     0.501954  0.575300  0.421647  11.00000  0.37995
H092  2     0.503946  0.522694  0.405478  11.00000  0.37995
H093  2     0.577963  0.599339  0.438833  12.00000  0.37995
H094  2     0.549128  0.629037  0.405512  12.00000  0.37995
H095  2     0.590453  0.616712  0.383919  12.00000  0.37995
N020  3     0.556126  0.472565  0.543982  11.00000  0.15363
C061  1     0.569612  0.466510  0.589472  11.00000  0.17402
C062  1     0.594105  0.430420  0.592658  11.00000  0.16743
O022  4     0.607693  0.422437  0.629402  11.00000  0.16211
C063  1     0.540265  0.460637  0.623233  11.00000  0.22240
C064  1     0.525911  0.500306  0.642148  11.00000  0.26622
C065  1     0.487268  0.501835  0.637673  11.00000  0.32600
N021  3     0.474089  0.543308  0.642453  11.00000  0.34880
C066  1     0.441167  0.552881  0.643097  11.00000  0.36222
N022  3     0.417209  0.524560  0.639334  11.00000  0.33677
N023  3     0.432074  0.591877  0.647639  11.00000  0.40212
H096  2     0.533982  0.457028  0.538016  11.00000  0.15363
H097  2     0.584093  0.493761  0.598556  11.00000  0.17402
H098  2     0.519500  0.444611  0.606217  11.00000  0.22240
H099  2     0.549511  0.442531  0.651707  11.00000  0.22240
H100  2     0.536638  0.526058  0.623538  11.00000  0.26622
H101  2     0.532424  0.502753  0.677977  11.00000  0.26622
H102  2     0.479938  0.490427  0.604352  11.00000  0.32600
H103  2     0.475800  0.482934  0.664045  11.00000  0.32600
H104  2     0.491329  0.566522  0.645504  11.00000  0.34880
H105  2     0.392077  0.532726  0.639911  11.00000  0.33677
H106  2     0.423876  0.494556  0.635843  11.00000  0.33677
H107  2     0.406840  0.599676  0.648181  11.00000  0.40212
H108  2     0.450336  0.613989  0.650588  11.00000  0.40212
N024  3     0.600031  0.408552  0.555303  11.00000  0.14818
C067  1     0.623273  0.373685  0.556422  11.00000  0.13995
C068  1     0.660486  0.387264  0.554015  11.00000  0.13514
O023  4     0.682374  0.369923  0.577269  11.00000  0.13615
C069  1     0.615074  0.344262  0.517643  11.00000  0.14692
C070  1     0.636222  0.304961  0.515474  11.00000  0.18086
C071  1     0.632800  0.280799  0.559303  11.00000  0.19922
C072  1     0.625699  0.279086  0.474967  11.00000  0.17035
H109  2     0.589304  0.414179  0.524423  11.00000  0.14818
H110  2     0.619927  0.357414  0.588353  11.00000  0.13995
H111  2     0.619084  0.360442  0.485780  11.00000  0.14692
H112  2     0.588384  0.335209  0.520728  11.00000  0.14692
H113  2     0.663040  0.313005  0.511271  11.00000  0.18086
H114  2     0.606033  0.272908  0.564422  11.00000  0.19922
H115  2     0.641867  0.299547  0.587404  11.00000  0.19922
H116  2     0.648098  0.253028  0.556761  11.00000  0.19922
H117  2     0.629786  0.296611  0.444019  11.00000  0.17035
H118  2     0.598779  0.271165  0.478221  11.00000  0.17035
H119  2     0.641177  0.251376  0.474323  11.00000  0.17035
N025  3     0.668046  0.418125  0.526253  11.00000  0.12538
C073  1     0.702858  0.433264  0.521542  11.00000  0.12387
C074  1     0.714224  0.459536  0.561405  11.00000  0.12906
O024  4     0.744898  0.467550  0.565913  11.00000  0.14920
C075  1     0.706408  0.457701  0.477543  11.00000  0.37995
C076  1     0.700610  0.432775  0.435477  11.00000  0.37995
C077  1     0.726688  0.409377  0.417206  11.00000  0.37995
C078  1     0.669068  0.432744  0.414155  11.00000  0.37995
C079  1     0.721426  0.386378  0.378496  11.00000  0.37995
C080  1     0.663653  0.409836  0.375445  11.00000  0.37995
C081  1     0.689858  0.386622  0.357581  11.00000  0.37995
H120  2     0.650627  0.433172  0.506695  11.00000  0.12538
H121  2     0.719843  0.406869  0.520491  11.00000  0.12387
H122  2     0.732153  0.470272  0.475780  11.00000  0.37995
H123  2     0.687712  0.482414  0.477747  11.00000  0.37995
H124  2     0.751437  0.409102  0.433375  11.00000  0.37995
H125  2     0.648456  0.450850  0.427952  11.00000  0.37995
H126  2     0.741987  0.368241  0.364632  11.00000  0.37995
H127  2     0.638929  0.410081  0.359174  11.00000  0.37995
H128  2     0.685695  0.368700  0.327345  11.00000  0.37995
N026  3     0.690522  0.473116  0.590183  11.00000  0.14945
C082  1     0.700049  0.498318  0.628996  11.00000  0.17592
C083  1     0.716854  0.472963  0.666181  11.00000  0.18732
O025  4     0.742012  0.486359  0.686689  11.00000  0.19935
C084  1     0.668633  0.520522  0.647571  11.00000  0.18719
O026  4     0.651853  0.541962  0.611539  11.00000  0.23139
C085  1     0.678901  0.550893  0.684214  11.00000  0.18896
H129  2     0.665007  0.467244  0.587506  11.00000  0.14945
H130  2     0.718821  0.521165  0.618793  11.00000  0.17592
H131  2     0.651138  0.498287  0.662384  11.00000  0.18719
H132  2     0.636835  0.523275  0.594895  11.00000  0.23139
H133  2     0.696447  0.573526  0.670316  11.00000  0.18896
H134  2     0.691441  0.534928  0.712044  11.00000  0.18896
H135  2     0.656348  0.566522  0.697197  11.00000  0.18896
N027  3     0.703420  0.436353  0.675164  11.00000  0.19023
C086  1     0.718080  0.410448  0.710484  11.00000  0.21176
C087  1     0.752331  0.393076  0.695265  11.00000  0.20378
O027  4     0.776671  0.390843  0.721874  11.00000  0.20796
C088  1     0.692795  0.376346  0.723467  11.00000  0.22924
C089  1     0.707149  0.350502  0.761906  11.00000  0.25381
C090  1     0.658443  0.394544  0.736585  11.00000  0.24000
H136  2     0.682834  0.424089  0.658893  11.00000  0.19023
H137  2     0.722959  0.428768  0.740687  11.00000  0.21176
H138  2     0.688913  0.356466  0.694248  11.00000  0.22924
H139  2     0.711567  0.370045  0.791261  11.00000  0.25381
H140  2     0.731132  0.336738  0.750924  11.00000  0.25381
H141  2     0.688785  0.326645  0.770415  11.00000  0.25381
H142  2     0.648277  0.411610  0.707874  11.00000  0.24000
H143  2     0.661967  0.414913  0.765466  11.00000  0.24000
H144  2     0.641024  0.369831  0.745602  11.00000  0.24000
N028  3     0.754859  0.380719  0.652385  11.00000  0.19428
C091  1     0.786325  0.363408  0.634284  11.00000  0.19745
C092  1     0.812811  0.395736  0.622115  11.00000  0.20644
O028  4     0.841494  0.384604  0.609674  11.00000  0.20594
C093  1     0.778484  0.337197  0.592522  11.00000  0.37995
C094  1     0.760708  0.297835  0.603471  11.00000  0.37995
C095  1     0.726356  0.290892  0.603437  11.00000  0.37995
C096  1     0.776441  0.259940  0.616047  11.00000  0.37995
N029  3     0.719715  0.251040  0.615233  11.00000  0.37995
C097  1     0.750032  0.231374  0.623131  11.00000  0.37995
C098  1     0.810283  0.247370  0.621979  11.00000  0.37995
C099  1     0.756392  0.191033  0.635877  11.00000  0.37995
C100  1     0.816566  0.207334  0.634623  11.00000  0.37995
C101  1     0.789748  0.179655  0.641436  11.00000  0.37995
H145  2     0.735474  0.382340  0.629335  11.00000  0.19428
H146  2     0.797614  0.344323  0.660622  11.00000  0.19745
H147  2     0.802416  0.329735  0.575574  11.00000  0.37995
H148  2     0.762036  0.354569  0.569777  11.00000  0.37995
H149  2     0.707174  0.313463  0.595302  11.00000  0.37995
H150  2     0.696064  0.238164  0.617708  11.00000  0.37995
H151  2     0.831124  0.268504  0.616793  11.00000  0.37995
H152  2     0.735780  0.169593  0.641165  11.00000  0.37995
H153  2     0.842592  0.197241  0.639334  11.00000  0.37995
H154  2     0.795546  0.148611  0.651300  11.00000  0.37995
N030  3     0.803949  0.435038  0.625403  11.00000  0.19580
C102  1     0.826910  0.469048  0.614725  11.00000  0.22139
C103  1     0.840217  0.465562  0.566015  11.00000  0.21847
O029  4     0.869333  0.453664  0.557168  11.00000  0.21505
C104  1     0.856409  0.471281  0.648825  11.00000  0.37995
C105  1     0.845350  0.481619  0.696587  11.00000  0.37995
C106  1     0.839067  0.521838  0.709298  11.00000  0.37995
C107  1     0.841264  0.451217  0.728992  11.00000  0.37995
C108  1     0.828979  0.531441  0.752924  11.00000  0.37995
C109  1     0.831175  0.460515  0.772686  11.00000  0.37995
C110  1     0.825071  0.500673  0.784482  11.00000  0.37995
O030  4     0.815033  0.510062  0.827938  11.00000  0.37995
H155  2     0.780374  0.444611  0.635639  11.00000  0.19580
H156  2     0.812939  0.497798  0.617403  11.00000  0.22139
H157  2     0.874263  0.494831  0.637775  11.00000  0.37995
H158  2     0.869154  0.441675  0.649842  11.00000  0.37995
H159  2     0.842107  0.546000  0.684689  11.00000  0.37995
H160  2     0.846065  0.419746  0.719908  11.00000  0.37995
H161  2     0.824177  0.562852  0.762245  11.00000  0.37995
H162  2     0.828085  0.436475  0.797464  11.00000  0.37995
H163  2     0.812913  0.539760  0.831396  11.00000  0.37995
N031  3     0.817715  0.476480  0.534524  11.00000  0.20138
C111  1     0.826041  0.475196  0.486356  11.00000  0.21657
C112  1     0.814548  0.514222  0.462493  11.00000  0.21379
O031  4     0.786121  0.528903  0.469882  11.00000  0.22405
C113  1     0.809644  0.437852  0.463950  11.00000  0.37995
S004  5     0.826322  0.389283  0.483204  11.00000  0.37995
O032  4     0.833448  0.531778  0.435409  11.00000  0.37995
H164  2     0.793528  0.486359  0.541270  11.00000  0.20138
H165  2     0.853651  0.472627  0.482526  11.00000  0.21657
H166  2     0.782341  0.438219  0.471170  11.00000  0.37995
H167  2     0.813910  0.439656  0.427477  11.00000  0.37995
H168  2     0.845810  0.395431  0.518220  11.00000  0.37995
HKLF 4

END
'''

if (__name__ == "__main__"):
  t0 = time.time()
  if verbose: print('test 1')
  run(olex_res_str_1, 'tst_olex2_vs_pydiscamb_TAAM_1_1.mtz',
    target_scores = (0.0009, 0.0005, 0.001))
  if verbose: print('test 2')
  run(olex_res_str_2, 'tst_olex2_vs_pydiscamb_TAAM_1_2.mtz',
    target_scores = (0.005, 0.0005, 0.001))
  if verbose: print('test 3')
  run(olex_res_str_3, 'tst_olex2_vs_pydiscamb_TAAM_1_3.mtz',
    target_scores = (0.0007, 0.0005, 0.006))
  print("OK. Time: %8.3f"%(time.time()-t0))
