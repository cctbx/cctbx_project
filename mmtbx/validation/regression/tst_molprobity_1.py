
from __future__ import absolute_import, division, print_function
from mmtbx.command_line import molprobity
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO
from six.moves import zip

# derivative of 1yjp, with alternate conformation Asn3 in different rotamer,
# plus Asn2/Gln4 split.  the unit cell 'b' edge has been increased to 6A to
# compensate.
pdb_raw_0 = """\
CRYST1   21.937    6.000   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.166667  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.047   4.634   6.066  1.00 16.37           N
ATOM      2  CA  GLY A   1      -9.040   4.191   4.677  1.00 16.17           C
ATOM      3  C   GLY A   1      -7.991   3.118   4.432  1.00 15.46           C
ATOM      4  O   GLY A   1      -7.507   2.493   5.371  1.00 16.71           O
ATOM      5  H1  GLY A   1      -9.744   5.171   6.203  1.00 16.37           H
ATOM      6  H2  GLY A   1      -9.108   3.927   6.603  1.00 16.37           H
ATOM      7  H3  GLY A   1      -8.294   5.075   6.240  1.00 16.37           H
ATOM      8  HA2 GLY A   1      -9.911   3.831   4.446  1.00 16.17           H
ATOM      9  HA3 GLY A   1      -8.853   4.945   4.096  1.00 16.17           H
ATOM     10  N  AASN A   2      -7.663   2.878   3.168  0.60 14.94           N
ATOM     11  CA AASN A   2      -6.533   2.012   2.845  0.60 14.28           C
ATOM     12  C  AASN A   2      -5.242   2.518   3.452  0.60 12.87           C
ATOM     13  O  AASN A   2      -5.000   3.723   3.505  0.60 12.55           O
ATOM     14  CB AASN A   2      -6.356   1.880   1.341  0.60 15.48           C
ATOM     15  CG AASN A   2      -7.500   1.158   0.693  0.60 13.90           C
ATOM     16  OD1AASN A   2      -8.052   1.622  -0.304  0.60 17.99           O
ATOM     17  ND2AASN A   2      -7.867   0.008   1.252  0.60 11.61           N
ATOM     18  H  AASN A   2      -8.073   3.200   2.484  0.60 15.11           H
ATOM     19  HA AASN A   2      -6.701   1.129   3.209  0.60 13.97           H
ATOM     20  HB2AASN A   2      -6.298   2.766   0.949  0.60 15.07           H
ATOM     21  HB3AASN A   2      -5.544   1.382   1.160  0.60 15.07           H
ATOM     22 HD21AASN A   2      -8.569  -0.402   0.970  0.60 11.52           H
ATOM     23 HD22AASN A   2      -7.404  -0.326   1.894  0.60 11.52           H
ATOM     24  N  BASN A   2      -7.626   2.899   3.175  0.40 15.11           N
ATOM     25  CA BASN A   2      -6.447   2.085   2.893  0.40 13.97           C
ATOM     26  C  BASN A   2      -5.209   2.731   3.488  0.40 13.04           C
ATOM     27  O  BASN A   2      -5.077   3.953   3.479  0.40 11.55           O
ATOM     28  CB BASN A   2      -6.250   1.906   1.394  0.40 15.07           C
ATOM     29  CG BASN A   2      -7.412   1.207   0.741  0.40 13.81           C
ATOM     30  OD1BASN A   2      -7.953   1.683  -0.257  0.40 17.17           O
ATOM     31  ND2BASN A   2      -7.805   0.067   1.298  0.40 11.52           N
ATOM     32  H  BASN A   2      -8.032   3.202   2.480  0.40 15.11           H
ATOM     33  HA BASN A   2      -6.556   1.209   3.295  0.40 13.97           H
ATOM     34  HB2BASN A   2      -6.153   2.778   0.981  0.40 15.07           H
ATOM     35  HB3BASN A   2      -5.454   1.374   1.240  0.40 15.07           H
ATOM     36 HD21BASN A   2      -8.484  -0.355   0.980  0.40 11.52           H
ATOM     37 HD22BASN A   2      -7.381  -0.249   1.976  0.40 11.52           H
ATOM     38  N  AASN A   3      -4.414   1.582   3.898  0.60 12.35           N
ATOM     39  CA AASN A   3      -3.182   1.908   4.599  0.60 11.68           C
ATOM     40  C  AASN A   3      -1.933   1.346   3.923  0.60 11.09           C
ATOM     41  O  AASN A   3      -1.858   0.142   3.653  0.60 10.18           O
ATOM     42  CB AASN A   3      -3.259   1.407   6.050  0.60 11.62           C
ATOM     43  CG AASN A   3      -2.011   1.746   6.852  0.60 13.00           C
ATOM     44  OD1AASN A   3      -1.704   2.921   7.070  0.60 15.44           O
ATOM     45  ND2AASN A   3      -1.287   0.720   7.295  0.60 12.87           N
ATOM     46  H  AASN A   3      -4.547   0.737   3.805  0.60 12.00           H
ATOM     47  HA AASN A   3      -3.078   2.871   4.634  0.60 11.53           H
ATOM     48  HB2AASN A   3      -4.019   1.821   6.488  0.60 12.90           H
ATOM     49  HB3AASN A   3      -3.362   0.442   6.047  0.60 12.90           H
ATOM     50 HD21AASN A   3      -0.573   0.863   7.753  0.60 12.97           H
ATOM     51 HD22AASN A   3      -1.533  -0.086   7.123  0.60 12.97           H
ATOM     52  N  BASN A   3      -4.310   1.908   4.015  0.40 12.00           N
ATOM     53  CA BASN A   3      -3.078   2.406   4.614  0.40 11.53           C
ATOM     54  C  BASN A   3      -1.844   1.661   4.111  0.40 11.70           C
ATOM     55  O  BASN A   3      -1.810   0.428   4.100  0.40 10.95           O
ATOM     56  CB BASN A   3      -3.152   2.335   6.142  0.40 12.90           C
ATOM     57  CG BASN A   3      -4.074   3.380   6.732  0.40 13.32           C
ATOM     58  OD1BASN A   3      -4.367   4.398   6.099  0.40 15.04           O
ATOM     59  ND2BASN A   3      -4.528   3.143   7.958  0.40 12.97           N
ATOM     60  H  BASN A   3      -4.390   1.052   4.038  0.40 12.00           H
ATOM     61  HA BASN A   3      -2.956   3.336   4.370  0.40 11.53           H
ATOM     62  HB2BASN A   3      -3.483   1.461   6.402  0.40 12.90           H
ATOM     63  HB3BASN A   3      -2.265   2.476   6.508  0.40 12.90           H
ATOM     64 HD21BASN A   3      -5.055   3.706   8.338  0.40 12.97           H
ATOM     65 HD22BASN A   3      -4.295   2.426   8.372  0.40 12.97           H
ATOM     66  N  AGLN A   4      -0.970   2.230   3.649  0.60 10.54           N
ATOM     67  CA AGLN A   4       0.385   1.832   3.248  0.60 10.42           C
ATOM     68  C  AGLN A   4       1.438   2.482   4.154  0.60 11.71           C
ATOM     69  O  AGLN A   4       1.592   3.742   4.128  0.60  8.92           O
ATOM     70  CB AGLN A   4       0.671   2.161   1.771  0.60  9.88           C
ATOM     71  CG AGLN A   4       1.921   1.446   1.228  0.60 10.02           C
ATOM     72  CD AGLN A   4       2.481   2.048  -0.057  0.60 12.86           C
ATOM     73  OE1AGLN A   4       2.716   3.260  -0.143  0.60 14.16           O
ATOM     74  NE2AGLN A   4       2.719   1.195  -1.059  0.60  9.04           N
ATOM     75  H  AGLN A   4      -1.078   3.082   3.689  0.60 10.83           H
ATOM     76  HA AGLN A   4       0.472   0.872   3.352  0.60 10.25           H
ATOM     77  HB2AGLN A   4      -0.088   1.885   1.234  0.60  9.74           H
ATOM     78  HB3AGLN A   4       0.810   3.117   1.682  0.60  9.74           H
ATOM     79  HG2AGLN A   4       2.620   1.485   1.899  0.60 10.05           H
ATOM     80  HG3AGLN A   4       1.694   0.521   1.045  0.60 10.05           H
ATOM     81 HE21AGLN A   4       2.557   0.357  -0.958  0.60  8.91           H
ATOM     82 HE22AGLN A   4       3.034   1.485  -1.805  0.60  8.91           H
ATOM     83  N  BGLN A   4      -0.835   2.425   3.699  0.40 10.83           N
ATOM     84  CA BGLN A   4       0.428   1.863   3.233  0.40 10.25           C
ATOM     85  C  BGLN A   4       1.600   2.421   4.037  0.40 10.29           C
ATOM     86  O  BGLN A   4       2.033   3.551   3.826  0.40 10.50           O
ATOM     87  CB BGLN A   4       0.633   2.137   1.738  0.40  9.74           C
ATOM     88  CG BGLN A   4       1.808   1.376   1.129  0.40 10.05           C
ATOM     89  CD BGLN A   4       2.359   2.037  -0.120  0.40 12.79           C
ATOM     90  OE1BGLN A   4       2.503   3.262  -0.175  0.40 15.08           O
ATOM     91  NE2BGLN A   4       2.674   1.228  -1.135  0.40  8.91           N
ATOM     92  H  BGLN A   4      -0.859   3.284   3.680  0.40 10.83           H
ATOM     93  HA BGLN A   4       0.422   0.902   3.357  0.40 10.25           H
ATOM     94  HB2BGLN A   4      -0.169   1.873   1.260  0.40  9.74           H
ATOM     95  HB3BGLN A   4       0.791   3.086   1.615  0.40  9.74           H
ATOM     96  HG2BGLN A   4       2.526   1.323   1.778  0.40 10.05           H
ATOM     97  HG3BGLN A   4       1.513   0.484   0.889  0.40 10.05           H
ATOM     98 HE21BGLN A   4       2.561   0.379  -1.059  0.40  8.91           H
ATOM     99 HE22BGLN A   4       2.990   1.556  -1.865  0.40  8.91           H
ATOM    100  N  AGLN A   5       2.125   1.651   4.991  0.60 10.59           N
ATOM    101  CA AGLN A   5       3.274   2.224   5.659  0.60 11.43           C
ATOM    102  C  AGLN A   5       4.582   1.732   5.073  0.60 11.24           C
ATOM    103  O  AGLN A   5       4.749   0.545   4.827  0.60 11.99           O
ATOM    104  CB AGLN A   5       3.229   2.076   7.170  0.60 12.07           C
ATOM    105  CG AGLN A   5       2.235   3.000   7.859  0.60 10.78           C
ATOM    106  CD AGLN A   5       1.562   2.322   9.034  0.60 12.94           C
ATOM    107  OE1AGLN A   5       1.005   1.233   8.899  0.60 10.72           O
ATOM    108  NE2AGLN A   5       1.621   2.959  10.197  0.60 12.32           N
ATOM    109  H  AGLN A   5       1.948   0.827   5.163  0.60 10.59           H
ATOM    110  HA AGLN A   5       3.279   3.179   5.509  0.60 11.43           H
ATOM    111  HB2AGLN A   5       2.987   1.162   7.384  0.60 12.07           H
ATOM    112  HB3AGLN A   5       4.109   2.273   7.527  0.60 12.07           H
ATOM    113  HG2AGLN A   5       2.702   3.784   8.188  0.60 10.78           H
ATOM    114  HG3AGLN A   5       1.548   3.260   7.225  0.60 10.78           H
ATOM    115 HE21AGLN A   5       2.025   3.716  10.253  0.60 12.32           H
ATOM    116 HE22AGLN A   5       1.255   2.615  10.895  0.60 12.32           H
ATOM    117  N  BGLN A   5       2.125   1.651   4.991  0.40 10.59           N
ATOM    118  CA BGLN A   5       3.274   2.224   5.659  0.40 11.43           C
ATOM    119  C  BGLN A   5       4.582   1.732   5.073  0.40 11.24           C
ATOM    120  O  BGLN A   5       4.749   0.545   4.827  0.40 11.99           O
ATOM    121  CB BGLN A   5       3.229   2.076   7.170  0.40 12.07           C
ATOM    122  CG BGLN A   5       2.235   3.000   7.859  0.40 10.78           C
ATOM    123  CD BGLN A   5       1.562   2.322   9.034  0.40 12.94           C
ATOM    124  OE1BGLN A   5       1.005   1.233   8.899  0.40 10.72           O
ATOM    125  NE2BGLN A   5       1.621   2.959  10.197  0.40 12.32           N
ATOM    126  H  BGLN A   5       1.863   0.873   5.246  0.40 10.59           H
ATOM    127  HA BGLN A   5       3.288   3.182   5.538  0.40 11.43           H
ATOM    128  HB2BGLN A   5       2.987   1.162   7.384  0.40 12.07           H
ATOM    129  HB3BGLN A   5       4.109   2.273   7.527  0.40 12.07           H
ATOM    130  HG2BGLN A   5       2.702   3.784   8.188  0.40 10.78           H
ATOM    131  HG3BGLN A   5       1.548   3.260   7.225  0.40 10.78           H
ATOM    132 HE21BGLN A   5       2.025   3.716  10.253  0.40 12.32           H
ATOM    133 HE22BGLN A   5       1.255   2.615  10.895  0.40 12.32           H
ATOM    134  N   ASN A   6       5.508   2.663   4.851  1.00 11.72           N
ATOM    135  CA  ASN A   6       6.825   2.322   4.325  1.00 12.12           C
ATOM    136  C   ASN A   6       7.854   2.763   5.330  1.00 13.15           C
ATOM    137  O   ASN A   6       8.221   3.937   5.380  1.00 13.93           O
ATOM    138  CB  ASN A   6       7.061   3.013   2.994  1.00 11.96           C
ATOM    139  CG  ASN A   6       5.963   2.732   2.005  1.00 12.58           C
ATOM    140  OD1 ASN A   6       5.799   1.604   1.549  1.00 14.01           O
ATOM    141  ND2 ASN A   6       5.192   3.751   1.679  1.00  9.96           N
ATOM    142  HA  ASN A   6       6.919   1.367   4.193  1.00 12.12           H
ATOM    143  HB2 ASN A   6       7.099   3.972   3.133  1.00 11.96           H
ATOM    144  HB3 ASN A   6       7.896   2.697   2.615  1.00 11.96           H
ATOM    145 HD21 ASN A   6       4.551   3.642   1.116  1.00  9.96           H
ATOM    146 HD22 ASN A   6       5.330   4.524   2.029  1.00  9.96           H
ATOM    147  H  AASN A   6       5.396   3.503   4.998  0.60 11.72           H
ATOM    148  H  BASN A   6       5.396   3.503   4.998  0.40 11.72           H
ATOM    149  N   TYR A   7       8.297   1.822   6.155  1.00 14.62           N
ATOM    150  CA  TYR A   7       9.162   2.146   7.291  1.00 15.04           C
ATOM    151  C   TYR A   7      10.611   2.329   6.888  1.00 15.56           C
ATOM    152  O   TYR A   7      11.046   1.810   5.854  1.00 15.52           O
ATOM    153  CB  TYR A   7       9.056   1.072   8.370  1.00 14.86           C
ATOM    154  CG  TYR A   7       7.657   0.941   8.898  1.00 14.38           C
ATOM    155  CD1 TYR A   7       6.767   0.030   8.334  1.00 15.46           C
ATOM    156  CD2 TYR A   7       7.206   1.753   9.930  1.00 14.37           C
ATOM    157  CE1 TYR A   7       5.476  -0.089   8.802  1.00 13.24           C
ATOM    158  CE2 TYR A   7       5.905   1.644  10.409  1.00 13.84           C
ATOM    159  CZ  TYR A   7       5.049   0.722   9.836  1.00 14.98           C
ATOM    160  OH  TYR A   7       3.767   0.593  10.299  1.00 14.28           O
ATOM    161  OXT TYR A   7      11.361   3.003   7.603  1.00 17.34           O
ATOM    162  H   TYR A   7       8.113   0.985   6.081  1.00 14.62           H
ATOM    163  HA  TYR A   7       8.860   2.981   7.681  1.00 15.04           H
ATOM    164  HB2 TYR A   7       9.319   0.217   7.995  1.00 14.86           H
ATOM    165  HB3 TYR A   7       9.638   1.306   9.110  1.00 14.86           H
ATOM    166  HD1 TYR A   7       7.053  -0.517   7.638  1.00 15.46           H
ATOM    167  HD2 TYR A   7       7.785   2.372  10.313  1.00 14.37           H
ATOM    168  HE1 TYR A   7       4.895  -0.708   8.422  1.00 13.24           H
ATOM    169  HE2 TYR A   7       5.614   2.185  11.107  1.00 13.84           H
ATOM    170  HH  TYR A   7       3.563   1.261  10.764  1.00 14.28           H
TER
HETATM  171  O   HOH S   8      -6.473   5.220   7.122  1.00 22.61           O
HETATM  172  O   HOH S   9      10.427   1.864   3.212  1.00 19.32           O
HETATM  173  O   HOH S  10     -11.288   1.762  -1.464  1.00 16.97           O
HETATM  174  O   HOH S  11      11.803   4.188   9.965  1.00 23.89           O
HETATM  175  O   HOH S  12      13.608   1.315   9.196  1.00 26.08           O
HETATM  176  O   HOH S  13      -2.736   3.452  10.015  1.00 38.68           O
HETATM  177  O   HOH S  14      -1.495   0.667  10.978  1.00 44.24           O
TER
END
"""

# test for corner cases (synthetic data okay)
def exercise_synthetic():
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_raw_0)
  xrs = pdb_in.xray_structure_simple()
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  flags = fc.resolution_filter(d_min=1.6).generate_r_free_flags()
  ls = fc.lone_set(other=flags)
  # case 1: no work set in high-res shell
  flags2 = ls.array(data=flex.bool(ls.size(), True))
  flags_all = flags.concatenate(other=flags2)
  mtz_out = fc.as_mtz_dataset(column_root_label="F")
  mtz_out.add_miller_array(flags_all, column_root_label="FreeR_flag")
  mtz_out.mtz_object().write("tst_molprobity_1.mtz")
  open("tst_molprobity_1.pdb", "w").write(pdb_raw_0)
  args = [
    "tst_molprobity_1.pdb",
    "tst_molprobity_1.mtz",
    "--kinemage",
    "--maps",
    "flags.clashscore=False",
    "flags.xtriage=True",
  ]
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  out = StringIO()
  result.show(out=out)
  # case 2: no test set in high-res shell
  flags2 = ls.array(data=flex.bool(ls.size(), False))
  flags_all = flags.concatenate(other=flags2)
  mtz_out = fc.as_mtz_dataset(column_root_label="F")
  mtz_out.add_miller_array(flags_all, column_root_label="FreeR_flag")
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  out = StringIO()
  result.show(out=out)
  # case 3: multi-MODEL structure
  # XXX This is not a very sophisticated test - it only ensures that the
  # program does not crash.  We need a test for expected output...
  hierarchy = pdb_in.construct_hierarchy()
  model2 = hierarchy.only_model().detached_copy()
  hierarchy.append_model(model2)
  hierarchy.models()[0].id = "1"
  hierarchy.models()[1].id = "2"
  open("tst_molprobity_multi_model.pdb", "w").write(hierarchy.as_pdb_string())
  args = [
    "tst_molprobity_multi_model.pdb",
    "tst_molprobity_1.mtz",
    "--kinemage",
    "--maps",
  ]
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  out = StringIO()
  result.show(out=out)
  # test rotamer distributions
  open("tst_molprobity_misc1.pdb", "w").write(pdb_raw_0)
  args = [
    "tst_molprobity_1.pdb",
    "rotamer_library=8000",
  ]
  out = StringIO()
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  result.show(outliers_only=False, out=out)

def exercise_cdl():
  pdb_raw = """
ATOM   1270  N   LEU A 199       6.903  55.119  -0.416  1.00 25.48           N
ATOM   1271  CA  LEU A 199       7.726  56.192  -0.941  1.00 25.93           C
ATOM   1272  C   LEU A 199       6.996  56.972  -2.047  1.00 26.39           C
ATOM   1273  O   LEU A 199       7.020  58.180  -2.064  1.00 25.38           O
ATOM   1274  CB  LEU A 199       9.033  55.633  -1.490  1.00 25.66           C
ATOM   1278  N   ARG A 200       6.361  56.258  -2.980  1.00 27.31           N
ATOM   1279  CA  ARG A 200       5.576  56.913  -3.993  1.00 28.53           C
ATOM   1280  C   ARG A 200       4.520  57.823  -3.397  1.00 27.54           C
ATOM   1281  O   ARG A 200       4.397  58.949  -3.851  1.00 27.50           O
ATOM   1282  CB  ARG A 200       4.933  55.879  -4.899  1.00 30.38           C
ATOM   1289  N   ALA A 201       3.790  57.365  -2.357  1.00 26.90           N
ATOM   1290  CA  ALA A 201       2.764  58.200  -1.713  1.00 26.49           C
ATOM   1291  C   ALA A 201       3.406  59.407  -1.045  1.00 26.59           C
ATOM   1292  O   ALA A 201       2.866  60.516  -1.082  1.00 26.58           O
ATOM   1293  CB  ALA A 201       1.959  57.412  -0.715  1.00 25.11           C
ATOM   1294  N   ARG A 202       4.566  59.205  -0.419  1.00 25.66           N
ATOM   1295  CA  ARG A 202       5.245  60.296   0.240  1.00 26.93           C
ATOM   1296  C   ARG A 202       5.676  61.346  -0.767  1.00 26.93           C
ATOM   1297  O   ARG A 202       5.555  62.541  -0.489  1.00 25.79           O
ATOM   1298  CB  ARG A 202       6.493  59.779   0.996  1.00 28.25           C
ATOM   1305  N   ILE A 203       6.154  60.912  -1.931  1.00 26.99           N
ATOM   1306  CA  ILE A 203       6.611  61.848  -2.965  1.00 27.49           C
ATOM   1307  C   ILE A 203       5.430  62.674  -3.480  1.00 28.29           C
ATOM   1308  O   ILE A 203       5.548  63.905  -3.624  1.00 27.82           O
ATOM   1309  CB  ILE A 203       7.322  61.125  -4.075  1.00 28.09           C
ATOM   1313  N   SER A 204       4.288  62.025  -3.678  1.00 27.96           N
ATOM   1314  CA  SER A 204       3.119  62.736  -4.184  1.00 28.04           C
ATOM   1315  C   SER A 204       2.683  63.793  -3.199  1.00 28.16           C
ATOM   1316  O   SER A 204       2.311  64.910  -3.605  1.00 28.25           O
ATOM   1317  CB  SER A 204       1.962  61.780  -4.504  1.00 27.64           C
"""
  pdb_raw_2 = """\
REMARK   3    GEOSTD + MON.LIB. + CDL v1.2
""" + pdb_raw
  open("tst_molprobity_cdl_1.pdb", "w").write(pdb_raw)
  open("tst_molprobity_cdl_2.pdb", "w").write(pdb_raw_2)
  files = ["tst_molprobity_cdl_1.pdb","tst_molprobity_cdl_2.pdb"]
  rmsds = [0.9019, 0.8769]
  for file_name, rmsd, cdl_expected in zip(files, rmsds, [False, True]):
    result = molprobity.run(args=[file_name, "flags.clashscore=False"],
      ignore_missing_modules=True,
      out=null_out()).validation
    assert approx_equal(result.rms_angles(), rmsd, eps=0.001), rmsd
    if cdl_expected :
      out = StringIO()
      result.show(out=out)
      assert ("Geometry Restraints Library: GeoStd + Monomer Library + CDL v1.2" in out.getvalue())
    else:
      out = StringIO()
      result.show(out=out)
      assert ("CDL" not in out.getvalue())

if (__name__ == "__main__"):
  exercise_cdl()
  exercise_synthetic()
  print("OK")
