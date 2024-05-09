
from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx import easy_run
import libtbx.load_env # import dependency
import os
from six.moves import zip

pdb_str = """
CRYST1   34.896   15.928   22.972  90.00  90.00  90.00 P 1
SCALE1      0.028657  0.000000  0.000000        0.00000
SCALE2      0.000000  0.062783  0.000000        0.00000
SCALE3      0.000000  0.000000  0.043531        0.00000
ATOM      1  N   GLY A   1       7.241  10.342  12.931  1.00 16.37           N
ATOM      2  CA  GLY A   1       7.248   9.899  11.542  1.00 16.17           C
ATOM      3  C   GLY A   1       8.297   8.826  11.297  1.00 15.46           C
ATOM      4  O   GLY A   1       8.781   8.201  12.236  1.00 16.71           O
ATOM      5  H1  GLY A   1       6.544  10.879  13.068  1.00 16.37           H
ATOM      6  H2  GLY A   1       7.180   9.635  13.468  1.00 16.37           H
ATOM      7  H3  GLY A   1       7.994  10.783  13.105  1.00 16.37           H
ATOM      8  HA2 GLY A   1       6.377   9.539  11.311  1.00 16.17           H
ATOM      9  HA3 GLY A   1       7.435  10.653  10.961  1.00 16.17           H
ATOM     10  N  AASN A   2       8.625   8.586  10.033  0.60 14.94           N
ATOM     11  CA AASN A   2       9.755   7.720   9.710  0.60 14.28           C
ATOM     12  C  AASN A   2      11.046   8.226  10.317  0.60 12.87           C
ATOM     13  O  AASN A   2      11.288   9.431  10.370  0.60 12.55           O
ATOM     14  CB AASN A   2       9.932   7.588   8.206  0.60 15.48           C
ATOM     15  CG AASN A   2       8.788   6.866   7.558  0.60 13.90           C
ATOM     16  OD1AASN A   2       8.236   7.330   6.561  0.60 17.99           O
ATOM     17  ND2AASN A   2       8.421   5.716   8.117  0.60 11.61           N
ATOM     18  H  AASN A   2       8.215   8.908   9.349  0.60 15.11           H
ATOM     19  HA AASN A   2       9.587   6.837  10.074  0.60 13.97           H
ATOM     20  HB2AASN A   2       9.990   8.474   7.814  0.60 15.07           H
ATOM     21  HB3AASN A   2      10.744   7.090   8.025  0.60 15.07           H
ATOM     22 HD21AASN A   2       7.719   5.306   7.835  0.60 11.52           H
ATOM     23 HD22AASN A   2       8.884   5.382   8.759  0.60 11.52           H
ATOM     24  N  BASN A   2       8.662   8.607  10.040  0.40 15.11           N
ATOM     25  CA BASN A   2       9.841   7.793   9.758  0.40 13.97           C
ATOM     26  C  BASN A   2      11.079   8.439  10.353  0.40 13.04           C
ATOM     27  O  BASN A   2      11.211   9.661  10.344  0.40 11.55           O
ATOM     28  CB BASN A   2      10.038   7.614   8.259  0.40 15.07           C
ATOM     29  CG BASN A   2       8.876   6.915   7.606  0.40 13.81           C
ATOM     30  OD1BASN A   2       8.335   7.391   6.608  0.40 17.17           O
ATOM     31  ND2BASN A   2       8.483   5.775   8.163  0.40 11.52           N
ATOM     32  H  BASN A   2       8.256   8.910   9.345  0.40 15.11           H
ATOM     33  HA BASN A   2       9.732   6.917  10.160  0.40 13.97           H
ATOM     34  HB2BASN A   2      10.135   8.486   7.846  0.40 15.07           H
ATOM     35  HB3BASN A   2      10.834   7.082   8.105  0.40 15.07           H
ATOM     36 HD21BASN A   2       7.804   5.353   7.845  0.40 11.52           H
ATOM     37 HD22BASN A   2       8.907   5.459   8.841  0.40 11.52           H
ATOM     38  N  AASN A   3      11.874   7.290  10.763  0.60 12.35           N
ATOM     39  CA AASN A   3      13.106   7.616  11.464  0.60 11.68           C
ATOM     40  C  AASN A   3      14.355   7.054  10.788  0.60 11.09           C
ATOM     41  O  AASN A   3      14.430   5.850  10.518  0.60 10.18           O
ATOM     42  CB AASN A   3      13.029   7.115  12.915  0.60 11.62           C
ATOM     43  CG AASN A   3      14.277   7.454  13.717  0.60 13.00           C
ATOM     44  OD1AASN A   3      14.584   8.629  13.935  0.60 15.44           O
ATOM     45  ND2AASN A   3      15.001   6.428  14.160  0.60 12.87           N
ATOM     46  H  AASN A   3      11.741   6.445  10.670  0.60 12.00           H
ATOM     47  HA AASN A   3      13.210   8.579  11.499  0.60 11.53           H
ATOM     48  HB2AASN A   3      12.269   7.529  13.353  0.60 12.90           H
ATOM     49  HB3AASN A   3      12.926   6.150  12.912  0.60 12.90           H
ATOM     50 HD21AASN A   3      15.715   6.571  14.618  0.60 12.97           H
ATOM     51 HD22AASN A   3      14.755   5.622  13.988  0.60 12.97           H
ATOM     52  N  BASN A   3      11.978   7.616  10.880  0.40 12.00           N
ATOM     53  CA BASN A   3      13.210   8.114  11.479  0.40 11.53           C
ATOM     54  C  BASN A   3      14.444   7.369  10.976  0.40 11.70           C
ATOM     55  O  BASN A   3      14.478   6.136  10.965  0.40 10.95           O
ATOM     56  CB BASN A   3      13.136   8.043  13.007  0.40 12.90           C
ATOM     57  CG BASN A   3      12.214   9.088  13.597  0.40 13.32           C
ATOM     58  OD1BASN A   3      11.921  10.106  12.964  0.40 15.04           O
ATOM     59  ND2BASN A   3      11.760   8.851  14.823  0.40 12.97           N
ATOM     60  H  BASN A   3      11.898   6.760  10.903  0.40 12.00           H
ATOM     61  HA BASN A   3      13.332   9.044  11.235  0.40 11.53           H
ATOM     62  HB2BASN A   3      12.805   7.169  13.267  0.40 12.90           H
ATOM     63  HB3BASN A   3      14.023   8.184  13.373  0.40 12.90           H
ATOM     64 HD21BASN A   3      11.233   9.414  15.203  0.40 12.97           H
ATOM     65 HD22BASN A   3      11.993   8.134  15.237  0.40 12.97           H
ATOM     66  N  AGLN A   4      15.318   7.938  10.514  0.60 10.54           N
ATOM     67  CA AGLN A   4      16.673   7.540  10.113  0.60 10.42           C
ATOM     68  C  AGLN A   4      17.726   8.190  11.019  0.60 11.71           C
ATOM     69  O  AGLN A   4      17.880   9.450  10.993  0.60  8.92           O
ATOM     70  CB AGLN A   4      16.959   7.869   8.636  0.60  9.88           C
ATOM     71  CG AGLN A   4      18.209   7.154   8.093  0.60 10.02           C
ATOM     72  CD AGLN A   4      18.769   7.756   6.808  0.60 12.86           C
ATOM     73  OE1AGLN A   4      19.004   8.968   6.722  0.60 14.16           O
ATOM     74  NE2AGLN A   4      19.007   6.903   5.806  0.60  9.04           N
ATOM     75  H  AGLN A   4      15.210   8.790  10.554  0.60 10.83           H
ATOM     76  HA AGLN A   4      16.760   6.580  10.217  0.60 10.25           H
ATOM     77  HB2AGLN A   4      16.200   7.593   8.099  0.60  9.74           H
ATOM     78  HB3AGLN A   4      17.098   8.825   8.547  0.60  9.74           H
ATOM     79  HG2AGLN A   4      18.908   7.193   8.764  0.60 10.05           H
ATOM     80  HG3AGLN A   4      17.982   6.229   7.910  0.60 10.05           H
ATOM     81 HE21AGLN A   4      18.845   6.065   5.907  0.60  8.91           H
ATOM     82 HE22AGLN A   4      19.322   7.193   5.060  0.60  8.91           H
ATOM     83  N  BGLN A   4      15.453   8.133  10.564  0.40 10.83           N
ATOM     84  CA BGLN A   4      16.716   7.571  10.098  0.40 10.25           C
ATOM     85  C  BGLN A   4      17.888   8.129  10.902  0.40 10.29           C
ATOM     86  O  BGLN A   4      18.321   9.259  10.691  0.40 10.50           O
ATOM     87  CB BGLN A   4      16.921   7.845   8.603  0.40  9.74           C
ATOM     88  CG BGLN A   4      18.096   7.084   7.994  0.40 10.05           C
ATOM     89  CD BGLN A   4      18.647   7.745   6.745  0.40 12.79           C
ATOM     90  OE1BGLN A   4      18.791   8.970   6.690  0.40 15.08           O
ATOM     91  NE2BGLN A   4      18.962   6.936   5.730  0.40  8.91           N
ATOM     92  H  BGLN A   4      15.429   8.992  10.545  0.40 10.83           H
ATOM     93  HA BGLN A   4      16.710   6.610  10.222  0.40 10.25           H
ATOM     94  HB2BGLN A   4      16.119   7.581   8.125  0.40  9.74           H
ATOM     95  HB3BGLN A   4      17.079   8.794   8.480  0.40  9.74           H
ATOM     96  HG2BGLN A   4      18.814   7.031   8.643  0.40 10.05           H
ATOM     97  HG3BGLN A   4      17.801   6.192   7.754  0.40 10.05           H
ATOM     98 HE21BGLN A   4      18.849   6.087   5.806  0.40  8.91           H
ATOM     99 HE22BGLN A   4      19.278   7.264   5.000  0.40  8.91           H
ATOM    100  N  AGLN A   5      18.413   7.359  11.856  0.60 10.59           N
ATOM    101  CA AGLN A   5      19.562   7.932  12.524  0.60 11.43           C
ATOM    102  C  AGLN A   5      20.870   7.440  11.938  0.60 11.24           C
ATOM    103  O  AGLN A   5      21.037   6.253  11.692  0.60 11.99           O
ATOM    104  CB AGLN A   5      19.517   7.784  14.035  0.60 12.07           C
ATOM    105  CG AGLN A   5      18.523   8.708  14.724  0.60 10.78           C
ATOM    106  CD AGLN A   5      17.850   8.030  15.899  0.60 12.94           C
ATOM    107  OE1AGLN A   5      17.293   6.941  15.764  0.60 10.72           O
ATOM    108  NE2AGLN A   5      17.909   8.667  17.062  0.60 12.32           N
ATOM    109  H  AGLN A   5      18.236   6.535  12.028  0.60 10.59           H
ATOM    110  HA AGLN A   5      19.567   8.887  12.374  0.60 11.43           H
ATOM    111  HB2AGLN A   5      19.275   6.870  14.249  0.60 12.07           H
ATOM    112  HB3AGLN A   5      20.397   7.981  14.392  0.60 12.07           H
ATOM    113  HG2AGLN A   5      18.990   9.492  15.053  0.60 10.78           H
ATOM    114  HG3AGLN A   5      17.836   8.968  14.090  0.60 10.78           H
ATOM    115 HE21AGLN A   5      18.313   9.424  17.118  0.60 12.32           H
ATOM    116 HE22AGLN A   5      17.543   8.323  17.760  0.60 12.32           H
ATOM    117  N  BGLN A   5      18.413   7.359  11.856  0.40 10.59           N
ATOM    118  CA BGLN A   5      19.562   7.932  12.524  0.40 11.43           C
ATOM    119  C  BGLN A   5      20.870   7.440  11.938  0.40 11.24           C
ATOM    120  O  BGLN A   5      21.037   6.253  11.692  0.40 11.99           O
ATOM    121  CB BGLN A   5      19.517   7.784  14.035  0.40 12.07           C
ATOM    122  CG BGLN A   5      18.523   8.708  14.724  0.40 10.78           C
ATOM    123  CD BGLN A   5      17.850   8.030  15.899  0.40 12.94           C
ATOM    124  OE1BGLN A   5      17.293   6.941  15.764  0.40 10.72           O
ATOM    125  NE2BGLN A   5      17.909   8.667  17.062  0.40 12.32           N
ATOM    126  H  BGLN A   5      18.151   6.581  12.111  0.40 10.59           H
ATOM    127  HA BGLN A   5      19.576   8.890  12.403  0.40 11.43           H
ATOM    128  HB2BGLN A   5      19.275   6.870  14.249  0.40 12.07           H
ATOM    129  HB3BGLN A   5      20.397   7.981  14.392  0.40 12.07           H
ATOM    130  HG2BGLN A   5      18.990   9.492  15.053  0.40 10.78           H
ATOM    131  HG3BGLN A   5      17.836   8.968  14.090  0.40 10.78           H
ATOM    132 HE21BGLN A   5      18.313   9.424  17.118  0.40 12.32           H
ATOM    133 HE22BGLN A   5      17.543   8.323  17.760  0.40 12.32           H
ATOM    134  N   ASN A   6      21.796   8.371  11.716  1.00 11.72           N
ATOM    135  CA  ASN A   6      23.113   8.030  11.190  1.00 12.12           C
ATOM    136  C   ASN A   6      24.142   8.471  12.195  1.00 13.15           C
ATOM    137  O   ASN A   6      24.509   9.645  12.245  1.00 13.93           O
ATOM    138  CB  ASN A   6      23.349   8.721   9.859  1.00 11.96           C
ATOM    139  CG  ASN A   6      22.251   8.440   8.870  1.00 12.58           C
ATOM    140  OD1 ASN A   6      22.087   7.312   8.414  1.00 14.01           O
ATOM    141  ND2 ASN A   6      21.480   9.459   8.544  1.00  9.96           N
ATOM    142  HA  ASN A   6      23.207   7.075  11.058  1.00 12.12           H
ATOM    143  HB2 ASN A   6      23.387   9.680   9.998  1.00 11.96           H
ATOM    144  HB3 ASN A   6      24.184   8.405   9.480  1.00 11.96           H
ATOM    145 HD21 ASN A   6      20.839   9.350   7.981  1.00  9.96           H
ATOM    146 HD22 ASN A   6      21.618  10.232   8.894  1.00  9.96           H
ATOM    147  H  AASN A   6      21.684   9.211  11.863  0.60 11.72           H
ATOM    148  H  BASN A   6      21.684   9.211  11.863  0.40 11.72           H
ATOM    149  N   TYR A   7      24.585   7.530  13.020  1.00 14.62           N
ATOM    150  CA  TYR A   7      25.450   7.854  14.156  1.00 15.04           C
ATOM    151  C   TYR A   7      26.899   8.037  13.753  1.00 15.56           C
ATOM    152  O   TYR A   7      27.334   7.518  12.719  1.00 15.52           O
ATOM    153  CB  TYR A   7      25.344   6.780  15.235  1.00 14.86           C
ATOM    154  CG  TYR A   7      23.945   6.649  15.763  1.00 14.38           C
ATOM    155  CD1 TYR A   7      23.055   5.738  15.199  1.00 15.46           C
ATOM    156  CD2 TYR A   7      23.494   7.461  16.795  1.00 14.37           C
ATOM    157  CE1 TYR A   7      21.764   5.619  15.667  1.00 13.24           C
ATOM    158  CE2 TYR A   7      22.193   7.352  17.274  1.00 13.84           C
ATOM    159  CZ  TYR A   7      21.337   6.430  16.701  1.00 14.98           C
ATOM    160  OH  TYR A   7      20.055   6.301  17.164  1.00 14.28           O
ATOM    161  OXT TYR A   7      27.649   8.711  14.468  1.00 17.34           O
ATOM    162  H   TYR A   7      24.401   6.693  12.946  1.00 14.62           H
ATOM    163  HA  TYR A   7      25.148   8.689  14.546  1.00 15.04           H
ATOM    164  HB2 TYR A   7      25.607   5.925  14.860  1.00 14.86           H
ATOM    165  HB3 TYR A   7      25.926   7.014  15.975  1.00 14.86           H
ATOM    166  HD1 TYR A   7      23.341   5.191  14.503  1.00 15.46           H
ATOM    167  HD2 TYR A   7      24.073   8.080  17.178  1.00 14.37           H
ATOM    168  HE1 TYR A   7      21.183   5.000  15.287  1.00 13.24           H
ATOM    169  HE2 TYR A   7      21.902   7.893  17.972  1.00 13.84           H
ATOM    170  HH  TYR A   7      19.851   6.969  17.629  1.00 14.28           H
TER
HETATM  171  O   HOH S   8       9.815  10.928  13.987  1.00 22.61           O
HETATM  172  O   HOH S   9      26.715   7.572  10.077  1.00 19.32           O
HETATM  173  O   HOH S  10       5.000   7.470   5.401  1.00 16.97           O
HETATM  174  O   HOH S  11      28.091   9.896  16.830  1.00 23.89           O
HETATM  175  O   HOH S  12      29.896   7.023  16.061  1.00 26.08           O
HETATM  176  O   HOH S  13      13.552   9.160  16.880  1.00 38.68           O
HETATM  177  O   HOH S  14      14.793   6.375  17.843  1.00 44.24           O
TER
END
"""

def exercise():
  pdb_file = "tmp_ringer.pdb"
  f = open(pdb_file, "w")
  f.write(pdb_str)
  f.close()
  mtz_file = "tmp_ringer.mtz"
  cmd = " ".join([
    "phenix.fmodel",
    pdb_file,
    "high_resolution=2.0",
    "type=real",
    "r_free_flags_fraction=0.1",
    "random_seed=12345",
    "label=F",
    "output.file_name=%s" % mtz_file,
  ])
  print(cmd)
  assert not easy_run.call(cmd)
  result = easy_run.fully_buffered(
    "phenix.maps \"%s\" \"%s\" output.prefix=tmp_ringer" %
    (pdb_file, mtz_file)).raise_if_errors()
  assert (result.return_code == 0)
  result = easy_run.fully_buffered(
    "mmtbx.ringer \"%s\" tmp_ringer_map_coeffs.mtz" % pdb_file).raise_if_errors()
  with open("tmp_ringer_ringer.csv") as f:
    _lines1 = f.read().splitlines()
  lines1 = []
  for line in _lines1 :
    if ("2mFo-DFc" in line):
      lines1.append(line)
  os.remove("tmp_ringer_ringer.csv")
  assert (result.return_code == 0)
  # Now with ccp4 map as input
  result2 = easy_run.fully_buffered(
    "phenix.mtz2map tmp_ringer_map_coeffs.mtz")
  assert (result2.return_code == 0)
  result3 = easy_run.fully_buffered(
    "mmtbx.ringer \"%s\" tmp_ringer_map_coeffs_2mFo-DFc.ccp4" % pdb_file)
  assert (result3.return_code == 0) , "DL: crash is expected due to new sanity check. Will replace input map."
  with open("tmp_ringer_ringer.csv") as f:
    lines2 = f.read().splitlines()
  assert len(lines1) == len(lines2)
  for line1, line2 in zip(lines1, lines2):
    fields1 = line1.split(",")
    fields2 = line2.split(",")
    rho1 = flex.double([ float(x) for x in fields1[4:] ])
    rho2 = flex.double([ float(x) for x in fields2[4:] ])
    cc = flex.linear_correlation(x=rho1, y=rho2).coefficient()
    assert (cc >= 0.99), cc

if (__name__ == "__main__"):
  exercise()
  print("OK")
