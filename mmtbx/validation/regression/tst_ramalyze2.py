from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.validation import ramalyze

# ------------------------------------------------------------------------------
# Make sure that pdb file with several models produces results for all models
# ------------------------------------------------------------------------------

def exercise():
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp)
  pdb_hierarchy = model.get_hierarchy()

  r = ramalyze.ramalyze(
    pdb_hierarchy = pdb_hierarchy,
    outliers_only = False)

  assert (len(r.results) == 3), 'Supposed to fail until fixed. Ramalyze results not available for all models.'

pdb_str = """
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000
MODEL        1
ATOM      1  N   ASN A   1       1.329   0.000   0.000  1.00  1.00           N
ATOM      2  CA  ASN A   1       2.093  -0.001  -1.242  1.00 64.21           C
ATOM      3  C   ASN A   1       1.973  -1.345  -1.954  1.00 21.54           C
ATOM      4  O   ASN A   1       2.071  -1.423  -3.178  1.00 42.13           O
ATOM      5  CB  ASN A   1       3.565   0.309  -0.960  1.00 52.42           C
ATOM      6  CG  ASN A   1       4.305   0.774  -2.199  1.00 64.34           C
ATOM      7  OD1 ASN A   1       4.331   0.081  -3.217  1.00 14.30           O
ATOM      8  ND2 ASN A   1       4.913   1.952  -2.118  1.00 64.45           N
ATOM      9  H1  ASN A   1       1.808   0.001   0.855  1.00 15.23           H
ATOM     10  HA  ASN A   1       1.688   0.769  -1.881  1.00 62.21           H
ATOM     11  HB2 ASN A   1       3.625   1.089  -0.215  1.00  4.30           H
ATOM     12  HB3 ASN A   1       4.049  -0.580  -0.585  1.00 41.42           H
ATOM     13 HD21 ASN A   1       4.851   2.448  -1.276  1.00  0.23           H
ATOM     14 HD22 ASN A   1       5.399   2.276  -2.905  1.00 40.32           H
ATOM     15  N   VAL A   2       1.759  -2.403  -1.177  1.00  2.11           N
ATOM     16  CA  VAL A   2       1.623  -3.744  -1.732  1.00 65.20           C
ATOM     17  C   VAL A   2       0.180  -4.031  -2.128  1.00 74.53           C
ATOM     18  O   VAL A   2      -0.536  -4.748  -1.429  1.00 30.45           O
ATOM     19  CB  VAL A   2       2.092  -4.816  -0.730  1.00 64.25           C
ATOM     20  CG1 VAL A   2       2.069  -6.195  -1.373  1.00 50.44           C
ATOM     21  CG2 VAL A   2       3.483  -4.486  -0.209  1.00 13.30           C
ATOM     22  H   VAL A   2       1.690  -2.278  -0.207  1.00 43.35           H
ATOM     23  HA  VAL A   2       2.247  -3.806  -2.612  1.00 23.21           H
ATOM     24  HB  VAL A   2       1.409  -4.822   0.106  1.00 41.00           H
ATOM     25 HG11 VAL A   2       2.853  -6.803  -0.948  1.00 21.35           H
ATOM     26 HG12 VAL A   2       1.111  -6.660  -1.192  1.00 63.34           H
ATOM     27 HG13 VAL A   2       2.227  -6.098  -2.437  1.00 44.14           H
ATOM     28 HG21 VAL A   2       3.405  -4.080   0.789  1.00 40.14           H
ATOM     29 HG22 VAL A   2       4.082  -5.384  -0.186  1.00 63.31           H
ATOM     30 HG23 VAL A   2       3.948  -3.760  -0.859  1.00  1.55           H
ATOM     31  N   ASP A   3      -0.242  -3.468  -3.255  1.00 10.11           N
ATOM     32  CA  ASP A   3      -1.601  -3.665  -3.746  1.00 41.42           C
ATOM     33  C   ASP A   3      -1.619  -4.643  -4.916  1.00 74.14           C
ATOM     34  O   ASP A   3      -1.696  -4.238  -6.077  1.00 21.55           O
ATOM     35  CB  ASP A   3      -2.210  -2.328  -4.173  1.00 71.21           C
ATOM     36  CG  ASP A   3      -2.287  -1.335  -3.030  1.00 22.40           C
ATOM     37  OD1 ASP A   3      -1.324  -0.559  -2.852  1.00 33.41           O
ATOM     38  OD2 ASP A   3      -3.310  -1.334  -2.314  1.00  4.32           O
ATOM     39  H   ASP A   3       0.376  -2.907  -3.769  1.00 51.30           H
ATOM     40  HA  ASP A   3      -2.189  -4.076  -2.940  1.00 22.31           H
ATOM     41  HB2 ASP A   3      -1.605  -1.898  -4.958  1.00 61.42           H
ATOM     42  HB3 ASP A   3      -3.209  -2.498  -4.546  1.00 71.35           H
TER      43      ASP A   3
ENDMDL
MODEL        2
ATOM      1  N   ASN A   1       1.728  -3.986  -1.323  1.00 51.14           N
ATOM      2  CA  ASN A   1       2.250  -2.656  -1.616  1.00 71.03           C
ATOM      3  C   ASN A   1       1.152  -1.749  -2.162  1.00 53.34           C
ATOM      4  O   ASN A   1       0.899  -1.718  -3.367  1.00 12.41           O
ATOM      5  CB  ASN A   1       3.399  -2.747  -2.622  1.00 42.32           C
ATOM      6  CG  ASN A   1       4.579  -3.531  -2.082  1.00 65.14           C
ATOM      7  OD1 ASN A   1       5.209  -3.132  -1.102  1.00 55.44           O
ATOM      8  ND2 ASN A   1       4.886  -4.654  -2.722  1.00 35.14           N
ATOM      9  H1  ASN A   1       1.721  -4.663  -2.032  1.00 31.13           H
ATOM     10  HA  ASN A   1       2.623  -2.236  -0.694  1.00 54.13           H
ATOM     11  HB2 ASN A   1       3.047  -3.236  -3.519  1.00 11.31           H
ATOM     12  HB3 ASN A   1       3.734  -1.750  -2.868  1.00 63.43           H
ATOM     13 HD21 ASN A   1       4.340  -4.910  -3.495  1.00 24.31           H
ATOM     14 HD22 ASN A   1       5.644  -5.181  -2.393  1.00 73.24           H
ATOM     15  N   VAL A   2       0.502  -1.010  -1.268  1.00 21.51           N
ATOM     16  CA  VAL A   2      -0.568  -0.101  -1.660  1.00  1.25           C
ATOM     17  C   VAL A   2      -1.748  -0.863  -2.251  1.00 51.04           C
ATOM     18  O   VAL A   2      -2.193  -0.573  -3.362  1.00 73.21           O
ATOM     19  CB  VAL A   2      -0.075   0.937  -2.685  1.00 42.21           C
ATOM     20  CG1 VAL A   2      -1.101   2.047  -2.857  1.00 63.03           C
ATOM     21  CG2 VAL A   2       1.272   1.505  -2.262  1.00 55.34           C
ATOM     22  H   VAL A   2       0.749  -1.079  -0.322  1.00 40.21           H
ATOM     23  HA  VAL A   2      -0.898   0.426  -0.776  1.00 74.01           H
ATOM     24  HB  VAL A   2       0.050   0.442  -3.637  1.00 31.22           H
ATOM     25 HG11 VAL A   2      -0.656   2.993  -2.587  1.00  5.23           H
ATOM     26 HG12 VAL A   2      -1.425   2.081  -3.887  1.00 42.51           H
ATOM     27 HG13 VAL A   2      -1.950   1.854  -2.218  1.00 60.04           H
ATOM     28 HG21 VAL A   2       2.045   0.776  -2.452  1.00 73.22           H
ATOM     29 HG22 VAL A   2       1.477   2.403  -2.827  1.00 45.03           H
ATOM     30 HG23 VAL A   2       1.249   1.741  -1.209  1.00 33.20           H
ATOM     31  N   ASP A   3      -2.251  -1.838  -1.502  1.00 61.14           N
ATOM     32  CA  ASP A   3      -3.382  -2.642  -1.952  1.00 25.01           C
ATOM     33  C   ASP A   3      -3.129  -3.202  -3.348  1.00 61.32           C
ATOM     34  O   ASP A   3      -3.762  -2.788  -4.319  1.00 13.41           O
ATOM     35  CB  ASP A   3      -4.662  -1.806  -1.948  1.00 74.23           C
ATOM     36  CG  ASP A   3      -5.117  -1.450  -0.547  1.00 13.13           C
ATOM     37  OD1 ASP A   3      -4.608  -2.061   0.416  1.00 22.02           O
ATOM     38  OD2 ASP A   3      -5.982  -0.559  -0.412  1.00 71.45           O
ATOM     39  H   ASP A   3      -1.852  -2.021  -0.626  1.00 35.24           H
ATOM     40  HA  ASP A   3      -3.498  -3.465  -1.263  1.00 62.12           H
ATOM     41  HB2 ASP A   3      -4.488  -0.890  -2.494  1.00 34.02           H
ATOM     42  HB3 ASP A   3      -5.450  -2.365  -2.432  1.00 51.44           H
TER      43      ASP A   3
ENDMDL
MODEL        3
ATOM      1  N   ASN A   1       0.315  -4.452  -3.331  1.00 42.01           N
ATOM      2  CA  ASN A   1       0.480  -3.854  -2.011  1.00 52.12           C
ATOM      3  C   ASN A   1       0.359  -2.335  -2.083  1.00 54.35           C
ATOM      4  O   ASN A   1       0.991  -1.690  -2.920  1.00 11.31           O
ATOM      5  CB  ASN A   1       1.836  -4.241  -1.419  1.00  2.14           C
ATOM      6  CG  ASN A   1       1.801  -4.332   0.095  1.00 41.02           C
ATOM      7  OD1 ASN A   1       1.394  -3.390   0.775  1.00 22.22           O
ATOM      8  ND2 ASN A   1       2.229  -5.470   0.629  1.00 42.11           N
ATOM      9  H1  ASN A   1       1.058  -4.413  -3.969  1.00 51.24           H
ATOM     10  HA  ASN A   1      -0.304  -4.236  -1.374  1.00 64.13           H
ATOM     11  HB2 ASN A   1       2.132  -5.203  -1.811  1.00 23.01           H
ATOM     12  HB3 ASN A   1       2.570  -3.501  -1.700  1.00  1.32           H
ATOM     13 HD21 ASN A   1       2.539  -6.177   0.025  1.00 24.23           H
ATOM     14 HD22 ASN A   1       2.217  -5.556   1.605  1.00 72.10           H
ATOM     15  N   VAL A   2      -0.458  -1.769  -1.200  1.00 71.03           N
ATOM     16  CA  VAL A   2      -0.661  -0.326  -1.162  1.00 54.31           C
ATOM     17  C   VAL A   2      -1.033   0.213  -2.539  1.00 12.52           C
ATOM     18  O   VAL A   2      -0.225   0.866  -3.200  1.00 73.12           O
ATOM     19  CB  VAL A   2       0.598   0.405  -0.659  1.00 53.34           C
ATOM     20  CG1 VAL A   2       0.336   1.899  -0.539  1.00 60.15           C
ATOM     21  CG2 VAL A   2       1.053  -0.173   0.672  1.00 63.12           C
ATOM     22  H   VAL A   2      -0.935  -2.336  -0.559  1.00  4.15           H
ATOM     23  HA  VAL A   2      -1.470  -0.120  -0.475  1.00 63.43           H
ATOM     24  HB  VAL A   2       1.387   0.258  -1.381  1.00 64.43           H
ATOM     25 HG11 VAL A   2       1.120   2.443  -1.046  1.00 11.14           H
ATOM     26 HG12 VAL A   2      -0.617   2.135  -0.990  1.00  1.10           H
ATOM     27 HG13 VAL A   2       0.321   2.179   0.504  1.00 74.33           H
ATOM     28 HG21 VAL A   2       0.189  -0.445   1.260  1.00 64.44           H
ATOM     29 HG22 VAL A   2       1.659  -1.050   0.495  1.00 23.13           H
ATOM     30 HG23 VAL A   2       1.635   0.564   1.205  1.00 14.11           H
ATOM     31  N   ASP A   3      -2.261  -0.063  -2.964  1.00 42.14           N
ATOM     32  CA  ASP A   3      -2.741   0.396  -4.263  1.00 34.24           C
ATOM     33  C   ASP A   3      -4.199   0.839  -4.177  1.00  5.10           C
ATOM     34  O   ASP A   3      -4.996   0.250  -3.447  1.00  4.24           O
ATOM     35  CB  ASP A   3      -2.593  -0.713  -5.306  1.00 10.31           C
ATOM     36  CG  ASP A   3      -1.146  -0.959  -5.687  1.00 75.20           C
ATOM     37  OD1 ASP A   3      -0.715  -2.131  -5.652  1.00 42.22           O
ATOM     38  OD2 ASP A   3      -0.446   0.019  -6.021  1.00 14.21           O
ATOM     39  H   ASP A   3      -2.859  -0.588  -2.391  1.00 44.52           H
ATOM     40  HA  ASP A   3      -2.138   1.240  -4.560  1.00 13.30           H
ATOM     41  HB2 ASP A   3      -3.002  -1.630  -4.908  1.00 31.01           H
ATOM     42  HB3 ASP A   3      -3.139  -0.436  -6.196  1.00 11.03           H
TER      43      ASP A   3
ENDMDL
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK. Time: %8.3f"%(time.time()-t0))
