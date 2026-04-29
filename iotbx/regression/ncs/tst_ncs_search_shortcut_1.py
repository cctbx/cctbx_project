from __future__ import absolute_import, division, print_function
import iotbx.ncs
import iotbx.ncs as ncs
from iotbx import pdb
from libtbx.test_utils import show_diff

pdb_str_00="""\
CRYST1   34.917   22.246   44.017  90.00  90.00  90.00 P 1
SCALE1      0.028639  0.000000  0.000000        0.00000
SCALE2      0.000000  0.044952  0.000000        0.00000
SCALE3      0.000000  0.000000  0.022718        0.00000
ATOM      1  N   ALA A   1      27.344  16.348  30.784  1.00 10.00           N
ATOM      2  CA  ALA A   1      26.429  15.281  31.335  1.00 10.00           C
ATOM      3  C   ALA A   1      26.610  14.025  30.603  1.00 10.00           C
ATOM      4  O   ALA A   1      26.479  13.979  29.356  1.00 10.00           O
ATOM      5  CB  ALA A   1      24.874  15.800  31.300  1.00 10.00           C
ATOM      1  N   ALA A   2      26.812  12.925  31.345  1.00 10.00           N
ATOM      2  CA  ALA A   2      27.084  11.577  30.797  1.00 10.00           C
ATOM      3  C   ALA A   2      25.856  10.737  30.707  1.00 10.00           C
ATOM      4  O   ALA A   2      25.741   9.860  29.891  1.00 10.00           O
ATOM      5  CB  ALA A   2      28.151  10.950  31.721  1.00 10.00           C
ATOM      1  N   ALA A   3      25.009  10.973  31.714  1.00 10.00           N
ATOM      2  CA  ALA A   3      23.621  10.543  31.560  1.00 10.00           C
ATOM      3  C   ALA A   3      23.023  11.008  30.214  1.00 10.00           C
ATOM      4  O   ALA A   3      22.786  10.233  29.249  1.00 10.00           O
ATOM      5  CB  ALA A   3      22.760  11.040  32.654  1.00 10.00           C
ATOM      1  N   ALA A   4      22.798  12.304  30.175  1.00 10.00           N
ATOM      2  CA  ALA A   4      22.329  13.084  28.981  1.00 10.00           C
ATOM      4  O   ALA A   4      22.533  12.805  26.670  1.00 10.00           O
ATOM      5  CB  ALA A   4      22.372  14.607  29.318  1.00 10.00           C
ATOM      1  N   ALA A   5      24.448  12.622  27.823  1.00 10.00           N
ATOM      2  CA  ALA A   5      25.228  12.407  26.573  1.00 10.00           C
ATOM      3  C   ALA A   5      25.222  10.947  26.143  1.00 10.00           C
ATOM      4  O   ALA A   5      25.386  10.664  24.983  1.00 10.00           O
ATOM      5  CB  ALA A   5      26.634  12.906  26.746  1.00 10.00           C
ATOM      1  N   ALA A   6      24.976  10.048  27.071  1.00 10.00           N
ATOM      2  CA  ALA A   6      24.857   8.614  26.805  1.00 10.00           C
ATOM      3  C   ALA A   6      23.537   8.349  26.054  1.00 10.00           C
ATOM      4  O   ALA A   6      23.439   7.570  25.057  1.00 10.00           O
ATOM      5  CB  ALA A   6      24.874   7.845  28.114  1.00 10.00           C
ATOM      1  N   ALA A   7      22.542   9.039  26.580  1.00 10.00           N
ATOM      2  CA  ALA A   7      21.228   8.903  25.942  1.00 10.00           C
ATOM      3  C   ALA A   7      21.329   9.698  24.628  1.00 10.00           C
ATOM      4  O   ALA A   7      20.707   9.383  23.632  1.00 10.00           O
ATOM      5  CB  ALA A   7      20.146   9.465  26.862  1.00 10.00           C
ATOM      1  N   ALA A   8      22.181  10.696  24.613  1.00 10.00           N
ATOM      2  CA  ALA A   8      22.526  11.372  23.378  1.00 10.00           C
ATOM      3  C   ALA A   8      23.351  10.555  22.448  1.00 10.00           C
ATOM      4  O   ALA A   8      23.618  10.883  21.252  1.00 10.00           O
ATOM      5  CB  ALA A   8      23.168  12.697  23.693  1.00 10.00           C
ATOM      1  N   ALA A   9      23.864   9.423  22.961  1.00 10.00           N
ATOM      2  CA  ALA A   9      24.785   8.541  22.264  1.00 10.00           C
ATOM      3  C   ALA A   9      24.057   7.451  21.484  1.00 10.00           C
ATOM      4  O   ALA A   9      24.127   7.381  20.257  1.00 10.00           O
ATOM      5  CB  ALA A   9      25.815   7.975  23.249  1.00 10.00           C
ATOM      1  N   ALA A  10      23.518   6.548  22.264  1.00 10.00           N
ATOM      2  CA  ALA A  10      22.629   5.525  21.690  1.00 10.00           C
ATOM      3  C   ALA A  10      21.549   6.308  21.009  1.00 10.00           C
ATOM      4  O   ALA A  10      21.114   5.933  19.930  1.00 10.00           O
ATOM      5  CB  ALA A  10      22.057   4.714  22.784  1.00 10.00           C
ATOM      1  N   ALA A  11      21.120   7.452  21.541  1.00 10.00           N
ATOM      2  CA  ALA A  11      20.186   8.260  20.874  1.00 10.00           C
ATOM      3  C   ALA A  11      20.978   9.215  19.937  1.00 10.00           C
ATOM      4  O   ALA A  11      20.386  10.177  19.507  1.00 10.00           O
ATOM      5  CB  ALA A  11      19.295   9.031  21.867  1.00 10.00           C
ATOM      1  N   ALA A  12      22.222   8.932  19.598  1.00 10.00           N
ATOM      2  CA  ALA A  12      22.896   9.709  18.563  1.00 10.00           C
ATOM      3  C   ALA A  12      22.924   8.925  17.308  1.00 10.00           C
ATOM      4  O   ALA A  12      22.982   9.445  16.193  1.00 10.00           O
ATOM      5  CB  ALA A  12      24.294  10.138  18.994  1.00 10.00           C
ATOM      1  N   ALA A  13      22.951   7.633  17.508  1.00 10.00           N
ATOM      2  CA  ALA A  13      22.709   6.629  16.554  1.00 10.00           C
ATOM      3  C   ALA A  13      21.275   6.673  16.206  1.00 10.00           C
ATOM      4  O   ALA A  13      20.870   6.521  15.092  1.00 10.00           O
ATOM      5  CB  ALA A  13      23.077   5.254  17.025  1.00 10.00           C
ATOM      1  N   ALA A  14      20.471   6.929  17.226  1.00 10.00           N
ATOM      2  CA  ALA A  14      19.039   6.992  17.025  1.00 10.00           C
ATOM      3  C   ALA A  14      18.676   8.380  16.528  1.00 10.00           C
ATOM      4  O   ALA A  14      17.748   8.556  15.761  1.00 10.00           O
ATOM      5  CB  ALA A  14      18.240   6.715  18.272  1.00 10.00           C
ATOM      1  N   ALA A  15      19.381   9.390  17.055  1.00 10.00           N
ATOM      2  CA  ALA A  15      19.204  10.743  16.669  1.00 10.00           C
ATOM      3  C   ALA A  15      19.407  10.807  15.174  1.00 10.00           C
ATOM      4  O   ALA A  15      18.402  10.987  14.424  1.00 10.00           O
ATOM      5  CB  ALA A  15      20.190  11.665  17.493  1.00 10.00           C
ATOM      1  N   ALA A  16      20.702  10.653  14.831  1.00 10.00           N
ATOM      2  CA  ALA A  16      21.206  10.546  13.480  1.00 10.00           C
ATOM      3  C   ALA A  16      20.484   9.612  12.585  1.00 10.00           C
ATOM      4  O   ALA A  16      20.380   9.918  11.386  1.00 10.00           O
ATOM      5  CB  ALA A  16      22.631  10.174  13.475  1.00 10.00           C
ATOM      1  N   ALA A  17      20.064   8.475  13.175  1.00 10.00           N
ATOM      2  CA  ALA A  17      19.355   7.473  12.426  1.00 10.00           C
ATOM      3  C   ALA A  17      17.924   7.807  12.064  1.00 10.00           C
ATOM      4  O   ALA A  17      17.535   7.721  10.871  1.00 10.00           O
ATOM      5  CB  ALA A  17      19.359   6.123  13.216  1.00 10.00           C
ATOM      1  N   ALA A  18      17.152   8.115  13.031  1.00 10.00           N
ATOM      2  CA  ALA A  18      15.835   8.594  12.861  1.00 10.00           C
ATOM      3  C   ALA A  18      15.811   9.835  11.861  1.00 10.00           C
ATOM      4  O   ALA A  18      15.020   9.889  10.868  1.00 10.00           O
ATOM      5  CB  ALA A  18      15.272   8.918  14.234  1.00 10.00           C
ATOM      1  N   ALA A  19      16.661  10.845  12.100  1.00 10.00           N
ATOM      2  CA  ALA A  19      16.435  12.061  11.275  1.00 10.00           C
ATOM      3  C   ALA A  19      17.004  11.815   9.833  1.00 10.00           C
ATOM      4  O   ALA A  19      16.334  12.117   8.857  1.00 10.00           O
ATOM      5  CB  ALA A  19      17.059  13.242  11.866  1.00 10.00           C
ATOM      1  N   ALA A  20      18.191  11.200   9.841  1.00 10.00           N
ATOM      2  CA  ALA A  20      19.091  11.247   8.697  1.00 10.00           C
ATOM      3  C   ALA A  20      19.549   9.835   8.231  1.00 10.00           C
ATOM      4  O   ALA A  20      20.670   9.692   7.663  1.00 10.00           O
ATOM      5  CB  ALA A  20      20.326  12.105   9.035  1.00 10.00           C
ATOM      1  N   ALA A  21      18.654   8.850   8.523  1.00 10.00           N
ATOM      2  CA  ALA A  21      18.827   7.437   8.168  1.00 10.00           C
ATOM      3  C   ALA A  21      17.565   6.607   8.282  1.00 10.00           C
ATOM      4  O   ALA A  21      16.485   6.992   7.820  1.00 10.00           O
ATOM      5  CB  ALA A  21      19.888   6.838   8.983  1.00 10.00           C
TER
ATOM      1  N   ALA B   1      16.348  17.420  35.897  1.00 50.00           N
ATOM      2  CA  ALA B   1      16.783  16.083  36.351  1.00 50.00           C
ATOM      3  C   ALA B   1      16.794  15.172  35.139  1.00 50.00           C
ATOM      4  O   ALA B   1      16.167  15.477  34.133  1.00 50.00           O
ATOM      5  CB  ALA B   1      15.785  15.534  37.468  1.00 50.00           C
ATOM      1  N   ALA B   2      17.491  14.058  35.255  1.00 50.00           N
ATOM      2  CA  ALA B   2      17.790  13.267  34.127  1.00 50.00           C
ATOM      3  C   ALA B   2      16.716  12.232  33.688  1.00 50.00           C
ATOM      4  O   ALA B   2      16.676  11.869  32.543  1.00 50.00           O
ATOM      5  CB  ALA B   2      19.125  12.656  34.415  1.00 50.00           C
ATOM      1  N   ALA B   3      15.904  11.687  34.605  1.00 50.00           N
ATOM      2  CA  ALA B   3      14.798  10.901  34.173  1.00 50.00           C
ATOM      3  C   ALA B   3      13.740  11.723  33.536  1.00 50.00           C
ATOM      4  O   ALA B   3      13.398  11.501  32.356  1.00 50.00           O
ATOM      5  CB  ALA B   3      14.148  10.176  35.403  1.00 50.00           C
ATOM      1  N   ALA B   4      13.239  12.708  34.247  1.00 50.00           N
ATOM      2  CA  ALA B   4      12.158  13.487  33.709  1.00 50.00           C
ATOM      3  C   ALA B   4      12.674  14.248  32.495  1.00 50.00           C
ATOM      4  O   ALA B   4      11.935  14.376  31.526  1.00 50.00           O
ATOM      5  CB  ALA B   4      11.553  14.432  34.712  1.00 50.00           C
ATOM      1  N   ALA B   5      13.947  14.627  32.479  1.00 50.00           N
ATOM      2  CA  ALA B   5      14.416  15.490  31.405  1.00 50.00           C
ATOM      3  C   ALA B   5      14.960  14.730  30.186  1.00 50.00           C
ATOM      4  O   ALA B   5      14.575  14.940  29.054  1.00 50.00           O
ATOM      5  CB  ALA B   5      15.464  16.431  31.928  1.00 50.00           C
ATOM      1  N   ALA B   6      15.867  13.827  30.546  1.00 50.00           N
ATOM      2  CA  ALA B   6      16.575  12.918  29.615  1.00 50.00           C
ATOM      3  C   ALA B   6      15.465  12.002  28.975  1.00 50.00           C
ATOM      4  O   ALA B   6      15.450  11.709  27.742  1.00 50.00           O
ATOM      5  CB  ALA B   6      17.632  12.157  30.362  1.00 50.00           C
ATOM      1  N   ALA B   7      14.542  11.597  29.783  1.00 50.00           N
ATOM      2  CA  ALA B   7      13.529  10.701  29.277  1.00 50.00           C
ATOM      3  C   ALA B   7      12.175  11.364  28.835  1.00 50.00           C
ATOM      4  O   ALA B   7      11.466  10.770  27.969  1.00 50.00           O
ATOM      5  CB  ALA B   7      13.161   9.644  30.376  1.00 50.00           C
ATOM      1  N   ALA B   8      11.753  12.455  29.452  1.00 50.00           N
ATOM      2  CA  ALA B   8      10.536  13.193  28.972  1.00 50.00           C
ATOM      3  C   ALA B   8      10.919  13.923  27.670  1.00 50.00           C
ATOM      4  O   ALA B   8      10.171  14.036  26.729  1.00 50.00           O
ATOM      5  CB  ALA B   8      10.032  14.139  30.014  1.00 50.00           C
ATOM      1  N   ALA B   9      12.185  14.247  27.579  1.00 50.00           N
ATOM      2  CA  ALA B   9      12.754  14.849  26.385  1.00 50.00           C
ATOM      3  C   ALA B   9      12.892  13.859  25.320  1.00 50.00           C
ATOM      4  O   ALA B   9      12.234  13.980  24.290  1.00 50.00           O
ATOM      5  CB  ALA B   9      14.108  15.448  26.695  1.00 50.00           C
ATOM      1  N   ALA B  10      13.655  12.794  25.566  1.00 50.00           N
ATOM      2  CA  ALA B  10      13.831  11.803  24.529  1.00 50.00           C
ATOM      3  C   ALA B  10      12.551  10.987  24.319  1.00 50.00           C
ATOM      4  O   ALA B  10      12.514  10.237  23.390  1.00 50.00           O
ATOM      5  CB  ALA B  10      15.024  10.750  24.992  1.00 50.00           C
ATOM      1  N   ALA B  11      11.558  11.184  25.126  1.00 50.00           N
ATOM      2  CA  ALA B  11      10.334  10.457  24.931  1.00 50.00           C
ATOM      3  C   ALA B  11       9.326  11.284  24.168  1.00 50.00           C
ATOM      4  O   ALA B  11       8.566  10.707  23.476  1.00 50.00           O
ATOM      5  CB  ALA B  11       9.644  10.042  26.251  1.00 50.00           C
ATOM      1  N   ALA B  12       9.277  12.611  24.334  1.00 50.00           N
ATOM      2  CA  ALA B  12       8.354  13.375  23.644  1.00 50.00           C
ATOM      3  C   ALA B  12       9.019  13.546  22.264  1.00 50.00           C
ATOM      4  O   ALA B  12       8.400  13.891  21.317  1.00 50.00           O
ATOM      5  CB  ALA B  12       8.056  14.678  24.287  1.00 50.00           C
ATOM      1  N   ALA B  13      10.333  13.339  22.264  1.00 50.00           N
ATOM      2  CA  ALA B  13      11.239  13.471  21.127  1.00 50.00           C
ATOM      3  C   ALA B  13      11.096  12.161  20.325  1.00 50.00           C
ATOM      4  O   ALA B  13      11.145  12.175  19.123  1.00 50.00           O
ATOM      5  CB  ALA B  13      12.584  13.665  21.596  1.00 50.00           C
ATOM      1  N   ALA B  14      11.051  11.078  21.086  1.00 50.00           N
ATOM      2  CA  ALA B  14      10.953   9.771  20.454  1.00 50.00           C
ATOM      3  C   ALA B  14       9.550   9.463  20.117  1.00 50.00           C
ATOM      4  O   ALA B  14       9.233   8.571  19.367  1.00 50.00           O
ATOM      1  N   ALA B  15       8.669  10.215  20.743  1.00 50.00           N
ATOM      2  CA  ALA B  15       7.282  10.010  20.486  1.00 50.00           C
ATOM      3  C   ALA B  15       6.825  10.982  19.376  1.00 50.00           C
ATOM      4  O   ALA B  15       5.855  10.783  18.619  1.00 50.00           O
ATOM      5  CB  ALA B  15       6.367  10.306  21.797  1.00 50.00           C
ATOM      1  N   ALA B  16       7.511  12.143  19.430  1.00 50.00           N
ATOM      2  CA  ALA B  16       7.233  13.302  18.551  1.00 50.00           C
ATOM      3  C   ALA B  16       7.912  13.082  17.205  1.00 50.00           C
ATOM      4  O   ALA B  16       7.492  13.573  16.111  1.00 50.00           O
ATOM      5  CB  ALA B  16       7.762  14.594  19.165  1.00 50.00           C
ATOM      1  N   ALA B  17       9.071  12.427  17.269  1.00 50.00           N
ATOM      2  CA  ALA B  17       9.595  11.771  16.091  1.00 50.00           C
ATOM      3  C   ALA B  17       8.883  10.519  15.763  1.00 50.00           C
ATOM      4  O   ALA B  17       8.890  10.193  14.597  1.00 50.00           O
ATOM      5  CB  ALA B  17      11.046  11.518  16.265  1.00 50.00           C
ATOM      1  N   ALA B  18       8.315   9.809  16.722  1.00 50.00           N
ATOM      2  CA  ALA B  18       7.515   8.647  16.448  1.00 50.00           C
ATOM      3  C   ALA B  18       6.253   9.063  15.707  1.00 50.00           C
ATOM      4  O   ALA B  18       5.559   8.173  15.198  1.00 50.00           O
ATOM      5  CB  ALA B  18       7.129   7.915  17.695  1.00 50.00           C
ATOM      1  N   ALA B  19       5.866  10.332  15.772  1.00 50.00           N
ATOM      2  CA  ALA B  19       4.686  10.808  15.089  1.00 50.00           C
ATOM      3  C   ALA B  19       5.011  11.578  13.803  1.00 50.00           C
ATOM      4  O   ALA B  19       4.291  11.514  12.837  1.00 50.00           O
ATOM      5  CB  ALA B  19       3.854  11.710  15.960  1.00 50.00           C
ATOM      1  N   ALA B  20       6.176  12.195  13.822  1.00 50.00           N
ATOM      2  CA  ALA B  20       6.614  13.121  12.789  1.00 50.00           C
ATOM      3  C   ALA B  20       7.933  12.759  12.098  1.00 50.00           C
ATOM      4  O   ALA B  20       8.620  13.613  11.585  1.00 50.00           O
ATOM      5  CB  ALA B  20       6.823  14.498  13.449  1.00 50.00           C
ATOM      1  N   ALA B  21       8.284  11.511  12.050  1.00 50.00           N
ATOM      2  CA  ALA B  21       9.513  11.117  11.323  1.00 50.00           C
ATOM      3  C   ALA B  21       9.313   9.628  11.029  1.00 50.00           C
ATOM      4  O   ALA B  21       9.731   8.751  11.795  1.00 50.00           O
ATOM      5  CB  ALA B  21      10.799  11.332  12.178  1.00 50.00           C
TER
"""

pdb_str_01 = """\
CRYST1   40.408   25.281   31.173  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1       5.043   9.706  12.193  1.00  1.00           N
ATOM      2  CA  GLY A   1       5.000   9.301  10.742  1.00  1.00           C
ATOM      3  C   GLY A   1       6.037   8.234  10.510  1.00  1.00           C
ATOM      4  O   GLY A   1       6.529   7.615  11.472  1.00  1.00           O
ATOM      5  N   ASN A   2       6.396   8.017   9.246  1.00  1.00           N
ATOM      6  CA  ASN A   2       7.530   7.132   8.922  1.00  1.00           C
ATOM      7  C   ASN A   2       8.811   7.631   9.518  1.00  1.00           C
ATOM      8  O   ASN A   2       9.074   8.836   9.517  1.00  1.00           O
ATOM      9  CB  ASN A   2       7.706   6.975   7.432  1.00  1.00           C
ATOM     10  CG  ASN A   2       6.468   6.436   6.783  1.00  1.00           C
ATOM     11  OD1 ASN A   2       6.027   5.321   7.107  1.00  1.00           O
ATOM     12  ND2 ASN A   2       5.848   7.249   5.922  1.00  1.00           N
ATOM     13  N   ASN A   3       9.614   6.684   9.996  1.00  1.00           N
ATOM     14  CA  ASN A   3      10.859   6.998  10.680  1.00  1.00           C
ATOM     15  C   ASN A   3      12.097   6.426   9.986  1.00  1.00           C
ATOM     16  O   ASN A   3      12.180   5.213   9.739  1.00  1.00           O
ATOM     17  CB  ASN A   3      10.793   6.472  12.133  1.00  1.00           C
ATOM     18  CG  ASN A   3      12.046   6.833  12.952  1.00  1.00           C
ATOM     19  OD1 ASN A   3      12.350   8.019  13.163  1.00  1.00           O
ATOM     20  ND2 ASN A   3      12.781   5.809  13.397  1.00  1.00           N
ATOM     21  N  AGLN A   4      13.047   7.322   9.689  0.50  1.00           N
ATOM     22  CA AGLN A   4      14.436   6.982   9.290  0.50  1.00           C
ATOM     23  C  AGLN A   4      15.487   7.700  10.179  0.50  1.00           C
ATOM     24  O  AGLN A   4      15.599   8.937  10.206  0.50  1.00           O
ATOM     25  CB AGLN A   4      14.708   7.242   7.802  0.50  1.00           C
ATOM     26  CG AGLN A   4      15.996   6.552   7.304  0.50  1.00           C
ATOM     27  CD AGLN A   4      16.556   7.138   6.002  0.50  1.00           C
ATOM     28  OE1AGLN A   4      16.796   8.362   5.901  0.50  1.00           O
ATOM     29  NE2AGLN A   4      16.802   6.255   5.000  0.50  1.00           N
ATOM     30  N  BGLN A   4      13.147   7.322   9.689  0.50  1.00           N
ATOM     31  CA BGLN A   4      14.536   6.982   9.290  0.50  1.00           C
ATOM     32  C  BGLN A   4      15.587   7.700  10.179  0.50  1.00           C
ATOM     33  O  BGLN A   4      15.699   8.937  10.206  0.50  1.00           O
ATOM     34  CB BGLN A   4      14.808   7.242   7.802  0.50  1.00           C
ATOM     35  CG BGLN A   4      14.009   6.298   6.878  0.50  1.00           C
ATOM     36  CD BGLN A   4      14.201   6.577   5.382  0.50  1.00           C
ATOM     37  OE1BGLN A   4      14.536   7.714   4.980  0.50  1.00           O
ATOM     38  NE2BGLN A   4      14.021   5.524   4.544  0.50  1.00           N
ATOM     39  N   GLN A   5      16.206   6.915  10.962  1.00  1.00           N
ATOM     40  CA  GLN A   5      17.322   7.455  11.731  1.00  1.00           C
ATOM     41  C   GLN A   5      18.646   6.862  11.263  1.00  1.00           C
ATOM     42  O   GLN A   5      18.820   5.640  11.145  1.00  1.00           O
ATOM     43  CB  GLN A   5      17.108   7.277  13.238  1.00  1.00           C
ATOM     44  CG  GLN A   5      15.881   8.044  13.738  1.00  1.00           C
ATOM     45  CD  GLN A   5      15.396   7.508  15.045  1.00  1.00           C
ATOM     46  OE1 GLN A   5      14.826   6.419  15.093  1.00  1.00           O
ATOM     47  NE2 GLN A   5      15.601   8.281  16.130  1.00  1.00           N
ATOM     48  N   ASN A   6      19.566   7.758  10.947  1.00  1.00           N
ATOM     49  CA  ASN A   6      20.883   7.404  10.409  1.00  1.00           C
ATOM     50  C   ASN A   6      21.906   7.855  11.415  1.00  1.00           C
ATOM     51  O   ASN A   6      22.271   9.037  11.465  1.00  1.00           O
ATOM     52  CB  ASN A   6      21.117   8.110   9.084  1.00  1.00           C
ATOM     53  CG  ASN A   6      20.013   7.829   8.094  1.00  1.00           C
ATOM     54  OD1 ASN A   6      19.850   6.698   7.642  1.00  1.00           O
ATOM     55  ND2 ASN A   6      19.247   8.841   7.770  1.00  1.00           N
ATOM     56  N   TYR A   7      22.344   6.911  12.238  1.00  1.00           N
ATOM     57  CA  TYR A   7      23.211   7.238  13.390  1.00  1.00           C
ATOM     58  C   TYR A   7      24.655   7.425  12.976  1.00  1.00           C
ATOM     59  O   TYR A   7      25.093   6.905  11.946  1.00  1.00           O
ATOM     60  CB  TYR A   7      23.113   6.159  14.460  1.00  1.00           C
ATOM     61  CG  TYR A   7      21.717   6.023  14.993  1.00  1.00           C
ATOM     62  CD1 TYR A   7      20.823   5.115  14.418  1.00  1.00           C
ATOM     63  CD2 TYR A   7      21.262   6.850  16.011  1.00  1.00           C
ATOM     64  CE1 TYR A   7      19.532   5.000  14.887  1.00  1.00           C
ATOM     65  CE2 TYR A   7      19.956   6.743  16.507  1.00  1.00           C
ATOM     66  CZ  TYR A   7      19.099   5.823  15.922  1.00  1.00           C
ATOM     67  OH  TYR A   7      17.818   5.683  16.382  1.00  1.00           O
ATOM     68  OXT TYR A   7      25.410   8.093  13.703  1.00  1.00           O
TER      69      TYR A   7
ATOM     70  N   GLY B   1      15.117  17.344  24.100  1.00  2.00           N
ATOM     71  CA  GLY B   1      14.893  17.193  22.618  1.00  2.00           C
ATOM     72  C   GLY B   1      16.041  16.412  22.035  1.00  2.00           C
ATOM     73  O   GLY B   1      16.791  15.749  22.776  1.00  2.00           O
ATOM     74  N   ASN B   2      16.207  16.491  20.716  1.00  2.00           N
ATOM     75  CA  ASN B   2      17.402  15.920  20.067  1.00  2.00           C
ATOM     76  C   ASN B   2      18.662  16.557  20.568  1.00  2.00           C
ATOM     77  O   ASN B   2      18.711  17.772  20.773  1.00  2.00           O
ATOM     78  CB  ASN B   2      17.340  16.059  18.566  1.00  2.00           C
ATOM     79  CG  ASN B   2      16.119  15.401  18.000  1.00  2.00           C
ATOM     80  OD1 ASN B   2      15.938  14.182  18.151  1.00  2.00           O
ATOM     81  ND2 ASN B   2      15.229  16.208  17.414  1.00  2.00           N
ATOM     82  N   ASN B   3      19.686  15.722  20.731  1.00  2.00           N
ATOM     83  CA  ASN B   3      20.958  16.158  21.285  1.00  2.00           C
ATOM     84  C   ASN B   3      22.136  15.974  20.326  1.00  2.00           C
ATOM     85  O   ASN B   3      22.381  14.862  19.832  1.00  2.00           O
ATOM     86  CB  ASN B   3      21.237  15.389  22.598  1.00  2.00           C
ATOM     87  CG  ASN B   3      22.532  15.849  23.291  1.00  2.00           C
ATOM     88  OD1 ASN B   3      22.661  17.018  23.692  1.00  2.00           O
ATOM     89  ND2 ASN B   3      23.498  14.932  23.416  1.00  2.00           N
ATOM     90  N  AGLN B   4      22.853  17.079  20.087  0.50  2.00           N
ATOM     91  CA AGLN B   4      24.189  17.098  19.439  0.50  2.00           C
ATOM     92  C  AGLN B   4      25.240  17.850  20.300  0.50  2.00           C
ATOM     93  O  AGLN B   4      25.142  19.061  20.559  0.50  2.00           O
ATOM     94  CB AGLN B   4      24.150  17.658  18.011  0.50  2.00           C
ATOM     95  CG AGLN B   4      25.431  17.336  17.211  0.50  2.00           C
ATOM     96  CD AGLN B   4      25.647  18.237  15.989  0.50  2.00           C
ATOM     97  OE1AGLN B   4      25.653  19.483  16.103  0.50  2.00           O
ATOM     98  NE2AGLN B   4      25.863  17.606  14.806  0.50  2.00           N
ATOM     99  N  BGLN B   4      22.953  17.079  20.087  0.50  2.00           N
ATOM    100  CA BGLN B   4      24.289  17.098  19.439  0.50  2.00           C
ATOM    101  C  BGLN B   4      25.340  17.850  20.300  0.50  2.00           C
ATOM    102  O  BGLN B   4      25.242  19.061  20.559  0.50  2.00           O
ATOM    103  CB BGLN B   4      24.250  17.658  18.011  0.50  2.00           C
ATOM    104  CG BGLN B   4      23.476  16.745  17.035  0.50  2.00           C
ATOM    105  CD BGLN B   4      23.355  17.310  15.615  0.50  2.00           C
ATOM    106  OE1BGLN B   4      23.416  18.542  15.406  0.50  2.00           O
ATOM    107  NE2BGLN B   4      23.215  16.401  14.615  0.50  2.00           N
ATOM    108  N   GLN B   5      26.207  17.103  20.803  1.00  2.00           N
ATOM    109  CA  GLN B   5      27.331  17.717  21.503  1.00  2.00           C
ATOM    110  C   GLN B   5      28.635  17.490  20.747  1.00  2.00           C
ATOM    111  O   GLN B   5      28.992  16.367  20.363  1.00  2.00           O
ATOM    112  CB  GLN B   5      27.415  17.244  22.958  1.00  2.00           C
ATOM    113  CG  GLN B   5      26.181  17.652  23.767  1.00  2.00           C
ATOM    114  CD  GLN B   5      26.029  16.814  24.994  1.00  2.00           C
ATOM    115  OE1 GLN B   5      25.671  15.641  24.900  1.00  2.00           O
ATOM    116  NE2 GLN B   5      26.284  17.416  26.173  1.00  2.00           N
ATOM    117  N   ASN B   6      29.319  18.593  20.493  1.00  2.00           N
ATOM    118  CA  ASN B   6      30.564  18.608  19.718  1.00  2.00           C
ATOM    119  C   ASN B   6      31.653  19.076  20.643  1.00  2.00           C
ATOM    120  O   ASN B   6      31.814  20.281  20.878  1.00  2.00           O
ATOM    121  CB  ASN B   6      30.440  19.563  18.542  1.00  2.00           C
ATOM    122  CG  ASN B   6      29.245  19.239  17.678  1.00  2.00           C
ATOM    123  OD1 ASN B   6      29.202  18.193  17.036  1.00  2.00           O
ATOM    124  ND2 ASN B   6      28.273  20.117  17.673  1.00  2.00           N
ATOM    125  N   TYR B   7      32.383  18.113  21.191  1.00  2.00           N
ATOM    126  CA  TYR B   7      33.367  18.405  22.254  1.00  2.00           C
ATOM    127  C   TYR B   7      34.664  18.946  21.690  1.00  2.00           C
ATOM    128  O   TYR B   7      34.999  18.709  20.526  1.00  2.00           O
ATOM    129  CB  TYR B   7      33.643  17.162  23.089  1.00  2.00           C
ATOM    130  CG  TYR B   7      32.405  16.659  23.771  1.00  2.00           C
ATOM    131  CD1 TYR B   7      31.593  15.703  23.155  1.00  2.00           C
ATOM    132  CD2 TYR B   7      31.999  17.192  24.987  1.00  2.00           C
ATOM    133  CE1 TYR B   7      30.442  15.252  23.765  1.00  2.00           C
ATOM    134  CE2 TYR B   7      30.837  16.742  25.628  1.00  2.00           C
ATOM    135  CZ  TYR B   7      30.061  15.782  24.994  1.00  2.00           C
ATOM    136  OH  TYR B   7      28.923  15.311  25.589  1.00  2.00           O
ATOM    137  OXT TYR B   7      35.408  19.618  22.425  1.00  2.00           O
TER     138      TYR B   7
END
"""

def exercise_00():
  """ Verify that shortcut1 is not failing"""
  pdb_inp = iotbx.pdb.input(lines=pdb_str_00,source_info=None)
  hierarchy = pdb_inp.construct_hierarchy()
  ncs_params = iotbx.ncs.input.get_default_params()
  ncs_params.ncs_search.try_shortcuts = True
  ncs_params.ncs_search.exclude_selection="water"

  ncs_inp = iotbx.ncs.input(
      hierarchy = hierarchy,
      params = ncs_params.ncs_search)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  ncs_inp.show(format='phil')
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1, len(ncs_groups[0].copies)
  assert ncs_groups[0].master_iselection.size() == 103, ncs_groups[0].master_iselection.size() # all chain
  assert ncs_groups[0].copies[0].iselection.size() == 103, ncs_groups[0].copies[0].iselection.size() # all chain

def exercise_01():
  """ Alternative conformations should remain.
  """
  pdb_inp = iotbx.pdb.input(lines=pdb_str_01,source_info=None)
  hierarchy = pdb_inp.construct_hierarchy()
  # print('n_atoms:', hierarchy.atoms_size())
  ncs_params = iotbx.ncs.input.get_default_params()
  ncs_params.ncs_search.try_shortcuts = True
  ncs_params.ncs_search.exclude_selection="water"

  ncs_inp = iotbx.ncs.input(
      hierarchy = hierarchy,
      params = ncs_params.ncs_search)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  ncs_phil_str = ncs_inp.show(format='phil')
  assert not show_diff(ncs_phil_str, """
NCS phil parameters:
-----------------------
ncs_group {
  reference = chain 'A'
  selection = chain 'B'
}""")
  assert len(ncs_groups) == 1
  assert len(ncs_groups[0].copies) == 1, len(ncs_groups[0].copies)
  assert ncs_groups[0].master_iselection.size() == 68, ncs_groups[0].master_iselection.size() # all chain
  assert ncs_groups[0].copies[0].iselection.size() == 68, ncs_groups[0].copies[0].iselection.size() # all chain

if(__name__=='__main__'):
  exercise_00()
  exercise_01()
  print("OK")
