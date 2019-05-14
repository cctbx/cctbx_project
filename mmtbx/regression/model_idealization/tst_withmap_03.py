from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os.path
import time

cryst_str = """\n
CRYST1   34.917   22.246   44.017  90.00  90.00  90.00 P 1
SCALE1      0.028639  0.000000  0.000000        0.00000
SCALE2      0.000000  0.044952  0.000000        0.00000
SCALE3      0.000000  0.000000  0.022718        0.00000
"""
# taken from phenix_regression/refinement/ncs/tst_ncs_0.py
pdb_str = """\
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
ATOM      3  C   ALA A   4      23.116  12.816  27.721  1.00 10.00           C
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
ATOM      5  CB  ALA B  14      11.461   8.697  21.413  1.00 50.00           C
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

def exercise_03(prefix="tst_mi_map_test_03"):
  """
  Run with reference map.
  Check if working with NCS in the model with symmetry.
  """
  # with cryst
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(cryst_str)
  pdb_file.write(pdb_str)
  pdb_file.close()
  cmd = " ".join([
      "phenix.model_idealization",
      "%s_start.pdb" % prefix,
      "use_map_for_reference=True",
      "loop_idealization.number_of_ccd_trials=1",
      "number_of_refinement_cycles=1",
      "n_macro=1",
      "debug=True",
      ">%s.log" % prefix])
  print(cmd)
  assert not easy_run.call(cmd)
  assert os.path.isfile("%s_start.pdb_all_idealized.pdb" % prefix)
  res_log = open("%s.log" % prefix, "r")
  log_lines = res_log.readlines()
  # NCS constraints with map are not implemented yet
  for l in [
      # "Using ncs\n",
      "  Minimizing... (NCS)\n",
      "Using map as reference\n",
      # "Ramachandran outliers:      0.00      0.00      0.00      0.00      0.00\n",
      "All done.\n"]:
    assert l in log_lines, "'%s' not in log file." % l
  res_log.close()

if (__name__ == "__main__"):
  t0 = time.time()
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping: probe not configured")
  else:
    exercise_03()
  print("Time: %.2f" % (time.time() - t0))
  print("OK")
