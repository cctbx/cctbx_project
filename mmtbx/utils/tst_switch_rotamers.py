from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.utils

pdb_str = """\
ATOM      1  N   TYR A  58       8.659  20.073  11.185  1.00  7.73           N
ATOM      2  CA  TYR A  58       9.250  19.144  10.233  1.00  8.65           C
ATOM      3  C   TYR A  58       9.039  17.721  10.706  1.00  9.84           C
ATOM      4  O   TYR A  58       9.023  17.464  11.919  1.00  8.58           O
ATOM      5  CB  TYR A  58      10.740  19.429  10.045  1.00 20.00           C
ATOM      6  CG  TYR A  58      11.155  20.825  10.457  1.00 20.00           C
ATOM      7  CD1 TYR A  58      10.204  21.809  10.700  1.00 20.00           C
ATOM      8  CD2 TYR A  58      12.494  21.155  10.603  1.00 20.00           C
ATOM      9  CE1 TYR A  58      10.579  23.084  11.076  1.00 20.00           C
ATOM     10  CE2 TYR A  58      12.879  22.429  10.980  1.00 20.00           C
ATOM     11  CZ  TYR A  58      11.918  23.388  11.215  1.00 20.00           C
ATOM     12  OH  TYR A  58      12.294  24.659  11.588  1.00 20.00           O
ATOM     21  N   SER A  59       8.887  16.797   9.756  1.00  9.65           N
ATOM     22  CA  SER A  59       8.720  15.382  10.080  1.00  5.80           C
ATOM     23  C   SER A  59       9.606  14.539   9.169  1.00 10.35           C
ATOM     24  O   SER A  59       9.972  14.972   8.075  1.00 10.56           O
ATOM     25  CB  SER A  59       7.257  14.957   9.943  1.00 20.00           C
ATOM     26  OG  SER A  59       6.811  15.089   8.604  1.00 20.00           O
ATOM     32  N   ILE A  60       9.964  13.344   9.625  1.00  9.39           N
ATOM     33  CA  ILE A  60      11.067  12.594   9.017  1.00 11.89           C
ATOM     34  C   ILE A  60      10.635  11.189   8.608  1.00  9.81           C
ATOM     35  O   ILE A  60      10.120  10.434   9.430  1.00  8.97           O
ATOM     36  CB  ILE A  60      12.262  12.515   9.987  1.00 20.00           C
ATOM     37  CG1 ILE A  60      12.852  13.907  10.220  1.00 20.00           C
ATOM     38  CG2 ILE A  60      13.321  11.565   9.451  1.00 20.00           C
ATOM     39  CD1 ILE A  60      13.852  13.963  11.353  1.00 20.00           C
ATOM     51  N   VAL A  61      10.858  10.832   7.348  1.00  7.72           N
ATOM     52  CA  VAL A  61      10.510   9.497   6.870  1.00  9.11           C
ATOM     53  C   VAL A  61      11.692   8.812   6.178  1.00 10.61           C
ATOM     54  O   VAL A  61      11.963   9.081   5.006  1.00 11.05           O
ATOM     55  CB  VAL A  61       9.349   9.563   5.856  1.00 10.49           C
ATOM     56  CG1 VAL A  61       9.023   8.166   5.330  1.00  8.40           C
ATOM     57  CG2 VAL A  61       8.127  10.228   6.485  1.00  9.13           C
ATOM     67  N   VAL A  62      12.408   7.927   6.895  1.00 10.62           N
ATOM     68  CA  VAL A  62      13.520   7.240   6.220  1.00  7.26           C
ATOM     69  C   VAL A  62      13.019   6.486   5.000  1.00 10.75           C
ATOM     70  O   VAL A  62      11.946   5.890   5.048  1.00 11.44           O
ATOM     71  CB  VAL A  62      14.239   6.259   7.164  1.00 20.00           C
ATOM     72  CG1 VAL A  62      15.297   5.470   6.407  1.00 20.00           C
ATOM     73  CG2 VAL A  62      14.857   7.004   8.337  1.00 20.00           C
"""

max_distant = """\
ATOM      1  N   TYR A  58       8.659  20.073  11.185  1.00  7.73           N
ATOM      2  CA  TYR A  58       9.250  19.144  10.233  1.00  8.65           C
ATOM      3  C   TYR A  58       9.039  17.721  10.706  1.00  9.84           C
ATOM      4  O   TYR A  58       9.023  17.464  11.919  1.00  8.58           O
ATOM      5  CB  TYR A  58      10.740  19.429  10.045  1.00 20.00           C
ATOM      6  CG  TYR A  58      11.415  18.526   9.036  1.00 20.00           C
ATOM      7  CD1 TYR A  58      11.229  18.716   7.671  1.00 20.00           C
ATOM      8  CD2 TYR A  58      12.237  17.488   9.448  1.00 20.00           C
ATOM      9  CE1 TYR A  58      11.844  17.894   6.747  1.00 20.00           C
ATOM     10  CE2 TYR A  58      12.857  16.660   8.530  1.00 20.00           C
ATOM     11  CZ  TYR A  58      12.657  16.868   7.183  1.00 20.00           C
ATOM     12  OH  TYR A  58      13.270  16.046   6.263  1.00 20.00           O
ATOM     21  N   SER A  59       8.887  16.797   9.756  1.00  9.65           N
ATOM     22  CA  SER A  59       8.720  15.382  10.080  1.00  5.80           C
ATOM     23  C   SER A  59       9.606  14.539   9.169  1.00 10.35           C
ATOM     24  O   SER A  59       9.972  14.972   8.075  1.00 10.56           O
ATOM     25  CB  SER A  59       7.257  14.957   9.943  1.00 20.00           C
ATOM     26  OG  SER A  59       6.438  15.638  10.879  1.00 20.00           O
ATOM     32  N   ILE A  60       9.964  13.344   9.625  1.00  9.39           N
ATOM     33  CA  ILE A  60      11.067  12.594   9.017  1.00 11.89           C
ATOM     34  C   ILE A  60      10.635  11.189   8.608  1.00  9.81           C
ATOM     35  O   ILE A  60      10.120  10.434   9.430  1.00  8.97           O
ATOM     36  CB  ILE A  60      12.262  12.515   9.987  1.00 20.00           C
ATOM     37  CG1 ILE A  60      11.838  11.859  11.302  1.00 20.00           C
ATOM     38  CG2 ILE A  60      12.836  13.900  10.238  1.00 20.00           C
ATOM     39  CD1 ILE A  60      11.874  10.347  11.268  1.00 20.00           C
ATOM     51  N   VAL A  61      10.858  10.832   7.348  1.00  7.72           N
ATOM     52  CA  VAL A  61      10.510   9.497   6.870  1.00  9.11           C
ATOM     53  C   VAL A  61      11.692   8.812   6.178  1.00 10.61           C
ATOM     54  O   VAL A  61      11.963   9.081   5.006  1.00 11.05           O
ATOM     55  CB  VAL A  61       9.349   9.563   5.856  1.00 10.49           C
ATOM     56  CG1 VAL A  61       8.097  10.135   6.519  1.00  8.40           C
ATOM     57  CG2 VAL A  61       9.757  10.373   4.628  1.00  9.13           C
ATOM     67  N   VAL A  62      12.408   7.927   6.895  1.00 10.62           N
ATOM     68  CA  VAL A  62      13.520   7.240   6.220  1.00  7.26           C
ATOM     69  C   VAL A  62      13.019   6.486   5.000  1.00 10.75           C
ATOM     70  O   VAL A  62      11.946   5.890   5.048  1.00 11.44           O
ATOM     71  CB  VAL A  62      14.239   6.259   7.164  1.00 20.00           C
ATOM     72  CG1 VAL A  62      14.849   7.004   8.342  1.00 20.00           C
ATOM     73  CG2 VAL A  62      13.279   5.181   7.645  1.00 20.00           C
"""

min_distant="""\
ATOM      1  N   TYR A  58       8.659  20.073  11.185  1.00  7.73           N
ATOM      2  CA  TYR A  58       9.250  19.144  10.233  1.00  8.65           C
ATOM      3  C   TYR A  58       9.039  17.721  10.706  1.00  9.84           C
ATOM      4  O   TYR A  58       9.023  17.464  11.919  1.00  8.58           O
ATOM      5  CB  TYR A  58      10.740  19.429  10.045  1.00 20.00           C
ATOM      6  CG  TYR A  58      11.032  20.755   9.375  1.00 20.00           C
ATOM      7  CD1 TYR A  58      10.946  20.890   7.995  1.00 20.00           C
ATOM      8  CD2 TYR A  58      11.391  21.866  10.124  1.00 20.00           C
ATOM      9  CE1 TYR A  58      11.212  22.098   7.381  1.00 20.00           C
ATOM     10  CE2 TYR A  58      11.658  23.080   9.517  1.00 20.00           C
ATOM     11  CZ  TYR A  58      11.567  23.190   8.147  1.00 20.00           C
ATOM     12  OH  TYR A  58      11.834  24.395   7.536  1.00 20.00           O
ATOM     21  N   SER A  59       8.887  16.797   9.756  1.00  9.65           N
ATOM     22  CA  SER A  59       8.720  15.382  10.080  1.00  5.80           C
ATOM     23  C   SER A  59       9.606  14.539   9.169  1.00 10.35           C
ATOM     24  O   SER A  59       9.972  14.972   8.075  1.00 10.56           O
ATOM     25  CB  SER A  59       7.257  14.957   9.943  1.00 20.00           C
ATOM     26  OG  SER A  59       7.098  13.577  10.226  1.00 20.00           O
ATOM     32  N   ILE A  60       9.964  13.344   9.625  1.00  9.39           N
ATOM     33  CA  ILE A  60      11.067  12.594   9.017  1.00 11.89           C
ATOM     34  C   ILE A  60      10.635  11.189   8.608  1.00  9.81           C
ATOM     35  O   ILE A  60      10.120  10.434   9.430  1.00  8.97           O
ATOM     36  CB  ILE A  60      12.262  12.515   9.987  1.00 20.00           C
ATOM     37  CG1 ILE A  60      12.803  13.915  10.280  1.00 20.00           C
ATOM     38  CG2 ILE A  60      13.354  11.625   9.415  1.00 20.00           C
ATOM     39  CD1 ILE A  60      12.132  14.595  11.453  1.00 20.00           C
ATOM     51  N   VAL A  61      10.858  10.832   7.348  1.00  7.72           N
ATOM     52  CA  VAL A  61      10.510   9.497   6.870  1.00  9.11           C
ATOM     53  C   VAL A  61      11.692   8.812   6.178  1.00 10.61           C
ATOM     54  O   VAL A  61      11.963   9.081   5.006  1.00 11.05           O
ATOM     55  CB  VAL A  61       9.349   9.563   5.856  1.00 10.49           C
ATOM     56  CG1 VAL A  61       9.784  10.315   4.599  1.00  8.40           C
ATOM     57  CG2 VAL A  61       8.845   8.160   5.528  1.00  9.13           C
ATOM     67  N   VAL A  62      12.408   7.927   6.895  1.00 10.62           N
ATOM     68  CA  VAL A  62      13.520   7.240   6.220  1.00  7.26           C
ATOM     69  C   VAL A  62      13.019   6.486   5.000  1.00 10.75           C
ATOM     70  O   VAL A  62      11.946   5.890   5.048  1.00 11.44           O
ATOM     71  CB  VAL A  62      14.239   6.259   7.164  1.00 20.00           C
ATOM     72  CG1 VAL A  62      13.305   5.129   7.570  1.00 20.00           C
ATOM     73  CG2 VAL A  62      15.496   5.710   6.506  1.00 20.00           C
"""

exact_match="""\
ATOM      1  N   TYR A  58       8.659  20.073  11.185  1.00  7.73           N
ATOM      2  CA  TYR A  58       9.250  19.144  10.233  1.00  8.65           C
ATOM      3  C   TYR A  58       9.039  17.721  10.706  1.00  9.84           C
ATOM      4  O   TYR A  58       9.023  17.464  11.919  1.00  8.58           O
ATOM      5  CB  TYR A  58      10.740  19.429  10.045  1.00 20.00           C
ATOM      6  CG  TYR A  58      11.032  20.755   9.375  1.00 20.00           C
ATOM      7  CD1 TYR A  58      10.946  20.890   7.995  1.00 20.00           C
ATOM      8  CD2 TYR A  58      11.391  21.866  10.124  1.00 20.00           C
ATOM      9  CE1 TYR A  58      11.212  22.098   7.381  1.00 20.00           C
ATOM     10  CE2 TYR A  58      11.658  23.080   9.517  1.00 20.00           C
ATOM     11  CZ  TYR A  58      11.567  23.190   8.147  1.00 20.00           C
ATOM     12  OH  TYR A  58      11.834  24.395   7.536  1.00 20.00           O
ATOM     21  N   SER A  59       8.887  16.797   9.756  1.00  9.65           N
ATOM     22  CA  SER A  59       8.720  15.382  10.080  1.00  5.80           C
ATOM     23  C   SER A  59       9.606  14.539   9.169  1.00 10.35           C
ATOM     24  O   SER A  59       9.972  14.972   8.075  1.00 10.56           O
ATOM     25  CB  SER A  59       7.257  14.957   9.943  1.00 20.00           C
ATOM     26  OG  SER A  59       6.829  15.030   8.594  1.00 20.00           O
ATOM     32  N   ILE A  60       9.964  13.344   9.625  1.00  9.39           N
ATOM     33  CA  ILE A  60      11.067  12.594   9.017  1.00 11.89           C
ATOM     34  C   ILE A  60      10.635  11.189   8.608  1.00  9.81           C
ATOM     35  O   ILE A  60      10.120  10.434   9.430  1.00  8.97           O
ATOM     36  CB  ILE A  60      12.262  12.515   9.987  1.00 20.00           C
ATOM     37  CG1 ILE A  60      12.808  13.915  10.275  1.00 20.00           C
ATOM     38  CG2 ILE A  60      13.351  11.619   9.419  1.00 20.00           C
ATOM     39  CD1 ILE A  60      13.808  13.957  11.408  1.00 20.00           C
ATOM     51  N   VAL A  61      10.858  10.832   7.348  1.00  7.72           N
ATOM     52  CA  VAL A  61      10.510   9.497   6.870  1.00  9.11           C
ATOM     53  C   VAL A  61      11.692   8.812   6.178  1.00 10.61           C
ATOM     54  O   VAL A  61      11.963   9.081   5.006  1.00 11.05           O
ATOM     55  CB  VAL A  61       9.349   9.563   5.856  1.00 10.49           C
ATOM     56  CG1 VAL A  61       9.064   8.175   5.283  1.00  8.40           C
ATOM     57  CG2 VAL A  61       8.107  10.171   6.504  1.00  9.13           C
ATOM     67  N   VAL A  62      12.408   7.927   6.895  1.00 10.62           N
ATOM     68  CA  VAL A  62      13.520   7.240   6.220  1.00  7.26           C
ATOM     69  C   VAL A  62      13.019   6.486   5.000  1.00 10.75           C
ATOM     70  O   VAL A  62      11.946   5.890   5.048  1.00 11.44           O
ATOM     71  CB  VAL A  62      14.239   6.259   7.164  1.00 20.00           C
ATOM     72  CG1 VAL A  62      15.313   5.487   6.412  1.00 20.00           C
ATOM     73  CG2 VAL A  62      14.836   7.001   8.350  1.00 20.00           C
"""

fix_outliers="""\
ATOM      1  N   TYR A  58       8.659  20.073  11.185  1.00  7.73           N
ATOM      2  CA  TYR A  58       9.250  19.144  10.233  1.00  8.65           C
ATOM      3  C   TYR A  58       9.039  17.721  10.706  1.00  9.84           C
ATOM      4  O   TYR A  58       9.023  17.464  11.919  1.00  8.58           O
ATOM      5  CB  TYR A  58      10.740  19.429  10.045  1.00 20.00           C
ATOM      6  CG  TYR A  58      11.032  20.755   9.375  1.00 20.00           C
ATOM      7  CD1 TYR A  58      10.946  20.890   7.995  1.00 20.00           C
ATOM      8  CD2 TYR A  58      11.391  21.866  10.124  1.00 20.00           C
ATOM      9  CE1 TYR A  58      11.212  22.098   7.381  1.00 20.00           C
ATOM     10  CE2 TYR A  58      11.658  23.080   9.517  1.00 20.00           C
ATOM     11  CZ  TYR A  58      11.567  23.190   8.147  1.00 20.00           C
ATOM     12  OH  TYR A  58      11.834  24.395   7.536  1.00 20.00           O
ATOM     21  N   SER A  59       8.887  16.797   9.756  1.00  9.65           N
ATOM     22  CA  SER A  59       8.720  15.382  10.080  1.00  5.80           C
ATOM     23  C   SER A  59       9.606  14.539   9.169  1.00 10.35           C
ATOM     24  O   SER A  59       9.972  14.972   8.075  1.00 10.56           O
ATOM     25  CB  SER A  59       7.257  14.957   9.943  1.00 20.00           C
ATOM     26  OG  SER A  59       6.811  15.089   8.604  1.00 20.00           O
ATOM     32  N   ILE A  60       9.964  13.344   9.625  1.00  9.39           N
ATOM     33  CA  ILE A  60      11.067  12.594   9.017  1.00 11.89           C
ATOM     34  C   ILE A  60      10.635  11.189   8.608  1.00  9.81           C
ATOM     35  O   ILE A  60      10.120  10.434   9.430  1.00  8.97           O
ATOM     36  CB  ILE A  60      12.262  12.515   9.987  1.00 20.00           C
ATOM     37  CG1 ILE A  60      12.852  13.907  10.220  1.00 20.00           C
ATOM     38  CG2 ILE A  60      13.321  11.565   9.451  1.00 20.00           C
ATOM     39  CD1 ILE A  60      13.852  13.963  11.353  1.00 20.00           C
ATOM     51  N   VAL A  61      10.858  10.832   7.348  1.00  7.72           N
ATOM     52  CA  VAL A  61      10.510   9.497   6.870  1.00  9.11           C
ATOM     53  C   VAL A  61      11.692   8.812   6.178  1.00 10.61           C
ATOM     54  O   VAL A  61      11.963   9.081   5.006  1.00 11.05           O
ATOM     55  CB  VAL A  61       9.349   9.563   5.856  1.00 10.49           C
ATOM     56  CG1 VAL A  61       9.023   8.166   5.330  1.00  8.40           C
ATOM     57  CG2 VAL A  61       8.127  10.228   6.485  1.00  9.13           C
ATOM     67  N   VAL A  62      12.408   7.927   6.895  1.00 10.62           N
ATOM     68  CA  VAL A  62      13.520   7.240   6.220  1.00  7.26           C
ATOM     69  C   VAL A  62      13.019   6.486   5.000  1.00 10.75           C
ATOM     70  O   VAL A  62      11.946   5.890   5.048  1.00 11.44           O
ATOM     71  CB  VAL A  62      14.239   6.259   7.164  1.00 20.00           C
ATOM     72  CG1 VAL A  62      15.297   5.470   6.407  1.00 20.00           C
ATOM     73  CG2 VAL A  62      14.857   7.004   8.337  1.00 20.00           C
"""

selection="""\
ATOM      1  N   TYR A  58       8.659  20.073  11.185  1.00  7.73           N
ATOM      2  CA  TYR A  58       9.250  19.144  10.233  1.00  8.65           C
ATOM      3  C   TYR A  58       9.039  17.721  10.706  1.00  9.84           C
ATOM      4  O   TYR A  58       9.023  17.464  11.919  1.00  8.58           O
ATOM      5  CB  TYR A  58      10.740  19.429  10.045  1.00 20.00           C
ATOM      6  CG  TYR A  58      11.155  20.825  10.457  1.00 20.00           C
ATOM      7  CD1 TYR A  58      10.204  21.809  10.700  1.00 20.00           C
ATOM      8  CD2 TYR A  58      12.494  21.155  10.603  1.00 20.00           C
ATOM      9  CE1 TYR A  58      10.579  23.084  11.076  1.00 20.00           C
ATOM     10  CE2 TYR A  58      12.879  22.429  10.980  1.00 20.00           C
ATOM     11  CZ  TYR A  58      11.918  23.388  11.215  1.00 20.00           C
ATOM     12  OH  TYR A  58      12.294  24.659  11.588  1.00 20.00           O
ATOM     21  N   SER A  59       8.887  16.797   9.756  1.00  9.65           N
ATOM     22  CA  SER A  59       8.720  15.382  10.080  1.00  5.80           C
ATOM     23  C   SER A  59       9.606  14.539   9.169  1.00 10.35           C
ATOM     24  O   SER A  59       9.972  14.972   8.075  1.00 10.56           O
ATOM     25  CB  SER A  59       7.257  14.957   9.943  1.00 20.00           C
ATOM     26  OG  SER A  59       6.811  15.089   8.604  1.00 20.00           O
ATOM     32  N   ILE A  60       9.964  13.344   9.625  1.00  9.39           N
ATOM     33  CA  ILE A  60      11.067  12.594   9.017  1.00 11.89           C
ATOM     34  C   ILE A  60      10.635  11.189   8.608  1.00  9.81           C
ATOM     35  O   ILE A  60      10.120  10.434   9.430  1.00  8.97           O
ATOM     36  CB  ILE A  60      12.262  12.515   9.987  1.00 20.00           C
ATOM     37  CG1 ILE A  60      12.852  13.907  10.220  1.00 20.00           C
ATOM     38  CG2 ILE A  60      13.321  11.565   9.451  1.00 20.00           C
ATOM     39  CD1 ILE A  60      13.852  13.963  11.353  1.00 20.00           C
ATOM     51  N   VAL A  61      10.858  10.832   7.348  1.00  7.72           N
ATOM     52  CA  VAL A  61      10.510   9.497   6.870  1.00  9.11           C
ATOM     53  C   VAL A  61      11.692   8.812   6.178  1.00 10.61           C
ATOM     54  O   VAL A  61      11.963   9.081   5.006  1.00 11.05           O
ATOM     55  CB  VAL A  61       9.349   9.563   5.856  1.00 10.49           C
ATOM     56  CG1 VAL A  61       9.023   8.166   5.330  1.00  8.40           C
ATOM     57  CG2 VAL A  61       8.127  10.228   6.485  1.00  9.13           C
ATOM     67  N   VAL A  62      12.408   7.927   6.895  1.00 10.62           N
ATOM     68  CA  VAL A  62      13.520   7.240   6.220  1.00  7.26           C
ATOM     69  C   VAL A  62      13.019   6.486   5.000  1.00 10.75           C
ATOM     70  O   VAL A  62      11.946   5.890   5.048  1.00 11.44           O
ATOM     71  CB  VAL A  62      14.239   6.259   7.164  1.00 20.00           C
ATOM     72  CG1 VAL A  62      15.297   5.470   6.407  1.00 20.00           C
ATOM     73  CG2 VAL A  62      14.857   7.004   8.337  1.00 20.00           C
"""

def core(mode, result, t1, t2):
  prefix = "exercise_%s"%mode
  ph = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  s0 = ph.atoms().extract_xyz()
  ph.write_pdb_file(file_name="%s_in.pdb"%prefix)
  ph = mmtbx.utils.switch_rotamers(pdb_hierarchy=ph, mode=mode)
  ph.write_pdb_file(file_name="%s_out.pdb"%prefix)
  s1 = iotbx.pdb.input(source_info=None,lines=result).atoms().extract_xyz()
  s2 = ph.atoms().extract_xyz()
  d = flex.sqrt((s1 - s2).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, t1, 1.e-3)
  d = flex.sqrt((s2 - s0).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, t2, 0.1)

def exercise_fix_outliers(prefix="exercise_fix_outliers"):
  ph = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  sel = ph.atom_selection_cache().selection(string =
    "resname TYR and not (name O or name CA or name N or name C or name CB)")
  s0 = ph.atoms().extract_xyz()
  ph.write_pdb_file(file_name="%s_in.pdb"%prefix)
  ph = mmtbx.utils.switch_rotamers(pdb_hierarchy=ph, mode="fix_outliers")
  ph.write_pdb_file(file_name="%s_out.pdb"%prefix)
  s1 =iotbx.pdb.input(source_info=None,lines=fix_outliers).atoms().extract_xyz()
  s2 = ph.atoms().extract_xyz()
  # assert fixed do not move
  d = flex.sqrt((s1.select(~sel) - s2.select(~sel)).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [0,0,0])
  d = flex.sqrt((s1.select(~sel) - s0.select(~sel)).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [0,0,0])
  d = flex.sqrt((s2.select(~sel) - s0.select(~sel)).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [0,0,0])
  #
  d = flex.sqrt((s1.select(sel)-s2.select(sel)).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [0,0,0], 1.e-3)
  d = flex.sqrt((s1.select(sel)-s0.select(sel)).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [1.1, 4.0, 2.6], 0.1)
  d = flex.sqrt((s2.select(sel)-s0.select(sel)).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [1.1, 4.0, 2.6], 0.1)

def exercise_selection(prefix="exercise_selection"):
  ph = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  sel = ph.atom_selection_cache().selection(string = "not resname TYR")
  s0 = ph.atoms().extract_xyz()
  ph.write_pdb_file(file_name="%s_in.pdb"%prefix)
  ph = mmtbx.utils.switch_rotamers(pdb_hierarchy=ph, mode="fix_outliers",
    selection = sel)
  ph.write_pdb_file(file_name="%s_out.pdb"%prefix)
  s1 =iotbx.pdb.input(source_info=None,lines=selection).atoms().extract_xyz()
  s2 = ph.atoms().extract_xyz()
  # assert fixed do not move
  d = flex.sqrt((s1 - s2).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [0,0,0])
  d = flex.sqrt((s1 - s0).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [0,0,0])
  d = flex.sqrt((s2 - s0).dot()).min_max_mean().as_tuple()
  assert approx_equal(d, [0,0,0])

if (__name__ == "__main__"):
  core(mode="max_distant", result=max_distant, t1=[0,0,0], t2=[0, 10.2, 1.5])
  core(mode="min_distant", result=min_distant, t1=[0,0,0], t2=[0,  4.1, 0.8])
  core(mode="exact_match", result=exact_match, t1=[0,0,0], t2=[0,  4.1, 0.4])
  exercise_fix_outliers()
  exercise_selection()
