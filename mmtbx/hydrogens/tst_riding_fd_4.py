from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.model
from cctbx.array_family import flex
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from six.moves import zip
from six.moves import range

#-----------------------------------------------------------------------------
# This finite difference test checks transformation of riding H gradients
# for nucleic acids
#-----------------------------------------------------------------------------

def exercise(pdb_str, eps, idealize):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)

  model = mmtbx.model.manager(
            model_input = pdb_inp,
            log         = null_out())
  model.process(make_restraints=True)
  geometry_restraints = model.restraints_manager.geometry
  xray_structure = model.get_xray_structure()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  riding_h_manager.idealize_riding_h_positions(xray_structure=xray_structure)
  sites_cart = xray_structure.sites_cart()

  g_analytical = geometry_restraints.energies_sites(
    sites_cart = sites_cart, compute_gradients = True).gradients
  hd_selection = xray_structure.hd_selection()

  g_analytical_reduced = riding_h_manager.gradients_reduced_cpp(
    gradients    = g_analytical,
    sites_cart   = sites_cart,
    hd_selection = hd_selection)
  #
  ex = [eps,0,0]
  ey = [0,eps,0]
  ez = [0,0,eps]
  g_fd = flex.vec3_double()
  for i_site in range(sites_cart.size()):
    g_fd_i = []
    for e in [ex,ey,ez]:
      ts = []
      for sign in [-1,1]:
        sites_cart_ = sites_cart.deep_copy()
        xray_structure_ = xray_structure.deep_copy_scatterers()
        sites_cart_[i_site] = [
          sites_cart_[i_site][j]+e[j]*sign for j in range(3)]
        xray_structure_.set_sites_cart(sites_cart_)
        # after shift, recalculate H position
        riding_h_manager.idealize_riding_h_positions(
          xray_structure=xray_structure_)
        sites_cart_ = xray_structure_.sites_cart()
        ts.append(geometry_restraints.energies_sites(
          sites_cart = sites_cart_,
          compute_gradients = False).target)
      g_fd_i.append((ts[1]-ts[0])/(2*eps))
    g_fd.append(g_fd_i)

  g_fd_reduced = g_fd.select(~hd_selection)

  for g1, g2 in zip(g_analytical_reduced, g_fd_reduced):
    assert approx_equal(g1,g2, 1.e-4)


# DNA and RNA nucleic acids
pdb_str = """
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
SCALE1      0.033333  0.000000  0.000000        0.00000
SCALE2      0.000000  0.033333  0.000000        0.00000
SCALE3      0.000000  0.000000  0.033333        0.00000
ATOM      1  P    DA A   1       4.321   6.421   9.951  1.00 76.88           P
ATOM      2  OP1  DA A   1       3.975   5.004  10.200  1.00 78.59           O
ATOM      3  OP2  DA A   1       3.621   7.495  10.691  1.00 75.20           O
ATOM      4  O5'  DA A   1       5.895   6.618  10.165  1.00 78.76           O
ATOM      5  C5'  DA A   1       6.704   5.513  10.551  1.00 77.93           C
ATOM      6  C4'  DA A   1       8.179   5.858  10.453  1.00 78.66           C
ATOM      7  O4'  DA A   1       8.521   6.124   9.066  1.00 79.49           O
ATOM      8  C3'  DA A   1       8.598   7.103  11.220  1.00 78.23           C
ATOM      9  O3'  DA A   1       8.961   6.778  12.559  1.00 77.17           O
ATOM     10  C2'  DA A   1       9.776   7.622  10.408  1.00 79.34           C
ATOM     11  C1'  DA A   1       9.411   7.226   8.980  1.00 78.22           C
ATOM     12  N9   DA A   1       8.758   8.292   8.225  1.00 20.00           N
ATOM     13  C8   DA A   1       7.414   8.503   8.091  1.00 20.00           C
ATOM     14  N7   DA A   1       7.112   9.543   7.350  1.00 20.00           N
ATOM     15  C5   DA A   1       8.344  10.050   6.971  1.00 20.00           C
ATOM     16  C6   DA A   1       8.713  11.153   6.176  1.00 20.00           C
ATOM     17  N6   DA A   1       7.829  11.974   5.600  1.00 20.00           N
ATOM     18  N1   DA A   1      10.031  11.381   5.996  1.00 20.00           N
ATOM     19  C2   DA A   1      10.913  10.557   6.575  1.00 20.00           C
ATOM     20  N3   DA A   1      10.687   9.491   7.341  1.00 20.00           N
ATOM     21  C4   DA A   1       9.369   9.291   7.502  1.00 20.00           C
ATOM     22  H5'  DA A   1       6.514   4.760   9.969  1.00 77.93           H
ATOM     23 H5''  DA A   1       6.494   5.269  11.466  1.00 77.93           H
ATOM     24  H4'  DA A   1       8.701   5.104  10.767  1.00 78.66           H
ATOM     25  H3'  DA A   1       7.879   7.754  11.217  1.00 78.23           H
ATOM     26  H2'  DA A   1       9.849   8.586  10.488  1.00 79.34           H
ATOM     27 H2''  DA A   1      10.599   7.190  10.683  1.00 79.34           H
ATOM     28  H1'  DA A   1      10.214   6.951   8.511  1.00 78.22           H
ATOM     29  H8   DA A   1       6.772   7.960   8.488  1.00 20.00           H
ATOM     30  H61  DA A   1       8.103  12.635   5.123  1.00 20.00           H
ATOM     31  H62  DA A   1       6.986  11.841   5.706  1.00 20.00           H
ATOM     32  H2   DA A   1      11.808  10.757   6.421  1.00 20.00           H
ATOM     33  P    DC A   2       8.575   7.768  13.765  1.00 76.88           P
ATOM     34  OP1  DC A   2       8.992   7.120  15.028  1.00 78.59           O
ATOM     35  OP2  DC A   2       7.165   8.176  13.578  1.00 75.20           O
ATOM     36  O5'  DC A   2       9.494   9.054  13.516  1.00 78.76           O
ATOM     37  C5'  DC A   2      10.911   8.922  13.497  1.00 77.93           C
ATOM     38  C4'  DC A   2      11.571  10.212  13.045  1.00 78.66           C
ATOM     39  O4'  DC A   2      11.202  10.488  11.669  1.00 79.49           O
ATOM     40  C3'  DC A   2      11.163  11.450  13.828  1.00 78.23           C
ATOM     41  O3'  DC A   2      11.993  11.627  14.972  1.00 77.17           O
ATOM     42  C2'  DC A   2      11.332  12.566  12.806  1.00 79.34           C
ATOM     43  C1'  DC A   2      10.997  11.881  11.483  1.00 78.22           C
ATOM     44  N1   DC A   2       9.591  12.099  11.036  1.00 20.00           N
ATOM     45  C2   DC A   2       9.254  13.276  10.359  1.00 20.00           C
ATOM     46  O2   DC A   2      10.133  14.120  10.143  1.00 20.00           O
ATOM     47  N3   DC A   2       7.972  13.459   9.959  1.00 20.00           N
ATOM     48  C4   DC A   2       7.054  12.524  10.212  1.00 20.00           C
ATOM     49  N4   DC A   2       5.802  12.749   9.798  1.00 20.00           N
ATOM     50  C5   DC A   2       7.378  11.318  10.900  1.00 20.00           C
ATOM     51  C6   DC A   2       8.646  11.149  11.290  1.00 20.00           C
ATOM     52  H5'  DC A   2      11.156   8.208  12.888  1.00 77.93           H
ATOM     53 H5''  DC A   2      11.222   8.701  14.389  1.00 77.93           H
ATOM     54  H4'  DC A   2      12.534  10.109  13.099  1.00 78.66           H
ATOM     55  H3'  DC A   2      10.233  11.382  14.095  1.00 78.23           H
ATOM     56  H2'  DC A   2      10.711  13.290  12.984  1.00 79.34           H
ATOM     57 H2''  DC A   2      12.247  12.890  12.802  1.00 79.34           H
ATOM     58  H1'  DC A   2      11.602  12.202  10.797  1.00 78.22           H
ATOM     59  H41  DC A   2       5.188  12.166   9.949  1.00 20.00           H
ATOM     60  H42  DC A   2       5.612  13.476   9.380  1.00 20.00           H
ATOM     61  H5   DC A   2       6.731  10.673  11.071  1.00 20.00           H
ATOM     62  H6   DC A   2       8.886  10.372  11.740  1.00 20.00           H
ATOM     63  P    DG A   3      11.380  12.226  16.332  1.00 76.88           P
ATOM     64  OP1  DG A   3      12.439  12.169  17.364  1.00 78.59           O
ATOM     65  OP2  DG A   3      10.081  11.560  16.572  1.00 75.20           O
ATOM     66  O5'  DG A   3      11.087  13.759  15.979  1.00 78.76           O
ATOM     67  C5'  DG A   3      12.153  14.613  15.582  1.00 77.93           C
ATOM     68  C4'  DG A   3      11.625  15.948  15.088  1.00 78.66           C
ATOM     69  O4'  DG A   3      10.821  15.738  13.897  1.00 79.49           O
ATOM     70  C3'  DG A   3      10.718  16.679  16.065  1.00 78.23           C
ATOM     71  O3'  DG A   3      11.476  17.512  16.939  1.00 77.17           O
ATOM     72  C2'  DG A   3       9.805  17.477  15.146  1.00 79.34           C
ATOM     73  C1'  DG A   3       9.678  16.578  13.920  1.00 78.22           C
ATOM     74  N9   DG A   3       8.483  15.738  13.931  1.00 20.00           N
ATOM     75  C8   DG A   3       8.377  14.458  14.419  1.00 20.00           C
ATOM     76  N7   DG A   3       7.181  13.952  14.294  1.00 20.00           N
ATOM     77  C5   DG A   3       6.447  14.961  13.684  1.00 20.00           C
ATOM     78  C6   DG A   3       5.086  14.991  13.297  1.00 20.00           C
ATOM     79  O6   DG A   3       4.233  14.102  13.422  1.00 20.00           O
ATOM     80  N1   DG A   3       4.744  16.208  12.711  1.00 20.00           N
ATOM     81  C2   DG A   3       5.609  17.261  12.523  1.00 20.00           C
ATOM     82  N2   DG A   3       5.096  18.354  11.940  1.00 20.00           N
ATOM     83  N3   DG A   3       6.886  17.245  12.881  1.00 20.00           N
ATOM     84  C4   DG A   3       7.234  16.067  13.454  1.00 20.00           C
ATOM     85  H5'  DG A   3      12.655  14.186  14.870  1.00 77.93           H
ATOM     86 H5''  DG A   3      12.740  14.763  16.340  1.00 77.93           H
ATOM     87  H4'  DG A   3      12.376  16.521  14.865  1.00 78.66           H
ATOM     88  H3'  DG A   3      10.199  16.041  16.579  1.00 78.23           H
ATOM     89  H2'  DG A   3       8.938  17.616  15.561  1.00 79.34           H
ATOM     90 H2''  DG A   3      10.213  18.325  14.909  1.00 79.34           H
ATOM     91  H1'  DG A   3       9.672  17.129  13.121  1.00 78.22           H
ATOM     92  H8   DG A   3       9.088  13.997  14.801  1.00 20.00           H
ATOM     93  H1   DG A   3       3.932  16.308  12.447  1.00 20.00           H
ATOM     94  H21  DG A   3       5.594  19.042  11.803  1.00 20.00           H
ATOM     95  H22  DG A   3       4.269  18.367  11.703  1.00 20.00           H
ATOM     96  P    DT A   4      11.038  17.682  18.476  1.00 76.88           P
ATOM     97  OP1  DT A   4      12.063  18.517  19.142  1.00 78.59           O
ATOM     98  OP2  DT A   4      10.727  16.337  19.007  1.00 75.20           O
ATOM     99  O5'  DT A   4       9.670  18.508  18.398  1.00 78.76           O
ATOM    100  C5'  DT A   4       9.648  19.790  17.782  1.00 77.93           C
ATOM    101  C4'  DT A   4       8.226  20.307  17.653  1.00 78.66           C
ATOM    102  O4'  DT A   4       7.471  19.434  16.774  1.00 79.49           O
ATOM    103  C3'  DT A   4       7.437  20.345  18.951  1.00 78.23           C
ATOM    104  O3'  DT A   4       7.650  21.574  19.642  1.00 77.17           O
ATOM    105  C2'  DT A   4       5.998  20.186  18.477  1.00 79.34           C
ATOM    106  C1'  DT A   4       6.134  19.307  17.235  1.00 78.22           C
ATOM    107  N1   DT A   4       5.854  17.865  17.489  1.00 20.00           N
ATOM    108  C2   DT A   4       4.559  17.407  17.408  1.00 20.00           C
ATOM    109  O2   DT A   4       3.613  18.124  17.136  1.00 20.00           O
ATOM    110  N3   DT A   4       4.408  16.069  17.658  1.00 20.00           N
ATOM    111  C4   DT A   4       5.402  15.162  17.975  1.00 20.00           C
ATOM    112  O4   DT A   4       5.166  13.976  18.183  1.00 20.00           O
ATOM    113  C5   DT A   4       6.738  15.708  18.044  1.00 20.00           C
ATOM    114  C7   DT A   4       7.899  14.820  18.380  1.00 20.00           C
ATOM    115  C6   DT A   4       6.898  17.018  17.800  1.00 20.00           C
ATOM    116  H5'  DT A   4      10.044  19.726  16.899  1.00 77.93           H
ATOM    117 H5''  DT A   4      10.164  20.411  18.319  1.00 77.93           H
ATOM    118  H4'  DT A   4       8.247  21.198  17.270  1.00 78.66           H
ATOM    119  H3'  DT A   4       7.687  19.597  19.515  1.00 78.23           H
ATOM    120 HO3'  DT A   4       7.975  21.569  20.416  1.00 77.17           H
ATOM    121  H2'  DT A   4       5.463  19.742  19.153  1.00 79.34           H
ATOM    122 H2''  DT A   4       5.617  21.047  18.247  1.00 79.34           H
ATOM    123  H1'  DT A   4       5.530  19.632  16.549  1.00 78.22           H
ATOM    124  H3   DT A   4       3.607  15.760  17.614  1.00 20.00           H
ATOM    125  H71  DT A   4       8.506  14.781  17.623  1.00 20.00           H
ATOM    126  H72  DT A   4       8.367  15.177  19.151  1.00 20.00           H
ATOM    127  H73  DT A   4       7.577  13.927  18.582  1.00 20.00           H
ATOM    128  H6   DT A   4       7.756  17.372  17.844  1.00 20.00           H
TER
ATOM    129  P     A B   1      18.553   9.499  13.673  1.00 76.88           P
ATOM    130  OP1   A B   1      18.244   8.077  13.965  1.00 78.59           O
ATOM    131  OP2   A B   1      18.063  10.559  14.590  1.00 75.20           O
ATOM    132  O5'   A B   1      20.135   9.645  13.551  1.00 78.76           O
ATOM    133  C5'   A B   1      20.984   8.516  13.690  1.00 77.93           C
ATOM    134  C4'   A B   1      22.397   8.832  13.268  1.00 78.66           C
ATOM    135  O4'   A B   1      22.420   9.197  11.863  1.00 79.49           O
ATOM    136  C3'   A B   1      23.053  10.012  13.969  1.00 78.23           C
ATOM    137  O3'   A B   1      23.578   9.678  15.242  1.00 77.17           O
ATOM    138  C2'   A B   1      24.112  10.447  12.965  1.00 79.34           C
ATOM    139  O2'   A B   1      25.261   9.616  13.036  1.00 79.34           O
ATOM    140  C1'   A B   1      23.396  10.192  11.639  1.00 78.22           C
ATOM    141  N9    A B   1      22.725  11.403  11.128  1.00 20.00           N
ATOM    142  C8    A B   1      21.387  11.706  11.176  1.00 20.00           C
ATOM    143  N7    A B   1      21.088  12.863  10.637  1.00 20.00           N
ATOM    144  C5    A B   1      22.312  13.355  10.204  1.00 20.00           C
ATOM    145  C6    A B   1      22.676  14.544   9.548  1.00 20.00           C
ATOM    146  N6    A B   1      21.805  15.494   9.199  1.00 20.00           N
ATOM    147  N1    A B   1      23.983  14.727   9.259  1.00 20.00           N
ATOM    148  C2    A B   1      24.857  13.775   9.609  1.00 20.00           C
ATOM    149  N3    A B   1      24.636  12.617  10.227  1.00 20.00           N
ATOM    150  C4    A B   1      23.329  12.466  10.500  1.00 20.00           C
ATOM    151  H5'   A B   1      20.643   7.794  13.139  1.00 77.93           H
ATOM    152 H5''   A B   1      20.986   8.233  14.618  1.00 77.93           H
ATOM    153  H4'   A B   1      22.948   8.043  13.397  1.00 78.66           H
ATOM    154  H3'   A B   1      22.401  10.723  14.070  1.00 78.23           H
ATOM    155  H2'   A B   1      24.339  11.384  13.071  1.00 79.34           H
ATOM    156 HO2'   A B   1      25.651   9.749  13.767  1.00 79.34           H
ATOM    157  H1'   A B   1      24.035   9.875  10.982  1.00 78.22           H
ATOM    158  H8    A B   1      20.752  11.143  11.556  1.00 20.00           H
ATOM    159  H61   A B   1      22.081  16.203   8.797  1.00 20.00           H
ATOM    160  H62   A B   1      20.969  15.398   9.376  1.00 20.00           H
ATOM    161  H2    A B   1      25.744  13.946   9.388  1.00 20.00           H
ATOM    162  P     C B   2      23.325  10.645  16.500  1.00 76.88           P
ATOM    163  OP1   C B   2      23.900   9.992  17.704  1.00 78.59           O
ATOM    164  OP2   C B   2      21.895  11.044  16.495  1.00 75.20           O
ATOM    165  O5'   C B   2      24.199  11.937  16.177  1.00 78.76           O
ATOM    166  C5'   C B   2      25.609  11.848  16.041  1.00 77.93           C
ATOM    167  C4'   C B   2      26.197  13.135  15.518  1.00 78.66           C
ATOM    168  O4'   C B   2      25.672  13.412  14.193  1.00 79.49           O
ATOM    169  C3'   C B   2      25.871  14.389  16.315  1.00 78.23           C
ATOM    170  O3'   C B   2      26.696  14.550  17.457  1.00 77.17           O
ATOM    171  C2'   C B   2      26.036  15.493  15.280  1.00 79.34           C
ATOM    172  O2'   C B   2      27.405  15.824  15.094  1.00 79.34           O
ATOM    173  C1'   C B   2      25.517  14.805  14.017  1.00 78.22           C
ATOM    174  N1    C B   2      24.086  15.094  13.772  1.00 20.00           N
ATOM    175  C2    C B   2      23.736  16.306  13.170  1.00 20.00           C
ATOM    176  O2    C B   2      24.634  17.102  12.859  1.00 20.00           O
ATOM    177  N3    C B   2      22.430  16.578  12.942  1.00 20.00           N
ATOM    178  C4    C B   2      21.494  15.695  13.292  1.00 20.00           C
ATOM    179  N4    C B   2      20.218  16.005  13.048  1.00 20.00           N
ATOM    180  C5    C B   2      21.823  14.452  13.907  1.00 20.00           C
ATOM    181  C6    C B   2      23.119  14.195  14.126  1.00 20.00           C
ATOM    182  H5'   C B   2      25.824  11.130  15.426  1.00 77.93           H
ATOM    183 H5''   C B   2      25.999  11.651  16.907  1.00 77.93           H
ATOM    184  H4'   C B   2      27.161  13.040  15.462  1.00 78.66           H
ATOM    185  H3'   C B   2      24.943  14.354  16.595  1.00 78.23           H
ATOM    186  H2'   C B   2      25.505  16.273  15.503  1.00 79.34           H
ATOM    187 HO2'   C B   2      27.694  16.194  15.791  1.00 79.34           H
ATOM    188  H1'   C B   2      26.041  15.093  13.253  1.00 78.22           H
ATOM    189  H41   C B   2      19.593  15.456  13.265  1.00 20.00           H
ATOM    190  H42   C B   2      20.022  16.754  12.674  1.00 20.00           H
ATOM    191  H5    C B   2      21.163  13.843  14.147  1.00 20.00           H
ATOM    192  H6    C B   2      23.364  13.393  14.526  1.00 20.00           H
ATOM    193  P     G B   3      26.058  15.013  18.858  1.00 76.88           P
ATOM    194  OP1   G B   3      27.123  14.941  19.891  1.00 78.59           O
ATOM    195  OP2   G B   3      24.791  14.266  19.059  1.00 75.20           O
ATOM    196  O5'   G B   3      25.696  16.547  18.630  1.00 78.76           O
ATOM    197  C5'   G B   3      26.713  17.507  18.380  1.00 77.93           C
ATOM    198  C4'   G B   3      26.131  18.843  17.993  1.00 78.66           C
ATOM    199  O4'   G B   3      25.392  18.715  16.750  1.00 79.49           O
ATOM    200  C3'   G B   3      25.121  19.437  18.964  1.00 78.23           C
ATOM    201  O3'   G B   3      25.730  20.093  20.063  1.00 77.17           O
ATOM    202  C2'   G B   3      24.308  20.364  18.072  1.00 79.34           C
ATOM    203  O2'   G B   3      24.988  21.590  17.846  1.00 79.34           O
ATOM    204  C1'   G B   3      24.268  19.570  16.766  1.00 78.22           C
ATOM    205  N9    G B   3      23.049  18.747  16.654  1.00 20.00           N
ATOM    206  C8    G B   3      22.948  17.389  16.836  1.00 20.00           C
ATOM    207  N7    G B   3      21.735  16.938  16.671  1.00 20.00           N
ATOM    208  C5    G B   3      20.988  18.066  16.360  1.00 20.00           C
ATOM    209  C6    G B   3      19.605  18.202  16.074  1.00 20.00           C
ATOM    210  O6    G B   3      18.735  17.324  16.037  1.00 20.00           O
ATOM    211  N1    G B   3      19.265  19.526  15.812  1.00 20.00           N
ATOM    212  C2    G B   3      20.142  20.583  15.823  1.00 20.00           C
ATOM    213  N2    G B   3      19.620  21.786  15.544  1.00 20.00           N
ATOM    214  N3    G B   3      21.433  20.469  16.088  1.00 20.00           N
ATOM    215  C4    G B   3      21.785  19.192  16.345  1.00 20.00           C
ATOM    216  H5'   G B   3      27.279  17.188  17.659  1.00 77.93           H
ATOM    217 H5''   G B   3      27.250  17.615  19.180  1.00 77.93           H
ATOM    218  H4'   G B   3      26.854  19.476  17.865  1.00 78.66           H
ATOM    219  H3'   G B   3      24.546  18.729  19.295  1.00 78.23           H
ATOM    220  H2'   G B   3      23.416  20.507  18.427  1.00 79.34           H
ATOM    221 HO2'   G B   3      25.624  21.654  18.390  1.00 79.34           H
ATOM    222  H1'   G B   3      24.316  20.182  16.014  1.00 78.22           H
ATOM    223  H8    G B   3      23.671  16.847  17.056  1.00 20.00           H
ATOM    224  H1    G B   3      18.442  19.695  15.629  1.00 20.00           H
ATOM    225  H21   G B   3      18.776  21.865  15.394  1.00 20.00           H
ATOM    226  H22   G B   3      20.128  22.479  15.515  1.00 20.00           H
ATOM    227  P     U B   4      25.038  20.052  21.513  1.00 76.88           P
ATOM    228  OP1   U B   4      25.946  20.728  22.475  1.00 78.59           O
ATOM    229  OP2   U B   4      24.599  18.658  21.776  1.00 75.20           O
ATOM    230  O5'   U B   4      23.737  20.955  21.346  1.00 78.76           O
ATOM    231  C5'   U B   4      23.842  22.322  20.974  1.00 77.93           C
ATOM    232  C4'   U B   4      22.492  22.914  20.655  1.00 78.66           C
ATOM    233  O4'   U B   4      21.908  22.220  19.522  1.00 79.49           O
ATOM    234  C3'   U B   4      21.438  22.796  21.745  1.00 78.23           C
ATOM    235  O3'   U B   4      21.564  23.790  22.747  1.00 77.17           O
ATOM    236  C2'   U B   4      20.136  22.882  20.959  1.00 79.34           C
ATOM    237  O2'   U B   4      19.818  24.228  20.636  1.00 79.34           O
ATOM    238  C1'   U B   4      20.506  22.137  19.676  1.00 78.22           C
ATOM    239  N1    U B   4      20.119  20.709  19.729  1.00 20.00           N
ATOM    240  C2    U B   4      18.803  20.391  19.462  1.00 20.00           C
ATOM    241  O2    U B   4      17.962  21.229  19.188  1.00 20.00           O
ATOM    242  N3    U B   4      18.508  19.052  19.527  1.00 20.00           N
ATOM    243  C4    U B   4      19.379  18.023  19.826  1.00 20.00           C
ATOM    244  O4    U B   4      18.961  16.865  19.845  1.00 20.00           O
ATOM    245  C5    U B   4      20.725  18.433  20.091  1.00 20.00           C
ATOM    246  C6    U B   4      21.037  19.732  20.034  1.00 20.00           C
ATOM    247  H5'   U B   4      24.413  22.395  20.194  1.00 77.93           H
ATOM    248 H5''   U B   4      24.240  22.819  21.706  1.00 77.93           H
ATOM    249  H4'   U B   4      22.605  23.850  20.427  1.00 78.66           H
ATOM    250  H3'   U B   4      21.503  21.920  22.156  1.00 78.23           H
ATOM    251 HO3'   U B   4      20.932  24.324  22.892  1.00 77.17           H
ATOM    252  H2'   U B   4      19.409  22.442  21.427  1.00 79.34           H
ATOM    253 HO2'   U B   4      19.611  24.636  21.340  1.00 79.34           H
ATOM    254  H1'   U B   4      20.072  22.562  18.920  1.00 78.22           H
ATOM    255  H5    U B   4      21.375  17.802  20.301  1.00 20.00           H
ATOM    256  H3    U B   4      17.694  18.830  19.364  1.00 20.00           H
ATOM    257  H6    U B   4      21.915  19.985  20.208  1.00 20.00           H
TER
END
"""

def run():
  for idealize in [True,False]:
    exercise(pdb_str=pdb_str, eps=1.e-4, idealize=idealize)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time:", round(time.time()-t0, 2), "seconds")
