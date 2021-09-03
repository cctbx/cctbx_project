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
# for all amino acids, one by one
#-----------------------------------------------------------------------------
#
def exercise(pdb_str, eps):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)

  model = mmtbx.model.manager(
            model_input = pdb_inp,
            log         = null_out())
  model.process(make_restraints=True)
  geometry_restraints = model.restraints_manager.geometry
  xray_structure = model.get_xray_structure()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  riding_h_manager.idealize_riding_h_positions(xray_structure = xray_structure)

  sites_cart = xray_structure.sites_cart()

  g_analytical = geometry_restraints.energies_sites(
    sites_cart = sites_cart,
    compute_gradients = True).gradients

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


#Alanine
pdb_str_01 = """
CRYST1   11.863   12.347   12.550  90.00  90.00  90.00 P 1
SCALE1      0.084296  0.000000  0.000000        0.00000
SCALE2      0.000000  0.080991  0.000000        0.00000
SCALE3      0.000000  0.000000  0.079681        0.00000
ATOM      1  N   ALA B  84       7.017   5.016   7.287  1.00 10.00           N
ATOM      2  CA  ALA B  84       5.784   5.706   7.445  1.00 10.00           C
ATOM      3  C   ALA B  84       5.129   5.694   6.011  1.00 10.00           C
ATOM      4  O   ALA B  84       5.714   5.812   4.900  1.00 10.00           O
ATOM      5  CB  ALA B  84       6.261   7.426   7.444  1.00 10.00           C
ATOM      6  HA  ALA B  84       5.269   5.787   7.982  1.00 10.00           H
ATOM      7  HB1 ALA B  84       5.337   8.016   7.304  1.00 10.00           H
ATOM      8  HB2 ALA B  84       6.549   7.313   8.459  1.00 10.00           H
ATOM      9  HB3 ALA B  84       6.616   7.667   6.795  1.00 10.00           H
"""
# Arginine:
pdb_str_02 = """
CRYST1   14.446   16.451   11.913  90.00  90.00  90.00 P 1
SCALE1      0.069223  0.000000  0.000000        0.00000
SCALE2      0.000000  0.060787  0.000000        0.00000
SCALE3      0.000000  0.000000  0.083942        0.00000
ATOM      1  N   ARG A  50       6.642  11.452   6.691  1.00 10.00           N
ATOM      2  CA  ARG A  50       7.570  10.240   6.682  1.00 10.00           C
ATOM      3  C   ARG A  50       9.005  10.280   6.707  1.00 10.00           C
ATOM      4  O   ARG A  50       9.602  11.108   6.810  1.00 10.00           O
ATOM      5  CB  ARG A  50       6.962   9.371   5.075  1.00 10.00           C
ATOM      6  CG  ARG A  50       5.382   9.275   4.806  1.00 10.00           C
ATOM      7  CD  ARG A  50       4.830   7.928   6.084  1.00 10.00           C
ATOM      8  NE  ARG A  50       5.434   6.841   5.374  1.00 10.00           N
ATOM      9  CZ  ARG A  50       6.450   6.116   6.035  1.00 10.00           C
ATOM     10  NH1 ARG A  50       7.392   6.694   6.803  1.00 10.00           N
ATOM     11  NH2 ARG A  50       6.898   5.026   5.305  1.00 10.00           N
ATOM     12  HA  ARG A  50       6.817   9.405   7.203  1.00 10.00           H
ATOM     13  HB2 ARG A  50       7.451   9.961   4.409  1.00 10.00           H
ATOM     14  HB3 ARG A  50       7.391   8.588   5.297  1.00 10.00           H
ATOM     15  HG2 ARG A  50       5.016  10.089   5.174  1.00 10.00           H
ATOM     16  HG3 ARG A  50       5.385   9.031   4.141  1.00 10.00           H
ATOM     17  HD2 ARG A  50       5.268   8.147   6.643  1.00 10.00           H
ATOM     18  HD3 ARG A  50       3.909   8.226   5.899  1.00 10.00           H
ATOM     19  HE  ARG A  50       5.208   6.478   4.962  1.00 10.00           H
ATOM     20 HH11 ARG A  50       6.970   7.301   7.302  1.00 10.00           H
ATOM     21 HH12 ARG A  50       7.866   6.356   7.069  1.00 10.00           H
ATOM     22 HH21 ARG A  50       6.604   4.626   5.005  1.00 10.00           H
ATOM     23 HH22 ARG A  50       7.522   4.396   5.734  1.00 10.00           H
"""
#Asparagine
pdb_str_03 = """
CRYST1   13.889   13.738   13.126  90.00  90.00  90.00 P 1
SCALE1      0.071999  0.000000  0.000000        0.00000
SCALE2      0.000000  0.072791  0.000000        0.00000
SCALE3      0.000000  0.000000  0.076185        0.00000
ATOM      1  N   ASN A  35       5.108   6.598   6.668  1.00 10.00           N
ATOM      2  CA  ASN A  35       6.071   7.634   7.159  1.00 10.00           C
ATOM      3  C   ASN A  35       6.113   8.702   6.163  1.00 10.00           C
ATOM      4  O   ASN A  35       5.556   8.558   4.817  1.00 10.00           O
ATOM      5  CB  ASN A  35       7.439   6.520   6.538  1.00 10.00           C
ATOM      6  CG  ASN A  35       7.597   5.465   7.487  1.00 10.00           C
ATOM      7  OD1 ASN A  35       6.743   4.980   8.032  1.00 10.00           O
ATOM      8  ND2 ASN A  35       8.841   4.827   7.318  1.00 10.00           N
ATOM      9  HA  ASN A  35       6.387   7.612   7.973  1.00 10.00           H
ATOM     10  HB2 ASN A  35       7.452   6.464   5.548  1.00 10.00           H
ATOM     11  HB3 ASN A  35       8.153   7.444   6.485  1.00 10.00           H
ATOM     12 HD21 ASN A  35       9.022   4.464   7.804  1.00 10.00           H
ATOM     13 HD22 ASN A  35       9.521   5.284   7.192  1.00 10.00           H
"""
# Aspartic Acid
pdb_str_04 = """
CRYST1   14.178   13.514   11.978  90.00  90.00  90.00 P 1
SCALE1      0.070532  0.000000  0.000000        0.00000
SCALE2      0.000000  0.073997  0.000000        0.00000
SCALE3      0.000000  0.000000  0.083486        0.00000
ATOM      1  N   ASP A 121       4.844   6.380   4.967  1.00 10.00           N
ATOM      2  CA  ASP A 121       6.237   6.469   5.750  1.00 10.00           C
ATOM      3  C   ASP A 121       5.959   7.525   7.138  1.00 10.00           C
ATOM      4  O   ASP A 121       5.683   8.452   6.797  1.00 10.00           O
ATOM      5  CB  ASP A 121       7.582   6.234   5.277  1.00 10.00           C
ATOM      6  CG  ASP A 121       8.550   6.163   6.136  1.00 10.00           C
ATOM      7  OD1 ASP A 121       9.035   5.174   6.337  1.00 10.00           O
ATOM      8  OD2 ASP A 121       9.115   7.253   6.644  1.00 10.00           O
ATOM      9  HA  ASP A 121       6.173   5.442   6.266  1.00 10.00           H
ATOM     10  HB2 ASP A 121       7.509   5.748   4.630  1.00 10.00           H
ATOM     11  HB3 ASP A 121       7.714   7.340   4.588  1.00 10.00           H
"""
# Cysteine
pdb_str_05 = """
CRYST1   12.081   11.399   14.639  90.00  90.00  90.00 P 1
SCALE1      0.082775  0.000000  0.000000        0.00000
SCALE2      0.000000  0.087727  0.000000        0.00000
SCALE3      0.000000  0.000000  0.068311        0.00000
ATOM      1  N   CYS A  96       4.823   5.221   7.409  1.00 10.00           N
ATOM      2  CA  CYS A  96       6.273   5.028   7.456  1.00 10.00           C
ATOM      3  C   CYS A  96       6.932   5.026   9.008  1.00 10.00           C
ATOM      4  O   CYS A  96       7.163   5.942   9.454  1.00 10.00           O
ATOM      5  CB  CYS A  96       6.905   6.223   6.813  1.00 10.00           C
ATOM      6  SG  CYS A  96       6.265   6.407   5.142  1.00 10.00           S
ATOM      7  HA  CYS A  96       6.714   4.104   6.960  1.00 10.00           H
ATOM      8  HB2 CYS A  96       6.849   7.034   6.999  1.00 10.00           H
ATOM      9  HB3 CYS A  96       8.130   6.106   6.575  1.00 10.00           H
ATOM     10  HG  CYS A  96       6.882   7.238   4.630  1.00 10.00           H
"""
# Glutamine
pdb_str_06 = """
CRYST1   14.547   13.375   15.374  90.00  90.00  90.00 P 1
SCALE1      0.068743  0.000000  0.000000        0.00000
SCALE2      0.000000  0.074766  0.000000        0.00000
SCALE3      0.000000  0.000000  0.065045        0.00000
ATOM      1  N   GLN A  70       6.344   7.129   5.035  1.00 10.00           N
ATOM      2  CA  GLN A  70       6.644   6.208   6.336  1.00 10.00           C
ATOM      3  C   GLN A  70       5.352   6.178   6.970  1.00 10.00           C
ATOM      4  O   GLN A  70       5.126   4.898   7.373  1.00 10.00           O
ATOM      5  CB  GLN A  70       7.877   6.998   7.230  1.00 10.00           C
ATOM      6  CG  GLN A  70       7.974   6.681   8.575  1.00 10.00           C
ATOM      7  CD  GLN A  70       9.249   7.237   9.063  1.00 10.00           C
ATOM      8  OE1 GLN A  70       9.559   8.209   8.935  1.00 10.00           O
ATOM      9  NE2 GLN A  70       9.506   6.762  10.402  1.00 10.00           N
ATOM     10  HA  GLN A  70       7.249   5.584   6.195  1.00 10.00           H
ATOM     11  HB2 GLN A  70       8.663   7.126   6.476  1.00 10.00           H
ATOM     12  HB3 GLN A  70       7.331   8.179   7.116  1.00 10.00           H
ATOM     13  HG2 GLN A  70       7.234   6.385   8.986  1.00 10.00           H
ATOM     14  HG3 GLN A  70       8.549   5.809   8.279  1.00 10.00           H
ATOM     15 HE21 GLN A  70       9.114   6.119  10.657  1.00 10.00           H
ATOM     16 HE22 GLN A  70      10.235   7.102  10.671  1.00 10.00           H
"""
# Glutamic Acid
pdb_str_07 = """
CRYST1   15.127   14.767   13.161  90.00  90.00  90.00 P 1
SCALE1      0.066107  0.000000  0.000000        0.00000
SCALE2      0.000000  0.067719  0.000000        0.00000
SCALE3      0.000000  0.000000  0.075982        0.00000
ATOM      1  N   GLU A  24       9.022   8.042   4.852  1.00 10.00           N
ATOM      2  CA  GLU A  24       8.611   8.384   6.596  1.00 10.00           C
ATOM      3  C   GLU A  24       9.851   8.805   7.028  1.00 10.00           C
ATOM      4  O   GLU A  24      10.208   9.672   7.827  1.00 10.00           O
ATOM      5  CB  GLU A  24       8.163   7.285   7.084  1.00 10.00           C
ATOM      6  CG  GLU A  24       6.610   6.769   6.720  1.00 10.00           C
ATOM      7  CD  GLU A  24       6.069   5.633   7.216  1.00 10.00           C
ATOM      8  OE1 GLU A  24       5.180   5.275   7.007  1.00 10.00           O
ATOM      9  OE2 GLU A  24       6.826   5.071   8.182  1.00 10.00           O
ATOM     10  HA  GLU A  24       8.239   8.942   6.618  1.00 10.00           H
ATOM     11  HB2 GLU A  24       8.834   6.577   7.107  1.00 10.00           H
ATOM     12  HB3 GLU A  24       8.074   7.246   7.880  1.00 10.00           H
ATOM     13  HG2 GLU A  24       5.925   7.518   6.600  1.00 10.00           H
ATOM     14  HG3 GLU A  24       6.944   6.683   5.520  1.00 10.00           H
"""
# Glycine
pdb_str_08 = """
CRYST1   13.147   11.167   10.709  90.00  90.00  90.00 P 1
SCALE1      0.076063  0.000000  0.000000        0.00000
SCALE2      0.000000  0.089550  0.000000        0.00000
SCALE3      0.000000  0.000000  0.093379        0.00000
ATOM      1  N   GLY B  93       8.304   6.036   5.544  1.00 10.00           N
ATOM      2  CA  GLY B  93       7.319   4.905   5.595  1.00 10.00           C
ATOM      3  C   GLY B  93       5.882   5.364   4.847  1.00 10.00           C
ATOM      4  O   GLY B  93       5.143   5.776   5.738  1.00 10.00           O
ATOM      5  HA2 GLY B  93       7.870   4.339   5.372  1.00 10.00           H
ATOM      6  HA3 GLY B  93       7.123   4.618   6.529  1.00 10.00           H
"""

# Histidine
pdb_str_09 = """
CRYST1   12.286   14.857   15.512  90.00  90.00  90.00 P 1
SCALE1      0.081393  0.000000  0.000000        0.00000
SCALE2      0.000000  0.067308  0.000000        0.00000
SCALE3      0.000000  0.000000  0.064466        0.00000
ATOM      1  N   HIS A  34       7.292   7.714   9.065  1.00 10.00           N
ATOM      2  CA  HIS A  34       6.363   8.592   8.255  1.00 10.00           C
ATOM      3  C   HIS A  34       6.024   9.792   9.401  1.00 10.00           C
ATOM      4  O   HIS A  34       6.246   9.507  10.701  1.00 10.00           O
ATOM      5  CB  HIS A  34       4.951   8.082   8.105  1.00 10.00           C
ATOM      6  CG  HIS A  34       5.323   6.839   7.073  1.00 10.00           C
ATOM      7  ND1 HIS A  34       5.021   7.108   5.763  1.00 10.00           N
ATOM      8  CD2 HIS A  34       5.317   5.623   7.295  1.00 10.00           C
ATOM      9  CE1 HIS A  34       5.333   5.771   5.081  1.00 10.00           C
ATOM     10  NE2 HIS A  34       5.601   4.927   6.077  1.00 10.00           N
ATOM     11  HA  HIS A  34       6.539   9.000   7.667  1.00 10.00           H
ATOM     12  HB2 HIS A  34       4.722   7.760   8.940  1.00 10.00           H
ATOM     13  HB3 HIS A  34       4.532   8.544   7.770  1.00 10.00           H
ATOM     14  HD2 HIS A  34       5.433   5.063   7.987  1.00 10.00           H
ATOM     15  HE1 HIS A  34       5.162   5.891   3.986  1.00 10.00           H
"""
# Isoleucine
pdb_str_10 = """
CRYST1   14.040   13.236   13.685  90.00  90.00  90.00 P 1
SCALE1      0.071225  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075552  0.000000        0.00000
SCALE3      0.000000  0.000000  0.073073        0.00000
ATOM      1  N   ILE B  81       6.774   6.654   8.528  1.00 10.00           N
ATOM      2  CA  ILE B  81       6.248   6.915   7.265  1.00 10.00           C
ATOM      3  C   ILE B  81       5.064   5.915   6.912  1.00 10.00           C
ATOM      4  O   ILE B  81       4.891   4.882   7.591  1.00 10.00           O
ATOM      5  CB  ILE B  81       7.374   7.060   6.201  1.00 10.00           C
ATOM      6  CG1 ILE B  81       7.928   5.644   6.013  1.00 10.00           C
ATOM      7  CG2 ILE B  81       8.450   8.146   6.766  1.00 10.00           C
ATOM      8  CD1 ILE B  81       9.014   5.797   4.814  1.00 10.00           C
ATOM      9  HA  ILE B  81       5.805   7.933   7.605  1.00 10.00           H
ATOM     10  HB  ILE B  81       6.968   7.172   5.320  1.00 10.00           H
ATOM     11 HG12 ILE B  81       8.488   5.662   7.117  1.00 10.00           H
ATOM     12 HG13 ILE B  81       7.287   4.921   5.901  1.00 10.00           H
ATOM     13 HG21 ILE B  81       8.841   8.468   5.740  1.00 10.00           H
ATOM     14 HG22 ILE B  81       7.916   9.169   6.877  1.00 10.00           H
ATOM     15 HG23 ILE B  81       8.650   8.115   7.619  1.00 10.00           H
ATOM     16 HD11 ILE B  81       9.263   4.839   4.994  1.00 10.00           H
ATOM     17 HD12 ILE B  81       8.679   5.953   4.216  1.00 10.00           H
ATOM     18 HD13 ILE B  81       9.551   6.406   5.137  1.00 10.00           H
"""

# Leucine
pdb_str_11 = """
CRYST1   13.550   15.138   12.702  90.00  90.00  90.00 P 1
SCALE1      0.073801  0.000000  0.000000        0.00000
SCALE2      0.000000  0.066059  0.000000        0.00000
SCALE3      0.000000  0.000000  0.078728        0.00000
ATOM      1  N   LEU B  94       6.273   6.593   4.935  1.00 10.00           N
ATOM      2  CA  LEU B  94       7.641   7.224   5.628  1.00 10.00           C
ATOM      3  C   LEU B  94       8.422   6.084   6.113  1.00 10.00           C
ATOM      4  O   LEU B  94       8.638   5.157   5.672  1.00 10.00           O
ATOM      5  CB  LEU B  94       7.066   8.117   6.861  1.00 10.00           C
ATOM      6  CG  LEU B  94       6.206   9.526   6.450  1.00 10.00           C
ATOM      7  CD1 LEU B  94       5.918  10.237   7.727  1.00 10.00           C
ATOM      8  CD2 LEU B  94       4.891   9.251   5.577  1.00 10.00           C
ATOM      9  HA  LEU B  94       7.928   7.858   4.837  1.00 10.00           H
ATOM     10  HB2 LEU B  94       6.600   7.347   7.185  1.00 10.00           H
ATOM     11  HB3 LEU B  94       7.833   8.235   7.219  1.00 10.00           H
ATOM     12  HG  LEU B  94       7.048   9.715   5.954  1.00 10.00           H
ATOM     13 HD11 LEU B  94       5.564  11.132   7.460  1.00 10.00           H
ATOM     14 HD12 LEU B  94       5.689   9.764   8.196  1.00 10.00           H
ATOM     15 HD13 LEU B  94       6.901  10.482   8.316  1.00 10.00           H
ATOM     16 HD21 LEU B  94       4.622   9.894   5.701  1.00 10.00           H
ATOM     17 HD22 LEU B  94       5.223   8.782   4.700  1.00 10.00           H
ATOM     18 HD23 LEU B  94       4.327   8.453   6.346  1.00 10.00           H
"""

# Lysine
pdb_str_12 = """
CRYST1   14.470   17.497   11.416  90.00  90.00  90.00 P 1
SCALE1      0.069109  0.000000  0.000000        0.00000
SCALE2      0.000000  0.057153  0.000000        0.00000
SCALE3      0.000000  0.000000  0.087596        0.00000
ATOM      1  N   LYS A  76       4.981  10.159   4.860  1.00 10.00           N
ATOM      2  CA  LYS A  76       5.876  10.071   6.050  1.00 10.00           C
ATOM      3  C   LYS A  76       6.304  11.829   6.325  1.00 10.00           C
ATOM      4  O   LYS A  76       5.161  12.466   6.576  1.00 10.00           O
ATOM      5  CB  LYS A  76       7.337   9.388   5.906  1.00 10.00           C
ATOM      6  CG  LYS A  76       7.118   7.831   5.788  1.00 10.00           C
ATOM      7  CD  LYS A  76       8.453   7.428   5.736  1.00 10.00           C
ATOM      8  CE  LYS A  76       8.134   5.835   5.433  1.00 10.00           C
ATOM      9  NZ  LYS A  76       9.426   4.919   5.392  1.00 10.00           N
ATOM     10  HA  LYS A  76       5.397  10.002   6.803  1.00 10.00           H
ATOM     11  HB2 LYS A  76       7.655   9.920   5.112  1.00 10.00           H
ATOM     12  HB3 LYS A  76       7.624   9.427   6.543  1.00 10.00           H
ATOM     13  HG2 LYS A  76       6.592   7.548   6.676  1.00 10.00           H
ATOM     14  HG3 LYS A  76       6.469   7.536   5.090  1.00 10.00           H
ATOM     15  HD2 LYS A  76       8.768   7.405   4.630  1.00 10.00           H
ATOM     16  HD3 LYS A  76       9.025   7.220   6.546  1.00 10.00           H
ATOM     17  HE2 LYS A  76       7.580   5.564   6.304  1.00 10.00           H
ATOM     18  HE3 LYS A  76       7.440   5.632   4.831  1.00 10.00           H
ATOM     19  HZ1 LYS A  76       9.379   4.186   5.442  1.00 10.00           H
ATOM     20  HZ2 LYS A  76       9.725   5.254   4.549  1.00 10.00           H
ATOM     21  HZ3 LYS A  76       9.965   5.006   6.093  1.00 10.00           H
"""

# Methionine
pdb_str_13 = """
CRYST1   13.775   12.351   15.645  90.00  90.00  90.00 P 1
SCALE1      0.072595  0.000000  0.000000        0.00000
SCALE2      0.000000  0.080965  0.000000        0.00000
SCALE3      0.000000  0.000000  0.063918        0.00000
ATOM      1  N   MET B  37       7.603   5.375   6.446  1.00 10.00           N
ATOM      2  CA  MET B  37       6.723   6.349   6.607  1.00 10.00           C
ATOM      3  C   MET B  37       6.323   7.146   5.241  1.00 10.00           C
ATOM      4  O   MET B  37       5.062   7.131   4.970  1.00 10.00           O
ATOM      5  CB  MET B  37       6.932   7.378   7.483  1.00 10.00           C
ATOM      6  CG  MET B  37       7.258   6.900   9.095  1.00 10.00           C
ATOM      7  SD  MET B  37       8.785   5.488   9.096  1.00 10.00           S
ATOM      8  CE  MET B  37       8.590   5.149  10.517  1.00 10.00           C
ATOM      9  HA  MET B  37       5.885   6.031   6.839  1.00 10.00           H
ATOM     10  HB2 MET B  37       7.862   7.670   7.351  1.00 10.00           H
ATOM     11  HB3 MET B  37       6.364   7.936   7.854  1.00 10.00           H
ATOM     12  HG2 MET B  37       7.509   7.571   9.764  1.00 10.00           H
ATOM     13  HG3 MET B  37       6.614   6.240   9.246  1.00 10.00           H
ATOM     14  HE1 MET B  37       9.624   4.355  10.601  1.00 10.00           H
ATOM     15  HE2 MET B  37       9.070   5.811  11.280  1.00 10.00           H
ATOM     16  HE3 MET B  37       8.003   4.499  11.109  1.00 10.00           H
"""

# Phenylanaline
pdb_str_14 = """
CRYST1   11.802   16.578   13.976  90.00  90.00  90.00 P 1
SCALE1      0.084731  0.000000  0.000000        0.00000
SCALE2      0.000000  0.060321  0.000000        0.00000
SCALE3      0.000000  0.000000  0.071551        0.00000
ATOM      1  N   PHE A  63       6.412   9.770   7.572  1.00 10.00           N
ATOM      2  CA  PHE A  63       6.289   9.157   6.198  1.00 10.00           C
ATOM      3  C   PHE A  63       6.292  10.223   5.110  1.00 10.00           C
ATOM      4  O   PHE A  63       6.451  11.437   5.241  1.00 10.00           O
ATOM      5  CB  PHE A  63       5.580   8.319   6.063  1.00 10.00           C
ATOM      6  CG  PHE A  63       5.743   6.958   6.849  1.00 10.00           C
ATOM      7  CD1 PHE A  63       6.569   6.117   6.768  1.00 10.00           C
ATOM      8  CD2 PHE A  63       5.050   7.055   8.045  1.00 10.00           C
ATOM      9  CE1 PHE A  63       6.802   5.054   7.620  1.00 10.00           C
ATOM     10  CE2 PHE A  63       5.331   6.153   8.840  1.00 10.00           C
ATOM     11  CZ  PHE A  63       6.270   5.103   8.840  1.00 10.00           C
ATOM     12  HA  PHE A  63       7.224   9.004   5.959  1.00 10.00           H
ATOM     13  HB2 PHE A  63       4.619   8.277   5.936  1.00 10.00           H
ATOM     14  HB3 PHE A  63       5.713   7.612   5.105  1.00 10.00           H
ATOM     15  HD1 PHE A  63       7.301   6.004   5.919  1.00 10.00           H
ATOM     16  HD2 PHE A  63       4.436   7.641   8.314  1.00 10.00           H
ATOM     17  HE1 PHE A  63       7.224   4.307   7.397  1.00 10.00           H
ATOM     18  HE2 PHE A  63       4.727   6.172   9.862  1.00 10.00           H
ATOM     19  HZ  PHE A  63       6.045   4.444   9.418  1.00 10.00           H
"""

# Proline
pdb_str_15 = """
CRYST1   12.293   14.006   12.486  90.00  90.00  90.00 P 1
SCALE1      0.081347  0.000000  0.000000        0.00000
SCALE2      0.000000  0.071398  0.000000        0.00000
SCALE3      0.000000  0.000000  0.080090        0.00000
ATOM      1  N   PRO A   4       4.920   7.325   6.411  1.00 10.00           N
ATOM      2  CA  PRO A   4       6.162   6.834   5.666  1.00 10.00           C
ATOM      3  C   PRO A   4       6.405   5.341   5.894  1.00 10.00           C
ATOM      4  O   PRO A   4       7.098   4.853   4.957  1.00 10.00           O
ATOM      5  CB  PRO A   4       7.380   7.894   6.618  1.00 10.00           C
ATOM      6  CG  PRO A   4       6.640   8.868   7.204  1.00 10.00           C
ATOM      7  CD  PRO A   4       5.415   8.476   7.560  1.00 10.00           C
ATOM      8  HA  PRO A   4       6.250   7.067   4.776  1.00 10.00           H
ATOM      9  HB2 PRO A   4       7.921   7.337   7.346  1.00 10.00           H
ATOM     10  HB3 PRO A   4       8.016   8.080   5.830  1.00 10.00           H
ATOM     11  HG2 PRO A   4       7.167   9.579   7.723  1.00 10.00           H
ATOM     12  HG3 PRO A   4       6.492   9.649   6.528  1.00 10.00           H
ATOM     13  HD2 PRO A   4       5.306   8.082   8.255  1.00 10.00           H
ATOM     14  HD3 PRO A   4       4.496   9.094   7.222  1.00 10.00           H
"""

# Serine
pdb_str_16 = """
CRYST1   12.893   12.708   12.721  90.00  90.00  90.00 P 1
SCALE1      0.077561  0.000000  0.000000        0.00000
SCALE2      0.000000  0.078691  0.000000        0.00000
SCALE3      0.000000  0.000000  0.078610        0.00000
ATOM      1  N   SER A  73       5.347   5.750   5.075  1.00 10.00           N
ATOM      2  CA  SER A  73       5.894   5.594   6.421  1.00 10.00           C
ATOM      3  C   SER A  73       5.946   7.282   6.996  1.00 10.00           C
ATOM      4  O   SER A  73       5.026   7.627   7.452  1.00 10.00           O
ATOM      5  CB  SER A  73       7.501   5.019   6.235  1.00 10.00           C
ATOM      6  OG  SER A  73       8.063   4.929   7.844  1.00 10.00           O
ATOM      7  HA  SER A  73       5.300   5.295   7.021  1.00 10.00           H
ATOM      8  HB2 SER A  73       7.399   4.067   6.115  1.00 10.00           H
ATOM      9  HB3 SER A  73       8.072   5.776   5.734  1.00 10.00           H
ATOM     10  HG  SER A  73       8.499   4.647   7.768  1.00 10.00           H
"""

# Threonine
pdb_str_17 = """
CRYST1   11.909   12.199   14.459  90.00  90.00  90.00 P 1
SCALE1      0.083970  0.000000  0.000000        0.00000
SCALE2      0.000000  0.081974  0.000000        0.00000
SCALE3      0.000000  0.000000  0.069161        0.00000
ATOM      1  N   THR A  68       6.967   7.035   7.598  1.00 10.00           N
ATOM      2  CA  THR A  68       5.580   6.835   7.409  1.00 10.00           C
ATOM      3  C   THR A  68       4.981   6.180   8.674  1.00 10.00           C
ATOM      4  O   THR A  68       6.079   5.829   9.469  1.00 10.00           O
ATOM      5  CB  THR A  68       5.071   6.452   6.355  1.00 10.00           C
ATOM      6  OG1 THR A  68       5.819   5.087   6.249  1.00 10.00           O
ATOM      7  CG2 THR A  68       5.178   7.219   5.045  1.00 10.00           C
ATOM      8  HA  THR A  68       4.978   7.946   7.649  1.00 10.00           H
ATOM      9  HB  THR A  68       3.853   6.156   6.182  1.00 10.00           H
ATOM     10  HG1 THR A  68       6.272   5.006   6.437  1.00 10.00           H
ATOM     11 HG21 THR A  68       5.183   6.722   4.241  1.00 10.00           H
ATOM     12 HG22 THR A  68       4.943   8.119   5.225  1.00 10.00           H
ATOM     13 HG23 THR A  68       6.305   7.117   4.778  1.00 10.00           H
"""

# Triptophane
pdb_str_18 = """
CRYST1   12.502   12.982   18.312  90.00  90.00  90.00 P 1
SCALE1      0.079987  0.000000  0.000000        0.00000
SCALE2      0.000000  0.077030  0.000000        0.00000
SCALE3      0.000000  0.000000  0.054609        0.00000
ATOM      1  N   TRP A  49       7.019   6.036   6.855  1.00 10.00           N
ATOM      2  CA  TRP A  49       6.009   5.654   7.261  1.00 10.00           C
ATOM      3  C   TRP A  49       5.132   5.145   5.877  1.00 10.00           C
ATOM      4  O   TRP A  49       5.267   5.983   5.027  1.00 10.00           O
ATOM      5  CB  TRP A  49       5.058   6.563   8.316  1.00 10.00           C
ATOM      6  CG  TRP A  49       6.065   6.967   9.207  1.00 10.00           C
ATOM      7  CD1 TRP A  49       6.863   7.994   9.455  1.00 10.00           C
ATOM      8  CD2 TRP A  49       5.959   6.505  10.743  1.00 10.00           C
ATOM      9  NE1 TRP A  49       7.435   8.004  10.831  1.00 10.00           N
ATOM     10  CE2 TRP A  49       6.803   7.200  11.335  1.00 10.00           C
ATOM     11  CE3 TRP A  49       5.119   5.187  11.318  1.00 10.00           C
ATOM     12  CZ2 TRP A  49       7.129   6.620  12.976  1.00 10.00           C
ATOM     13  CZ3 TRP A  49       5.299   5.155  12.457  1.00 10.00           C
ATOM     14  CH2 TRP A  49       6.374   5.658  13.211  1.00 10.00           C
ATOM     15  HA  TRP A  49       6.139   4.842   7.608  1.00 10.00           H
ATOM     16  HB2 TRP A  49       5.070   7.319   7.481  1.00 10.00           H
ATOM     17  HB3 TRP A  49       4.305   6.336   8.251  1.00 10.00           H
ATOM     18  HD1 TRP A  49       7.226   8.582   8.796  1.00 10.00           H
ATOM     19  HE1 TRP A  49       8.120   8.555  10.929  1.00 10.00           H
ATOM     20  HE3 TRP A  49       4.311   5.074  10.564  1.00 10.00           H
ATOM     21  HZ2 TRP A  49       7.806   7.092  13.214  1.00 10.00           H
ATOM     22  HZ3 TRP A  49       4.914   4.129  12.933  1.00 10.00           H
ATOM     23  HH2 TRP A  49       6.353   5.513  14.064  1.00 10.00           H
"""

# Tyrosine
pdb_str_19 = """
CRYST1   17.955   13.272   13.095  90.00  90.00  90.00 P 1
SCALE1      0.055695  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075347  0.000000        0.00000
SCALE3      0.000000  0.000000  0.076365        0.00000
ATOM      1  N   TYR A 139      10.057   7.968   5.049  1.00 10.00           N
ATOM      2  CA  TYR A 139      10.657   7.531   6.379  1.00 10.00           C
ATOM      3  C   TYR A 139      12.203   7.725   6.416  1.00 10.00           C
ATOM      4  O   TYR A 139      12.999   8.272   7.373  1.00 10.00           O
ATOM      5  CB  TYR A 139      10.644   6.145   6.711  1.00 10.00           C
ATOM      6  CG  TYR A 139       9.159   5.899   6.690  1.00 10.00           C
ATOM      7  CD1 TYR A 139       8.503   5.230   5.513  1.00 10.00           C
ATOM      8  CD2 TYR A 139       8.317   6.046   8.121  1.00 10.00           C
ATOM      9  CE1 TYR A 139       6.876   4.938   5.643  1.00 10.00           C
ATOM     10  CE2 TYR A 139       7.209   5.706   8.077  1.00 10.00           C
ATOM     11  CZ  TYR A 139       6.420   5.365   6.855  1.00 10.00           C
ATOM     12  OH  TYR A 139       5.027   4.949   7.088  1.00 10.00           O
ATOM     13  HA  TYR A 139      10.303   8.174   6.785  1.00 10.00           H
ATOM     14  HB2 TYR A 139      10.989   5.647   5.882  1.00 10.00           H
ATOM     15  HB3 TYR A 139      10.828   5.883   7.586  1.00 10.00           H
ATOM     16  HD1 TYR A 139       8.618   5.356   4.741  1.00 10.00           H
ATOM     17  HD2 TYR A 139       8.841   6.121   8.546  1.00 10.00           H
ATOM     18  HE1 TYR A 139       6.432   5.030   4.892  1.00 10.00           H
ATOM     19  HE2 TYR A 139       6.780   5.619   9.066  1.00 10.00           H
ATOM     20  HH  TYR A 139       4.693   5.141   7.840  1.00 10.00           H
"""

# Valine
pdb_str_20 = """
CRYST1   12.396   13.122   13.130  90.00  90.00  90.00 P 1
SCALE1      0.080671  0.000000  0.000000        0.00000
SCALE2      0.000000  0.076208  0.000000        0.00000
SCALE3      0.000000  0.000000  0.076161        0.00000
ATOM      1  N   VAL B  78       4.820   5.520   7.634  1.00 10.00           N
ATOM      2  CA  VAL B  78       6.105   5.869   7.138  1.00 10.00           C
ATOM      3  C   VAL B  78       6.717   7.176   8.084  1.00 10.00           C
ATOM      4  O   VAL B  78       5.750   8.012   8.022  1.00 10.00           O
ATOM      5  CB  VAL B  78       6.029   6.374   5.701  1.00 10.00           C
ATOM      6  CG1 VAL B  78       5.673   4.961   4.874  1.00 10.00           C
ATOM      7  CG2 VAL B  78       7.349   6.810   5.233  1.00 10.00           C
ATOM      8  HA  VAL B  78       7.125   5.144   7.411  1.00 10.00           H
ATOM      9  HB  VAL B  78       5.404   6.962   5.498  1.00 10.00           H
ATOM     10 HG11 VAL B  78       5.583   5.210   4.024  1.00 10.00           H
ATOM     11 HG12 VAL B  78       6.062   4.120   5.126  1.00 10.00           H
ATOM     12 HG13 VAL B  78       4.907   4.924   5.317  1.00 10.00           H
ATOM     13 HG21 VAL B  78       7.299   7.038   4.409  1.00 10.00           H
ATOM     14 HG22 VAL B  78       7.542   7.634   5.823  1.00 10.00           H
ATOM     15 HG23 VAL B  78       7.889   5.963   5.346  1.00 10.00           H
"""

pdb_list = [pdb_str_01, pdb_str_02, pdb_str_03, pdb_str_04, pdb_str_05,
  pdb_str_06, pdb_str_07, pdb_str_08, pdb_str_09, pdb_str_10, pdb_str_11,
  pdb_str_12, pdb_str_13, pdb_str_14, pdb_str_15, pdb_str_16, pdb_str_17,
  pdb_str_18, pdb_str_19, pdb_str_20]
#
pdb_list_name = ['pdb_str_01', 'pdb_str_02', 'pdb_str_03', 'pdb_str_04', 'pdb_str_05',
  'pdb_str_06', 'pdb_str_07', 'pdb_str_08', 'pdb_str_09', 'pdb_str_10', 'pdb_str_11',
  'pdb_str_12', 'pdb_str_13', 'pdb_str_14', 'pdb_str_15', 'pdb_str_16', 'pdb_str_17',
  'pdb_str_18', 'pdb_str_19', 'pdb_str_20']

#pdb_list = [pdb_str_02]
#pdb_list_name = ['pdb_str_02']

def run():
  #for use_ideal_bonds_angles in [True, False]:
  for pdb_str, str_name in zip(pdb_list,pdb_list_name):
      #print 'pdb_string:', str_name, 'idealize =', idealize
    exercise(pdb_str=pdb_str, eps=1.e-4)


if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time:", round(time.time()-t0, 2), "seconds")
