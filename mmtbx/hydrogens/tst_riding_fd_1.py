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
# Finite difference test for modified gradients of riding H
# Every type of H geometry is tested:
# planar H
# tetragonal H
# two tetragonal H
# X-H2 group (ARG etc)
# propeller group
# H with rotational degree of freedom (SER, etc)
#-----------------------------------------------------------------------------

def exercise(pdb_str, eps, use_ideal_bonds_angles):
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

# alg2a shaked
pdb_str_00 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      6  CG  TYR A   7       8.528   6.445   8.572  1.00 15.00           C
ATOM      8  CD2 TYR A   7       8.348   7.014   9.545  1.00 15.00           C
ATOM     10  CE2 TYR A   7       7.068   7.273  10.230  1.00 15.00           C
ATOM     17  HD2 TYR A   7       9.064   7.637  10.081  1.00 15.00           H
TER
"""
# alg2a ideal geometry
pdb_str_01 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      1  CG  TYR A   7       8.544   6.391   8.506  1.00 15.00           C
ATOM      2  CD2 TYR A   7       8.334   7.144   9.654  1.00 15.00           C
ATOM      3  CE2 TYR A   7       7.078   7.251  10.220  1.00 15.00           C
ATOM      4  HD2 TYR A   7       9.052   7.584  10.049  1.00 15.00           H
TER
"""
# alg2b shaked
pdb_str_02 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      2  CA  TYR A   7       9.887   7.460   6.606  1.00 15.00           C
ATOM      5  CB  TYR A   7       9.830   6.306   7.775  1.00 15.00           C
ATOM      6  CG  TYR A   7       8.528   6.445   8.572  1.00 15.00           C
ATOM     14  HB2 TYR A   7       9.943   5.616   7.588  1.00 15.00           H
END
"""
# alg2b ideal geometry
pdb_str_03 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      2  CA  TYR A   7       9.962   6.614   6.452  1.00 15.00           C
ATOM      5  CB  TYR A   7       9.920   6.991   7.935  1.00 15.00           C
ATOM      6  CG  TYR A   7       8.575   6.781   8.602  1.00 15.00           C
ATOM     14  HB2 TYR A   7      10.571   6.452   8.411  1.00 15.00           H
"""

# alg3 shaked
pdb_str_04 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      1  N   TYR A   7       9.035   7.089   5.637  1.00 15.00           N
ATOM      2  CA  TYR A   7       9.887   7.460   6.606  1.00 15.00           C
ATOM      3  C   TYR A   7      11.285   7.169   6.250  1.00 15.00           C
ATOM      5  CB  TYR A   7       9.830   6.306   7.775  1.00 15.00           C
ATOM     13  HA  TYR A   7       9.831   8.340   7.271  1.00 15.00           H
"""

# alg3 ideal geometry
pdb_str_05 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      1  N   TYR A   7       8.981   7.170   5.658  1.00 15.00           N
ATOM      2  CA  TYR A   7       9.949   7.360   6.732  1.00 15.00           C
ATOM      3  C   TYR A   7      11.374   7.278   6.195  1.00 15.00           C
ATOM      4  CB  TYR A   7       9.740   6.317   7.832  1.00 15.00           C
ATOM      5  HA  TYR A   7       9.823   8.239   7.122  1.00 15.00           H
"""

# alg1a shaked
pdb_str_06 = """
CRYST1   14.446   16.451   11.913  90.00  90.00  90.00 P 1
SCALE1      0.069223  0.000000  0.000000        0.00000
SCALE2      0.000000  0.060787  0.000000        0.00000
SCALE3      0.000000  0.000000  0.083942        0.00000
ATOM      8  NE  ARG A  50       5.682   6.869   5.406  1.00 10.00           N
ATOM      9  CZ  ARG A  50       6.378   6.261   5.888  1.00 10.00           C
ATOM     10  NH1 ARG A  50       7.214   6.770   6.820  1.00 10.00           N
ATOM     11  NH2 ARG A  50       6.924   5.048   5.452  1.00 10.00           N
ATOM     21 HH12 ARG A  50       7.949   6.540   7.393  1.00 10.00           H
"""

# alg1a idealized
pdb_str_07 = """
CRYST1   14.446   16.451   11.913  90.00  90.00  90.00 P 1
SCALE1      0.069223  0.000000  0.000000        0.00000
SCALE2      0.000000  0.060787  0.000000        0.00000
SCALE3      0.000000  0.000000  0.083942        0.00000
ATOM      1  NE  ARG A  50       5.543   6.877   5.312  1.00 10.00           N
ATOM      2  CZ  ARG A  50       6.559   6.257   5.905  1.00 10.00           C
ATOM      3  NH1 ARG A  50       7.213   6.851   6.894  1.00 10.00           N
ATOM      4  NH2 ARG A  50       6.920   5.044   5.510  1.00 10.00           N
ATOM      6 HH12 ARG A  50       7.870   6.449   7.277  1.00 10.00           H
"""

#alg1b_flat shaked
pdb_str_08 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      9  CE1 TYR A   7       6.338   5.396   8.586  1.00 15.00           C
ATOM     10  CE2 TYR A   7       7.068   7.273  10.230  1.00 15.00           C
ATOM     11  CZ  TYR A   7       6.090   6.501   9.765  1.00 15.00           C
ATOM     12  OH  TYR A   7       4.857   6.463  10.553  1.00 15.00           O
ATOM     20  HH  TYR A   7       5.208   7.016  11.256  1.00 15.00           H
"""

#alg1b_flat idealized
pdb_str_09 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      1  CE1 TYR A   7       6.239   5.565   8.724  1.00 15.00           C
ATOM      2  CE2 TYR A   7       7.170   7.259  10.131  1.00 15.00           C
ATOM      3  CZ  TYR A   7       6.130   6.416   9.803  1.00 15.00           C
ATOM      4  OH  TYR A   7       4.978   6.424  10.555  1.00 15.00           O
ATOM      5  HH  TYR A   7       5.045   6.985  11.177  1.00 15.00           H
"""

#alg1b shaked (HH out of plane)
pdb_str_10 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      9  CE1 TYR A   7       6.231   6.395   8.363  1.00 15.00           C
ATOM     10  CE2 TYR A   7       7.119   6.748  10.501  1.00 15.00           C
ATOM     11  CZ  TYR A   7       5.926   6.249   9.983  1.00 15.00           C
ATOM     12  OH  TYR A   7       4.831   6.083  10.481  1.00 15.00           O
ATOM     20  HH  TYR A   7       4.708   6.641  11.301  1.00 15.00           H
"""

# alg1b idealized (HH out of plane)
pdb_str_11 = """
CRYST1   16.660   13.261   16.215  90.00  90.00  90.00 P 1
SCALE1      0.060024  0.000000  0.000000        0.00000
SCALE2      0.000000  0.075409  0.000000        0.00000
SCALE3      0.000000  0.000000  0.061671        0.00000
ATOM      9  CE1 TYR A   7       6.193   6.410   8.461  1.00 15.00           C
ATOM     10  CE2 TYR A   7       7.242   6.572  10.604  1.00 15.00           C
ATOM     11  CZ  TYR A   7       6.111   6.396   9.837  1.00 15.00           C
ATOM     12  OH  TYR A   7       4.893   6.206  10.448  1.00 15.00           O
ATOM     20  HH  TYR A   7       4.878   6.603  11.188  1.00 15.00           H
"""

#alg1a shaked with both H atoms
pdb_str_06_bla = """
CRYST1   14.446   16.451   11.913  90.00  90.00  90.00 P 1
SCALE1      0.069223  0.000000  0.000000        0.00000
SCALE2      0.000000  0.060787  0.000000        0.00000
SCALE3      0.000000  0.000000  0.083942        0.00000
ATOM      8  NE  ARG A  50       5.682   6.869   5.406  1.00 10.00           N
ATOM      9  CZ  ARG A  50       6.378   6.261   5.888  1.00 10.00           C
ATOM     10  NH1 ARG A  50       7.214   6.770   6.820  1.00 10.00           N
ATOM     11  NH2 ARG A  50       6.924   5.048   5.452  1.00 10.00           N
ATOM     20 HH11 ARG A  50       6.939   7.628   7.092  1.00 10.00           H
ATOM     21 HH12 ARG A  50       7.949   6.540   7.393  1.00 10.00           H
"""

pdb_list = [pdb_str_00, pdb_str_01, pdb_str_02, pdb_str_03,
  pdb_str_04, pdb_str_05, pdb_str_06, pdb_str_07, pdb_str_08, pdb_str_09,
  pdb_str_10, pdb_str_11]

pdb_list_name = ['pdb_str_00', 'pdb_str_01', 'pdb_str_02', 'pdb_str_03',
  'pdb_str_04', 'pdb_str_05', 'pdb_str_06', 'pdb_str_07', 'pdb_str_08',
  'pdb_str_09', 'pdb_str_10', 'pdb_str_11']

def run():
  for use_ideal_bonds_angles in [True, False]:
    for pdb_str, str_name in zip(pdb_list,pdb_list_name):
      #print 'pdb_string:', str_name, 'idealize =', idealize
      exercise(pdb_str=pdb_str, eps=1.e-4,
        use_ideal_bonds_angles=use_ideal_bonds_angles)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time:", round(time.time()-t0, 2), "seconds")
