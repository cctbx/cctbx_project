from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
import mmtbx.refinement.geometry_minimization
import scitbx.lbfgs
import cctbx.geometry_restraints
from mmtbx.hydrogens import riding
from libtbx.utils import null_out
from six.moves import zip

#-----------------------------------------------------------------------------
# This test checks the parameterization of hydrogen atoms for all H geometries
# for each fragment, the coordinates are minimized before computation of
# new H atom positions.
#-----------------------------------------------------------------------------


def exercise(pdb_str, use_ideal_bonds_angles):
# --------------------------------------------------------------
#          code to switch off CDL
# --------------------------------------------------------------
  #params_line = grand_master_phil_str
  #params = iotbx.phil.parse(
  #    input_string=params_line, process_includes=True).extract()
  #params.pdb_interpretation.restraints_library.cdl=False
# ---------------------------------------------------------------

  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  pdb_hierarchy = model.get_hierarchy()
  geometry_restraints = model.get_restraints_manager().geometry
  xray_structure = model.get_xray_structure()

  sites_cart = model.get_sites_cart()

  grf = cctbx.geometry_restraints.flags.flags(default=True)
  minimized = mmtbx.refinement.geometry_minimization.lbfgs(
      sites_cart                         = sites_cart,
      correct_special_position_tolerance = 1.0,
      geometry_restraints_manager        = geometry_restraints,
      geometry_restraints_flags          = grf,
      lbfgs_termination_params           = scitbx.lbfgs.termination_parameters(
        max_iterations=500))
  xray_structure.set_sites_cart(sites_cart)
  pdb_hierarchy.adopt_xray_structure(xray_structure)
  atoms = pdb_hierarchy.atoms()
  sites_cart = xray_structure.sites_cart()

  riding_h_manager = riding.manager(
    pdb_hierarchy       = pdb_hierarchy,
    geometry_restraints = geometry_restraints,
    use_ideal_bonds_angles = use_ideal_bonds_angles)

  h_parameterization = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart = sites_cart,
    threshold  = 0.05)
  h_distances   = diagnostics.h_distances

  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_parameterization) - h_parameterization.count(None)

  assert (number_h_para == number_h), 'Not all H atoms are parameterized'

  for ih in h_distances:
    labels = atoms[ih].fetch_labels()
    if use_ideal_bonds_angles:
      assert (h_distances[ih] < 0.03), \
        'distance too large: %s  atom: %s (%s) residue: %s ' \
        % (h_parameterization[ih].htype, atoms[ih].name, ih, labels.resseq.strip())
    else:
      assert (h_distances[ih] < 1e-7), \
        'distance too large: %s  atom: %s (%s) residue: %s  distance %s' \
        % (h_parameterization[ih].htype, atoms[ih].name, ih, labels.resseq.strip(), h_distances[ih])

# alg2a - flat two neighbors
pdb_str_00 = """
CRYST1   27.805   30.931   25.453  90.00  90.00  90.00 P 1
SCALE1      0.035965  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032330  0.000000        0.00000
SCALE3      0.000000  0.000000  0.039288        0.00000
ATOM      6  C   ASN A   1      15.183  17.286  10.105  1.00  9.67           C
ANISOU    6  C   ASN A   1     1289   1143   1241   -122   -106    160       C
ATOM     15  N   SER A   2      15.981  16.452  10.758  1.00  9.28           N
ANISOU   15  N   SER A   2     1274    902   1350    -21    -87    -74       N
ATOM     17  CA  SER A   2      15.404  15.325  11.479  1.00  9.47           C
ANISOU   17  CA  SER A   2     1420    764   1413     89    320   -137       C
ATOM     22  H   SER A   2      16.836  16.507  10.826  1.00 11.14           H
"""

# 2tetra
pdb_str_01 = """
CRYST1   27.805   30.931   25.453  90.00  90.00  90.00 P 1
SCALE1      0.035965  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032330  0.000000        0.00000
SCALE3      0.000000  0.000000  0.039288        0.00000
ATOM      1  CB  SER A   3      10.251  13.658  12.560  1.00  6.41           C
ANISOU    1  CB  SER A   3      638    591   1207   -162    148    106       C
ATOM      2  OG  SER A   3       9.057  14.277  12.222  1.00  7.60           O
ANISOU    2  OG  SER A   3      680    955   1251   -244   -164    339       O
ATOM      3  HB2 SER A   3      10.809  13.573  11.771  1.00  6.41           H
ATOM      4  HB3 SER A   3      10.073  12.759  12.879  1.00  6.41           H
ATOM     12  CA  SER A   3      11.103  14.453  13.643  1.00  5.95           C
ANISOU   12  CA  SER A   3      582    671   1007   -168    135    158       C
"""

# 3tetra
pdb_str_02 = """
CRYST1   27.805   30.931   25.453  90.00  90.00  90.00 P 1
SCALE1      0.035965  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032330  0.000000        0.00000
SCALE3      0.000000  0.000000  0.039288        0.00000
ATOM     16  CA  SER A   4       7.228  18.726  13.388  1.00  5.95           C
ANISOU   16  CA  SER A   4      582    671   1007   -168    135    158       C
ATOM     18  HA  SER A   4       6.548  18.900  14.058  1.00  5.95           H
ATOM     20  N   SER A   4       7.555  19.994  12.853  1.00  5.10           N
ANISOU   20  N   SER A   4      500    632    808   -107     58    104       N
ATOM     21  C   SER A   4       8.298  17.964  14.019  1.00  5.79           C
ANISOU   21  C   SER A   4      661    588    952   -116    206    135       C
ATOM     23  CB  SER A   4       6.376  17.931  12.305  1.00  6.41           C
ANISOU   23  CB  SER A   4      638    591   1207   -162    148    106       C
"""

# alg1a
pdb_str_03 = """
CRYST1   27.805   30.931   25.453  90.00  90.00  90.00 P 1
SCALE1      0.035965  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032330  0.000000        0.00000
SCALE3      0.000000  0.000000  0.039288        0.00000
ATOM      8  NE  ARG A   5       6.123   6.987   5.466  1.00 10.00           N
ATOM      9  CZ  ARG A   5       6.527   5.903   6.122  1.00 10.00           C
ATOM     10  NH1 ARG A   5       6.329   5.797   7.429  1.00 10.00           N
ATOM     20 HH11 ARG A   5       5.938   6.429   7.862  1.00 10.00           H
ATOM     21 HH12 ARG A   5       6.593   5.093   7.846  1.00 10.00           H
"""

# alg1b
pdb_str_04 = """
CRYST1   27.805   30.931   25.453  90.00  90.00  90.00 P 1
SCALE1      0.035965  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032330  0.000000        0.00000
SCALE3      0.000000  0.000000  0.039288        0.00000
ATOM     32  OG  SER A   6       8.321  21.302  15.330  1.00  7.60           O
ANISOU   32  OG  SER A   6      680    955   1251   -244   -164    339       O
ATOM     33  HG  SER A   6       8.020  20.964  14.622  1.00  7.60           H
ATOM     35  CA  SER A   6      10.267  21.478  16.751  1.00  5.95           C
ANISOU   35  CA  SER A   6      582    671   1007   -168    135    158       C
ATOM     38  CB  SER A   6       9.515  20.683  15.668  1.00  6.41           C
ANISOU   38  CB  SER A   6      638    591   1207   -162    148    106       C
"""

# propeller
pdb_str_05 = """
CRYST1   27.805   30.931   25.453  90.00  90.00  90.00 P 1
SCALE1      0.035965  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032330  0.000000        0.00000
SCALE3      0.000000  0.000000  0.039288        0.00000
ATOM     83  CA AVAL A   7      17.913   8.982  14.720  0.44 10.37           C
ANISOU   83  CA AVAL A   7     1783    902   1255     67   -148    -14       C
ATOM     84  CB AVAL A   7      16.666   8.901  13.831  0.44 11.26           C
ANISOU   84  CB AVAL A   7     1571    963   1743     96   -176     62       C
ATOM     85  CG1AVAL A   7      16.663  10.050  12.901  0.44 11.06           C
ANISOU   85  CG1AVAL A   7     1704   1072   1425    125   -105     37       C
ATOM     88 HG11AVAL A   7      16.650  10.872  13.416  0.44 13.27           H
ATOM     89 HG12AVAL A   7      15.875   9.999  12.339  0.44 13.27           H
ATOM     90 HG13AVAL A   7      17.463  10.015  12.354  0.44 13.27           H
"""


#type_list_known = ['flat_2neigbs', '2tetra', '2tetra', '3neigbs', 'alg1a',
#  'alg1a', 'alg1b', 'prop', 'prop', 'prop']

pdb_list = [pdb_str_00, pdb_str_01, pdb_str_02, pdb_str_03,
  pdb_str_04, pdb_str_05]

pdb_list_name = ['pdb_str_00', 'pdb_str_01', 'pdb_str_02', 'pdb_str_03',
  'pdb_str_04', 'pdb_str_05']

def run():
  for use_ideal_bonds_angles in [True, False]:
    print('use_ideal_bonds_angles =', use_ideal_bonds_angles)
    for pdb_str, str_name in zip(pdb_list,pdb_list_name):
      exercise(pdb_str=pdb_str, use_ideal_bonds_angles=use_ideal_bonds_angles)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
