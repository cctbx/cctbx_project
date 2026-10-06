from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from cctbx import crystal
from mmtbx.hydrogens import connectivity
from libtbx.utils import null_out

pdb_str = """\
CRYST1   17.955   13.272   13.095  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A 139      10.241   7.920   5.000  1.00 10.00           N
ATOM      2  CA  TYR A 139      10.853   7.555   6.271  1.00 10.00           C
ATOM      3  C   TYR A 139      12.362   7.771   6.227  1.00 10.00           C
ATOM      4  O   TYR A 139      12.955   8.272   7.181  1.00 10.00           O
ATOM      5  CB  TYR A 139      10.540   6.098   6.617  1.00 10.00           C
ATOM      6  CG  TYR A 139       9.063   5.805   6.749  1.00 10.00           C
ATOM      7  CD1 TYR A 139       8.316   5.391   5.654  1.00 10.00           C
ATOM      8  CD2 TYR A 139       8.414   5.943   7.969  1.00 10.00           C
ATOM      9  CE1 TYR A 139       6.966   5.122   5.770  1.00 10.00           C
ATOM     10  CE2 TYR A 139       7.064   5.676   8.095  1.00 10.00           C
ATOM     11  CZ  TYR A 139       6.345   5.266   6.993  1.00 10.00           C
ATOM     12  OH  TYR A 139       5.000   5.000   7.113  1.00 10.00           O
ATOM     13  HA  TYR A 139      10.480   8.127   6.960  1.00 10.00           H
ATOM     14  HB2 TYR A 139      10.915   5.524   5.931  1.00 10.00           H
ATOM     15  HB3 TYR A 139      10.982   5.870   7.450  1.00 10.00           H
ATOM     16  HD1 TYR A 139       8.732   5.293   4.828  1.00 10.00           H
ATOM     17  HD2 TYR A 139       8.896   6.220   8.714  1.00 10.00           H
ATOM     18  HE1 TYR A 139       6.479   4.845   5.028  1.00 10.00           H
ATOM     19  HE2 TYR A 139       6.643   5.772   8.919  1.00 10.00           H
ATOM     20  HH  TYR A 139       4.759   5.128   7.907  1.00 10.00           H
TER
END
"""

#----------------------------------------------------
# This test checks for residue Tyr (pdb_str above):
# - if all bonds involving H atoms are recognized
# - if all angles involving H atoms are recognized
# - if 3rd neighbors of HH atom are correctly found
#----------------------------------------------------

def exercise():
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  restraints_manager = model.get_restraints_manager()
  angle_proxies = restraints_manager.geometry.get_all_angle_proxies()

  connectivity_manager = connectivity.determine_connectivity(
    pdb_hierarchy       = model.get_hierarchy(),
    geometry_restraints = restraints_manager.geometry)
  h_connectivity = connectivity_manager.h_connectivity

# get bonds stored in connectivity
  bond_list = {}
  angle_list = {}
  for neighbors in h_connectivity:
    if (neighbors is None): continue
    ih = neighbors.ih
    a0 = neighbors.a0
    i_a0 = a0['iseq']
    a1 = neighbors.a1
    i_a1 = a1['iseq']
    bond_list[ih] = [i_a0, a0['dist_ideal']]
    selected_atoms = tuple(sorted([ih, i_a0, i_a1]))
    angle_list[selected_atoms] = a1['angle_ideal']
    if neighbors.a2:
      a2 = neighbors.a2
      selected_atoms2 = tuple(sorted([ih, i_a0, a2['iseq']]))
      angle_list[selected_atoms2] = a2['angle_ideal']
    if neighbors.a3:
      a3 = neighbors.a3
      selected_atoms3 = tuple(sorted([ih, i_a0, a3['iseq']]))
      angle_list[selected_atoms3] = a3['angle_ideal']
    if neighbors.h1:
      h1 = neighbors.h1
      selected_atoms4 = tuple(sorted([ih, i_a0, h1['iseq']]))
      angle_list[selected_atoms4] =h1['angle_ideal']
    if neighbors.b1:
      i_b1 = neighbors.b1['iseq']
      third_nb_dict = {ih: i_b1}

  bond_ctrl = {}
  for i in model.xh_connectivity_table():
    bond_ctrl[i[1]]=[i[0],i[3]]

# List of angle restraints
  angles = [
    (4, 1, 12),

    (0, 1, 12),
    (2, 1, 12),
    (13, 4, 14),
    (5, 4, 14),
    (5, 4, 13),
    (1, 4, 13),
    (1, 4, 14),
    (8, 6, 15),
    (5, 6, 15),
    (9, 7, 16),
    (5, 7, 16),
    (10, 8, 17),
    (6, 8, 17),
    (10, 11, 19),
    (7, 9, 18),
    (10, 9, 18)]

  angle_ctrl = {}
  for ap in angle_proxies:
    if(ap.i_seqs in angles):
      angle_ctrl[tuple(sorted(list(ap.i_seqs)))] = ap.angle_ideal

# HH needs also third neighbors:
  third_nb_ctrl = {19: 8}

  assert (bond_list == bond_ctrl), '1-2 neighbors and distance_ideal are wrong'
  assert (angle_list == angle_ctrl), '1-3 neighbors and angle_ideal are wrong'
  assert (third_nb_dict == third_nb_ctrl), '1-4 neighbors are wrong'

#----------------------------------------------------
# ARG 63 altloc B of 6b8f. HH22 is at fractional x = -0.001, right on the ASU
# boundary. phenix.fit_h processes the model under the PDB cell and the fmodel
# setup then puts the MTZ cell (7e-6 A smaller) on the same geometry manager
# without reprocessing; that moves HH22 across the boundary, so NH2-HH22 becomes
# a symmetry bond. This checks that the dihedral loop skips such an H instead
# of dereferencing None.
#----------------------------------------------------

pdb_str_arg63 = """\
CRYST1  180.040  180.040  180.040  90.00  90.00  90.00 F 4 3 2
ATOM      1  CA BARG A  63       5.193  32.311  31.420  0.47  6.93           C
ATOM      2  CB BARG A  63       3.750  31.902  31.178  0.47  8.28           C
ATOM      3  CG BARG A  63       3.238  30.969  32.226  0.47 10.10           C
ATOM      4  CD BARG A  63       1.820  30.536  31.904  0.47  9.82           C
ATOM      5  NE BARG A  63       1.703  29.805  30.654  0.47 10.85           N
ATOM      6  CZ BARG A  63       0.543  29.497  30.085  0.47 12.24           C
ATOM      7  NH1BARG A  63      -0.592  29.880  30.645  0.47 12.44           N
ATOM      8  NH2BARG A  63       0.549  28.828  28.932  0.47 14.99           N
ATOM      9  HA BARG A  63       5.246  32.760  32.278  0.47  8.32           H
ATOM     10  HB2BARG A  63       3.191  32.695  31.183  0.47  9.94           H
ATOM     11  HB3BARG A  63       3.686  31.455  30.320  0.47  9.94           H
ATOM     12  HG2BARG A  63       3.800  30.179  32.258  0.47 12.13           H
ATOM     13  HG3BARG A  63       3.234  31.418  33.086  0.47 12.13           H
ATOM     14  HD2BARG A  63       1.499  29.961  32.616  0.47 11.79           H
ATOM     15  HD3BARG A  63       1.259  31.325  31.841  0.47 11.79           H
ATOM     16  HE BARG A  63       2.419  29.473  30.313  0.47 13.03           H
ATOM     17 HH11BARG A  63      -0.580  30.315  31.387  0.47 14.94           H
ATOM     18 HH12BARG A  63      -1.343  29.683  30.274  0.47 14.94           H
ATOM     19 HH21BARG A  63       1.295  28.589  28.577  0.47 17.99           H
ATOM     20 HH22BARG A  63      -0.194  28.625  28.549  0.47 17.99           H
TER
END
"""

def exercise_symmetry_bond():
  pdb_inp = iotbx.pdb.input(lines=pdb_str_arg63.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  # Shrink the cell by 7e-6 A on the processed model, as the fmodel setup does.
  model.set_crystal_symmetry(crystal.symmetry(
    unit_cell          = (180.0399932861, 180.0399932861, 180.0399932861,
                          90, 90, 90),
    space_group_symbol = "F 4 3 2"))
  geometry = model.get_restraints_manager().geometry
  bond_proxies_simple, asu = geometry.get_all_bond_proxies(
    sites_cart = model.get_sites_cart())
  # Without a symmetry bond the test would pass for the wrong reason.
  assert (asu.size() > 0), 'no symmetry bond, test does not exercise the bug'

  connectivity_manager = connectivity.determine_connectivity(
    pdb_hierarchy       = model.get_hierarchy(),
    geometry_restraints = geometry)
  h_connectivity = connectivity_manager.h_connectivity

  # add_slipped must give the H of the symmetry bond its entry back
  hd_sel = model.get_hd_selection()
  without_entry = [model.get_hierarchy().atoms()[ih].name.strip()
    for ih in range(len(hd_sel))
      if hd_sel[ih] and h_connectivity[ih] is None]
  assert (not without_entry), 'H atoms without connectivity: %s' % without_entry

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  exercise_symmetry_bond()
  print("OK. Time: %8.3f"%(time.time()-t0))
