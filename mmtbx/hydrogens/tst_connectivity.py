from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
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

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK. Time: %8.3f"%(time.time()-t0))
