from __future__ import division
import time

import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.model
from mmtbx import monomer_library
from cctbx import geometry_restraints
from mmtbx.hydrogens import riding

#from mmtbx import hydrogens

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
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = None,
    raw_records    = pdb_str,
    force_symmetry = True)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  geometry_restraints = processed_pdb_file.geometry_restraints_manager(
    show_energies = False)

# necessary for comparison
  xray_structure = processed_pdb_file.xray_structure()
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry_restraints,
    normalization = False)
  angle_proxies = restraints_manager.geometry.get_all_angle_proxies()

  riding_h_manager = riding.manager(
    pdb_hierarchy           = pdb_hierarchy,
    geometry_restraints = geometry_restraints)
  h_connectivity = riding_h_manager.h_connectivity

  bond_list = {}
  angle_list = {}
  third_nb_list = {}
  for ih in h_connectivity.keys():
    a0 = (h_connectivity[ih][0])
    bond_list[ih]=[a0.iseq, a0.dist_ideal]
    for atom in h_connectivity[ih][1]+h_connectivity[ih][2]:
      helper = tuple(sorted([ih, a0.iseq, atom.iseq]))
      angle_list[helper]=atom.angle_ideal
    if(len(h_connectivity[ih])==4):
      third_nb_list[ih]=[]
      for third in h_connectivity[ih][3]:
        third_nb_list[ih].append(third.iseq)

#-----------------------------------------------------------------------------
# This is useful to keep for debugging: human readable output of connectivity
#-----------------------------------------------------------------------------
#  for ih in connectivity.keys():
#    if(len(connectivity[ih])==3):
#      string = (" ".join([names[p.iseq] for p in connectivity[ih][2]]))
#    else:
#      string = 'n/a'
#    print  names[ih],': ', names[(connectivity[ih][0][0]).iseq], \
#      ',', (" ".join([names[p.iseq] for p in connectivity[ih][1]])), \
#      ',', string
#-----------------------------------------------------------------------------

# determine bonds from pdb_str
  model = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure     = xray_structure,
    pdb_hierarchy      = pdb_hierarchy)
  bond_ctrl = {}
  for i in model.xh_connectivity_table():
    bond_ctrl[i[1]]=[i[0],i[3]]

# List of angle restraints
  angles = [
    (4, 1, 12),
    (2, 1, 12),
    (0, 1, 12),
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
      angle_ctrl[tuple(sorted(list(ap.i_seqs)))]=ap.angle_ideal

# HH needs also third neighbors:
  third_nb_ctrl = {19: [8, 9]}

  assert (bond_list == bond_ctrl), '1-2 neighbors and distance_ideal are wrong'
  assert (angle_list == angle_ctrl), '1-3 neighbors and angle_ideal are wrong'
  assert (third_nb_list == third_nb_ctrl), '1-4 neighbors are wrong'

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print "OK. Time: %8.3f"%(time.time()-t0)
