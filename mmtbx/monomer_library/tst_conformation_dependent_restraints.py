from mmtbx.monomer_library import conformation_dependent_restraints
import iotbx.pdb
from cctbx.array_family import flex
from cctbx import geometry_restraints
import sys

# For tests, add assertions. Look at some tests -- tst_geometry_restraints.py

def read_pdb():

    pdbstring = """\
ATOM      0  CA  GLY A   3       5.804  -2.100   7.324  1.00  1.36           C
ATOM      1  C   GLY A   3       4.651  -1.149   7.578  1.00  1.01           C
ATOM      2  O   GLY A   3       3.598  -1.553   8.071  1.00  1.38           O
ATOM      3  N   GLY A   3       6.706  -1.622   6.294  1.00  1.11           N
ATOM      4  CA  PHE A   4       3.819   1.134   7.419  1.00  0.89           C
ATOM      5  CB  PHE A   4       4.397   2.380   8.094  1.00  1.13           C
ATOM      6  C   PHE A   4       3.185   1.509   6.084  1.00  0.94           C
ATOM      7  N   PHE A   4       4.852   0.121   7.242  1.00  0.88           N
ATOM      8  O   PHE A   4       2.361   2.421   6.010  1.00  1.47           O
ATOM      9  CA  LEU A   5       3.055   1.059   3.693  1.00  0.87           C
ATOM     10  CB  LEU A   5       3.965   0.435   2.634  1.00  1.13           C
ATOM     11  C   LEU A   5       1.634   0.527   3.541  1.00  0.87           C
ATOM     12  N   LEU A   5       3.576   0.800   5.030  1.00  0.92           N
ATOM     13  O   LEU A   5       1.246  -0.440   4.196  1.00  1.23           O
"""

    pdb_inp = iotbx.pdb.input(
        lines=flex.split_lines(pdbstring), source_info=None)

    sites_cart = pdb_inp.atoms().extract_xyz()

# TRANS    phi      1 C      2 N      2 CA     2 C        60.00  20.0 3
# TRANS    psi      1 N      1 CA     1 C      2 N       160.00  30.0 2

    dihedral_proxies = geometry_restraints.shared_dihedral_proxy()

    # residue 1
    psi = geometry_restraints.dihedral_proxy(
        i_seqs=[3, 0, 1, 7],
        angle_ideal=160.0,
        weight=1/30.0**2,
        periodicity=3
        )
    dihedral_proxies.append(psi)

    # residue 2
    phi = geometry_restraints.dihedral_proxy(
        i_seqs=[1, 7, 4, 6],
        angle_ideal=60.0,
        weight=1/20.0**2,
        periodicity=3
        )
    dihedral_proxies.append(phi)

    psi = geometry_restraints.dihedral_proxy(
        i_seqs=[7, 4, 6, 8],
        angle_ideal=160.0,
        weight=1/30.0**2,
        periodicity=3
        )
    dihedral_proxies.append(psi)

    # residue 3
    phi = geometry_restraints.dihedral_proxy(
        i_seqs=[6, 12, 9, 11],
        angle_ideal=60.0,
        weight=1/20.0**2,
        periodicity=3
        )
    dihedral_proxies.append(phi)

    angle_proxies = geometry_restraints.shared_angle_proxy()

    ## Residue 1
    # a3
    a = geometry_restraints.angle_proxy(
        i_seqs=[3, 0, 1],
        angle_ideal=0,
        weight=1
        )
    angle_proxies.append(a)

    # a7
    a = geometry_restraints.angle_proxy(
        i_seqs=[2, 1, 7],
        angle_ideal=0,
        weight=1
        )
    angle_proxies.append(a)

    ## Residue 2
    # a1
    a = geometry_restraints.angle_proxy(
        i_seqs=[1, 7, 4],
        angle_ideal=0,
        weight=1
        )
    angle_proxies.append(a)

    # a3
    a = geometry_restraints.angle_proxy(
        i_seqs=[7, 4, 6],
        angle_ideal=0,
        weight=1
        )
    angle_proxies.append(a)

    # a7
    a = geometry_restraints.angle_proxy(
        i_seqs=[8, 6, 12],
        angle_ideal=0,
        weight=1
        )
    angle_proxies.append(a)

    ## Residue 3
    # a1
    a = geometry_restraints.angle_proxy(
        i_seqs=[6, 12, 9],
        angle_ideal=0,
        weight=1
        )
    angle_proxies.append(a)

    # a3
    a = geometry_restraints.angle_proxy(
        i_seqs=[12, 9, 11],
        angle_ideal=0,
        weight=1
        )
    angle_proxies.append(a)

    # compute dihedral
    #dihedral = geometry_restraints.dihedral(
    #    sites_cart=sites_cart,
    #    proxy=dihedral_proxies[0])

    # Shows real dihedral value
    #print dihedral.angle_model, dihedral.delta
    cfd_list = []

    cfd = conformation_dependent_restraints.conformation_dependent_restraints(
        residue_name='GLY',
        next_residue_name='PHE',
        conformation_proxies=None,
        i_phi_proxy=None, # index into dihedral_proxies
        i_psi_proxy=0,
        i_dynamic_angles=[None, None, 0, None, None, None, 1], # indexes into angles in angle_proxies
        i_dynamic_dihedrals=None
        )
    cfd_list.append(cfd)

    cfd = conformation_dependent_restraints.conformation_dependent_restraints(
        residue_name='PHE',
        next_residue_name='LEU',
        conformation_proxies=None,
        i_phi_proxy=1, # index into dihedral_proxies
        i_psi_proxy=2,
        i_dynamic_angles=[2, None, 3, None, None, None, 4], # indexes into angles in angle_proxies
        i_dynamic_dihedrals=None
        )
    cfd_list.append(cfd)

    cfd = conformation_dependent_restraints.conformation_dependent_restraints(
        residue_name='LEU',
        next_residue_name=None,
        conformation_proxies=None,
        i_phi_proxy=3, # index into dihedral_proxies
        i_psi_proxy=None,
        i_dynamic_angles=[5, None, 6, None, None, None, None], # indexes into angles in angle_proxies
        i_dynamic_dihedrals=None
        )
    cfd_list.append(cfd)

    for x in range(1, 4):
        print
        print 'Starting cycle', x
        print
        for cfd in cfd_list:
            cfd.update_restraints(sites_cart, dihedral_proxies, angle_proxies)

def pdb_interpretation_run(args):
  from mmtbx.monomer_library.pdb_interpretation import run
  all_processed_pdb_files = run(args=args, return_all_processed_pdb_files=True)
  for processed_pdb_file in all_processed_pdb_files:
    all_proxies = processed_pdb_file.all_chain_proxies
    print "update_restraints:", \
      len(all_proxies.conformation_dependent_restraints_list)
    for x in all_proxies.conformation_dependent_restraints_list:
      x.update_restraints(sites_cart=all_proxies.sites_cart,
        dihedral_proxies=all_proxies.geometry_proxy_registries.dihedral.proxies,
        angle_proxies=all_proxies.geometry_proxy_registries.angle.proxies)

def exercise(args):
  if (not conformation_dependent_restraints.is_available):
    print "Skipping mmtbx.monomer_library.tst_conformation_dependent_restraints"
    return
  read_pdb()
  pdb_interpretation_run(args=args)
  print "OK"

if (__name__ == "__main__"):
  exercise(args=sys.argv[1:])
