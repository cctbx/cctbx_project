from cctbx.array_family import flex
import math, time, sys, os
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
import mmtbx.model
from libtbx import introspection
import libtbx.load_env
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cStringIO import StringIO
from mmtbx import model_statistics

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
                   relative_path="phenix_regression/pdb/enk.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
################
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  aal= processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  mol = mmtbx.model.manager(
             restraints_manager     = restraints_manager,
             xray_structure         = xray_structure,
             atom_attributes_list   = aal)
  mol.xray_structure.scattering_type_registry(table = "wk1995")
################



  mol_copy = mol.deep_copy()
  assert mol.number_of_ordered_solvent_molecules() == 9
  mol.write_pdb_file(out = open("test_model_out.pdb","w"))
  mol.remove_atoms(atom_type = "H")
  mol.write_pdb_file(out = open("test_model_out_nohydrogens.pdb","w"))
  mol.remove_solvent()
  assert mol.number_of_ordered_solvent_molecules() == 0
  mol.write_pdb_file(out = open("test_model_out_nosolvent.pdb","w"))
  mol.remove_atoms(atom_type = "O")
  mol.write_pdb_file(out = open("test_model_out_noO.pdb","w"))
  mol.show_geometry_statistics()
  mol.show_geometry_statistics()

################
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  aal= processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  mol_other = mmtbx.model.manager(
             restraints_manager     = restraints_manager,
             xray_structure         = xray_structure,
             atom_attributes_list   = aal)
  mol_other.xray_structure.scattering_type_registry(table = "wk1995")
################

  mol_other.show_geometry_statistics()
  print
  mol.show_geometry_statistics()

#####
  class iso: pass
  iso.sphere_radius = 5.0
  iso.distance_power = 1.57
  iso.average_power = 0.58
  iso.wilson_b_weight_auto = False
  iso.wilson_b_weight = None
  iso.plain_pairs_radius = 5.0
  iso.refine_ap_and_dp = False
  iso.use_u_local_only = False
#####

  mol.show_adp_statistics()
  print
  mol.show_adp_statistics()

  rm = mol.restraints_manager

  mol_copy.write_pdb_file(out = open("XXX.pdb","w"))
  mol_copy.remove_atoms(leave_only_labels = ["CA", "N", "O"])
  mol_copy.write_pdb_file(out = open("XXXr.pdb","w"))


def exercise_2():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
                   relative_path="phenix_regression/pdb/adp_out_stat.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  aal= processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  mol = mmtbx.model.manager(
             restraints_manager     = restraints_manager,
             xray_structure         = xray_structure,
             atom_attributes_list   = aal)
  mol.xray_structure.scattering_type_registry(table = "wk1995")

  out = StringIO()
  adp_stat = mol.show_adp_statistics(out = out)
  expected_result = \
  """|-ADP statistics-------------------------------------------------------|
| Atom    | Number of   | Isotropic or equivalent| Anisotropy lmin/max |
| type    |iso    aniso | min     max     mean   | min   max    mean   |
| - - - - |- - - - - - -| - - - - - - - - - - - -| - - - - - - - - - - |
| all     : 14     7      1.31    18.27   13.11    0.02  0.73   0.20   |
| all(noH): 9      5      1.45    18.27   15.14    0.03  0.73   0.18   |
| Sol.    : 1      1      1.45    2.00    1.73     0.73  0.73   0.73   |
| Mac.    : 8      4      15.00   18.27   17.38    0.03  0.05   0.04   |
| Hyd.    : 5      2      1.31    17.15   9.04     0.02  0.50   0.26   |
| - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  |
|    Distribution of isotropic (or equivalent) ADP for non-H atoms:    |
| Bin#      value range     #atoms | Bin#      value range     #atoms  |
|   0:     1.450 -   3.132:    2   |   5:     9.860 -  11.542:    0    |
|   1:     3.132 -   4.814:    0   |   6:    11.542 -  13.224:    0    |
|   2:     4.814 -   6.496:    0   |   7:    13.224 -  14.906:    0    |
|   3:     6.496 -   8.178:    0   |   8:    14.906 -  16.588:    1    |
|   4:     8.178 -   9.860:    0   |   9:    16.588 -  18.270:   11    |
|                            =>continue=>                              |
| - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  |
|                     Distribution of anisotropy:                      |
| Bin#      value range     #atoms | Bin#      value range     #atoms  |
|   0:     0.028 -   0.098:    4   |   5:     0.377 -   0.447:    0    |
|   1:     0.098 -   0.168:    0   |   6:     0.447 -   0.517:    0    |
|   2:     0.168 -   0.238:    0   |   7:     0.517 -   0.587:    0    |
|   3:     0.238 -   0.307:    0   |   8:     0.587 -   0.657:    0    |
|   4:     0.307 -   0.377:    0   |   9:     0.657 -   0.726:    1    |
|                            =>continue=>                              |
|----------------------------------------------------------------------|
"""
  assert out.getvalue() == expected_result

def run():
  exercise()
  exercise_2()

if (__name__ == "__main__"):
  run()
