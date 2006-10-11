from cctbx.array_family import flex
import math, time, sys, os
from scitbx.python_utils.misc import adopt_init_args
from libtbx.test_utils import approx_equal
import mmtbx.model
from libtbx import introspection
import libtbx.load_env
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation



def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
                   relative_path="regression/pdb/enk.pdb", test=os.path.isfile)
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
  restraints_manager_ini = mmtbx.restraints.manager(
                                  geometry      = geometry.select(selection),
                                  normalization = False)
  aal= processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  mol = mmtbx.model.manager(
             restraints_manager     = restraints_manager,
             restraints_manager_ini = restraints_manager_ini,
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
  mol.geometry_statistics(show = True)
  mol.geometry_statistics(other = mol, show = True)

################
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  restraints_manager_ini = mmtbx.restraints.manager(
                                  geometry      = geometry.select(selection),
                                  normalization = False)
  aal= processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  mol_other = mmtbx.model.manager(
             restraints_manager     = restraints_manager,
             restraints_manager_ini = restraints_manager_ini,
             xray_structure         = xray_structure,
             atom_attributes_list   = aal)
  mol_other.xray_structure.scattering_type_registry(table = "wk1995")
################

  mol_other.geometry_statistics(show = True)
  print
  mol.geometry_statistics(other = mol_other, show = True)

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

  mol.adp_statistics(iso_restraints = iso, show = True, other = None, wilson_b = None)
  print
  mol.adp_statistics(iso_restraints = iso, show = True, other = mol, wilson_b = None)

  rm = mol.restraints_manager

  mol_copy.write_pdb_file(out = open("XXX.pdb","w"))
  mol_copy.remove_atoms(leave_only_labels = ["CA", "N", "O"])
  mol_copy.write_pdb_file(out = open("XXXr.pdb","w"))


def run():
  exercise()

if (__name__ == "__main__"):
  run()
