from cctbx.array_family import flex
import os
import mmtbx.model
import libtbx.load_env
from libtbx import easy_pickle
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cStringIO import StringIO
from mmtbx import utils
from libtbx.utils import format_cpu_times

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
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  mol.xray_structure.scattering_type_registry(table = "wk1995")
################



  mol_copy = mol.deep_copy()
  assert mol.number_of_ordered_solvent_molecules() == 9
  mol.write_pdb_file(out = open("test_model_out.pdb","w"))
  mol.write_pdb_file(out = open("test_model_out_nohydrogens.pdb","w"))
  mol = mol.remove_solvent()
  assert mol.number_of_ordered_solvent_molecules() == 0
  mol.write_pdb_file(out = open("test_model_out_nosolvent.pdb","w"))
  mol.write_pdb_file(out = open("test_model_out_noO.pdb","w"))
  mol.show_geometry_statistics(ignore_hd = True)
  mol.show_geometry_statistics(ignore_hd = True)

################
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  mol_other = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  mol_other.xray_structure.scattering_type_registry(table = "wk1995")
################

  mol_other.show_geometry_statistics(ignore_hd = True)
  print
  mol.show_geometry_statistics(ignore_hd = True)

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
  mol_copy.write_pdb_file(out = open("XXXr.pdb","w"))


def exercise_2():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
                   relative_path="phenix_regression/pdb/adp_out_stat.pdb", test=os.path.isfile)
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.nonbonded_weight = 16
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       params                    = params,
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
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
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
  # XXX phenix GUI support (see wxtbx.adp_statistics)
  stats = mol.adp_statistics()
  tables = stats.format_tables()
  stats_pkl = easy_pickle.dumps(stats)
  stats2 = easy_pickle.loads(stats_pkl)
  tables2 = stats.format_tables()
  assert (tables2 == tables)
  t1,t2,t3 = tables2
  assert (len(t1) == 5) and (len(t2) == 10) and (len(t3) == 10)
  assert (t1[0][1:3] == ["14","7"]) and (t1[-1][1:3] == ["5","2"])
  assert (t2[1][1] == "3.132 - 4.814")
  p1, p2 = stats.format_plots()
  y1, yrange1 = p1
  assert (len(y1) == 10)

def exercise_3():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/arg_h_hohh1_sh.pdb",
    test=os.path.isfile)
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.nonbonded_weight = 16
  processed_pdb_files_srv = utils.process_pdb_file_srv(
    pdb_interpretation_params=params)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  out = StringIO()
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    log = out)
  #
  mol.idealize_h()
  assert out.getvalue().splitlines()[0] == \
  "X-H deviation from ideal before regularization (bond): mean= 0.154 max= 0.496"
  assert out.getvalue().splitlines()[1] == \
  "X-H deviation from ideal after  regularization (bond): mean= 0.000 max= 0.000"

def exercise_4():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lysozyme_noH.pdb",
    test=os.path.isfile)
  processed_pdb_files_srv = utils.process_pdb_file_srv()
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  out = StringIO()
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    log = out)
  #
  result = mol.isolated_atoms_selection()
  solvent_sel = mol.solvent_selection()
  assert solvent_sel.count(True)+1 == result.count(True)

def run():
  exercise()
  exercise_2()
  exercise_3()
  exercise_4()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
