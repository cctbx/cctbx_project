 # -*- coding: utf-8; py-indent-offset: 2 -*-

from __future__ import absolute_import, division, print_function
from mmtbx.ions.environment import ChemicalEnvironment
import mmtbx.ions.identify
from mmtbx import ions
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
from mmtbx.ions.environment import chem_carboxy, chem_amide, chem_backbone, \
     chem_water, chem_phosphate, \
     chem_nitrogen_primary, chem_nitrogen_secondary, \
     chem_chloride, chem_oxygen, chem_nitrogen, chem_sulfur
import libtbx.load_env
from collections import OrderedDict, Counter
import os
import sys
from six.moves import zip
from six.moves import range


def exercise():
  if not libtbx.env.has_module("phenix_regression"):
    print("Skipping {}".format(os.path.split(__file__)[1]))
    return

  models = OrderedDict([
    ("2qng", [
      Counter({chem_oxygen: 7, chem_carboxy: 2, chem_water: 2,
               chem_backbone: 3}),
      Counter({chem_oxygen: 6, chem_carboxy: 3, chem_water: 1,
               chem_backbone: 2}),
      ]),
    ("3rva", [
      Counter({chem_oxygen: 6, chem_carboxy: 4, chem_water: 2}),
      Counter({chem_nitrogen: 1, chem_oxygen: 4, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 1}),
      Counter({chem_nitrogen: 4, chem_nitrogen_primary: 1,
               chem_nitrogen_secondary: 3, chem_backbone: 3}),
      ]),
    ("1mjh", [
      Counter({chem_oxygen: 6, chem_water: 3, chem_phosphate: 3}),
      Counter({chem_oxygen: 6, chem_water: 3, chem_phosphate: 3}),
      ]),
    ("4e1h", [
      Counter({chem_oxygen: 6, chem_carboxy: 4}),
      Counter({chem_oxygen: 6, chem_carboxy: 3}),
      Counter({chem_oxygen: 6, chem_carboxy: 3}),
      ]),
    ("2xuz", [
      Counter({chem_oxygen: 6}),
      ]),
    ("3zli", [
      Counter({chem_nitrogen: 2, chem_oxygen: 4, chem_nitrogen_secondary: 2,
               chem_carboxy: 1, chem_water: 1}),
      Counter({chem_sulfur: 4}),
      Counter({chem_nitrogen: 2, chem_oxygen: 4, chem_nitrogen_secondary: 2,
               chem_carboxy: 1, chem_water: 1}),
      Counter({chem_sulfur: 4}),
      ]),
    ("3e0f", [
      Counter({chem_nitrogen: 2, chem_oxygen: 4, chem_nitrogen_secondary: 2,
               chem_carboxy: 2, chem_phosphate: 2}),
      Counter({chem_nitrogen: 2, chem_oxygen: 2, chem_nitrogen_secondary: 2,
               chem_carboxy: 1, chem_phosphate: 1}),
      Counter({chem_nitrogen: 2, chem_oxygen: 3, chem_nitrogen_secondary: 2,
               chem_carboxy: 2, chem_phosphate: 1}),
      ]),
    ("3dkq", [
      Counter({chem_nitrogen: 4, chem_oxygen: 1, chem_nitrogen_secondary: 4,
               chem_carboxy: 1}),
      Counter({chem_nitrogen: 2, chem_oxygen: 1, chem_nitrogen_secondary: 2,
               chem_carboxy: 1}),
      Counter({chem_nitrogen: 4, chem_oxygen: 1, chem_nitrogen_secondary: 4,
               chem_carboxy: 1}),
      ]),
    ("2o8q", [
      Counter({chem_nitrogen: 3, chem_oxygen: 3, chem_nitrogen_secondary: 3,
               chem_water: 3}),
      Counter({chem_nitrogen: 3, chem_oxygen: 3, chem_nitrogen_secondary: 3,
               chem_water: 3}),
      ]),
    ("1tgg", [
      Counter({chem_oxygen: 5, chem_chloride: 1, chem_carboxy: 4,
               chem_water: 1}),
      Counter({chem_oxygen: 3, chem_chloride: 2, chem_carboxy: 3}),
      Counter({chem_oxygen: 4, chem_chloride: 2, chem_carboxy: 4}),
      ]),
    ("3zu8", [
      Counter({chem_oxygen: 7, chem_carboxy: 3, chem_water: 1,
                      chem_backbone: 2}),
      Counter({chem_nitrogen: 4, chem_oxygen: 1, chem_nitrogen_primary: 1,
               chem_nitrogen_secondary: 3, chem_carboxy: 1, chem_backbone: 3}),
      ]),
    ("1ofs", [
      Counter({chem_nitrogen: 1, chem_oxygen: 4, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 1}),
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      ]),
    ("3ul2", [
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_backbone: 1,
               chem_water: 2}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      ]),
    ("3snm", [
      Counter({chem_oxygen: 5, chem_amide: 1, chem_carboxy: 3,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 3, chem_nitrogen_secondary: 1,
               chem_carboxy: 3}),
      ]),
    ("3qlq", [
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      Counter({chem_oxygen: 7, chem_amide: 1, chem_carboxy: 3, chem_water: 2,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 5, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 2}),
      ]),
    ("2gdf", [
      Counter({chem_nitrogen: 1, chem_oxygen: 4, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 1}),
      Counter({chem_oxygen: 6, chem_amide: 1, chem_carboxy: 3, chem_water: 1,
               chem_backbone: 1}),
      Counter({chem_nitrogen: 1, chem_oxygen: 4, chem_nitrogen_secondary: 1,
               chem_carboxy: 3, chem_water: 1}),
      Counter({chem_oxygen: 6, chem_amide: 1, chem_carboxy: 3, chem_water: 1,
               chem_backbone: 1}),
      ]),
    ("1q8h", [
      Counter({chem_oxygen: 7, chem_carboxy: 6, chem_water: 1}),
      Counter({chem_oxygen: 7, chem_carboxy: 4, chem_water: 3}),
      Counter({chem_oxygen: 8, chem_carboxy: 6, chem_water: 2}),
      ]),
  ])

  for model, expected_environments in models.items():
    pdb_path = libtbx.env.find_in_repositories(
      relative_path = os.path.join(
        "phenix_regression", "mmtbx", "ions", model + ".pdb"),
      test = os.path.isfile
      )

    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv = mon_lib_srv,
      ener_lib = ener_lib,
      file_name = pdb_path,
      raw_records = None,
      force_symmetry = True,
      log = libtbx.utils.null_out()
      )

    geometry = \
      processed_pdb_file.geometry_restraints_manager(show_energies = False)
    xray_structure = processed_pdb_file.xray_structure()
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    connectivity = geometry.shell_sym_tables[0].full_simple_connectivity()

    manager = mmtbx.ions.identify.manager(
      fmodel = None,
      pdb_hierarchy = pdb_hierarchy,
      xray_structure = xray_structure,
      connectivity = connectivity)

    elements = set(ions.DEFAULT_IONS + ions.TRANSITION_METALS)
    elements.difference_update(["CL"])

    metals = [i_seq for i_seq, atom in enumerate(manager.pdb_atoms)
             if atom.fetch_labels().resname.strip().upper() in elements]
    assert len(metals) == len(expected_environments)

    for index, metal, expected_environment in \
      zip(range(len(metals)), metals, expected_environments):
      env = ChemicalEnvironment(
        metal,
        manager.find_nearby_atoms(metal, filter_by_two_fofc = False),
        manager
        )
      if env.chemistry != expected_environment:
        print("Problem detecting chemistry environment in", model, index)
        print("Found:    ", env.chemistry)
        print("Should be:", expected_environment)
        sys.exit()

  print("OK")

if __name__ == "__main__":
  exercise()
