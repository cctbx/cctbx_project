 # -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division

import os
import sys

from collections import OrderedDict, Counter

import libtbx
from mmtbx.ions.environment import get_environments
from mmtbx import ions
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library

def exercise () :
  if not libtbx.env.has_module("phenix_regression"):
    print "Skipping {}".format(os.path.split(__file__)[1])
    return

  models = OrderedDict([
    ("3rva", [Counter(O = 6, Carboxy = 4, HOH = 2),
              Counter(N = 1, O = 4, X2N = 1, Carboxy = 3, HOH = 1),
              Counter(N = 4, XN = 1, X2N = 3, Backbone = 3)]),
    ("1mjh", [Counter(O = 6, HOH = 3, PO = 3),
              Counter(O = 6, HOH = 3, PO = 3)]),
    ("4e1h", [Counter(O = 6, Carboxy = 4),
              Counter(O = 6, Carboxy = 3),
              Counter(O = 6, Carboxy = 3)]),
    ("2xuz", [Counter(O = 6)]),
    ("3zli", [Counter(N = 2, O = 4, X2N = 2, Carboxy = 1, HOH = 1),
              Counter(S = 4),
              Counter(N = 2, O = 4, X2N = 2, Carboxy = 1, HOH = 1),
              Counter(S = 4)]),
    ("3e0f", [Counter(N = 2, O = 4, X2N = 2, Carboxy = 2, PO = 2),
              Counter(N = 2, O = 2, X2N = 2, Carboxy = 1, PO = 1),
              Counter(N = 2, O = 3, X2N = 2, Carboxy = 2, PO = 1)]),
    ("3dkq", [Counter(N = 4, O = 1, X2N = 4, Carboxy = 1),
              Counter(N = 2, O = 1, X2N = 2, Carboxy = 1),
              Counter(N = 4, O = 1, X2N = 4, Carboxy = 1)]),
    ("2o8q", [Counter(N = 3, O = 3, X2N = 3, HOH = 3),
              Counter(N = 3, O = 3, X2N = 3, HOH = 3)]),
    ("1tgg", [Counter(O = 5, CL = 1, Carboxy = 4, HOH = 1),
              Counter(O = 3, CL = 2, Carboxy = 3),
              Counter(O = 4, CL = 2, Carboxy = 4)]),
    ("3zu8", [Counter(O = 7, Carboxy = 3, HOH = 1, Backbone = 2),
              Counter(N = 4, O = 1, XN = 1, X2N = 3, Carboxy = 1,
                      Backbone = 3)]),
    ("1ofs", [Counter(N = 1, O = 4, X2N = 1, Carboxy = 3, HOH = 1),
              Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2),
              Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1)]),
    ("3ul2", [Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2),
              Counter(O = 7, Amide = 1, Carboxy = 3, Backbone = 1, HOH = 2),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2),
              Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2),
              Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2)]),
    ("3snm", [Counter(O = 5, Amide = 1, Carboxy = 3, Backbone = 1),
              Counter(N = 1, O = 3, X2N = 1, Carboxy = 3)]),
    ("3qlq", [Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2),
              Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2),
              Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1),
              Counter(O = 7, Amide = 1, Carboxy = 3, HOH = 2, Backbone = 1),
              Counter(N = 1, O = 5, X2N = 1, Carboxy = 3, HOH = 2)]),
    ("2gdf", [Counter(N = 1, O = 4, X2N = 1, Carboxy = 3, HOH = 1),
              Counter(O = 6, Amide = 1, Carboxy = 3, HOH = 1, Backbone = 1),
              Counter(N = 1, O = 4, X2N = 1, Carboxy = 3, HOH = 1),
              Counter(O = 6, Amide = 1, Carboxy = 3, HOH = 1, Backbone = 1)]),
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

    manager = ions.Manager(
      fmodel = None,
      pdb_hierarchy = pdb_hierarchy,
      xray_structure = xray_structure,
      connectivity = connectivity)

    elements = set(ions.DEFAULT_IONS + ions.TRANSITION_METALS)
    elements.difference_update(["CL"])

    metals = [i_seq for i_seq, atom in enumerate(manager.pdb_atoms)
             if atom.fetch_labels().resname.strip().upper() in elements]
    assert len(metals) == len(expected_environments)

    for index, metal, expected_environment in zip(xrange(100), metals,
                                                   expected_environments):
      contacts = manager.find_nearby_atoms(metal, filter_by_two_fofc = False)
      environments = Counter([env for contact in contacts
                              for env in get_environments(contact, manager)])
      if environments != expected_environment:
        print "Problem detecting environments in", model, index
        print "Found environments:", environments
        print "Should be:", expected_environment
        sys.exit()

  print "OK"

if __name__ == "__main__":
  exercise()
