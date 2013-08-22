 # -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division

import os
import sys

import libtbx
from mmtbx.ions.environment import get_environments
from mmtbx.ions import Manager
from iotbx import file_reader
from cctbx.eltbx import chemical_elements
from collections import OrderedDict, Counter

def exercise () :
  if not libtbx.env.has_module("phenix_regression"):
    print "Skipping {}".format(os.path.split(__file__)[1])
    return

  models = OrderedDict([
    ("3rva", [Counter(O = 6, HOH = 2),
              Counter(O = 4, X2N = 1, HOH = 1, N = 1),
              Counter(N = 4, X2N = 4, Backbone = 3)]),
    ("1mjh", [Counter(O = 6, HOH = 3),
              Counter(O = 6, HOH = 3)]),
    ("4e1h", [Counter(O = 6),
              Counter(O = 6),
              Counter(O = 6)]),
    ("2xuz", [Counter(O = 6)]),
    ("3zli", [Counter(O = 4, N = 2, HOH = 1),
              Counter(S = 4),
              Counter(O = 4, N = 2, HOH = 1),
              Counter(S = 4)]),
    ("3e0f", [Counter(O = 4, N = 2), Counter(O = 2, N = 2),
              Counter(O = 3, N = 2)]),
    ("3dkq", [Counter(N = 4, O = 1), Counter(N = 2, O = 1),
              Counter(N = 4, O = 1)]),
    ("2o8q", [Counter(O = 3, HOH = 3, N = 3),
              Counter(O = 3, HOH = 3, N = 3)]),
    ("1tgg", [Counter(O = 6, HOH = 2),
              Counter(O = 5, HOH = 2),
              Counter(O = 6, HOH = 2)]),
    ("3zu8", [Counter(O = 7, Backbone = 2, HOH = 1),
              Counter(N = 4, Backbone = 3)]),
    ("1ofs", [Counter(O = 4, X2N = 1, HOH = 1, N = 1),
              Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1),
              Counter(O = 7, HOH = 2, Backbone = 1)]),
    ("3ul2", [Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1),
              Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1),
              Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1),
              Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1)]),
    ("3snm", [Counter(O = 5, Backbone = 1),
              Counter(O = 3, X2N = 1, N = 1)]),
    ("3qlq", [Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1),
              Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1),
              Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 7, HOH = 2, Backbone = 1),
              Counter(O = 5, HOH = 2, X2N = 1, N = 1)]),
    ("2gdf", [Counter(O = 4, X2N = 1, HOH = 1, N = 1),
              Counter(O = 6, Backbone = 1, HOH = 1),
              Counter(O = 4, X2N = 1, HOH = 1, N = 1),
              Counter(O = 6, Backbone = 1, HOH = 1)]),
  ])

  for model, expected_environments in models.items():
    pdb_path = libtbx.env.find_in_repositories(
      relative_path = os.path.join("phenix_regression", "mmtbx",
                                   "geometry", model + ".pdb"),
      test = os.path.isfile)
    pdb_in = file_reader.any_file(pdb_path, force_type = "pdb")
    manager = Manager(
      fmodel = None,
      pdb_hierarchy = pdb_in.file_object.construct_hierarchy(),
      xray_structure = pdb_in.file_object.xray_structure_simple(),
      connectivity = pdb_in.file_object.extract_connectivity())

    elements = chemical_elements.proper_upper_list()
    elements = [element for element in elements
                if element.strip().upper() not in ["O", "C", "H", "N"]]

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
        sys.exit(1)

  print "OK"

if __name__ == "__main__":
  exercise()
