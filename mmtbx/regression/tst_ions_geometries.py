 # -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division

import os
import sys

import libtbx
from mmtbx.ions.geometry import find_coordination_geometry
from mmtbx.ions import Manager
from iotbx import file_reader
from cctbx.eltbx import chemical_elements
from collections import OrderedDict

def exercise () :
  if not libtbx.env.has_module("phenix_regression"):
    print "Skipping {}".format(os.path.split(__file__)[1])
    return

  models = OrderedDict([
    ("3rva", [["octahedron"], ["trigonal_bipyramid"], ["square_plane"]]),
    ("1mjh", [["octahedron"], ["octahedron"]]),
    ("4e1h", [["octahedron"], ["octahedron"], ["octahedron"]]),
    ("2xuz", [["trigonal_prism"]]),
    ("3zli", [["octahedron"], ["tetrahedron"],
              ["octahedron"], ["tetrahedron"],]),
    ("3e0f", [["octahedron"], ["trigonal_pyramid"], ["square_pyramid"]]),
    ("3dkq", [["square_pyramid"], ["three_legs"], ["square_pyramid"]]),
    ("2o8q", [["octahedron"], ["octahedron"]]),
    ("1tgg", [["octahedron"], ["square_pyramid"], ["see_saw"]]),
    ("3zu8", [["pentagonal_bipyramid"], []]),
    ("1ofs", [["square_pyramid"], ["ring_pop"], ["octahedron"],
              ["ring_pop"]]),
    ("3ul2", [["ring_pop"], ["octahedron"],
              ["ring_pop"], ["octahedron"],
              ["ring_pop"], ["octahedron"],
              ["ring_pop"], ["octahedron"],]),
    ("3snm", [[], ["see_saw"]]),
    ("3qlq", [["ring_pop"], ["octahedron"], ["octahedron"], ["ring_pop"],
              ["octahedron"], ["ring_pop"], ["ring_pop"], ["octahedron"]]),
    ("2gdf", [["square_pyramid"], ["pentagonal_pyramid"],
              ["square_pyramid"], ["pentagonal_pyramid"]]),
  ])

  for model, geometries in models.items():
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
    for index, metal, expected_geometry in zip(xrange(100), metals, geometries):
      contacts = manager.find_nearby_atoms(metal, filter_by_two_fofc = False)
      found = find_coordination_geometry(contacts, minimizer_method = True)
      geometry_names = [i[0] for i in found]
      if geometry_names != expected_geometry:
        print "Problem detecting geometries in", model, index
        print "Found geometries:", geometry_names
        print "Should be:", expected_geometry
        sys.exit(1)

  print "OK"

if __name__ == "__main__":
  exercise()
