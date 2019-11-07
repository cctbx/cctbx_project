 # -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import absolute_import, division, print_function
from mmtbx.ions.geometry import find_coordination_geometry
import mmtbx.ions.identify
from mmtbx import ions
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
import libtbx.load_env
from collections import OrderedDict
import os
import sys
from six.moves import zip
from six.moves import range

def exercise():
  if not libtbx.env.has_module("phenix_regression"):
    print("Skipping {}".format(os.path.split(__file__)[1]))
    return

  models = OrderedDict([
    ("2qng", [["pentagonal_bipyramid"],
              ["octahedron"]]),
    ("3rva", [["square_pyramid_bidentate_miss"],
              ["trigonal_bipyramid"],
              ["square_plane"]]),
    ("1mjh", [["octahedron"],
              ["octahedron"]]),
    ("4e1h", [["octahedron"],
              ["octahedron"],
              ["octahedron"]]),
    ("2xuz", [["trigonal_prism"]]),
    ("3zli", [["octahedron"],
              ["tetrahedron"],
              ["octahedron"],
              ["tetrahedron"],]),
    ("3e0f", [["square_pyramid_bidentate_miss"],
              ["trigonal_pyramid"],
              ["square_pyramid"]]),
    ("3dkq", [["square_pyramid"],
              ["three_legs"],
              ["square_pyramid"]]),
    ("2o8q", [["octahedron"],
              ["octahedron"]]),
    ("1tgg", [["octahedron"],
              ["square_pyramid"],
              ["see_saw"]]),
    ("3zu8", [["pentagonal_bipyramid"],
              []]),
    ("1ofs", [["square_pyramid"],
              ["square_pyramid_bidentate"],
              ["octahedron"],
              ["square_pyramid_bidentate"]]),
    ("3ul2", [["square_pyramid_bidentate"],
              ["octahedron"],
              ["square_pyramid_bidentate"],
              ["octahedron"],
              ["square_pyramid_bidentate"],
              ["octahedron"],
              ["square_pyramid_bidentate"],
              ["octahedron"],]),
    ("3snm", [[],
              ["see_saw"]]),
    ("3qlq", [["square_pyramid_bidentate"],
              ["octahedron"],
              ["octahedron"],
              ["square_pyramid_bidentate"],
              ["octahedron"],
              ["square_pyramid_bidentate"],
              ["square_pyramid_bidentate"],
              ["octahedron"]]),
    ("2gdf", [["square_pyramid"],
              ["square_pyramid_bidentate_miss"],
              ["square_pyramid"],
              ["square_pyramid_bidentate_miss"]]),
    ("1q8h", [["trigonal_prism"],
              ["square_pyramid_bidentate"],
              ["square_pyramid_bidentate_miss"]]),
  ])

  for model, geometries in models.items():
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
    assert len(metals) == len(geometries)

    for index, metal, expected_geometry in \
      zip(range(len(metals)), metals, geometries):
      contacts = manager.find_nearby_atoms(metal, filter_by_two_fofc = False)
      found = find_coordination_geometry(contacts, minimizer_method = True)
      geometry_names = [i[0] for i in found]

      if geometry_names != expected_geometry:
        print("Problem detecting geometries in", model, index)
        print(manager.pdb_atoms[metal].id_str())
        print("Found geometries:", geometry_names)
        print("Should be:", expected_geometry)
        sys.exit()

  print("OK")

if __name__ == "__main__":
  exercise()
