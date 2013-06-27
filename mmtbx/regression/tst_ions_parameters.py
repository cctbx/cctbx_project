# -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division

def exercise () :
  from mmtbx.ions.parameters import server
  import iotbx.pdb.hierarchy
  import iotbx.pdb
  from cctbx.eltbx import chemical_elements
  s = server()

  # Assert that valence parameters exist for all common ions with their
  # coordinating atoms
  for elem in ["NA", "MG", "K", "CA", "MN", "FE", "NI", "CO", "CU", "ZN", "CD"]:
    params = s.get_metal_parameters(elem)
    assert (len(params.allowed_coordinating_atoms) > 0)
    assert (params.charge is not None)
    for elem2 in params.allowed_coordinating_atoms :
      atom1 = iotbx.pdb.hierarchy.atom()
      atom1.name = elem
      atom1.element = elem
      atom1.charge = "+" + str(params.charge)
      atom2 = iotbx.pdb.hierarchy.atom()
      atom2.name = elem2
      atom2.element = elem2
      atom2.charge = "{:+}".format(s.get_charge(elem2))
      assert (s.get_valence_params(atom1, atom2) != (None, None))

  # Make sure we don't crash on any ion residue names
  for elem in chemical_elements.proper_upper_list():
    s.get_element(elem)
    s.get_charge(elem)

  # Make sure we don't crash on any common residue names either
  from mmtbx import monomer_library
  from mmtbx.rotamer.rotamer_eval import mon_lib_query
  mon_lib_srv = monomer_library.server.server()
  common_residues = [getattr(iotbx.pdb, i) for i in dir(iotbx.pdb)
                     if i.startswith("common_residue_names") and
                     isinstance(getattr(iotbx.pdb, i), list)]
  for common in common_residues:
    for resn in common:
      mlq = mon_lib_query(resn, mon_lib_srv)
      if mlq is None:
        continue
      for atom in mlq.atom_dict().values():
        s._get_charge_params(resn, element = atom.type_symbol)

  # And we support all waters
  for resn in iotbx.pdb.common_residue_names_water:
    assert s.get_element(resn) == "O"
    assert s.get_charge(resn) == -2

  print "OK"

if (__name__ == "__main__") :
  exercise()
