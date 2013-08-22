 # -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division

from iotbx.pdb import common_residue_names_get_class as get_class
from libtbx.utils import Sorry
from mmtbx import ions

SUPPORTED_ENVIRONMENTS = set([
  "O", "N", "S", "CL",
  "Carboxy", "Amide",
  "Backbone",
  "HOH", "SO", "PO",
  "XN", "X2N", "X3N",
  "SS"
  ])

# XXX: Deal with alt loc environments?
def get_environments(contact, manager):
  def _non_hydrogen_neighbors(seq):
    return [other
            for other in manager.connectivity[seq]
            if manager.server.get_element(manager.pdb_atoms[other]) not in ["H"]]

  def _n_non_altlocs(i_seqs):
    """
    Count the number of neighbors, filtering out those that are alt-locs of
    other neighboring atoms.
    """
    return len(
      [i_seq
       for index, i_seq in enumerate(i_seqs)
       if all(not ions.same_atom_different_altloc(
           manager.pdb_atoms[i_seq],
           manager.pdb_atoms[j_seq])
           for j_seq in i_seqs[index + 1:])]
      )

  environments = set()

  # Add any of our included elements
  element = contact.element

  if element in SUPPORTED_ENVIRONMENTS:
    environments.add(element)

  i_seq = contact.atom.i_seq
  neighbor_elements = [manager.server.get_element(manager.pdb_atoms[i])
                       for i in _non_hydrogen_neighbors(i_seq)]

  # Check for waters, sulfates, and phosphates
  if element in ["O"]:
    if get_class(contact.resname()) in ["common_water"]:
      environments.add("HOH")
    elif "S" in neighbor_elements:
      environments.add("SO")
    elif "P" in neighbor_elements:
      environments.add("PO")

  n_neighbors = _n_non_altlocs(_non_hydrogen_neighbors(i_seq))

  # Check the degree of connectivity on nitrogens
  if element in ["N"]:
    if n_neighbors in [1]:
      environments.add("XN")
    elif n_neighbors in [2]:
      environments.add("X2N")
    elif n_neighbors in [3]:
      environments.add("X3N")
    elif n_neighbors not in [0]:
      raise Sorry("Nitrogen with more than three coordinating atoms: " +
                  contact.atom.id_str())

  # Check for disulfide bridges
  if element in ["S"]:
    if "S" in neighbor_elements:
      environments.add("SS")

  # Check for backbone nitrogens
  if contact.atom.name.strip().upper() in ["N", "C", "O", "CA"] and \
    get_class(contact.resname()) in ["common_amino_acid"]:
    environments.add("Backbone")

  # Check for carboxy / amide groups
  if element in ["O", "N"] and n_neighbors in [1]:
    for j_seq in _non_hydrogen_neighbors(i_seq):
      j_atom = manager.pdb_atoms[j_seq]
      if manager.server.get_element(j_atom) in ["C"]:
        for k_seq in manager.connectivity[j_seq]:
          k_neighbors = _non_hydrogen_neighbors(k_seq)
          k_atom = manager.pdb_atoms[k_seq]
          if k_seq != i_seq and _n_non_altlocs(k_neighbors) in [1]:
            k_element = manager.server.get_element(k_atom)
            if [element, k_element] in [["O", "O"]]:
              environments.add("Carboxy")
            elif [element, k_element] in [["N", "O"], ["O", "N"]]:
              environments.add("Amide")

  return environments
