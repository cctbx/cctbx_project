# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
Deals with examing the atoms around a site and recognizing distinct and useful
chemical environments.
"""
from __future__ import division

from collections import Counter

from iotbx.pdb import common_residue_names_get_class as get_class
from libtbx.utils import Sorry
from mmtbx.ions import same_atom_different_altloc

# Enums for the chemical environments supported by this module
chem_carboxy = 1
chem_amide = 2
chem_backbone = 3
chem_water = 4
chem_sulfate = 5
chem_phosphate = 6
chem_disulfide = 7
chem_nitrogen_primary = 8
chem_nitrogen_secondary = 9
chem_nitrogen_tertiary = 10
chem_chlorine = 11
chem_oxygen = 12
chem_nitrogen = 13
chem_sulfur = 14

class Environment (object):
  def __init__(self, contacts, manager):
    # XXX: This constructor needs a better API...should we just take the i_seq
    # and generate all other information (Including the list of contacts) from
    # that and the manager?
    self.chemistry_env = self._get_chemical_environment(contacts, manager)
    # from mmtbx.ions import geometry
    # self.geometry = geometry.find_coordination_geometry(contacts, minimizer_method = True)
    # self.electron_density = None
    # self.anomalous = None
    # self.valences = None

  def _get_chemical_environment(self, contacts, manager):
    """
    Examines an atom contact object for chemical properties useful to ion
    picking. These properties include degree of connectivitiy, presence in
    carboxy and amide groups, disulfide bonds, and so on.

    Parameters
    ----------
    contacts: list of mmtbx.ions.atom_contact
        The contact objects to examine the chemical environments of. Generally
        created by Manager.find_nearby_atoms.
    manager: mmtbx.ions.Manager
        The ions manager that created contact. Used for auxillary information
        such as bond connectivity.

    Returns
    -------
    collections.Counter
        The chemical environments found for contacts.
    """
    def _non_hydrogen_neighbors(seq):
      """
      Parameters
      ----------
      seq: int
          The atom.i_seq to find the neighbors of.

      Returns
      -------
      list of int
          A list of neighboring atoms that does not include any hydrogens.
      """
      neighbors = []
      for other in manager.connectivity[seq]:
        if manager.server.get_element(manager.pdb_atoms[other]) not in ["H"]:
          neighbors.append(other)
      return neighbors

    def _n_non_altlocs(i_seqs):
      """
      Count the number of neighbors, filtering out those that are alt-locs of
      other neighboring atoms.
      """
      non_alt_locs = []
      for index, i_seq in enumerate(i_seqs):
        no_alt_loc = True
        for j_seq in i_seqs[index + 1:]:
          if same_atom_different_altloc(
              manager.pdb_atoms[i_seq], manager.pdb_atoms[j_seq]):
            no_alt_loc = False
            break
        if no_alt_loc:
          non_alt_locs.append(i_seq)
      return len(non_alt_locs)

    chem_env = Counter()

    # Filter down the list of contacts so it doesn't include alt locs
    non_alt_locs = []
    for index, contact in enumerate(contacts):
      no_alt_loc = True
      for other_contact in contacts[index + 1:]:
        if same_atom_different_altloc(contact.atom, other_contact.atom) and \
          contact.rt_mx == other_contact.rt_mx:
          no_alt_loc = False
          break
      if no_alt_loc:
        non_alt_locs.append(contact)

    for contact in non_alt_locs:
      # Add any of our included elements
      element = contact.element
      i_seq = contact.atom.i_seq
      neighbor_elements = [manager.server.get_element(manager.pdb_atoms[i])
                           for i in _non_hydrogen_neighbors(i_seq)]

      n_neighbors = _n_non_altlocs(_non_hydrogen_neighbors(i_seq))

      # Check for waters, sulfates, phosphates, chlorides, etc
      if element == "CL":
        chem_env[chem_chlorine] += 1
      elif element == "O":
        chem_env[chem_oxygen] += 1
        if get_class(contact.resname()) in ["common_water"]:
          chem_env[chem_water] += 1
        elif "S" in neighbor_elements:
          chem_env[chem_sulfate] += 1
        elif "P" in neighbor_elements:
          chem_env[chem_phosphate] += 1
      elif element == "N":
        chem_env[chem_nitrogen] += 1
        # Check the degree of connectivity on nitrogens
        if n_neighbors == 1:
          chem_env[chem_nitrogen_primary] += 1
        elif n_neighbors == 2:
          chem_env[chem_nitrogen_secondary] += 1
        elif n_neighbors == 3:
          chem_env[chem_nitrogen_tertiary] += 1
        elif n_neighbors != 0:
          raise Sorry("Nitrogen with more than three coordinating atoms: " +
                      contact.atom.id_str())
      elif element == "S":
        chem_env[chem_sulfur] += 1
        # Check for disulfide bridges
        if "S" in neighbor_elements:
          chem_env[chem_disulfide] += 1

      # Check for backbone nitrogens
      if contact.atom.name.strip().upper() in ["N", "C", "O", "CA"] and \
        get_class(contact.resname()) in ["common_amino_acid"]:
        chem_env[chem_backbone] += 1

      # Check for carboxy / amide groups
      if element in ["O", "N"] and n_neighbors == 1:
        for j_seq in _non_hydrogen_neighbors(i_seq):
          j_atom = manager.pdb_atoms[j_seq]
          if manager.server.get_element(j_atom) in ["C"]:
            for k_seq in manager.connectivity[j_seq]:
              k_neighbors = _non_hydrogen_neighbors(k_seq)
              k_atom = manager.pdb_atoms[k_seq]
              if k_seq != i_seq and _n_non_altlocs(k_neighbors) in [1]:
                k_element = manager.server.get_element(k_atom)
                if set([element, k_element]) == set(["O", "O"]):
                  chem_env[chem_carboxy] += 1
                elif set([element, k_element]) == set(["N", "O"]):
                  chem_env[chem_amide] += 1

    return chem_env
