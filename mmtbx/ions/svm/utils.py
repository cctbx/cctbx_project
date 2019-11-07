# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
Utility functions used within this module.
"""
from __future__ import absolute_import, division, print_function

try : # XXX required third-party dependencies
  import numpy as np
except ImportError :
  np = None

from libtbx import Auto
from mmtbx.ions import server
from mmtbx.ions import halides
from scitbx.matrix import col

def iterate_sites(pdb_hierarchy, split_sites=False, res_filter=None):
  """
  Returns a generator iterating over all atoms in pdb_hierarchy. Optionally
  skips sites with alternate conformations and can filter by residue name.

  Parameters
  ----------
  pdb_hierarchy: iotbx.pdb.hierarchy.root
  split_sites: bool, optional
      Indicates whether to iterate over sites with alternate conformations, by
      default they are not included.
  res_filter: list of str, optional
      List of residue names to include, by default, all residues are examined.

  Returns
  -------
  generator of iotbx.pdb.hierarchy.atom
  """

  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for atom_group in residue_group.atom_groups():
          resname = atom_group.resname.strip().upper()
          if atom_group.altloc.strip() != "" and not split_sites:
            continue
          if res_filter is None or resname in res_filter:
            atoms = atom_group.atoms()
            if len(atoms) == 1:
              atom = atoms[0]
              element = atom.element.strip().upper()
              if element in ["H", "D"]:
                continue
              yield atom

def _is_favorable_halide_environment(
    chem_env, scatter_env, assume_hydrogens_all_missing=Auto,
    ):
  assume_hydrogens_all_missing = True
  atom = chem_env.atom
  binds_amide_hydrogen = False
  near_cation = False
  near_lys = False
  near_hydroxyl = False
  xyz = col(atom.xyz)
  min_distance_to_cation = None
  min_distance_to_hydroxyl = min_distance_to_cation
  for contact in chem_env.contacts:
    other = contact.atom
    resname = contact.resname()
    atom_name = contact.atom_name()
    element = contact.element
    distance = abs(contact)
    # XXX need to figure out exactly what this should be - CL has a
    # fairly large radius though (1.67A according to ener_lib.cif)
    if distance < 2.5:
      return False
    if not element in ["C", "N", "H", "O", "S"]:
      charge = server.get_charge(element)
      if charge < 0 and distance <= 3.5:
        # Nearby anion that is too close
        return False
      if charge > 0 and distance <= 3.5:
        # Nearby cation
        near_cation = True
        if min_distance_to_cation is None or \
          distance < min_distance_to_cation:
          min_distance_to_cation = distance
    # Lysine sidechains (can't determine planarity)
    elif (atom_name in ["NZ"] and #, "NE", "NH1", "NH2"] and
          resname in ["LYS"] and
          distance <= 3.5):
      near_lys = True
      if min_distance_to_cation is None or \
        distance < min_distance_to_cation:
        min_distance_to_cation = distance
    # sidechain amide groups, no hydrogens (except Arg)
    # XXX this would be more reliable if we also calculate the expected
    # hydrogen positions and use the vector method below
    elif (atom_name in ["NZ", "NH1", "NH2", "ND2", "NE2"] and
          resname in ["ARG", "ASN", "GLN"] and
          (assume_hydrogens_all_missing or resname == "ARG") and
          distance <= 3.5):
      binds_amide_hydrogen = True
      if resname == "ARG" and (
          min_distance_to_cation is None or
          distance < min_distance_to_cation):
        min_distance_to_cation = distance
    # hydroxyl groups - note that the orientation of the hydrogen is usually
    # arbitrary and we can't determine precise bonding
    elif ((atom_name in ["OG1", "OG2", "OH1"]) and
          (resname in ["SER", "THR", "TYR"]) and
          (distance <= 3.5)):
      near_hydroxyl = True
      if distance < min_distance_to_hydroxyl:
        min_distance_to_hydroxyl = distance
    # Backbone amide, implicit H
    elif atom_name in ["N"] and assume_hydrogens_all_missing:
      binds_amide_hydrogen = True
      # xyz_n = col(contact.site_cart)
      # bonded_atoms = connectivity[j_seq]
      # ca_same = c_prev = None
      # for k_seq in bonded_atoms:
      #   other2 = pdb_atoms[k_seq]
      #   if other2.name.strip().upper() in ["CA"]:
      #     ca_same = col(get_site(k_seq))
      #   elif other2.name.strip().upper() in ["C"]:
      #     c_prev = col(get_site(k_seq))
      # if ca_same is not None and c_prev is not None:
      #   xyz_cca = (ca_same + c_prev) / 2
      #   vec_ncca = xyz_n - xyz_cca
      #   # 0.86 is the backbone N-H bond distance in geostd
      #   xyz_h = xyz_n + (vec_ncca.normalize() * 0.86)
      #   vec_nh = xyz_n - xyz_h
      #   vec_nx = xyz_n - xyz
      #   angle = abs(vec_nh.angle(vec_nx, deg=True))
      #   if abs(angle - 180) <= 20:
      #     binds_amide_hydrogen = True
  # now check again for negatively charged sidechain (etc.) atoms (e.g.
  # carboxyl groups), but with some leeway if a cation is also nearby.
  # backbone carbonyl atoms are also excluded.
  for contact in chem_env.contacts:
    if contact.altloc() not in ["", "A"]:
      continue
    resname = contact.resname()
    atom_name = contact.atom_name()
    distance = abs(contact)
    if ((distance < 3.2) and
        (min_distance_to_cation is not None and
         distance < (min_distance_to_cation + 0.2)) and
        halides.is_negatively_charged_oxygen(atom_name, resname)):
      return False
  return binds_amide_hydrogen or near_cation or near_lys

def filter_svm_outputs(chem_env, scatter_env, predictions):
  """
  Applies a simple set of filters to the accepted ions that might match a given
  chemical and scattering environment to help catch corner cases where the SVM
  might fail.

  Parameters
  ----------
  chem_env : mmtbx.ions.environment.ChemicalEnvironment
  scatter_env : mmtbx.ions.environment.ScatteringEnvironment
  elements : list of str

  Returns
  -------
  list of str
  """
  bvs_ratio = 0.5
  vecsum_cutoff = 0.6
  ok_elements = []
  for element, score in predictions:
    if element != "HOH":
      if scatter_env.fo_density[0] < 1:
        continue
      if element not in ["NA", "MG"]:
        if scatter_env.fofc_density[0] < 0:
          continue
      bvs, vecsum = chem_env.get_valence(element)
      # bvs <= 0.5 * lower or bvs >= 1.5 * upper
      if element not in ["F", "CL", "BR", "I"]:
        # Require cations have okay BVS values
        charges = server.get_charges(element)
        if bvs <= (1 - bvs_ratio) * min(charges):
          continue
        if bvs >= (1 + bvs_ratio) * max(charges):
          continue
      else:
        # Require halides be touching at least one atom with a positive charge
        if not _is_favorable_halide_environment(chem_env, scatter_env):
          continue
        # if len(chem_env.contacts) == 0 or \
        #   not any(server.get_charge(i.atom) > 0 for i in chem_env.contacts):
        #   print [(i.atom.id_str(), server.get_charge(i.atom)) for i in chem_env.contacts]
        #   continue
      if vecsum > vecsum_cutoff:
        continue
      if any(abs(i.vector) < 1.8 for i in chem_env.contacts):
        continue
    ok_elements.append((element, score))
  return ok_elements

def scale_to(matrix, source, target):
  """
  Given an upper and lower bound for each row of matrix, scales the values to be
  within the range specified by target.

  Parameters
  ----------
  matrix : numpy.array of float
      The matrix to be scaled.
  source : tuple of numpy.array of float
      The upper and lower bound on the values of each row in the original
      matrix.
  target : tuple of float
      The target range to scale to.

  Returns
  -------
  matrix : numpy.array of float
      The matrix with scaled values.

  Examples
  --------
  >>> from mmtbx.ions.svm.utils import scale_to
  >>> import numpy as np
  >>> matrix = np.array([[0, 1, 2],
                         [2, 3, 4],
                         [1, 2, 3]])
  >>> source = (np.array([2, 3, 4]),
                np.array([0, 1, 2]))
  >>> target = (0, 1)
  >>> scale_to(matrix, source, target)
  array([[ 1. ,  1. ,  1. ],
         [ 0. ,  0. ,  0. ],
         [ 0.5,  0.5,  0.5]])
  """
  matrix = np.array(matrix)
  keep_rows = source[0] != source[1]
  matrix = matrix[:, keep_rows]
  source = (source[0][keep_rows], source[1][keep_rows])
  return (matrix - source[0]) * (target[1] - target[0]) / \
    (source[1] - source[0]) + target[0]
