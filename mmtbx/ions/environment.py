# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
Deals with examing the atoms around a site and recognizing distinct and useful
chemical environments.
"""
from __future__ import division

from collections import Counter
from math import pi, cos, sin

from iotbx.pdb import common_residue_names_get_class as get_class
from libtbx.utils import Sorry, xfrange
from mmtbx import ions
from mmtbx.ions.parameters import MetalParameters, server
from scitbx.array_family import flex
from scitbx.math import gaussian_fit_1d_analytical
from scitbx.matrix import col

# Enums for the chemical environments supported by this module
N_SUPPORTED_ENVIRONMENTS = 14
chem_carboxy, \
  chem_amide, \
  chem_backbone, \
  chem_water, \
  chem_sulfate, \
  chem_phosphate, \
  chem_disulfide, \
  chem_nitrogen_primary, \
  chem_nitrogen_secondary, \
  chem_nitrogen_tertiary, \
  chem_chlorine, \
  chem_oxygen, \
  chem_nitrogen, \
  chem_sulfur = range(N_SUPPORTED_ENVIRONMENTS)

METAL_SERVER = server()

class ScatteringEnvironment (object):
  def __init__(self, i_seq, manager, fo_map, fofc_map):
    atom = manager.pdb_atoms[i_seq]
    self.d_min = manager.fmodel.f_obs().d_min()
    self.wavelength = manager.wavelength
    self.fp, self.fpp = manager.get_fp(i_seq), manager.get_fpp(i_seq)
    self.fo_density = _fit_gaussian(manager, atom.xyz, fo_map)
    self.fofc_density = _fit_gaussian(manager, atom.xyz, fofc_map)
    self.b_iso = manager.get_b_iso(i_seq)
    self.occ = atom.occ

class ChemicalEnvironment (object):
  def __init__(self, i_seq, contacts, manager):
    self.atom = manager.pdb_atoms[i_seq].fetch_labels()

    self.contacts = contacts
    self.contacts_no_alts = []

    # Filter down the list of contacts so it doesn't include alt locs
    for index, contact in enumerate(contacts):
      if contact.element in ["H", "D"]:
        continue

      no_alt_loc = True
      for other_contact in contacts[index + 1:]:
        if ions.same_atom_different_altloc(contact.atom, other_contact.atom) and \
          contact.rt_mx == other_contact.rt_mx:
          no_alt_loc = False
          break
      if no_alt_loc:
        self.contacts_no_alts.append(contact)

    self.chemistry = self._get_chemical_environment(
      self.contacts_no_alts, manager)
    self.geometry = ions.geometry.find_coordination_geometry(
      self.contacts_no_alts, minimizer_method = True)

    # We can either store the contacts or generate a list of valences for all of
    # the atom types we want to test. I'm choosing the former because it's more
    # flexible down the road and contacts is pickable anyways.
    # self.valences = None

  def get_valence(self, element, charge = None):
    """
    Calculates the BVS and VECSUM for a given element / charge identity.

    Parameters
    ----------
    element : str
        The element identity to calculate valences for.
    charge : int, optional
        The charge to calculate valences for.

    Returns
    -------
    float
        BVS
    float
        VECSUM

    Notes
    -----
    .. [1] Müller, P., Köpke, S. & Sheldrick, G. M. Is the bond-valence method
           able to identify metal atoms in protein structures? Acta
           Crystallogr. D. Biol. Crystallogr. 59, 32–7 (2003).
    """
    if charge is None:
      charge = METAL_SERVER.get_charge(element)
    ion_params = MetalParameters(element = element, charge = charge)
    vectors = METAL_SERVER.calculate_valences(ion_params, self.contacts)
    bvs = sum([abs(i) for i in vectors])
    vecsum = abs(sum(vectors, col((0, 0, 0))))
    return bvs, vecsum

  def _get_chemical_environment(self, contacts, manager):
    """
    Examines an atom contact object for chemical properties useful to ion
    picking. These properties include degree of connectivitiy, presence in
    carboxy and amide groups, disulfide bonds, and so on.

    Parameters
    ----------
    contacts : list of mmtbx.ions.atom_contact
        The contact objects to examine the chemical environments of. Generally
        created by Manager.find_nearby_atoms.
    manager : mmtbx.ions.Manager
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
      seq : int
          The atom.i_seq to find the neighbors of.

      Returns
      -------
      list of int
          A list of neighboring atoms that does not include any hydrogens.
      """
      neighbors = []
      for other in manager.connectivity[seq]:
        other_atom = manager.pdb_atoms[other]
        if manager.server.get_element(other_atom) not in ["H", "D"]:
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
          if ions.same_atom_different_altloc(
              manager.pdb_atoms[i_seq], manager.pdb_atoms[j_seq]):
            no_alt_loc = False
            break
        if no_alt_loc:
          non_alt_locs.append(i_seq)
      return len(non_alt_locs)

    chem_env = Counter()

    for contact in contacts:
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

def _get_points_within_radius(point, radius, radius_step = 0.2,
                              angle_step = pi / 5):
  """
  Generates a list of points and their associated radius in steps around a
  sphere.

  Parameters
  ----------
  point : tuple of float, float, float
      X, Y, Z, coordinates to center the sampling around.
  radius : float
      Max radius around the center to sample.
  radius_step : float, optional
      Steps along the radius to use when sampling.
  angle_step : float, optional
      Steps around each radii distance to use when sampling. Amount is in
      radians.

  Returns
  -------
  list of tuple of float, float, float
      List of points to be sampled.
  list of float
      List of radii corresponding to each point.
  """

  points = [point]
  radiuses = [0]
  for r in xfrange(radius_step, radius, radius_step):
    for theta in xfrange(-pi, pi, angle_step):
      for phi in xfrange(-pi, pi, angle_step):
        x = r * cos(theta) * sin(phi) + point[0]
        y = r * sin(theta) * sin(phi) + point[1]
        z = r * cos(phi) + point[2]
        points.append((x, y, z))
        radiuses.append(r)

  return points, radiuses

def _fit_gaussian(manager, site_cart, real_map, radius = 1.6):
  """
  Fit a gaussian function to the map around a site. Samples points in concentric
  spheres up to radius away from the site.

  f(x) = a * exp(-b * x ** 2)

  Parameters
  ----------
  manager : mmtbx.ions.Manager
  site_cart : tuple of float, float, float
      The site's cartesian coordinates to sample the density around.
  real_map : scitbx.array_family.flex
      Real space map of the electron density in the unit cell.
  radius : float, optional
      The max radius to use for sampling.

  Returns
  -------
  float
      Height of gaussian curve.
  float
      Spread of guassian curve.
  """
  points, radiuses = _get_points_within_radius(site_cart, radius)

  map_heights = \
    [real_map.tricubic_interpolation(manager.unit_cell.fractionalize(i))
     for i in points]

  # Gaussian functions can't have negative values, filter sampled points below
  # zero to allow us to find the analytical solution (radius = 2.0 is too big
  # for most atoms anyways)
  x, y = flex.double(), flex.double()
  for rad, height in zip(radiuses, map_heights):
    if height > 0:
      x.append(rad)
      y.append(height)

  try:
    fit = gaussian_fit_1d_analytical(x = x, y = y)
  except RuntimeError as err:
    print(err)
    return 0., 0.
  else:
    return fit.a, fit.b
