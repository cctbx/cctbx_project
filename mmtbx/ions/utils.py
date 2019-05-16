# -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import absolute_import, division, print_function

from collections import defaultdict
from math import cos, pi, sin, sqrt
import sys

from cctbx import crystal
from cctbx.eltbx import sasaki, chemical_elements
from iotbx.pdb import common_residue_names_water as WATER_RES_NAMES
from libtbx.utils import null_out, xfrange
from mmtbx.ions import server
from scitbx.array_family import flex
from scitbx.math import gaussian_fit_1d_analytical
from six.moves import zip
from six.moves import range

def anonymize_ions(pdb_hierarchy, log=sys.stdout):
  """
  Convert any elemental ions in the PDB hierarchy to water, resetting the
  occupancy and scaling the B-factor.  The atom segids will be set to the old
  resname.  NOTE: this does not change the corresponding scatterer in the xray
  structure, but a new xray structure can be obtained by calling
  hierarchy.extract_xray_structure(crystal_symmetry).

  Parameters
  ----------
  pdb_hierarchy : iotbx.pdb.hierarchy.root
  log : file, optional

  Returns
  -------
  iotbx.pdb.hierarchy.root
      New pdb hierarchy with its ions anonymized
  int
      Number of atoms that were anonymized.
  """
  ion_resnames = set(chemical_elements.proper_upper_list())
  for resname in server.params["_lib_charge.resname"]:
    if resname not in WATER_RES_NAMES:
      ion_resnames.add(resname)
  n_converted = 0
  pdb_hierarchy = pdb_hierarchy.deep_copy()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for atom_group in residue_group.atom_groups():
          if atom_group.resname.strip() in ion_resnames:
            atoms = atom_group.atoms()
            id_strs = []
            for atom in atoms:
              elem = atom.element.strip()
              if elem in ["H", "D"]:
                atomic_number = 1
              elif elem in ["HE"]:
                atomic_number = 2
              else:
                atomic_number = sasaki.table(elem).atomic_number()
              id_strs.append(atom.id_str())
              atom.segid = atom_group.resname
              atom.name = " O  "
              atom.element = "O"
              atom.charge = ""
              atom.occupancy = 1.0
              atom.b = atom.b * (10 / atomic_number)
            atom_group.resname = "HOH"
            for atom, id_str in zip(atoms, id_strs):
              print("%s --> %s, B-iso = %.2f" % (id_str, atom.id_str(),
                atom.b), file=log)
              n_converted += 1
  return pdb_hierarchy, n_converted

def sort_atoms_permutation(pdb_atoms, xray_structure):
  """
  Creates a list of atoms in pdb_atoms, sorted first by their atomic number,
  then by occupancy, and finally, by isotropic b-factor.

  Parameters
  ----------
  pdb_atoms : iotbx.pdb.hierarchy.af_shared_atom
  xray_structure : cctbx.xray.structure.structure

  Returns
  -------
  flex.size_t of int
      i_seqs of sorted atoms
  """
  assert pdb_atoms.size() == xray_structure.scatterers().size()
  pdb_atoms.reset_i_seq()
  atoms_sorted = sorted(
    pdb_atoms,
    key=lambda x:
    (sasaki.table(x.element.strip().upper()).atomic_number(), x.occ, x.b),
    reverse=True,
    )
  sele = flex.size_t([atom.i_seq for atom in atoms_sorted])
  return sele

def collect_ions(pdb_hierarchy):
  """
  Collects a list of all ions in pdb_hierarchy.

  Parameters
  ----------
  pdb_hierarchy : iotbx.pdb.hierarchy.root

  Returns
  -------
  list of iotbx.pdb.hierarchy.atom
  """
  elements = chemical_elements.proper_upper_list()
  ions = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for atom_group in residue_group.atom_groups():
          if atom_group.resname.strip() in elements:
            atoms = atom_group.atoms()
            assert len(atoms) == 1
            for atom in atoms:
              ions.append(atom)
  return ions

# TODO add test
def compare_ions(hierarchy, reference_hierarchy, reference_xrs,
    distance_cutoff=2.0, log=None, ignore_elements=(), only_elements=(),
    sel_str_base="segid ION"):
  """
  Compares two pdb structures to determine the number of ions that appear in the
  reference structure and are either matched or missing in the other structure.

  Parameters
  ----------
  hierarchy : iotbx.pdb.hierarchy.root
  reference_hierarchy : iotbx.pdb.hierarchy.root
  reference_xrs : ...
  distance_cutoff : float, optional
  log : file, optional
  ignore_element : iterable, optional
  only_elements : iterable, optional
  sel_str_base : str, optional

  Returns
  -------
  int
      Number of ions in reference_hierarchy that were also found in hierarchy.
  int
      Number of ions in reference_hierarchy that were not found in hierarchy.
  """
  if log is None:
    log = null_out()
  sel_cache = hierarchy.atom_selection_cache()
  sel_str = sel_str_base
  if len(only_elements) > 0:
    sel_str += " and (%s)" % " or ".join(
      ["element %s" % e for e in only_elements])
  elif len(ignore_elements) > 0:
    sel_str += " and not (%s)" % " or ".join(
      ["element %s" % e for e in ignore_elements])
  ion_isel = sel_cache.selection(sel_str).iselection()
  if len(ion_isel) == 0:
    return [], []
  pdb_atoms = hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  ions = pdb_atoms.select(ion_isel)
  asu_mappings = reference_xrs.asu_mappings(
    buffer_thickness=distance_cutoff+0.1)
  unit_cell = reference_xrs.unit_cell()
  sites_cart = ions.extract_xyz()
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  asu_mappings.process_sites_frac(sites_frac,
    min_distance_sym_equiv=reference_xrs.min_distance_sym_equiv())
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff)
  reference_atoms = reference_hierarchy.atoms()
  n_xray = reference_xrs.scatterers().size()
  ion_ref_i_seqs = []
  for k in range(len(ions)):
    ion_ref_i_seqs.append([])
  for pair in pair_generator:
    if ((pair.i_seq < n_xray and pair.j_seq < n_xray) or
        (pair.i_seq >= n_xray and pair.j_seq >= n_xray)):
      continue
    if pair.i_seq < n_xray:
      ion_seq, ref_seq = pair.j_seq, pair.i_seq
    else:
      ion_seq, ref_seq = pair.i_seq, pair.j_seq
    site_frac = sites_frac[ion_seq - n_xray]
    dxyz = sqrt(pair.dist_sq)
    j_seq = ion_seq - n_xray
    ion_ref_i_seqs[j_seq].append(ref_seq)
  # FIXME better filtering needed - right now we risk double-counting ions in
  # the reference model, although I haven't found a case of this yet
  matched = []
  missing = []
  for i_seq, ref_i_seqs in enumerate(ion_ref_i_seqs):
    ion = ions[i_seq]
    if len(ref_i_seqs) == 0:
      print("No match for %s" % ion.id_str(), file=log)
      missing.append(ion.id_str())
    else:
      ref_ions = []
      for i_seq in ref_i_seqs:
        ref_atom = reference_atoms[i_seq]
        if ion.element.upper() == ref_atom.element.upper():
          ref_ions.append(ref_atom.id_str())
      if len(ref_ions) >= 1:
        matched.append(ion.id_str())
        if len(ref_ions) > 1:
          print("Multiple matches for %s:" % ion.id_str(), file=log)
          for ref_ion in ref_ions:
            print("  %s" % ref_ion, file=log)
        else:
          print("Ion %s matches %s" % (ion.id_str(),
            ref_ions[0]), file=log)
      else:
        print("No match for %s" % ion.id_str(), file=log)
        missing.append(ion.id_str())
  return matched, missing

def _get_points_within_radius(point, radius, radius_step=0.2,
                              angle_step=pi / 5):
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

def fit_gaussian(unit_cell, site_cart, real_map, radius=1.6):
  """
  Fit a gaussian function to the map around a site. Samples points in concentric
  spheres up to radius away from the site.

  f(x) = a * exp(-b * x ** 2)

  Parameters
  ----------
  unit_cell : uctbx.unit_cell
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

  See Also
  --------
  scitbx.math.gaussian_fit_1d_analytical
  """
  points, radiuses = _get_points_within_radius(site_cart, radius)

  map_heights = \
    [real_map.tricubic_interpolation(unit_cell.fractionalize(i))
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
    fit = gaussian_fit_1d_analytical(x=x, y=y)
  except RuntimeError as err:
    print(err)
    return 0., 0.
  else:
    return fit.a, fit.b

def count_coordinating_residues(nearby_atoms, distance_cutoff=3.0):
  """
  Count the number of residues of each type involved in the coordination
  sphere.  This may yield additional clues to the identity of ions, e.g. only
  Zn will have 4 Cys residues.

  Parameters
  ----------
  nearby_atoms : list of mmtbx.ions.environment.atom_contact
  distance_cutoff : float, optional

  Returns
  -------
  dict of str, int
  """
  unique_residues = []
  residue_counts = defaultdict(int)
  for contact in nearby_atoms:
    if contact.distance() <= distance_cutoff:
      parent = contact.atom.parent()
      for residue in unique_residues:
        if residue == parent:
          break
      else:
        resname = parent.resname
        residue_counts[resname] += 1
        unique_residues.append(parent)
  return residue_counts
