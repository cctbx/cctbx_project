# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
A set of functions to act upon ChemicalEnvironment and ScatteringEnvironment and
produce a single class and vector of features for use with a classifier.

See Also
--------
mmtbx.ions.environment
mmtbx.ions.geometry
"""
from __future__ import division

from collections import Iterable
import numpy as np

from cctbx.eltbx import sasaki
from mmtbx import ions
from mmtbx.ions.environment import N_SUPPORTED_ENVIRONMENTS
from mmtbx.ions.geometry import SUPPORTED_GEOMETRY_NAMES, \
     find_coordination_geometry
from mmtbx.ions.parameters import MetalParameters
from mmtbx.ions.parameters import server as ions_server

METAL_SERVER = None
ALLOWED_IONS = [ions.WATER_RES_NAMES[0]] + ["MN", "ZN", "FE", "NI", "CA"]

def ion_class(chem_env):
  """
  Returns the class name associated with the ion, analogous to the chemical
  ID.

  Parameters
  ----------
  chem_env: mmtbx.ions.environment.ChemicalEnvironment
      The object to extract the class from.
  Returns
  -------
  str
      The class associated with the ion.
  """
  if hasattr(chem_env.atom, "segid") and chem_env.atom.segid.strip():
    return chem_env.atom.segid.strip().upper()
  return chem_env.atom.resname.strip().upper()

def ion_vector(chem_env, scatter_env, server = None):
  """
  Creates a vector containing all of the useful properties contained
  within ion. Merges together the vectors from ion_*_vector().

  Parameters
  ----------
  chem_env: mmtbx.ions.environment.ChemicalEnvironment
      A object containing information about the chemical environment at a site.
  scatter_env: mmtbx.ions.environment.ScatteringEnvironment, optional
      An object containing information about the scattering environment at a
      site.
  server: mmtbx.ions.parameters.server, optional
      A server used to find charges of elements.

  Returns
  -------
  numpy.array of float
      A vector containing quantitative properties for classification.

  See Also
  --------
  ion_model_vector()
  ion_electron_density_vector()
  ion_geometry_vector()
  ion_nearby_atoms_vector()
  ion_valence_vector()
  ion_anomalous_vector()
  """
  if server is None:
    global METAL_SERVER
    if METAL_SERVER is None:
      METAL_SERVER = ions_server()
    server = METAL_SERVER
  return np.concatenate((
    ion_model_vector(scatter_env),
    ion_electron_density_vector(scatter_env),
    ion_geometry_vector(chem_env),
    ion_nearby_atoms_vector(chem_env),
    ion_valence_vector(chem_env, server = server),
    ion_anomalous_vector(scatter_env),
    ))

def ion_model_vector(scatter_env):
  """
  Creates a vector containing information about the general properties of the
  model in which the site is found. Currently this only includes the minimum
  resolution of the data set.

  Parameters
  ----------
  scatter_env: mmtbx.ions.environment.ScatteringEnvironment
      An object containing information about the scattering environment at a
      site.

  Returns
  -------
  numpy.array of float
      A vector containing quantitative properties for classification.
  """
  return np.array([
    scatter_env.d_min,
    ], dtype = float)

def ion_electron_density_vector(scatter_env):
  """
  Creates a vector containing information about the electron density around
  the site. Currently this only includes the site's peak in the 2FoFc map. May
  be expanded in the future to include information about the volume of density
  around the site.

  Parameters
  ----------
  scatter_env: mmtbx.ions.environment.ScatteringEnvironment
      An object containing information about the scattering environment at a
      site.

  Returns
  -------
  numpy.array of float
      A vector containing quantitative properties for classification.
  """
  return np.array([
    scatter_env.electron_density[0],
    scatter_env.electron_density[1],
    ], dtype = float)

def ion_b_factor_occ_vector(scatter_env, elements = None):
  """
  Calculates the theoretical b-factors and occupancies for a set of ion
  identities at a site.

  Unimplemented, currently pending on other improvements to phenix's
  refinement code.

  Parameters
  ----------
  scatter_env: mmtbx.ions.environment.Environment
      The object to extract features from.
  elements: list of str, optional
      A list of elements to compare the experimental b-factors against. If
      unset, takes the list from mmtbx.ions.ALLOWED_IONS.

  Returns
  -------
  numpy.array of float
      A vector containing quantitative properties for classification.
  """

  if elements is None:
    elements = ["ZN"]

  raise NotImplementedError("b_factor_occ_vector not yet implemented")

def ion_geometry_vector(chem_env, geometry_names = None):
  """
  Creates a vector for a site's geometry. For each geometry in geometry_names
  the vector contains a 1 if that geometry is present at the site and 0
  otherwise.

  A single boolean was chosen after some trial and error with an SVM as
  differences in deviations < 15 degrees were not found to be significant in
  helping to diffentiate ions.

  Parameters
  ----------
  chem_env: mmtbx.ions.environment.ChemicalEnvironment
      A object containing information about the chemical environment at a site.
  geometry_names: list of str, optional
      A list of geometry names to check for. If unset, take names from
      mmtbx.ions.SUPPORTED_GEOMETRY_NAMES.

  Returns
  -------
  numpy.array of float
      A vector containing quantitative properties for classification.
  """
  if geometry_names is None:
    geometry_names = SUPPORTED_GEOMETRY_NAMES

  present_geometry_names = [i[0] for i in chem_env.geometry]
  return np.fromiter((i in present_geometry_names for i in geometry_names),
                     dtype = float)

def ion_nearby_atoms_vector(chem_env, environments = None):
  """
  Creates a vector for the identities of the ions surrounding a site. Returns
  a vector with a count of coordinating nitrogens, oxygens, sulfurs, and
  chlorines.

  Parameters
  ----------
  chem_env: mmtbx.ions.environment.ChemicalEnvironment
      A object containing information about the chemical environment at a site.
  environments: list of int, optional
      A list of environments to check for. If unset, takes values from
      mmtbx.ions.environment.N_SUPPORTED_ENVIRONMENTS.

  Returns
  -------
  numpy.array of float
      A vector containing quantitative properties for classification.
  """

  if environments is None:
    environments = range(N_SUPPORTED_ENVIRONMENTS)

  return np.fromiter((chem_env.chemistry[i] for i in environments), dtype = float)

def ion_valence_vector(chem_env, elements = None, server = None):
  """
  Calculate the BVS and VECSUM values for a variety of ion identities.

  Parameters
  ----------
  chem_env: mmtbx.ions.environment.ChemicalEnvironment
      A object containing information about the chemical environment at a site.
  elements: list of str, optional
      A list of elements to compare the experimental BVS and VECSUM
      against. If unset, takes the list from mmtbx.ions.ALLOWED_IONS.
  server: mmtbx.ions.parameters.server, optional
      A server used to find charges of elements.

  Returns
  -------
  numpy.array of float
      A vector containing quantitative properties for classification.
  """

  if server is None:
    global METAL_SERVER
    if METAL_SERVER is None:
      METAL_SERVER = ions_server()
    server = METAL_SERVER
  if elements is None:
    elements = [i for i in ALLOWED_IONS if i not in ions.WATER_RES_NAMES]

  ret = []

  for element in elements:
    ret.append(chem_env.get_valence(
      element = element,
      charge = server.get_charge(element)
      ))

  # Flatten the list
  return _flatten_list(ret)

def ion_anomalous_vector(scatter_env, elements = None, ratios = True):
  """
  Calculate the f'' / f''_expected for a variety of ion identities.

  Parameters
  ----------
  scatter_env: mmtbx.ions.environment.ScatteringEnvironment
      An object containing information about the scattering environment at a
      site.
  elements: list of str, optional
      A list of elements to compare the experimental f'' against. If unset,
      takes the list from mmtbx.ions.ALLOWED_IONS.
  ratios: bool, optional
      If False, instead of calculating ratios, just return a vector of the
      wavelength, f', and f''.

  Returns
  -------
  numpy.array of float
      A vector containing quantitative properties for classification.
  """

  if elements is None:
    elements = [i for i in ALLOWED_IONS if i not in ions.WATER_RES_NAMES]

  if scatter_env.fpp is None or scatter_env.wavelength is None:
    if ratios:
      return np.zeros(len(elements))
    else:
      return np.zeros(3)

  if ratios:
    ret = np.fromiter(
      (scatter_env.fpp /
       sasaki.table(element).at_angstrom(scatter_env.wavelength).fdp()
       for element in elements),
       float)
  else:
    ret = _flatten_list([
      scatter_env.wavelength, scatter_env.fpp, scatter_env.fp
      ])
  return ret

def _flatten_list(lst):
  """
  Turn a tree main out of lists into one flat numpy array. Converts all Nones
  into zeros and integers into floats in the process.

  Parameters
  ----------
  lst: list or list of list or list of list of list or ...
      A list to be flattened

  Returns
  -------
  numpy.array of float
      A flat list of values.
  """

  def _flatten(lst):
    """
    Returns a generator for each element in the flattened version of lst.
    """

    for item in lst:
      if isinstance(item, Iterable) and not isinstance(item, basestring):
        for sub in _flatten(item):
          yield sub
      else:
        yield item

  return np.fromiter(
    (float(i) if i is not None else 0. for i in _flatten(lst)),
    dtype = float
    )
