from __future__ import division

import mmtbx.geometry.primitive # import dependency
import mmtbx.geometry.shared_types # import dependency

import boost.python
ext = boost.python.import_ext( "mmtbx_geometry_asa_ext" )
from mmtbx_geometry_asa_ext import *

# Atom representation
def Regular(sphere, processor):
  """
  An object without an altloc identifier
  """

  processor.process_regular( sphere = sphere )


class Altloc(object):
  """
  An object assigned with an altloc identifier
  """

  def __init__(self, identifier):

    self.identifier = identifier


  def __call__(self, sphere, processor):

    processor.process_altloc( sphere = sphere, identifier = self.identifier )


class Description(object):
  """
  An internal format for fast calculation
  """

  def __init__(self, sphere, strategy):

    self.sphere = sphere
    self.strategy = strategy


  def accept(self, processor):

    self.strategy( self.sphere, processor = processor )


  @classmethod
  def from_parameters(cls, centre, radius, index, strategy):

    return cls(
      sphere = sphere( centre = centre, radius = radius, index = index ),
      strategy = strategy,
      )


# Indexing with altloc support
class Indexer(object):
  """
  Indexer that takes into account altloc
  """

  def __init__(self, factory):

    self.factory = factory
    self.regular = self.factory()

    self.altlocs = {}


  def add(self, altloc):

    self.altlocs[ altloc ] = self.factory()


  @classmethod
  def create(cls, factory, descriptions):

    indexer = cls( factory = factory )
    inserter = Inserter( indexer = indexer )

    for d in descriptions:
      d.accept( processor = inserter )

    return indexer


def get_linear_indexer_for(descriptions):

  return Indexer.create(
    factory = indexing.linear_spheres,
    descriptions = descriptions,
    )


def get_hash_indexer_for(descriptions):

  voxelizer = get_voxelizer_for( descriptions = descriptions )
  return Indexer.create(
    factory = lambda: indexing.hash_spheres( voxelizer = voxelizer ),
    descriptions = descriptions,
    )


def get_optimal_indexer_for(descriptions):

  return get_linear_indexer_for( descriptions = descriptions )


def get_voxelizer_for(descriptions, step = 7):

  lows = [ d.sphere.low for d in descriptions ]
  ( low_xs, low_ys, low_zs ) = zip( *lows )
  low = ( min( low_xs ), min( low_ys ), min( low_zs ) )

  return mmtbx.geometry.shared_types.voxelizer(
    base = low,
    step = ( step, step, step ),
    )


# Visitors
class Inserter(object):
  """
  Fills up an indexer with spheres
  """

  def __init__(self, indexer):

    self.indexer = indexer


  def process_regular(self, sphere):

    self.indexer.regular.add( object = sphere )


  def process_altloc(self, sphere, identifier):

    if identifier not in self.indexer.altlocs:
      self.indexer.add( altloc = identifier )

    self.indexer.altlocs[ identifier ].add( object = sphere )


class CompositeCheckerBuilder(object):
  """
  Finds neighbours for a sphere
  """

  def __init__(self, indexer, description):

    self.indexer = indexer
    self.checker = accessibility.pythagorean_checker()
    description.accept( processor = self )


  def process_regular(self, sphere):

    self.append_neighbours( indexer = self.indexer.regular, sphere = sphere )

    for indexer in self.indexer.altlocs.values():
      self.append_neighbours( indexer = indexer, sphere = sphere )


  def process_altloc(self, sphere, identifier):

    self.append_neighbours( indexer = self.indexer.regular, sphere = sphere )
    self.append_neighbours(
      indexer = self.indexer.altlocs[ identifier ],
      sphere = sphere,
      )


  def append_neighbours(self, indexer, sphere):

    self.checker.add(
      neighbours = accessibility.filter(
        range = indexer.close_to( object = sphere ),
        predicate = accessibility.overlap_equality_predicate( object = sphere )
        )
      )


class SeparateCheckerBuilder(object):
  """
  Finds neighbours for a sphere
  """

  def __init__(self, indexer, description):

    self.indexer = indexer
    self.regular = accessibility.pythagorean_checker()
    self.altlocs = {}
    description.accept( processor = self )


  def process_regular(self, sphere):

    self.append_neighbours(
      indexer = self.indexer.regular,
      sphere = sphere,
      checker = self.regular,
      )

    for ( identifier, indexer ) in self.indexer.altlocs.items():
      checker = accessibility.pythagorean_checker()
      self.append_neighbours(
        indexer = indexer,
        sphere = sphere,
        checker = checker,
        )

      if checker.neighbours():
        self.altlocs[ identifier ] = checker


  def process_altloc(self, sphere, identifier):

    self.append_neighbours(
      indexer = self.indexer.regular,
      sphere = sphere,
      checker = self.regular,
      )
    self.altlocs[ identifier ] = accessibility.pythagorean_checker()
    self.append_neighbours(
      indexer = self.indexer.altlocs[ identifier ],
      sphere = sphere,
      checker = self.altlocs[ identifier ],
      )


  @staticmethod
  def append_neighbours(indexer, sphere, checker):

    checker.add(
      neighbours = accessibility.filter(
        range = indexer.close_to( object = sphere ),
        predicate = accessibility.overlap_equality_predicate( object = sphere )
        )
      )


# Results
class AccessibleSurfaceResult(object):
  """
  Result of the calculation
  """

  def __init__(self, count, radius_sq):

    self.count = count
    self.radius_sq = radius_sq


  @property
  def surface(self):

    return self.count * self.radius_sq


class AccessibleSurfaceAreas(object):
  """
  Result of a series of calculations
  """

  def __init__(self, values, unit):

    self.values = values
    self.unit = unit


  @property
  def points(self):

    return [ v.count for v in self.values ]


  @property
  def areas(self):

    return [ self.unit * v.surface for v in self.values ]


# Ways for calculating ASA
def simple_surface_calculation(indexer, sampling, description):
  """
  Calculates ASA by not worrying about altlocs
  """

  builder = CompositeCheckerBuilder( indexer = indexer, description = description )
  overlapped = accessibility.filter(
    range = accessibility.transform(
      range = sampling.points,
      transformation = accessibility.transformation(
        centre = description.sphere.centre,
        radius = description.sphere.radius,
        ),
      ),
    predicate = builder.checker,
    )

  return AccessibleSurfaceResult(
    count = len( overlapped ),
    radius_sq = description.sphere.radius_sq,
    )


def altloc_averaged_calculation(indexer, sampling, description):
  """
  For atoms with altloc identifier, use empty altloc + atom with same altloc.
  For atoms with empty altloc, run a calculation for each know altloc and
  average the results.
  """

  builder = SeparateCheckerBuilder( indexer = indexer, description = description )

  overlapped = accessibility.filter(
    range = accessibility.transform(
      range = sampling.points,
      transformation = accessibility.transformation(
        centre = description.sphere.centre,
        radius = description.sphere.radius,
        ),
      ),
    predicate = builder.regular,
    )

  accessible_for = {}

  for ( identifier, checker ) in builder.altlocs.items():
    accessible_for[ identifier ] = [
      p for p in overlapped if checker( point = p )
      ]

  if not accessible_for:
    count = len( overlapped )

  else:
    count = sum( [ len( l ) for l in accessible_for.values() ] ) / len( accessible_for )

  return AccessibleSurfaceResult(
    count = count,
    radius_sq = description.sphere.radius_sq,
    )


# Module-level function
def calculate(
  atoms,
  calculation = simple_surface_calculation,
  indexer_selector = get_optimal_indexer_for,
  probe = 1.4,
  precision = 960,
  ):

  from cctbx.eltbx import van_der_waals_radii
  radius_for = van_der_waals_radii.vdw.table

  from mmtbx.geometry import altloc
  descriptions = [
    Description.from_parameters(
      centre = a.xyz,
      radius = radius_for[ a.determine_chemical_element_simple().strip().capitalize() ] + probe,
      index = i,
      strategy = altloc.from_atom( entity = a ),
      )
    for ( i, a ) in enumerate( atoms )
    ]

  indexer = indexer_selector( descriptions = descriptions )

  from mmtbx.geometry import sphere_surface_sampling
  sampling = sphere_surface_sampling.golden_spiral( count = precision )

  values = [
    calculation( indexer = indexer, sampling = sampling, description = d )
    for d in descriptions
    ]

  return AccessibleSurfaceAreas(
    values = values,
    unit = sampling.unit_area,
    )

