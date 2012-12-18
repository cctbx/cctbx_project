from __future__ import division

import mmtbx.geometry.primitive

import boost.python
ext = boost.python.import_ext( "mmtbx_geometry_asa_ext" )
from mmtbx_geometry_asa_ext import *

def get_linear_indexer_for(atoms):

  return Indexer( factory = indexing.linear_spheres )


def get_optimal_indexer_for(atoms):

  return Indexer( factory = indexing.linear_spheres )


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

    processor.process_atloc( sphere = sphere, identifier = self.identifier )


class Description(object):
  """
  An internal format for fast calculation
  """

  def __init__(self, sphere, strategy):

    self.sphere = sphere
    self.strategy = strategy


  def accept(self, processor):

    self.strategy( sphere = self.sphere, processor = processor )


  @classmethod
  def from_atom(cls, atom, index, params):

    radius = params.van_der_waals_radii[
      atom.determine_chemical_element_simple().strip().capitalize()
      ]

    altloc = atom.parent().altloc

    return cls(
      sphere = sphere( centre = atom.xyz, radius = radius + params.probe, index = index ),
      strategy = Regular if not altloc else Altloc( identifier = altloc ),
      )


class Indexer(object):
  """
  Indexer that takes into account altloc
  """

  def __init__(self, factory):

    self.factory = factory
    self.regular = self.factory()

    self.indexer_for = {}


  @property
  def altloc(self, identifier):

    return self.indexer_for[ identifier ]


  @property
  def altlocs(self):

    return self.indexer_for


  def add_altloc(self, identifier):

    self.indexer_for[ identifier ] = self.factory()


class Inserter(object):
  """
  Fills up an indexer with spheres
  """

  def __init__(self, indexer):

    self.indexer = indexer


  def process_regular(self, sphere):

    self.indexer.regular.add( object = sphere )


  def process_altloc(self, sphere, identifier):

    if identifier not in self.indexer.atlocs:
      self.indexer.add_altloc( identifier = identifier )

    self.inserter.altloc( identifier = identifier ).add( object = sphere )


class SimpleASACalculator(object):
  """
  Calculates ASA by not worrying about altlocs
  """

  def __init__(self, indexer, params):

    self.indexer = indexer
    self.sampling = params.sampling
    self.values = []


  def process(self, sphere, checker):

    overlapped = checker.filter(
      range = self.sampling.transformed(
        centre = sphere.centre,
        radius = sphere.radius
        )
      )

    self.values.append(
      self.sampling.unit_area * sphere.radius_sq * len( overlapped )
      )


  def process_regular(self, sphere):

    checker = containment.pythagorean_checker()
    checker.add_from_range(
      neighbours = self.indexer.regular.prefiltered_overlapping_with(
        object = sphere,
        prefilter = index_filter( index = sphere.index ),
        )
      )

    for identifier in self.indexer.altlocs:
      indexer = self.indexer.altloc( identifier = identifier )
      checker.add_from_range(
        neighbours = indexer.overlapping_with( object = sphere )
        )

    self.process( sphere = sphere, checker = checker )


  def process_altloc(self, sphere, identifier):

    checker = containment.pythagorean_checker()
    checker.add_from_range(
      neighbours = self.indexer.regular.overlapping_with( object = sphere )
      )

    indexer = self.indexer.altloc( identifier = identifier )
    checker.add_from_range(
      neighbours = indexer.prefiltered_overlapping_with(
        object = sphere,
        prefilter = index_filter( index = sphere.index ),
        )
      )

    self.process( sphere = sphere, checker = checker )


class CalcParams(object):
  """
  Parameters of the calculation
  """

  def __init__(
    self,
    probe = 1.4,
    precision = 960,
    indexer_for = get_optimal_indexer_for,
    calculator = SimpleASACalculator,
    ):

    self.probe = probe

    from mmtbx.geometry import sphere_surface_sampling
    self.sampling = sphere_surface_sampling.golden_spiral( count = precision )

    self.indexer_for = indexer_for
    self.calculator = calculator


  @property
  def van_der_waals_radii(self):

    from cctbx.eltbx import van_der_waals_radii
    return van_der_waals_radii.vdw.table


class Calculation(object):
  """
  The actual calculation for a single atom
  """

  def __init__(self, centre, neighbours, params):

    self.count = iterator_length(
      iterator = accessible_points(
        neighbours = neighbours,
        points = params.points.transformed(
          centre = centre.xyz,
          radius = centre.radius,
          ),
        )
      )
    self.centre = centre
    self.params = params


  @property
  def area(self):

    return self.params.unit * self.count * self.centre.radius_sq


class Averaged(object):
  """
  For atoms with altloc identifier, use empty altloc + atom with same altloc.
  For atoms with empty altloc, run a calculation for each know altloc and
  average the results.
  """

  @staticmethod
  def process_regular(centre, model, params):

    regulars = model.regular.overlapping( description = centre.geometric )

    alternatives = set( frozenset( neighbours ) for neighbours
      in model.disordereds( description = centre.geometric ) )

    results = []

    for neighbours in alternatives:
      area = Calculation(
        centre = centre.geometric,
        neighbours = ( model.regular.overlapping( description = centre.geometric )
          + disordered.overlapping( description = centre.geometric ) ),
        params = params,
        ) \
        .area
      results.append( area )

    assert results
    return sum( results ) / len( results )


  @staticmethod
  def process_disordered(centre, model, params):

    disordered = model.disordered( altloc = centre.altloc )

    return Calculation(
      centre = centre.geometric,
      neighbours = ( model.regular.overlapping( description = centre.geometric )
        + disordered.overlapping( description = centre.geometric ) ),
      params = params,
      ) \
      .area


def convert_to_spheres(atoms, params):

  table = params.van_der_waals_radii
  vdw_radii = [
    table[
      atom.determine_chemical_element_simple().strip().capitalize()
      ]
    for atom in atoms
    ]

  return [
    sphere( centre = atom.xyz, radius = radius + params.probe, index = index )
    for ( index, ( atom, radius ) ) in enumerate( zip( atoms, vdw_radii ) )
    ]


def single(indexer, sphere, sampling):

  checker = containment.pythagorean_checker()
  checker.add_from_range(
    neighbours = indexer.prefiltered_overlapping_with(
      object = sphere,
      prefilter = index_filter( index = sphere.index ) ),
    )
  overlapped = checker.filter(
    range = sampling.transformed( centre = sphere.centre, radius = sphere.radius )
    )

  return len( overlapped )


def calculate(atoms, params):

  descriptions = [
    Description.from_atom( atom = a, index = i, params = params )
    for ( i, a ) in enumerate( atoms )
    ]

  indexer = params.indexer_for( atoms = atoms )
  inserter = Inserter( indexer = indexer )

  for d in descriptions:
    d.accept( processor = inserter )

  calculator = params.calculator( indexer = indexer, params = params )

  for d in descriptions:
    d.accept( processor = calculator )

  return calculator.values

