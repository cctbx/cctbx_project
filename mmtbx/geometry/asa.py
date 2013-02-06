from __future__ import division

import mmtbx.geometry.primitive

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

    self.strategy( sphere = self.sphere, processor = processor )


  @classmethod
  def from_parameters(cls, centre, radius, altloc):
        
    return cls(
      sphere = sphere.create( centre = centre, radius = radius ),
      strategy = Regular if not altloc else Altloc( identifier = altloc ),
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
    

def get_linear_indexer_for(atoms):

  return Indexer( factory = indexing.linear_spheres )


def get_optimal_indexer_for(atoms):

  return Indexer( factory = indexing.linear_spheres )


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
  
  def __init__(self, indexer):
    
    self.indexer = indexer
    self.checker = accessibility.pythagorean_checker()
    

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
    
    self.checker.add_from_range(
      neighbours = indexer.overlapping_with( object = sphere ),
      )
    
    
class SeparateCheckerBuilder(object):
  """
  Finds neighbours for a sphere
  """
  
  def __init__(self, indexer):
    
    self.indexer = indexer
    self.regular = accessibility.pythagorean_checker()
    self.altlocs = {}
    

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
      
      if list( checker.neighbours() ):
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
    
    checker.add_from_range(
      neighbours = indexer.overlapping_with( object = sphere ),
      )


# Parameters
class CalcParams(object):
  """
  Parameters of the calculation
  """

  def __init__(self, probe = 1.4, precision = 960):

    self.probe = probe

    from mmtbx.geometry import sphere_surface_sampling
    self.sampling = sphere_surface_sampling.golden_spiral( count = precision )


  @property
  def van_der_waals_radii(self):

    from cctbx.eltbx import van_der_waals_radii
    return van_der_waals_radii.vdw.table


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
class SimpleSurfaceCalculator(object):
  """
  Calculates ASA by not worrying about altlocs
  """

  def __init__(self, indexer, params):

    self.indexer = indexer
    self.sampling = params.sampling
    self.values = []


  def process(self, description):

    builder = CompositeCheckerBuilder( indexer = self.indexer )
    description.accept( processor = builder )
    overlapped = builder.checker.filter(
      range = self.sampling.transformed(
        centre = description.sphere.centre,
        radius = description.sphere.radius,
        )
      )

    self.values.append(
      AccessibleSurfaceResult(
        count = len( overlapped ),
        radius_sq = description.sphere.radius_sq,
        )
      )
    
    
class AltlocAveragedCalculator(object):
  """
  For atoms with altloc identifier, use empty altloc + atom with same altloc.
  For atoms with empty altloc, run a calculation for each know altloc and
  average the results.
  """
  
  def __init__(self, indexer, params):
    
    self.indexer = indexer
    self.sampling = params.sampling
    self.values = []
    
    
  def process(self, description):
    
    builder = SeparateCheckerBuilder( indexer = self.indexer )
    description.accept( processor = builder )
    
    overlapped = builder.regular.filter(
      range = self.sampling.transformed(
        centre = description.sphere.centre,
        radius = description.sphere.radius,
        )
      )
    
    accessible_for = {}
    
    for ( identifier, checker ) in builder.altlocs.items():
      accessible_for[ identifier ] = [
        p for p in overlapped if checker.is_selected( point = p )
        ]
    
    if not accessible_for:
      count = len( overlapped )
    
    else:
      count = sum( [ len( l ) for l in accessible_for.values() ] ) / len( accessible_for )

    self.values.append(
      AccessibleSurfaceResult(
        count = count,
        radius_sq = description.sphere.radius_sq,
        )
      )


# Module-level function
def calculate(
  atoms,
  calculator = SimpleSurfaceCalculator,
  indexer_selector = get_optimal_indexer_for,
  params = CalcParams(),
  ):
  
  centres = [ a.xyz for a in atoms ]
  
  radius_for = params.van_der_waals_radii
  radii = [
    radius_for[ a.determine_chemical_element_simple().strip().capitalize() ] + params.probe
    for a in atoms
    ]
  altlocs = [ a.parent().altloc if a.parent() else None for a in atoms ]
  
  descriptions = [
    Description.from_parameters( centre = c, radius = r, altloc = a )
    for ( c, r, a ) in zip( centres, radii, altlocs )
    ]

  indexer = indexer_selector( atoms = atoms )
  inserter = Inserter( indexer = indexer )

  for d in descriptions:
    d.accept( processor = inserter )

  calculator = calculator( indexer = indexer, params = params )

  for d in descriptions:
    calculator.process( description = d )

  return AccessibleSurfaceAreas(
    values = calculator.values,
    unit = params.sampling.unit_area,
    )
