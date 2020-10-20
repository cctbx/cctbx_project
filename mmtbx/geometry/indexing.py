from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
ext = bp.import_ext( "mmtbx_geometry_indexing_ext" )
from mmtbx_geometry_indexing_ext import *

def filter(range, predicate):

  return ext.filter(
    range = range,
    predicate = code_predicate( callable = predicate ),
    )


def direct(obj):

  return obj


class overlap_predicate(object):
  """
  Simple predicate to evaluate actual distance
  """

  def __init__(self, centre, distance, getter = direct):

    ( self.x, self.y, self.z ) = centre
    self.distance_sq = distance ** 2
    self.getter = getter


  def __call__(self, obj):

    ( x, y, z ) = self.getter( obj )
    dist_sq = ( self.x - x ) ** 2 + ( self.y - y ) ** 2 + ( self.z - z ) ** 2
    return dist_sq <= self.distance_sq


class containment_predicate(object):
  """
  Simple predicate to evaluate whether the object is in a predefined set
  """

  def __init__(self, selection, getter, relation):

    self.selection = selection
    self.getter = getter
    self.relation = relation


  def __call__(self, obj):

    return self.relation( self.getter( obj ), self.selection )


  @staticmethod
  def include(obj, selection):

    return obj in selection


  @staticmethod
  def exclude(obj, selection):

    return obj not in selection


  @classmethod
  def inclusion(cls, selection, getter):

    return cls( selection = selection, getter = getter, relation = cls.include )


  @classmethod
  def exclusion(cls, selection, getter):

    return cls( selection = selection, getter = getter, relation = cls.exclude )


class identity_predicate(object):
  """
  Simple predicate to evaluate whether the object is identical to something
  """

  def __init__(self, element, getter):

    self.element = element
    self.getter = getter


  def __call__(self, obj):

    return self.element == self.element_of( obj )


class composite_predicate(object):
  """
  Layer several conditions
  """

  def __init__(self, conditions):

    self.conditions = conditions


  def __call__(self, obj):

    return all( cond( obj ) for cond in self.conditions )


class structure_indexer(object):
  """
  Indexes iotbx.pdb.hierarchy objects
  """

  def __init__(self, indexer, datagen, maxdist):

    self.indexer = indexer
    self.datagen = datagen
    self.maxdist = maxdist


  def add_atoms_in(self, obj):

    for atom in obj.atoms():
      xyz = atom.xyz
      self.indexer.add(
        object = ( self.datagen( obj, atom ), xyz ),
        position = xyz,
        )


  def interaction_ranges_with(self, obj, predicates = []):

    if not predicates:
      predgen = self.simple_distance_predicate

    else:
      predgen = self.composite_distance_predicate

    for atom in obj.atoms():
      centre = atom.xyz
      yield list(filter(
        range = self.indexer.close_to( centre = centre ),
        predicate = predgen( predicates = predicates, centre = centre ),
        ))


  def interaction_counts_with(self, obj, predicates = []):

    rangiter = self.interaction_ranges_with( obj = obj, predicates = predicates )

    return sum( len( r ) for r in rangiter )


  def interactions_with(self, obj, predicates = []):

    rangiter = self.interaction_ranges_with( obj = obj, predicates = predicates )

    for r in rangiter:
      for i in r:
        yield i


  def overlap_predicate(self, centre):

    return overlap_predicate(
      centre = centre,
      distance = self.maxdist,
      getter = self.coordinate_of,
      )


  @classmethod
  def include_predicate(cls, selection):

    return containment_predicate.Inclusion(
      selection = selection,
      getter = cls.data_of,
      )


  @classmethod
  def exclude_predicate(cls, selection):

    return containment_predicate.Exclusion(
      selection = selection,
      getter = cls.data_of,
      )


  @classmethod
  def identity_predicate(cls, element):

    return identity_predicate( element = element, getter = cls.data_of )


  @staticmethod
  def data_of(obj):

    return obj[0]


  @staticmethod
  def coordinate_of(obj):

    return obj[1]


  def composite_distance_predicate(self, predicates, centre):

    return composite_predicate(
      conditions = predicates + [ self.overlap_predicate( centre = centre ) ],
      )


  def simple_distance_predicate(self, predicates, centre):

    return self.overlap_predicate( centre = centre )


  @classmethod
  def hash(
    cls,
    base,
    maxdist,
    datagen = lambda o, a: None,
    subdivision = 1,
    safescale = 1.01,
    ):

    boxsize = safescale / subdivision * maxdist

    from mmtbx.geometry import shared_types
    voxelizer = shared_types.voxelizer(
      base = base,
      step = ( boxsize, boxsize, boxsize ),
      )

    return cls(
      indexer = hash( voxelizer = voxelizer, margin = subdivision ),
      datagen = datagen,
      maxdist = maxdist,
      )
