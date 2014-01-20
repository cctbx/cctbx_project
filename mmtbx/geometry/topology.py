from __future__ import division

from boost_adaptbx import graph

import math
import operator

class Atom(object):
  """
  Match results by identity
  """

  def __init__(self, label, coordinates):

    self.label = label
    self.coordinates = coordinates


class Molecule(object):
  """
  A graph-based description of a molecule
  """

  def __init__(self):

    self.graph = graph.adjacency_list()
    self.atom_for = {}


  def add(self, atom):

    descriptor = self.graph.add_vertex( label = atom.label )
    ( x2, y2, z2 ) = atom.coordinates
    assert descriptor not in self.atom_for

    for ( d, a ) in self.atom_for.items():
      ( x1, y1, z1 ) = a.coordinates
      self.graph.add_edge(
        vertex1 = d,
        vertex2 = descriptor,
        weight = math.sqrt(
          ( x1 - x2 ) ** 2  + ( y1 - y2 ) ** 2 + ( z1 - z2 ) ** 2
          ),
        )

    self.atom_for[ descriptor ] = atom


  def size(self):

    return len( self.atom_for )


class Match(object):
  """
  Geometry and label-based correspondence matching
  """

  def __init__(
    self,
    molecule1,
    molecule2,
    is_valid,
    maxsteps = 500,
    vertex_equality = operator.eq,
    edge_equality = operator.eq
    ):

    self.molecule1 = molecule1.atom_for
    self.molecule2 = molecule2.atom_for
    self.is_valid = is_valid
    self.best = []
    self.steps = 0
    self.maxsteps = maxsteps
    
    from graph import graph_structure_comparison as gsc
    gsc.mcgregor_common_subgraphs_unique(
      graph1 = molecule1.graph,
      graph2 = molecule2.graph,
      vertex_equality = vertex_equality,
      edge_equality = edge_equality,
      callback = self,
      )


  def remapped(self):

    return [
      ( self.molecule1[ l ], self.molecule2[ r ] ) for ( l, r ) in self.best
      ]


  def __call__(self, match):

    if len( match ) <= len( self.best ):
      self.steps += 1

    elif ( self.is_valid( ( self.molecule1[ m[0] ] for m in match ) )
      and self.is_valid( ( self.molecule2[ m[1] ] for m in match ) ) ):
      self.best = match
      self.steps = 0

    else:
      self.steps += 1

    return self.steps <= self.maxsteps


def is_valid_sidechain(molecule):

  return any( m.label == "CA" for m in molecule )


def sidechain_match(molecule1, molecule2, tolerance = 0.1):

  m = Match(
    molecule1 = molecule1,
    molecule2 = molecule2,
    is_valid = is_valid_sidechain,
    edge_equality = lambda l, r: abs( l - r ) <= tolerance,
    )

  return m.remapped()

