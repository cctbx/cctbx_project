from __future__ import division

from boost_adaptbx import graph

import math
import operator

class Atom(object):
  """
  Match results by identity
  """

  def __init__(self, **kwargs):

    for ( name, value ) in kwargs.items():
      setattr( self, name, value )


class Molecule(object):
  """
  A graph-based 3D-description of a molecule
  """

  def __init__(self, vertex_type = "vector", edge_type = "set"):

    self.graph = graph.adjacency_list(
      graph_type = "undirected",
      vertex_type = vertex_type,
      edge_type = edge_type,
      )
    self.atom_for = {}


  @property
  def descriptor_for(self):

    return dict( zip( self.atom_for.values(), self.atom_for.keys() ) )


  @property
  def atoms(self):

    return self.atom_for.values()


  def add(self, atom):

    descriptor = self.graph.add_vertex( label = atom )
    ( x2, y2, z2 ) = atom.xyz
    assert descriptor not in self.atom_for

    for ( d, a ) in self.atom_for.items():
      ( x1, y1, z1 ) = a.xyz
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


class Compound(object):
  """
  A graph-based 2D-representation of a compound
  """

  def __init__(self, graph):

    self.graph = graph
    self.descriptor_for = dict(
      ( self.graph.vertex_label( vertex = v ), v ) for v in self.graph.vertices()
      )


  @property
  def atom_for(self):

    return dict( zip( self.descriptor_for.values(), self.descriptor_for.keys() ) )


  @property
  def atoms(self):

    return self.descriptor_for.keys()


  @property
  def descriptors(self):

    return self.descriptor_for.values()


  def add_atom(self, atom):

    self.descriptor_for[ atom ] = self.graph.add_vertex( label = atom )


  def add_bond(self, left, right, order = 1):

    if left not in self.descriptor_for:
      self.add_atom( atom = left )

    if right not in self.descriptor_for:
      self.add_atom( atom = right )

    self.graph.add_edge(
      vertex1 = self.descriptor_for[ left ],
      vertex2 = self.descriptor_for[ right ],
      weight = order,
      )


  def distances_from(self, atom):

    if atom not in self.descriptor_for:
      raise ValueError, "Unknown atom: %s" % atom

    from graph import breadth_first_search as bfs

    vertex = self.descriptor_for[ atom ]
    visitor = bfs.distance_recording_visitor( start_vertex = vertex )
    bfs.breadth_first_search(
      graph = self.graph,
      vertex = vertex,
      visitor = visitor,
      )

    atom_for = self.atom_for
    return dict(
      ( atom_for[ v ], visitor.distance_for.get( v, None ) )
      for v in self.descriptors
      )


  def subset(self, atoms):

    if not all( a in self.descriptor_for for a in atoms ):
      raise ValueError, "Unknown atoms: %s" % atoms

    from graph import maximum_clique
    subgraph = maximum_clique.selected_subgraph(
      graph = self.graph,
      vertices = ( self.descriptor_for[ a ] for a in atoms ),
      )

    return self.__class__( graph = subgraph )


  def connected_segments(self):

    from graph import connected_component_algorithm as cca
    res = cca.connected_components( graph = self.graph )
    atom_for = self.atom_for

    return [ [ atom_for[ v ] for v in comp ] for comp in res ]


  @classmethod
  def create(cls, vertex_type = "vector", edge_type = "set"):

    return cls(
      graph = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = vertex_type,
        edge_type = edge_type,
        ),
      )


class RascalMatch(object):
  """
  Geometry and label-based correspondence matching using the Rascal algorithm
  """

  def __init__(
    self,
    molecule1,
    molecule2,
    prematched = [],
    vertex_equality = operator.eq,
    edge_equality = operator.eq
    ):

    self.molecule1 = molecule1.atom_for
    self.molecule2 = molecule2.atom_for
    self.best = [ [] ]
    self.largest = 0

    from graph import maximum_clique

    self.compat_graph = maximum_clique.compatibility_graph(
      first = molecule1.graph,
      second = molecule2.graph,
      vertex_equality = vertex_equality,
      edge_equality = edge_equality
      )

    self.prematched = prematched

    if self.prematched:
      desc_for_1 = molecule1.descriptor_for
      desc_for_2 = molecule2.descriptor_for
      pairs = set(
        [ ( desc_for_1[ l ], desc_for_2[ r ] ) for ( l, r ) in self.prematched ]
        )
      vertices = [ v for v in self.compat_graph.vertices()
        if self.compat_graph.vertex_label( vertex = v ) in pairs ]

      assert len( vertices ) == len( self.prematched )
      neighbours = set( self.compat_graph.adjacent_vertices( vertex = vertices[0] ) )

      for v in vertices[1:]:
        neighbours.intersection_update( self.compat_graph.adjacent_vertices( vertex = v ) )

      self.compat_graph = maximum_clique.selected_subgraph(
        graph = self.compat_graph,
        vertices = neighbours,
        )

    maximum_clique.rascal( graph = self.compat_graph, callable = self )


  def count(self):

    return len( self.best )


  def length(self):

    return len( self.prematched ) + self.largest


  def remapped(self):

    matches = [
      [ self.compat_graph.vertex_label( vertex = v ) for v in match ]
      for match in self.best
      ]
    return [
      self.prematched + [
        ( self.molecule1[ l ], self.molecule2[ r ] ) for ( l, r ) in pairs
        ]
      for pairs in matches
      ]


  def __call__(self, result):

    size = len( result )

    if self.largest < size:
      self.largest = size
      self.best = [ result ]

    elif self.largest == size:
      self.best.append( result )


class McGregorMatch(object):
  """
  Geometry and label-based correspondence matching using the McGregor algorithm
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


  def length(self):

    return len( self.best )


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

