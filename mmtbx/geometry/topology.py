from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph

import math
from six.moves import zip

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
    self.xyz_for = {}


  @property
  def descriptor_for(self):
    # TODO: six.moves.zip this file
    return dict( zip( self.atom_for.values(), self.atom_for.keys() ))


  @property
  def atoms(self):

    return list(self.atom_for.values())


  def add(self, atom, xyz):

    d2 = self.graph.add_vertex( label = atom )
    ( x2, y2, z2 ) = xyz
    assert d2 not in self.xyz_for

    for ( d1, ( x1, y1, z1 ) ) in self.xyz_for.items():
      self.graph.add_edge(
        vertex1 = d1,
        vertex2 = d2,
        weight = math.sqrt(
          ( x1 - x2 ) ** 2  + ( y1 - y2 ) ** 2 + ( z1 - z2 ) ** 2
          ),
        )

    self.atom_for[ d2 ] = atom
    self.xyz_for[ d2 ] = xyz


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

    return dict( zip( self.descriptor_for.values(), self.descriptor_for.keys() ))


  @property
  def atoms(self):

    return list(self.descriptor_for.keys())


  @property
  def descriptors(self):

    return list(self.descriptor_for.values())


  @property
  def bonds(self):

    atom_for = self.atom_for

    for edge in self.graph.edges():
      yield (
        atom_for[ self.graph.source( edge = edge) ],
        atom_for[ self.graph.target( edge = edge) ],
        )


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
      raise ValueError("Unknown atom: %s" % atom)

    from boost_adaptbx.graph import breadth_first_search as bfs

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
      raise ValueError("Unknown atoms: %s" % atoms)

    from boost_adaptbx.graph import maximum_clique
    subgraph = maximum_clique.selected_subgraph(
      graph = self.graph,
      vertices = ( self.descriptor_for[ a ] for a in atoms ),
      )

    return self.__class__( graph = subgraph )


  def connected_segments(self):

    from boost_adaptbx.graph import connected_component_algorithm as cca
    res = cca.connected_components( graph = self.graph )
    atom_for = self.atom_for

    return [ [ atom_for[ v ] for v in comp ] for comp in res ]


  def connected_segment_from(self, atom):

    if atom not in self.descriptor_for:
      raise ValueError("Unknown atom: %s" % atom)

    from boost_adaptbx.graph import breadth_first_search as bfs

    vertex = self.descriptor_for[ atom ]
    visitor = bfs.vertex_recording_visitor( start_vertex = vertex )
    bfs.breadth_first_search(
      graph = self.graph,
      vertex = vertex,
      visitor = visitor,
      )

    atom_for = self.atom_for

    return [ atom_for[ v ] for v in visitor.visited_vertices ]


  @classmethod
  def create(cls, vertex_type = "vector", edge_type = "set"):

    return cls(
      graph = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = vertex_type,
        edge_type = edge_type,
        ),
      )


  @classmethod
  def from_structure(cls, atoms, tolerance = 0.1, vertex_type = "vector", edge_type = "set"):

    compound = cls.create( vertex_type = vertex_type, edge_type = edge_type )

    for a in atoms:
      compound.add_atom( atom = a )

    import itertools
    from mmtbx.monomer_library import bondlength_defaults

    for ( a1, a2 ) in itertools.combinations( atoms, 2 ):
      ( x1, y1, z1 ) = a1.xyz
      ( x2, y2, z2 ) = a2.xyz

      dist2 = ( ( x1 - x2 ) ** 2 + ( y1 - y2 ) ** 2 + ( z1 - z2 ) ** 2 )
      expected = bondlength_defaults.get_default_bondlength( a1.element, a2.element )

      if dist2 < ( expected + tolerance ) ** 2:
        compound.add_bond( a1, a2 )

    return compound


class RascalMatch(object):
  """
  Geometry and label-based correspondence matching using the Rascal algorithm
  """

  def __init__(
    self,
    molecule1,
    molecule2,
    prematched = [],
    vertex_equality = None,
    edge_equality = None,
    ):

    if vertex_equality is None:
      import operator
      vertex_equality = operator.eq

    if edge_equality is None:
      import operator
      edge_equality = operator.eq

    self.molecule1 = molecule1.atom_for
    self.molecule2 = molecule2.atom_for
    self.best = [ [] ]
    self.largest = 0

    from boost_adaptbx.graph import maximum_clique

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
    vertex_equality = None,
    edge_equality = None,
    ):

    if vertex_equality is None:
      import operator
      vertex_equality = operator.eq

    if edge_equality is None:
      import operator
      edge_equality = operator.eq

    self.molecule1 = molecule1.atom_for
    self.molecule2 = molecule2.atom_for
    self.is_valid = is_valid
    self.best = []
    self.steps = 0
    self.maxsteps = maxsteps

    from boost_adaptbx.graph import graph_structure_comparison as gsc
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


class GreedyMatch(object):
  """
  Geometry and label-based correspondence matching using the Greedy algorithm
  """

  def __init__(
    self,
    molecule1,
    molecule2,
    prematched = [],
    vertex_equality = None,
    edge_equality = None,
    maxsol = 0,
    ):

    if vertex_equality is None:
      import operator
      vertex_equality = operator.eq

    if edge_equality is None:
      import operator
      edge_equality = operator.eq

    self.molecule1 = molecule1.atom_for
    self.molecule2 = molecule2.atom_for
    self.best = [ [] ]
    self.largest = 0

    from boost_adaptbx.graph import maximum_clique

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

    self.result = maximum_clique.greedy( graph = self.compat_graph, maxsol = maxsol )


  def count(self):

    return len( self.result )


  def length(self):

    assert self.result
    return len( self.prematched ) + len( self.result[0] )


  def remapped(self):

    matches = [
      [ self.compat_graph.vertex_label( vertex = v ) for v in match ]
      for match in self.result
      ]
    return [
      self.prematched + [
        ( self.molecule1[ l ], self.molecule2[ r ] ) for ( l, r ) in pairs
        ]
      for pairs in matches
      ]
