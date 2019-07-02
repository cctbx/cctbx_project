from __future__ import absolute_import, division, print_function

class Graph(object):
  """
  Basic graph implementation

  Graph interface:
      * Properties:
        - vertices: an iterable that support constant time membership tests
        - edges: an iterable returning a unique list of edges

      * Methods:
        - add_edge: insert an edge, and the vertices it connects
        - add vertex: insert a vertex
        - edges_from: return edges out from a vertex
  """

  def __init__(self):

    self.vertices = set()
    self.edges = set()


  def add_edge(self, edge):

    self.vertices.update( edge.vertices )
    self.edges.add( edge )


  def add_vertex(self, vertex):

    self.vertices.add( vertex )


  def edges_from(self, vertex):

    return ( e for e in self.edges if vertex in e.vertices )


class VertexIndexedGraph(object):
  """
  Represents a graph, with a quick way to get all edges from a vertex
  """

  def __init__(self):

    self._edges_from = dict()


  @property
  def vertices(self):

    return self._edges_from


  @property
  def edges(self):

    result = set()
    result.update( *list(self._edges_from.values()) )
    return result


  def add_edge(self, edge):

    for v in edge.vertices:
      self._edges_from.setdefault( v, set() ).add( edge )


  def add_vertex(self, vertex):

    self._edges_from.setdefault( vertex, set() )


  def edges_from(self, vertex):

    return self._edges_from.get( vertex, set() )


class Edge(object):
  """
  Edge in an undirected graph
  """

  def __init__(self, vertex1, vertex2):

    self.vertices = frozenset( ( vertex1, vertex2 ) )


  def other_vertex(self, current):

    assert current in self.vertices

    if len( self.vertices ) == 2:
      diff = self.vertices.difference( [ current ] )
      assert len( diff ) == 1
      return list( diff )[0]

    else:
      return current


  def __eq__(self, other):

    return self.vertices == other.vertices


  def __ne__(self, other):

    return not ( self == other )


  def __hash__(self):

    return hash( self.vertices )


def preorder_depth_first_iteration(vertex, graph):
  """
  Preorder depth-first traversal of a graph
  """

  assert vertex in graph.vertices
  yield vertex

  discovereds = set( [ vertex ] )
  stack = [ vertex ]

  while stack:
    current = stack[-1]

    for edge in graph.edges_from( vertex = current ):
      other = edge.other_vertex( current = current )

      if other in discovereds:
        continue

      yield other
      discovereds.add( other )
      stack.append( other )
      break

    else:
      stack.pop()


def connected_components(graph):

  unused = set( graph.vertices )

  while unused:
    vertex = unused.pop()
    component = list( preorder_depth_first_iteration( vertex = vertex, graph = graph ) )
    unused.difference_update( component )
    yield component

