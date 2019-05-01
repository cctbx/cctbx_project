from __future__ import absolute_import, division, print_function

from boost_adaptbx.graph import min_cut_max_flow
from boost_adaptbx.graph import utility

import operator
from functools import reduce


def no_heap_element_limit(heap):

  return False


class heap_element_limit(object):

  def __init__(self, maxcount):

    self.maxcount = maxcount
    self.iteration = 0


  def __call__(self, heap):

    if self.maxcount < len( heap ):
      heap.sort()

      while self.maxcount < len( heap ):
        heap.pop()

      import heapq
      heapq.heapify( heap )

    self.iteration += 1
    return self.maxcount <= self.iteration


def stirling_adaptive_tree_enumeration(
  graph,
  reference,
  sw_path_vertex_selector,
  bk_path_vertex_selector,
  maxiter = None,
  ):
  """
  Enumerates cuts in order of weight based on a tree built using the recursive
  formula for Stirling numbers of the second kind
  """

  assert maxiter is None or 0 < maxiter

  import heapq
  heap = []

  labelledgraph = construct_vertex_labelled_graph( graph = graph )
  current = stirling_tree_sw_node(
    graph = labelledgraph,
    source_vertex_label = frozenset( ( reference, ) ),
    sw_path_vertex_selector = sw_path_vertex_selector,
    bk_path_vertex_selector = bk_path_vertex_selector,
    )
  heapq.heappush( heap, ( current.weight, current ) )

  if maxiter is None:
    handler = no_heap_element_limit

  else:
    handler = heap_element_limit( maxcount = maxiter )

  while heap:
    ( weight, current ) = heapq.heappop( heap )
    node = current

    while not node.is_terminal_leaf():
      for child in node.immediate_offpath_children():
        heapq.heappush( heap, ( child.weight, child ) )

      node = node.immediate_onpath_child()

    yield (
      weight,
      reduce( operator.or_, current.source_set_labels(), set() ),
      reduce( operator.or_, current.sink_set_labels(), set() ),
      )

    if handler( heap = heap ):
      break


class stirling_tree_node(object):
  """
  An element in the Stirling tree
  """

  def is_terminal_leaf(self):

    return self.graph.num_vertices() == 2


  def available_source_set_vertices(self):

    remaining = self.source_set.copy()
    remaining.remove( self.source_vertex )
    return remaining


  def available_sink_set_vertices(self):

    remaining = self.sink_set.copy()
    remaining.remove( self.sink_vertex )
    return remaining


  def source_set_labels(self):

    for vertex in self.source_set:
      yield self.graph.vertex_label( vertex = vertex )


  def sink_set_labels(self):

    for vertex in self.sink_set:
      yield self.graph.vertex_label( vertex = vertex )


class stirling_tree_sw_node(stirling_tree_node):
  """
  A node in the Stirling-tree that uses Stoer-Wagner minimum cut

  Branches 3-way
  """

  def __init__(
    self,
    graph,
    source_vertex_label,
    sw_path_vertex_selector,
    bk_path_vertex_selector,
    ):

    assert 1 < graph.num_vertices()

    self.graph = graph
    self.source_vertex = find_vertex_by_label(
      graph = self.graph,
      label = source_vertex_label,
      )

    result = min_cut_max_flow.stoer_wagner_min_cut( graph = self.graph )
    self.weight = result.weight

    ( self.source_set, self.sink_set ) = result.cutsets

    if self.source_vertex not in self.source_set:
      assert self.source_vertex in self.sink_set
      ( self.source_set, self.sink_set ) = ( self.sink_set, self.source_set )

    assert self.source_vertex in self.source_set

    self.sw_path_vertex_selector = sw_path_vertex_selector
    self.bk_path_vertex_selector = bk_path_vertex_selector

    self.sink_vertex = self.sw_path_vertex_selector( vertices = self.sink_set )

    if self.is_terminal_leaf():
      self.bk_route = TerminalRoute

    else:
      self.bk_route = self.bk_path_vertex_selector( node = self )


  def immediate_onpath_child(self):

    data = self.bk_route.onpath_child( node = self )

    return stirling_tree_onpath_node(
      graph = data.graph,
      source_vertex_label = data.source_vertex_label,
      source_set_labels = data.source_set_labels,
      sink_vertex_label = data.sink_vertex_label,
      sink_set_labels = data.sink_set_labels,
      path_vertex_selector = self.bk_path_vertex_selector,
      )


  def immediate_offpath_children(self):

    if self.is_terminal_leaf():
      raise RuntimeError("No more children")

    ( contracted, label ) = copy_graph_and_merge_vertices(
      graph = self.graph,
      merge1 = self.source_vertex,
      merge2 = self.sink_vertex,
      )
    yield stirling_tree_sw_node(
      graph = contracted,
      source_vertex_label = label,
      sw_path_vertex_selector = self.sw_path_vertex_selector,
      bk_path_vertex_selector= self.bk_path_vertex_selector,
      )

    data = self.bk_route.offpath_child( node = self )

    yield stirling_tree_bk_node(
      graph = data.graph,
      source_vertex_label = data.source_vertex_label,
      sink_vertex_label = data.sink_vertex_label,
      path_vertex_selector = self.bk_path_vertex_selector,
      )


class stirling_tree_bk_node(stirling_tree_node):
  """
  A node in the Stirling-tree that uses Boykov-Kolmogorov minimum cut

  Branches 2-way
  """

  def __init__(
    self,
    graph,
    source_vertex_label,
    sink_vertex_label,
    path_vertex_selector,
    ):

    assert 1 < graph.num_vertices()
    self.graph = graph
    ( self.source_vertex, self.sink_vertex ) = find_two_vertices_by_label(
      graph = self.graph,
      label1 = source_vertex_label,
      label2 = sink_vertex_label,
      )

    result = min_cut_max_flow.boykov_kolmogorov_min_st_cut(
      graph = self.graph,
      source = self.source_vertex,
      sink = self.sink_vertex,
      )
    self.weight = result.weight

    ( self.source_set, self.sink_set ) = result.cutsets

    assert self.source_vertex in self.source_set
    assert self.sink_vertex in self.sink_set

    self.path_vertex_selector = path_vertex_selector

    if self.is_terminal_leaf():
      self.route = TerminalRoute

    else:
      self.route = self.path_vertex_selector( node = self )


  def immediate_onpath_child(self):

    data = self.route.onpath_child( node = self )

    return stirling_tree_onpath_node(
      graph = data.graph,
      source_vertex_label = data.source_vertex_label,
      source_set_labels = data.source_set_labels,
      sink_vertex_label = data.sink_vertex_label,
      sink_set_labels = data.sink_set_labels,
      path_vertex_selector = self.path_vertex_selector,
      )


  def immediate_offpath_children(self):

    data = self.route.offpath_child( node = self )

    yield stirling_tree_bk_node(
      graph = data.graph,
      source_vertex_label = data.source_vertex_label,
      sink_vertex_label = data.sink_vertex_label,
      path_vertex_selector = self.path_vertex_selector,
      )


class stirling_tree_onpath_node(stirling_tree_node):
  """
  A node in the Stirling tree en route to the minimum cut

  Branches 2-way
  """

  def __init__(
    self,
    graph,
    source_vertex_label,
    source_set_labels,
    sink_vertex_label,
    sink_set_labels,
    path_vertex_selector
    ):

    assert 1 < graph.num_vertices()
    assert source_vertex_label in source_set_labels
    assert sink_vertex_label in sink_set_labels

    self.graph = graph

    mapping = calculate_vertex_by_label_dict( graph = self.graph )

    self.source_vertex = mapping[ source_vertex_label ]
    self.sink_vertex = mapping[ sink_vertex_label ]
    self.source_set = set( mapping[ l ] for l in source_set_labels )
    self.sink_set = set( mapping[ l ] for l in sink_set_labels )

    self.path_vertex_selector = path_vertex_selector

    if self.is_terminal_leaf():
      self.route = TerminalRoute

    else:
      self.route = self.path_vertex_selector( node = self )


  def immediate_onpath_child(self):

    data = self.route.onpath_child( node = self )

    return stirling_tree_onpath_node(
      graph = data.graph,
      source_vertex_label = data.source_vertex_label,
      source_set_labels = data.source_set_labels,
      sink_vertex_label = data.sink_vertex_label,
      sink_set_labels = data.sink_set_labels,
      path_vertex_selector = self.path_vertex_selector,
      )


  def immediate_offpath_children(self):

    data = self.route.offpath_child( node = self )

    yield stirling_tree_bk_node(
      graph = data.graph,
      source_vertex_label = data.source_vertex_label,
      sink_vertex_label = data.sink_vertex_label,
      path_vertex_selector = self.path_vertex_selector,
      )


def arbitrary_selection_from_set(vertices):

  return next( iter( vertices ) )


class BKState(object):
  """
  A named tuple to hold all the data
  """

  def __init__(
    self,
    graph,
    source_vertex_label,
    source_set_labels,
    sink_vertex_label,
    sink_set_labels,
    ):

    self.graph = graph
    self.source_vertex_label = source_vertex_label
    self.source_set_labels = source_set_labels
    self.sink_vertex_label = sink_vertex_label
    self.sink_set_labels = sink_set_labels


class TerminalRoute(object):
  """
  A singleton object that signifies no more steps in the tree
  """

  @staticmethod
  def onpath_child(node):

    raise RuntimeError("No more children")


  @staticmethod
  def offpath_child(node):

    raise RuntimeError("No more children")


class BKRoute(object):
  """
  Route using the Boykov-Kolmogorov algorithm
  """

  def __init__(self, vertex, onpath, offpath):

    self.vertex = vertex
    self.onpath = onpath
    self.offpath = offpath


  def onpath_child(self, node):

    return self.get_child_state( node = node, method = self.onpath )


  def offpath_child(self, node):

    return self.get_child_state( node = node, method = self.offpath )


  def get_child_state(self, node, method):

    ( graph, source_vertex_label, sink_vertex_label ) = method(
      node = node,
      vertex = self.vertex,
      )

    source_set_labels = set(
      node.graph.vertex_label( vertex = v ) for v in node.source_set
      if v != node.source_vertex and v != self.vertex
      )
    source_set_labels.add( source_vertex_label )

    sink_set_labels = set(
      node.graph.vertex_label( vertex = v ) for v in node.sink_set
      if v != node.sink_vertex and v != self.vertex
      )
    sink_set_labels.add( sink_vertex_label )

    return BKState(
      graph = graph,
      source_vertex_label = source_vertex_label,
      source_set_labels = source_set_labels,
      sink_vertex_label = sink_vertex_label,
      sink_set_labels = sink_set_labels,
      )


  @staticmethod
  def merge_with_source(node, vertex):

    ( contracted, label ) = copy_graph_and_merge_vertices(
      graph = node.graph,
      merge1 = node.source_vertex,
      merge2 = vertex,
      )

    return ( contracted, label, node.graph.vertex_label( vertex = node.sink_vertex ) )


  @staticmethod
  def merge_with_sink(node, vertex):

    ( contracted, label ) = copy_graph_and_merge_vertices(
      graph = node.graph,
      merge1 = node.sink_vertex,
      merge2 = vertex,
      )

    return ( contracted, node.graph.vertex_label( vertex = node.source_vertex ), label )


  @classmethod
  def SourceSideChoice(cls, vertex):

    return cls(
      vertex = vertex,
      onpath = cls.merge_with_source,
      offpath = cls.merge_with_sink,
      )


  @classmethod
  def SinkSideChoice(cls, vertex):

    return cls(
      vertex = vertex,
      onpath = cls.merge_with_sink,
      offpath = cls.merge_with_source,
      )


class BKEqualSizeSelector(object):
  """
  Selects vertices to make the source and sink sets equal in size
  """

  def __init__(self, selector = arbitrary_selection_from_set):

    self.selector = selector


  def __call__(self, node):

    assert not node.is_terminal_leaf()

    if len( node.sink_set ) <= len( node.source_set ):
      selected = self.selector( vertices = node.available_source_set_vertices() )
      return BKRoute.SourceSideChoice( vertex = selected )

    else:
      selected = self.selector( vertices = node.available_sink_set_vertices() )
      return BKRoute.SinkSideChoice( vertex = selected )


class BKPreferentialSourceSideSelector(object):
  """
  Selects preferentially from the source set until it is depleted
  """

  def __init__(self, selector = arbitrary_selection_from_set):

    self.selector = selector


  def __call__(self, node):

    assert not node.is_terminal_leaf()

    remaining_source = node.available_source_set_vertices()

    if remaining_source:
      selected = self.selector( vertices = remaining_source )
      return BKRoute.SourceSideChoice( vertex = selected )

    else:
      remaining_sink = node.available_sink_set_vertices()
      selected = self.selector( vertices = remaining_sink )
      return BKRoute.SinkSideChoice( vertex = selected )


def construct_vertex_labelled_graph(graph):

  ( copy, mapping ) = utility.copy_graph_and_map_vertices( graph = graph )

  for v in graph.vertices():
    copy.set_vertex_label( vertex = mapping[ v ], label = frozenset( ( v, ) ) )

  return copy


def find_vertex_by_label(graph, label):

  for vertex in graph.vertices():
    if graph.vertex_label( vertex = vertex ) == label:
      return vertex

  raise KeyError("vertex not found")


def find_two_vertices_by_label(graph, label1, label2):

  vertex1 = None
  vertex2 = None

  for vertex in graph.vertices():
    label = graph.vertex_label( vertex = vertex )

    if vertex1 is None and label == label1:
      vertex1 = vertex

      if vertex2 is not None:
        return ( vertex1, vertex2 )

    elif vertex2 is None and label == label2:
      vertex2 = vertex

      if vertex1 is not None:
        return ( vertex1, vertex2 )

  raise KeyError("vertices not found")


def calculate_vertex_by_label_dict(graph):

  return dict(
    ( graph.vertex_label( vertex  = v ), v ) for v in graph.vertices()
    )


contract_undirected_graph = utility.undirected_graph_vertex_contraction(
  joint_vertex_label = operator.or_,
  joint_edge_weight = operator.add,
  )


def copy_graph_and_merge_vertices(graph, merge1, merge2):

  ( copy, mapping ) = utility.copy_graph_and_map_vertices( graph = graph )
  contract_undirected_graph(
    graph = copy,
    v1 = mapping[ merge1 ],
    v2 = mapping[ merge2 ],
    )

  return (
    copy,
    contract_undirected_graph.joint_vertex_label(
      graph.vertex_label( vertex = merge1 ),
      graph.vertex_label( vertex = merge2 ),
      ),
    )
