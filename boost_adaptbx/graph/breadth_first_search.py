from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
ext = bp.import_ext( "boost_adaptbx_graph_breadth_first_search_ext" )
from boost_adaptbx_graph_breadth_first_search_ext import *

class visitor(object):
  """
  Default visitor for the breadth-first search algorithm and its variants.

  See also:
    breadth_first_search
    breadth_first_visit

  Complete C++ documentation is available at:
    http://www.boost.org/libs/graph/doc/BFSVisitor.html
  # Copyright 2005 The Trustees of Indiana University.

  # Use, modification and distribution is subject to the Boost Software
  # License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
  # http:#www.boost.org/LICENSE_1_0.txt)

  #  Authors: Douglas Gregor
  #           Andrew Lumsdaine
  """

  def initialize_vertex(self, vertex, graph):
    """
    Invoked on each vertex in the graph when the algorithm initializes it.
    """
    pass


  def discover_vertex(self, vertex, graph):
    """
    Invoked on a vertex when it is first "discovered" by the algorithm.
    """
    pass


  def examine_vertex(self, vertex, graph):
    """
    Invoked on a vertex just before its outgoing edges will be
    examined.
    """
    pass


  def examine_edge(self, edge, graph):
    """
    Invoked on an edge as it is being examined.
    """
    pass


  def tree_edge(self, edge, graph):
    """
    Invoked on an edge when it is determined that it is an edge in
    the breadth-first search tree.
    """
    pass


  def non_tree_edge(self, edge, graph):
    """
    Invoked on an edge when it is determined that it is not an
    edge in the breadth-first search tree.
    """
    pass


  def gray_target(self, edge, graph):
    """
    Invoked on a vertex that is the target of the edge being
    examined. The vertex was marked"gray", meaning that it has been
    seen before in the breadth-first traversal but has not yet
    been examined.
    """
    pass


  def black_target(self, edge, graph):
    """
    Invoked on a vertex that is the target of the edge being
    examined. The edge was marked"black", meaning that it has
    already been examined in the breadth-first traversal.
    """
    pass


  def finish_vertex(self, vertex, graph):
    """
    Invoked on a vertex after all of its outgoing edges have been
    examined.
    """
    pass


class vertex_recording_visitor(visitor):

  def __init__(self, start_vertex):

    self.visited_vertices = [ start_vertex ]


  def tree_edge(self, edge, graph):

    self.visited_vertices.append( graph.target( edge ) )


class distance_recording_visitor(visitor):

  def __init__(self, start_vertex):

    self.distance_for = {}
    self.distance_for[ start_vertex ] = 0


  def tree_edge(self, edge, graph):

    source = graph.source( edge )
    target = graph.target( edge )
    self.distance_for[ target ] = self.distance_for[ source ] + 1
