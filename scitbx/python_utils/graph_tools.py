# A *VERY* lightweight graph class
# More efficient implementation can be found in the BGL I guess
#
# For my (PHZ) purposes, a simple, lightweight python
# implementation was good enough
#

import sys
from libtbx.utils import Sorry

class graph(object):
  def __init__(self):

    self.node_objects ={}
    self.edge_objects = {}

    self.o = {}

  def insert_node(self,
               name,
               edge_object,
               node_object=None):

    if self.node_objects.has_key( name ):
      if self.node_objects[ name ] == []: # only allow an update if it was blanco
        self.node_objects.update( {name:node_object} )
    else:
      self.node_objects.update( {name:node_object} )

    # The edge objects are stored a s adouble dictionairy
    if not self.edge_objects.has_key( name  ):
      # edge object is not there, place it fully!
      self.edge_objects.update( {name: {} } )
      self.o.update( {name : []} )
      for item in edge_object:
        self.edge_objects[name].update( { item : edge_object[item] } )
        if not item in self.o[name]:
          self.o[name].append( item )
    else:
      # edge object is there
      tmp_edges = self.edge_objects[ name ]
      for item in edge_object:
        if self.edge_objects[ name ].has_key( item ):
          if self.edge_objects[ name ][ item ] is None:
            self.edge_objects[ name ].update( {item :  edge_object[item] } )
        else:
          self.edge_objects[ name ].update( {item :  edge_object[item] } )
        if not item in self.o[name]:
          self.o[name].append( item )


  def remove_node(self,name):
    # remove it from the object list please
    if self.node_objects.has_key( name ):
      del self.node_objects[name]
    # take care of outgoing and incomming
    if self.o.has_key( name ):
      del self.o[name]

    if self.edge_objects.has_key( name ):
      del self.edge_objects[name]

    for item in self.edge_objects:
      if self.edge_objects[item].has_key( name ):
        del self.edge_objects[item][name]

    for item in self.o:
      if name in self.o[ item ]:
        del self.o[item][self.o[item].index( name )]


  def assert_is_clean(self):
    # check if the dictionairy have the same items
    clean=True
    for item in self.node_objects:
      if not self.o.has_key( item ):
        clean = False
    for item in self.o:
      if not self.node_objects.has_key( item ):
        clean = False

    if not clean:
      raise Sorry("The graph is not clean")
    assert( clean )
    # we now also check wether or not the graph is 'done'
    # by checking if outgoing nodes do not point into thin air

    clean = True

    for item in self.node_objects:
      clean = True
      if self.edge_objects.has_key( item ):
        con_list = self.o[ item ]
        if len(con_list) > 0:
          for trial_node in con_list:
            if not self.node_objects.has_key( trial_node ):
              clean=False
      else:
        clean=False

      if not clean:
        raise Sorry("The graph does not seem to be finished, \nsome connections point into thin air!")


  def is_contained_in(self,
                      new_graph):
    # This routine does !NOT! compare topologies
    # but relies on nodes being named in the same way
    # It compares the incommnig and outygoing connections
    # More efficient implementations might be possible
    #
    new_graph.assert_is_clean()

    is_contained_in = True

    for item in self.o:
      if len(self.o[item] ) > len( new_graph.o[item]):
        is_contained_in = False

      for this_one in self.o[ item ]:
        if not this_one in new_graph.o[ item ]:
          is_contained_in = False
    return( is_contained_in )


  def is_equivalent_to(self,new_graph):
    self_in_new = self.is_contained_in(new_graph)
    new_in_self = new_graph.is_contained_in(self)
    is_equivalent = True

    if not self_in_new:
      is_equivalent = False
    if not new_in_self:
      is_equivalent = False

    return is_equivalent


  def find_all_paths(self, start, end, path=[]):
    ## please check
    ## http://www.python.org/doc/essays/graphs.html
    path = path + [start]
    if start == end:
      return [path]
    if not self.o.has_key(start):
      return []
    paths = []
    for node in self.o[start]:
      if node not in path:
        newpaths = self.find_all_paths(node, end, path)
        for newpath in newpaths:
          paths.append(newpath)
    return paths


  def find_paths_of_length(self, start, length, path=[]):
    ## please check
    ## http://www.python.org/doc/essays/graphs.html
    path = path + [start]
    if len(path)==length:
      return [path]
    if not self.o.has_key(start):
      return []
    paths = []
    for node in self.o[start]:
      if node not in path:
        newpaths = self.find_paths_of_length(node, length, path)
        for newpath in newpaths:
          paths.append(newpath)
    return paths




  def find_shortest_path(self, start, end, path=[]):
    ## please check
    ## http://www.python.org/doc/essays/graphs.html
    ## This is not the most optimal way of doing it
    ## and a proper Dijkstra method might be implemented at a later stage
    ## assumed are distance between nodes of equal length
    ## the properties that can be stored there, can represent methods
    ## such as symmetry operators etc etc etc

    path = path + [start]
    if start == end:
      return path
    if not self.o.has_key(start):
      return None
    shortest = None
    for node in self.o[start]:
      if node not in path:
        newpath = self.find_shortest_path(node, end, path)
        if newpath:
          if not shortest or len(newpath) < len(shortest):
            shortest = newpath
    return shortest


  def show(self,out=None):
    if out is None:
      out=sys.stdout
    print >> out
    print >> out, "----------------------------"
    print >> out
    print >> out, "Outgoing edges"
    for node in self.o:
      print >> out, node, "---->", self.o[node]
    print >> out
    for node_1 in self.edge_objects:
      for node_2 in self.edge_objects[ node_1 ] :
        print node_1 , " ---> ", node_2 , " :: ",  \
              self.edge_objects[ node_1 ][ node_2 ]
    print >> out
    print >> out, "----------------------------"
    print >> out
