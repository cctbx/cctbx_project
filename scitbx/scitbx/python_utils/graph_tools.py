import sys, os

# This graph class is intended as follows
#
# self.object: a python dictionairy with a text label:data object dictionary
# self.in: incomming edges/nodes
# self.out: outgoing edges/nodes
#
# weights to edges are not supported in this particular class
# but can be if desired

class graph(object):
  def __init__(self):
    self.objects ={}
    self.i = {}
    self.o = {}


  def insert_node(self,
               name,
               in_names,
               out_names,
               node_object=None):
    if self.objects.has_key( name ):
      if self.objects[ name ] == []: # only allow an update if it was blanco
        self.objects.update( {name:node_object} )
    else:
      self.objects.update( {name:node_object} )

    if not self.i.has_key( name ):
      self.i.update( {name:in_names} )
    else: # this entry is allready here, we have to be carefull now
      for this_one in in_names:
        if not (this_one in self.i[name]):
          self.i[name].append( this_one )

    if not self.o.has_key( name ):
      self.o.update( {name:out_names} )
    else: # this entry is allready here, we have to be carefull now
      for this_one in out_names:
        if not (this_one in self.o[name]):
          self.o[name].append( this_one )

    # update the other nodes: from incomming
    for this_one in in_names:
      # the node this_one is allready in our list
      if self.o.has_key( this_one ):
        if not ( name in self.o[ this_one ] ):
          self.o[ this_one ].append( name )
      # the node is not in our list yet, lets update it please
      else:
        self.insert_node( name=this_one,
                          out_names=[ name ],
                          in_names= [])

    # update the other nodes: from outgoing
    for this_one in out_names:
      if self.i.has_key( this_one ):
        if not (name in self.i[ this_one ]):
          self.i[ this_one ].append( name )
      else:
        # node not present, please update stuff
        self.insert_node( name=this_one,
                          in_names=[ name ],
                          out_names = [])


  def remove_node(self,name):
    # remove it from the object list please
    if self.objects.has_key( name ):
      self.objects.pop( name )
    # take care of outgoing and incomming
    if self.o.has_key( name ):
      out_list = self.o.pop( name )
    if self.i.has_key( name ):
      in_list = self.i.pop( name )

    # the out_list has element that have name in their incomming list
    for item in out_list:
      self.i[ item ].pop( self.i[ item ].index( name ) )

    for item in in_list:
      self.o[ item ].pop( self.o[ item ].index( name ) )


  def assert_is_clean(self):
    # check if the dictionairy have the same items
    clean=True
    for item in self.objects:
      if not self.i.has_key( item ):
        clean = False
      if not self.o.has_key( item ):
        clean = False

    for item in self.i:
      if not self.objects.has_key( item ):
        clean = False
      if not self.o.has_key( item ):
        clean = False

    for item in self.o:
      if not self.objects.has_key( item ):
        clean = False
      if not self.i.has_key( item ):
        clean = False

    assert( clean )


  def is_contained_in(self,
                      new_graph):
    # check if self is a sub graph of new_graph
    new_graph.assert_is_clean()

    is_contained_in = True
    # check first the incomming nodes
    for item in self.i:
      if len( self.i[item] ) > len( new_graph.i[item]):
        is_contained_in = False

      for this_one in self.i[ item ]:
        if not this_one in new_graph.i[ item ]:
          is_contained_in = False

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


  def show(self,out=None):
    if out is None:
      out=sys.stdout
    print >> out
    print >> out, "----------------------------"
    print >> out
    print >> out, "Incomming edges"
    for node in self.i:
      print >> out, node, "<----", self.i[node]
    print >> out, "Outgoing edges"
    for node in self.o:
      print >> out, node, "---->", self.o[node]
    print >> out
    print >> out, "----------------------------"
    print >> out
