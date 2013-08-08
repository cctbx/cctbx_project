from __future__ import division

import boost.python
ext = boost.python.import_ext( "scitbx_suffixtree_shared_ext" )
from scitbx_suffixtree_shared_ext import *

def dump(module, root, word):
  
  for ( index, edge ) in enumerate( module.preorder_iteration( root = root ) ):
    print "%s: %s" % (
      index,
      "".join( str( word[ i ] ) for i in range( edge.start, edge.stop ) ),
      )
    
    
def tree_string(module, tree):
  
  result = []
  
  depth_for = calculate_edge_depth( module = module, root = tree.root )
  iter = module.preorder_iteration( root = tree.root )
  iter.next()
  
  for edge in iter: 
    parent = edge.parent()
    result.append(
      "-" * depth_for[ parent ]
      + "".join( str( tree.word[ i ] ) for i in range( edge.start, edge.stop ) )
      + " (%s-%s, %s)" % (
        edge.start,
        edge.stop,
        "root" if edge.is_root() else ( "leaf (index=%s)" % edge.label if edge.is_leaf() else "branch" ),
        )
      )
  
  return "\n".join( result )


def calculate_edge_depth(module, root):
  
  depth_for = { root: 0 }
  iter = module.preorder_iteration( root = root )
  iter.next() # discard root
  
  for edge in iter:
    parent = edge.parent()
    depth_for[ edge ] = depth_for[ parent ] + edge.stop - edge.start
  
  return depth_for


def label(edge, word, separator = "-"):
  
  return separator.join( str( word[ i ] ) for i in range( edge.start, edge.stop ) )