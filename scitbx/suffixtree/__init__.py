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


def edge_print(edge):

  if edge.is_root():
    return "(0-0 root)"

  else:
    return "(%s-%s, %s)" % (
      edge.start,
      edge.stop,
      "leaf (index=%s)" % edge.label if edge.is_leaf() else "branch",
      )


def tree_string(module, tree):

  depth_for = calculate_edge_depth( module = module, root = tree.root )
  it = iter( module.preorder_iteration( root = tree.root ) )
  it.next()
  result = [ edge_print( tree.root ) ]

  for edge in it:
    result.append(
      "-" * depth_for[ edge.parent ]
      + "".join( str( tree.word[ i ] ) for i in range( edge.start, edge.stop ) )
      + " %s" % edge_print( edge )
      )

  return "\n".join( result )


def calculate_edge_depth(module, root):

  depth_for = { root: 0 }
  it = iter( module.preorder_iteration( root = root ) )
  it.next() # discard root

  for edge in it:
    depth_for[ edge ] = depth_for[ edge.parent ] + edge.stop - edge.start

  return depth_for


def label(edge, word, separator = "-"):

  return separator.join( str( word[ i ] ) for i in range( edge.start, edge.stop ) )

