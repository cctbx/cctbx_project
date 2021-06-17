from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
from functools import reduce
from six.moves import range
ext = bp.import_ext( "scitbx_suffixtree_shared_ext" )
from scitbx_suffixtree_shared_ext import *

def dump(root, word):

  for ( index, edge ) in enumerate( root.preorder_iteration() ):
    print("%s: %s" % ( index, label( edge = edge, word = tree.word ) ))


def label(edge, word, separator = ""):

  return separator.join(
    str( word[ index ] ) for index in range( edge.start, edge.stop )
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


def tree_print(tree, stream = None):

  if stream is None:
    import sys
    stream = sys.stdout

  depth_for = calculate_edge_depth( root = tree.root )
  it = iter( tree.root.preorder_iteration() )
  next(it)
  stream.write( edge_print( tree.root ) )
  stream.write( "\n" )

  for edge in it:
    stream.write(
      "-" * depth_for[ edge.parent ]
      + label( edge = edge, word = tree.word )
      + " %s" % edge_print( edge )
      )
    stream.write( "\n" )


def calculate_edge_depth(root):

  depth_for = { root: 0 }
  it = iter( root.preorder_iteration() )
  next(it) # discard root

  for edge in it:
    depth_for[ edge ] = depth_for[ edge.parent ] + edge.stop - edge.start

  return depth_for


def calculate_leaf_indices(root):

  import operator
  leaf_indices_below = {}

  for edge in root.postorder_iteration():
    if edge.is_leaf():
      leaf_indices_below[ edge ] = [ edge.label ]

    else:
      leaf_indices_below[ edge ] = reduce(
        operator.add,
        [ leaf_indices_below[ c ] for c in edge.values() ],
        []
        )

  return leaf_indices_below
