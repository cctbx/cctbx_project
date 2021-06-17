from __future__ import absolute_import, division, print_function

from scitbx import suffixtree
from scitbx.suffixtree import single

import unittest
from six.moves import range

class TestWord(unittest.TestCase):

  def test_single_word(self):

    word = single.word()

    ld = word.length_descriptor()

    self.assertEqual( len( word ), 0 )
    self.assertEqual( ld(), 0 )
    self.assertEqual( list( word.substring( 0, ld() ) ), [] )

    c1 = object()
    c2 = object()
    word.append( c1 )
    self.assertEqual( len( word ), 1 )
    self.assertEqual( ld(), 1 )
    self.assertEqual( word[0], c1 )
    self.assertEqual( list( word ), [ c1 ] )
    word.append( c2 )
    self.assertEqual( len( word ), 2 )
    self.assertEqual( ld(), 2 )
    self.assertEqual( word[0], c1 )
    self.assertEqual( word[1], c2 )
    self.assertEqual( list( word ), [ c1, c2 ] )

    self.assertEqual( list( word.substring( 0, 0 ) ), [] )
    self.assertEqual( list( word.substring( 0, 1 ) ), [ c1 ] )
    self.assertEqual( list( word.substring( 0, 2 ) ), [ c1, c2 ] )

    self.assertEqual( list( word.substring( 1, 1 ) ), [] )
    self.assertEqual( list( word.substring( 1, 2 ) ), [ c2 ] )


class TestEdge(unittest.TestCase):

  def setUp(self):

    self.root = single.edge.root()
    self.branch = single.edge.branch( start = 1, stop = 2 )
    self.word = single.word()
    self.leaf = single.edge.leaf( start = 4, length = self.word.length_descriptor(), label = 6 )


  def test_root(self):

    self.assertEqual( self.root.start, 0 )
    self.assertRaises( RuntimeError, setattr, self.root, "start", 1 )
    self.assertEqual( self.root.stop, 0 )
    self.assertRaises( RuntimeError, getattr, self.root, "label" )
    self.assertRaises( RuntimeError, getattr, self.root, "parent" )
    self.assertRaises( RuntimeError, setattr, self.root, "parent", self.branch )
    self.assertRaises( RuntimeError, getattr, self.root, "suffix" )
    self.assertRaises( RuntimeError, setattr, self.root, "suffix", self.branch )

    self.assertTrue( self.root.is_root() )
    self.assertFalse( self.root.is_leaf() )

    self.assertTrue( self.root.is_empty() )
    self.assertEqual( self.root.keys(), [] )
    self.assertFalse( 1 in self.root )
    self.assertTrue( 1 not in self.root )
    self.assertRaises( KeyError, self.root.__getitem__, 1 )

    self.root[1] = self.branch
    self.root[1] = self.leaf
    self.root[2] = self.branch
    self.assertFalse( self.root.is_empty() )
    self.assertEqual( set( self.root.keys() ), set( [ 1, 2  ] ) )
    self.assertEqual( set( self.root.values() ), set( [ self.branch, self.leaf ] ) )
    self.assertTrue( 1 in self.root )
    self.assertTrue( 2 in self.root )
    self.assertEqual( self.root[1], self.leaf )
    self.assertEqual( self.root[2], self.branch )


  def test_branch(self):

    self.assertEqual( self.branch.start, 1 )
    self.branch.start = 3
    self.assertEqual( self.branch.start, 3 )
    self.assertEqual( self.branch.stop, 2 )
    self.assertRaises( RuntimeError, getattr, self.branch, "label" )

    self.assertEqual( self.branch.parent, None )
    self.branch.parent = self.root
    self.assertTrue( isinstance( self.branch.parent, single.edge ) )
    self.assertEqual( self.branch.parent, self.root )

    self.assertEqual( self.branch.suffix, None )
    self.branch.suffix = self.root
    self.assertTrue( isinstance( self.branch.parent, single.edge ) )
    self.assertEqual( self.branch.suffix, self.root )

    self.assertFalse( self.branch.is_root() )
    self.assertFalse( self.branch.is_leaf() )

    self.assertTrue( self.branch.is_empty() )
    self.assertEqual( self.branch.keys(), [] )
    self.assertFalse( 1 in self.branch )
    self.assertTrue( 1 not in self.branch )
    self.assertRaises( KeyError, self.branch.__getitem__, 1 )

    self.branch[1] = self.root
    self.branch[1] = self.leaf
    self.branch[2] = self.root
    self.assertFalse( self.branch.is_empty() )
    self.assertEqual( set( self.branch.keys() ), set( [ 1, 2  ] ) )
    self.assertEqual( set( self.branch.values() ), set( [ self.root, self.leaf ] ) )
    self.assertTrue( 1 in self.branch )
    self.assertTrue( 2 in self.branch )
    self.assertEqual( self.branch[1], self.leaf )
    self.assertEqual( self.branch[2], self.root )


  def test_leaf(self):

    self.assertEqual( self.leaf.start, 4 )
    self.leaf.start = 5
    self.assertEqual( self.leaf.start, 5 )
    ld = self.word.length_descriptor()
    self.assertEqual( self.leaf.stop, ld() )
    self.word.append( 1 )
    self.word.append( 2 )
    self.assertEqual( self.leaf.stop, ld() )

    self.assertEqual( self.leaf.label, 6 )

    self.assertEqual( self.leaf.parent, None )
    self.leaf.parent = self.root
    self.assertTrue( isinstance( self.leaf.parent, single.edge ) )
    self.assertEqual( self.leaf.parent, self.root )

    self.assertRaises( RuntimeError, getattr, self.leaf, "suffix" )
    self.assertRaises( RuntimeError, setattr, self.leaf, "suffix", self.branch )

    self.assertFalse( self.leaf.is_root() )
    self.assertTrue( self.leaf.is_leaf() )

    self.assertTrue( self.leaf.is_empty() )
    self.assertEqual( self.leaf.keys(), [] )
    self.assertEqual( self.leaf.values(), [] )
    self.assertFalse( 1 in self.leaf )
    self.assertTrue( 1 not in self.leaf )
    self.assertRaises( KeyError, self.leaf.__getitem__, 1 )
    self.assertRaises( RuntimeError, self.leaf.__setitem__, 1, self.root )


class TestPreOrder(unittest.TestCase):

  def setUp(self):

    self.root = single.edge.root()
    self.word = single.word()
    length = self.word.length_descriptor()
    self.edge_named = {
      "r_b1": single.edge.branch( start = 0, stop = 3 ),
      "r_b2": single.edge.branch( start = 1, stop = 3 ),
      "r_l1": single.edge.leaf( start = 4, length = length, label = 6 ),

      "r_b1_b1": single.edge.branch( start = 2, stop = 3 ),
      "r_b1_l1": single.edge.leaf( start = 4, length = length, label = 7 ),

      "r_b1_b1_l1": single.edge.leaf( start = 4, length = length, label = 8 ),
      "r_b1_b1_l2": single.edge.leaf( start = 4, length = length, label = 9 ),

      "r_b2_l1": single.edge.leaf( start = 4, length = length, label = 10 ),
      "r_b2_l2": single.edge.leaf( start = 4, length = length, label = 11 ),
      }
    self.root[ "r_b1" ] = self.edge_named[ "r_b1" ]
    self.root[ "r_b2" ] = self.edge_named[ "r_b2" ]
    self.root[ "r_l1" ] = self.edge_named[ "r_l1" ]

    self.edge_named[ "r_b1" ][ "r_b1_b1" ] = self.edge_named[ "r_b1_b1" ]
    self.edge_named[ "r_b1" ][ "r_b1_l1" ] = self.edge_named[ "r_b1_l1" ]

    self.edge_named[ "r_b1_b1" ][ "r_b1_b1_l1" ] = self.edge_named[ "r_b1_b1_l1" ]
    self.edge_named[ "r_b1_b1" ][ "r_b1_b1_l2" ] = self.edge_named[ "r_b1_b1_l2" ]

    self.edge_named[ "r_b2" ][ "r_b2_l1" ] = self.edge_named[ "r_b2_l1" ]
    self.edge_named[ "r_b2" ][ "r_b2_l2" ] = self.edge_named[ "r_b2_l2" ]

  def test_single_root(self):

    root = single.edge.root()
    res = list( root.preorder_iteration() )
    self.assertEqual( res, [ root ] )
    self.assertTrue( isinstance( res[0], single.edge ) )


  def test_single_branch(self):

    branch = single.edge.branch( start = 0, stop = 3 )
    res = list( branch.preorder_iteration() )
    self.assertEqual( res, [ branch ] )
    self.assertTrue( isinstance( res[0], single.edge ) )


  def test_single_leaf(self):

    word = single.word()
    leaf = single.edge.leaf( start = 4, length = word.length_descriptor(), label = 6 )
    res = list( leaf.preorder_iteration() )
    self.assertEqual( res, [ leaf ] )
    self.assertTrue( isinstance( res[0], single.edge ) )


  def check_leaf_named(self, result, leafname):

    self.assertEqual( result[0], self.edge_named[ leafname ] )
    result.pop( 0 )


  def check_r_b2(self, result):

    action_for = {
      "r_b2_l1": lambda result: self.check_leaf_named( result, "r_b2_l1" ),
      "r_b2_l2": lambda result: self.check_leaf_named( result, "r_b2_l2" ),
      }

    self.assertEqual( result[0], self.edge_named[ "r_b2" ] )
    result.pop( 0 )

    for key in self.edge_named[ "r_b2" ].keys():
      action_for[ key ]( result = result )


  def check_r_b1_b1(self, result):

    action_for = {
      "r_b1_b1_l1": lambda result: self.check_leaf_named( result, "r_b1_b1_l1" ),
      "r_b1_b1_l2": lambda result: self.check_leaf_named( result, "r_b1_b1_l2" ),
      }

    self.assertEqual( result[0], self.edge_named[ "r_b1_b1" ] )
    result.pop( 0 )

    for key in self.edge_named[ "r_b1_b1" ].keys():
      action_for[ key ]( result = result )


  def check_r_b1(self, result):

    action_for = {
      "r_b1_b1": self.check_r_b1_b1,
      "r_b1_l1": lambda result: self.check_leaf_named( result, "r_b1_l1" ),
      }

    self.assertEqual( result[0], self.edge_named[ "r_b1" ] )
    result.pop( 0 )

    for key in self.edge_named[ "r_b1" ].keys():
      action_for[ key ]( result = result )


  def check_root(self, result):

    action_for = {
      "r_b1": self.check_r_b1,
      "r_b2": self.check_r_b2,
      "r_l1": lambda result: self.check_leaf_named( result, "r_l1" ),
      }

    self.assertEqual( result[0], self.root )
    result.pop( 0 )

    for key in self.root.keys():
      action_for[ key ]( result = result )


  def test_r_b2(self):

    result = list( self.edge_named[ "r_b2" ].preorder_iteration() )
    self.assertEqual( len( result ), 3 )
    self.check_r_b2( result )
    self.assertEqual( result, [] )


  def test_r_b1_b1(self):

    result = list( self.edge_named[ "r_b1_b1" ].preorder_iteration() )
    self.assertEqual( len( result ), 3 )
    self.check_r_b1_b1( result )
    self.assertEqual( result, [] )


  def test_r_b1(self):

    result = list( self.edge_named[ "r_b1" ].preorder_iteration() )
    self.assertEqual( len( result ), 5 )
    self.check_r_b1( result )
    self.assertEqual( result, [] )


  def test_root(self):

    result = list( self.root.preorder_iteration() )
    self.assertEqual( len( result ), 10 )
    self.assertTrue( all( isinstance( e, single.edge ) for e in result ) )
    self.check_root( result )
    self.assertEqual( result, [] )


class TestPostOrder(unittest.TestCase):

  def setUp(self):

    self.root = single.edge.root()
    self.word = single.word()
    length = self.word.length_descriptor()
    self.edge_named = {
      "r_b1": single.edge.branch( start = 0, stop = 3 ),
      "r_b2": single.edge.branch( start = 1, stop = 3 ),
      "r_l1": single.edge.leaf( start = 4, length = length, label = 6 ),

      "r_b1_b1": single.edge.branch( start = 2, stop = 3 ),
      "r_b1_l1": single.edge.leaf( start = 4, length = length, label = 7 ),

      "r_b1_b1_l1": single.edge.leaf( start = 4, length = length, label = 8 ),
      "r_b1_b1_l2": single.edge.leaf( start = 4, length = length, label = 9 ),

      "r_b2_l1": single.edge.leaf( start = 4, length = length, label = 10 ),
      "r_b2_l2": single.edge.leaf( start = 4, length = length, label = 11 ),
      }
    self.root[ "r_b1" ] = self.edge_named[ "r_b1" ]
    self.root[ "r_b2" ] = self.edge_named[ "r_b2" ]
    self.root[ "r_l1" ] = self.edge_named[ "r_l1" ]

    self.edge_named[ "r_b1" ][ "r_b1_b1" ] = self.edge_named[ "r_b1_b1" ]
    self.edge_named[ "r_b1" ][ "r_b1_l1" ] = self.edge_named[ "r_b1_l1" ]

    self.edge_named[ "r_b1_b1" ][ "r_b1_b1_l1" ] = self.edge_named[ "r_b1_b1_l1" ]
    self.edge_named[ "r_b1_b1" ][ "r_b1_b1_l2" ] = self.edge_named[ "r_b1_b1_l2" ]

    self.edge_named[ "r_b2" ][ "r_b2_l1" ] = self.edge_named[ "r_b2_l1" ]
    self.edge_named[ "r_b2" ][ "r_b2_l2" ] = self.edge_named[ "r_b2_l2" ]

  def test_single_root(self):

    root = single.edge.root()
    res = list( root.postorder_iteration() )
    self.assertEqual( res, [ root ] )
    self.assertTrue( isinstance( res[0], single.edge ) )


  def test_single_branch(self):

    branch = single.edge.branch( start = 0, stop = 3 )
    res = list( branch.postorder_iteration() )
    self.assertEqual( res, [ branch ] )
    self.assertTrue( isinstance( res[0], single.edge ) )


  def test_single_leaf(self):

    word = single.word()
    leaf = single.edge.leaf( start = 4, length = word.length_descriptor(), label = 6 )
    res = list( leaf.postorder_iteration() )
    self.assertEqual( res, [ leaf ] )
    self.assertTrue( isinstance( res[0], single.edge ) )


  def check_leaf_named(self, result, leafname):

    self.assertEqual( result[0], self.edge_named[ leafname ] )
    result.pop( 0 )


  def check_r_b2(self, result):

    action_for = {
      "r_b2_l1": lambda result: self.check_leaf_named( result, "r_b2_l1" ),
      "r_b2_l2": lambda result: self.check_leaf_named( result, "r_b2_l2" ),
      }

    for key in self.edge_named[ "r_b2" ].keys():
      action_for[ key ]( result = result )

    self.assertEqual( result[0], self.edge_named[ "r_b2" ] )
    result.pop( 0 )


  def check_r_b1_b1(self, result):

    action_for = {
      "r_b1_b1_l1": lambda result: self.check_leaf_named( result, "r_b1_b1_l1" ),
      "r_b1_b1_l2": lambda result: self.check_leaf_named( result, "r_b1_b1_l2" ),
      }

    for key in self.edge_named[ "r_b1_b1" ].keys():
      action_for[ key ]( result = result )

    self.assertEqual( result[0], self.edge_named[ "r_b1_b1" ] )
    result.pop( 0 )


  def check_r_b1(self, result):

    action_for = {
      "r_b1_b1": self.check_r_b1_b1,
      "r_b1_l1": lambda result: self.check_leaf_named( result, "r_b1_l1" ),
      }

    for key in self.edge_named[ "r_b1" ].keys():
      action_for[ key ]( result = result )

    self.assertEqual( result[0], self.edge_named[ "r_b1" ] )
    result.pop( 0 )


  def check_root(self, result):

    action_for = {
      "r_b1": self.check_r_b1,
      "r_b2": self.check_r_b2,
      "r_l1": lambda result: self.check_leaf_named( result, "r_l1" ),
      }

    for key in self.root.keys():
      action_for[ key ]( result = result )

    self.assertEqual( result[0], self.root )
    result.pop( 0 )


  def test_r_b2(self):

    result = list( self.edge_named[ "r_b2" ].postorder_iteration() )
    self.assertEqual( len( result ), 3 )
    self.check_r_b2( result )
    self.assertEqual( result, [] )


  def test_r_b1_b1(self):

    result = list( self.edge_named[ "r_b1_b1" ].postorder_iteration() )
    self.assertEqual( len( result ), 3 )
    self.check_r_b1_b1( result )
    self.assertEqual( result, [] )


  def test_r_b1(self):

    result = list( self.edge_named[ "r_b1" ].postorder_iteration() )
    self.assertEqual( len( result ), 5 )
    self.check_r_b1( result )
    self.assertEqual( result, [] )


  def test_root(self):

    result = list( self.root.postorder_iteration() )
    self.assertEqual( len( result ), 10 )
    self.assertTrue( all( isinstance( e, single.edge ) for e in result ) )
    self.check_root( result )
    self.assertEqual( result, [] )


class TestTree(unittest.TestCase):

  def test_tree(self):

    tree = single.tree()
    self.assertEqual( tree.in_construction, False )
    w = tree.word
    self.assertTrue( isinstance( w, single.word ) )
    self.assertEqual( w.length_descriptor()(), 0 )
    r = tree.root
    self.assertTrue( isinstance( r, single.edge ) )
    self.assertEqual( r.keys(), [] )


  def test_ukkonen1(self):

    tree = single.tree()
    self.assertEqual( tree.in_construction, False )

    builder = single.ukkonen( tree )
    self.assertTrue( tree.in_construction )
    self.assertTrue( builder.is_attached )
    self.assertTrue( builder.is_valid )

    self.assertRaises( RuntimeError, single.ukkonen, tree )

    builder.detach()
    self.assertFalse( builder.is_attached )
    self.assertTrue( builder.is_valid )
    self.assertFalse( tree.in_construction )
    self.assertRaises( RuntimeError, builder.append, "a" )


  def test_ukkonen2(self):

    tree = single.tree()
    builder = single.ukkonen( tree )

    builder.append( glyph = "a" )
    self.assertTrue( builder.is_valid )

    builder.append( glyph = "n" )
    self.assertTrue( builder.is_valid )

    builder.detach()
    self.assertFalse( builder.is_attached )
    self.assertFalse( tree.in_construction )

    builder = single.ukkonen( tree )
    self.assertTrue( builder.is_valid )
    self.assertTrue( builder.is_attached )
    self.assertTrue( tree.in_construction )

    builder.append( glyph = "a" )
    self.assertFalse( builder.is_valid )
    self.assertRaises( RuntimeError, builder.detach )

    builder.append( glyph = "n" )
    self.assertFalse( builder.is_valid )
    self.assertRaises( RuntimeError, builder.detach )

    builder.append( glyph = "a" )
    self.assertFalse( builder.is_valid )
    self.assertRaises( RuntimeError, builder.detach )

    builder.append( glyph = "s" )
    self.assertTrue( builder.is_valid )

    builder.detach()
    self.assertFalse( builder.is_attached )
    self.assertFalse( tree.in_construction )

    builder = single.ukkonen( tree )
    self.assertTrue( builder.is_valid )
    self.assertTrue( builder.is_attached )
    self.assertTrue( tree.in_construction )

    builder.append( glyph = "$" )
    self.assertTrue( builder.is_valid )

    builder.detach()

    root = tree.root
    self.assertTrue( root.is_root() )
    self.assertEqual( set( root.keys() ), set( [ "a", "n", "s", "$" ] ) )

    b_a = root[ "a" ]
    self.assertFalse( b_a.is_root() )
    self.assertFalse( b_a.is_leaf() )
    self.assertEqual( b_a.start, 0 )
    self.assertEqual( b_a.stop, 1 )
    self.assertEqual( set( b_a.keys() ), set( [ "n", "s" ] ) )
    self.assertEqual( b_a.parent, root )

    b_a_n = b_a[ "n" ]
    self.assertFalse( b_a_n.is_root() )
    self.assertFalse( b_a_n.is_leaf() )
    self.assertEqual( b_a_n.start, 1 )
    self.assertEqual( b_a_n.stop, 3 )
    self.assertEqual( set( b_a_n.keys() ), set( [ "n", "s" ] ) )
    self.assertEqual( b_a_n.parent, b_a )

    b_a_n_n = b_a_n[ "n" ]
    self.assertTrue( b_a_n_n.is_leaf() )
    self.assertEqual( b_a_n_n.start, 3 )
    self.assertEqual( b_a_n_n.stop, 7 )
    self.assertEqual( b_a_n_n.label, 0 )
    self.assertEqual( b_a_n_n.parent, b_a_n )

    b_a_n_s = b_a_n[ "s" ]
    self.assertTrue( b_a_n_s.is_leaf() )
    self.assertEqual( b_a_n_s.start, 5 )
    self.assertEqual( b_a_n_s.stop, 7 )
    self.assertEqual( b_a_n_s.label, 2 )
    self.assertEqual( b_a_n_s.parent, b_a_n )

    b_a_s = b_a[ "s" ]
    self.assertTrue( b_a_s.is_leaf() )
    self.assertEqual( b_a_s.start, 5 )
    self.assertEqual( b_a_s.stop, 7 )
    self.assertEqual( b_a_s.label, 4 )
    self.assertEqual( b_a_s.parent, b_a )

    b_n = root[ "n" ]
    self.assertFalse( b_n.is_root() )
    self.assertFalse( b_n.is_leaf() )
    self.assertEqual( b_n.start, 1 )
    self.assertEqual( b_n.stop, 3 )
    self.assertEqual( set( b_n.keys() ), set( [ "n", "s" ] ) )
    self.assertEqual( b_n.parent, root )

    b_n_n = b_n[ "n" ]
    self.assertTrue( b_n_n.is_leaf() )
    self.assertEqual( b_n_n.start, 3 )
    self.assertEqual( b_n_n.stop, 7 )
    self.assertEqual( b_n_n.label, 1 )
    self.assertEqual( b_n_n.parent, b_n )

    b_n_s = b_n[ "s" ]
    self.assertTrue( b_n_s.is_leaf() )
    self.assertEqual( b_n_s.start, 5 )
    self.assertEqual( b_n_s.stop, 7 )
    self.assertEqual( b_n_s.label, 3 )
    self.assertEqual( b_n_s.parent, b_n )

    b_s = root[ "s" ]
    self.assertTrue( b_s.is_leaf() )
    self.assertEqual( b_s.start, 5 )
    self.assertEqual( b_s.stop, 7 )
    self.assertEqual( b_s.label, 5 )
    self.assertEqual( b_s.parent, root )

    b_dl = root[ "$" ]
    self.assertTrue( b_dl.is_leaf() )
    self.assertEqual( b_dl.start, 6 )
    self.assertEqual( b_dl.stop, 7 )
    self.assertEqual( b_dl.label, 6 )
    self.assertEqual( b_dl.parent, root )

    self.assertEqual( b_a.suffix, root )
    self.assertEqual( b_n.suffix, b_a )
    self.assertEqual( b_a_n.suffix, b_n )


class TestMatchingStatistics(unittest.TestCase):

  def test(self):

    tree = single.tree()
    result = list( single.matching_statistics( tree, list( "ABC" ) ) )
    self.assertEqual( result, [ ( 0, ( tree.root, 0 ) ) ] * 3 )

    builder = single.ukkonen( tree = tree )

    for c in "TTAGC$":
      builder.append( c )

    self.assertRaises( RuntimeError, single.matching_statistics, tree, list() );

    builder.detach()
    root = tree.root
    branch_t = root[ "T" ]
    leaf_0 = branch_t[ "T" ]
    leaf_1 = branch_t[ "A" ]
    leaf_2 = root[ "A" ]
    leaf_3 = root[ "G" ]
    leaf_4 = root[ "C" ]

    result = list( single.matching_statistics( tree, [] ) )
    self.assertEqual( result, [] )

    result = list( single.matching_statistics( tree, list( "QTTATTATTTAGCQWTTAGFK" ) ) )

    self.assertEqual(
      result,
      [
        ( 0, (   root, 0 ) ),
        ( 3, ( leaf_0, 3 ) ),
        ( 2, ( leaf_1, 3 ) ),
        ( 1, ( leaf_2, 3 ) ),
        ( 3, ( leaf_0, 3 ) ),
        ( 2, ( leaf_1, 3 ) ),
        ( 1, ( leaf_2, 3 ) ),
        ( 2, ( leaf_0, 2 ) ),
        ( 5, ( leaf_0, 5 ) ),
        ( 4, ( leaf_1, 5 ) ),
        ( 3, ( leaf_2, 5 ) ),
        ( 2, ( leaf_3, 5 ) ),
        ( 1, ( leaf_4, 5 ) ),
        ( 0, (   root, 0 ) ),
        ( 0, (   root, 0 ) ),
        ( 4, ( leaf_0, 4 ) ),
        ( 3, ( leaf_1, 4 ) ),
        ( 2, ( leaf_2, 4 ) ),
        ( 1, ( leaf_3, 4 ) ),
        ( 0, (   root, 0 ) ),
        ( 0, (   root, 0 ) ),
        ]
      )


class TestLeafIndices(unittest.TestCase):

  def test(self):

    tree = single.tree()
    result = list( single.matching_statistics( tree, list( "ABC" ) ) )
    self.assertEqual( result, [ ( 0, ( tree.root, 0 ) ) ] * 3 )

    builder = single.ukkonen( tree = tree )

    for c in "TTAGC$":
      builder.append( c )

    builder.detach()
    leaf_indices_below = suffixtree.calculate_leaf_indices( root = tree.root )

    root = tree.root
    branch_t = root[ "T" ]
    leaf_0 = branch_t[ "T" ]
    leaf_1 = branch_t[ "A" ]
    leaf_2 = root[ "A" ]
    leaf_3 = root[ "G" ]
    leaf_4 = root[ "C" ]
    leaf_5 = root[ "$" ]

    self.assertEqual( leaf_indices_below[ leaf_0 ], [ 0 ] )
    self.assertEqual( leaf_indices_below[ leaf_1 ], [ 1 ] )
    self.assertEqual( leaf_indices_below[ leaf_2 ], [ 2 ] )
    self.assertEqual( leaf_indices_below[ leaf_3 ], [ 3 ] )
    self.assertEqual( leaf_indices_below[ leaf_4 ], [ 4 ] )
    self.assertEqual( leaf_indices_below[ leaf_5 ], [ 5 ] )

    self.assertEqual( sorted( leaf_indices_below[ branch_t ] ), [ 0, 1 ] )
    self.assertEqual( sorted( leaf_indices_below[ root ] ), list(range( 6)) )


suite_word = unittest.TestLoader().loadTestsFromTestCase(
  TestWord
  )
suite_edge = unittest.TestLoader().loadTestsFromTestCase(
  TestEdge
  )
suite_preorder = unittest.TestLoader().loadTestsFromTestCase(
  TestPreOrder
  )
suite_postorder = unittest.TestLoader().loadTestsFromTestCase(
  TestPostOrder
  )
suite_tree = unittest.TestLoader().loadTestsFromTestCase(
  TestTree
  )
suite_msi = unittest.TestLoader().loadTestsFromTestCase(
  TestMatchingStatistics
  )
suite_leaf_indices = unittest.TestLoader().loadTestsFromTestCase(
  TestLeafIndices
  )


alltests = unittest.TestSuite(
  [
    suite_word,
    suite_edge,
    suite_preorder,
    suite_postorder,
    suite_tree,
    suite_msi,
    suite_leaf_indices,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )

