from __future__ import division

from scitbx.suffixtree import single

import unittest    

class TestWord(unittest.TestCase):
  
  def test_single_word(self):
    
    word = single.word()
    
    ld = word.length()
    
    self.assertEqual( ld(), 0 )
    self.assertEqual( list( word.substring( 0, ld() ) ), [] )
    
    c1 = object()
    c2 = object()
    word.append( c1 )
    self.assertEqual( ld(), 1 )
    self.assertEqual( word[0], c1 )
    word.append( c2 )
    self.assertEqual( ld(), 2 )
    self.assertEqual( word[0], c1 )
    self.assertEqual( word[1], c2 )
    
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
    self.leaf = single.edge.leaf( start = 4, length = self.word.length(), label = 6 )
    
    
  def test_root(self):
    
    wr = self.root.weakref()
    self.assertTrue( isinstance( wr, single.weak_edge ) )
    self.assertFalse( wr == self.root )
    self.assertTrue( wr != self.root )
    self.assertTrue( isinstance( wr(), single.edge ) )
    self.assertEqual( wr(), self.root )
    
    self.assertEqual( self.root.start, 0 )
    self.assertRaises( RuntimeError, setattr, self.root, "start", 1 )
    self.assertEqual( self.root.stop, 0 )
    self.assertRaises( RuntimeError, getattr, self.root, "label" )
    self.assertRaises( RuntimeError, getattr, self.root, "parent" )
    self.assertRaises( RuntimeError, setattr, self.root, "parent", self.branch.weakref() )
    self.assertRaises( RuntimeError, getattr, self.root, "suffix" )
    self.assertRaises( RuntimeError, setattr, self.root, "suffix", self.branch.weakref() )
    
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
    self.assertTrue( 1 in self.root )
    self.assertTrue( 2 in self.root )
    self.assertEqual( self.root[1], self.leaf )
    self.assertEqual( self.root[2], self.branch )
    
    
  def test_const_root(self):
    
    cr = single.const_edge.from_edge( self.root )
    
    self.assertEqual( self.root, cr )
    self.assertEqual( hash( self.root ), hash( cr ) )
    self.assertNotEqual( single.edge.root(), cr )
    
    wcr = cr.weakref()
    self.assertTrue( isinstance( wcr, single.const_weak_edge ) )
    self.assertFalse( wcr == cr )
    self.assertTrue( wcr != cr )
    self.assertTrue( isinstance( wcr(), single.const_edge ) )
    self.assertEqual( wcr(), cr )
    
    self.assertEqual( cr.start, 0 )
    self.assertEqual( cr.stop, 0 )
    self.assertRaises( RuntimeError, getattr, cr, "label" )
    self.assertRaises( RuntimeError, getattr, cr, "parent" )
    self.assertRaises( RuntimeError, getattr, cr, "suffix" )
    
    self.assertTrue( cr.is_root() )
    self.assertFalse( cr.is_leaf() )
    
    self.assertTrue( cr.is_empty() )
    self.assertEqual( cr.keys(), [] )
    self.assertFalse( 1 in cr )
    self.assertTrue( 1 not in cr )
    self.assertRaises( KeyError, cr.__getitem__, 1 )
    
  
  def test_branch(self):
    
    wb = self.branch.weakref()
    self.assertTrue( isinstance( wb, single.weak_edge ) )
    self.assertFalse( wb == self.branch )
    self.assertTrue( wb != self.branch )
    self.assertTrue( isinstance( wb(), single.edge ) )
    self.assertEqual( wb(), self.branch )
    
    self.assertEqual( self.branch.start, 1 )
    self.branch.start = 3
    self.assertEqual( self.branch.start, 3 )
    self.assertEqual( self.branch.stop, 2 )
    self.assertRaises( RuntimeError, getattr, self.branch, "label" )
    
    self.assertTrue( isinstance( self.branch.parent, single.weak_edge ) )
    self.assertEqual( self.branch.parent(), None )
    self.branch.parent = self.root.weakref()
    self.assertEqual( self.branch.parent(), self.root )
    
    self.assertTrue( isinstance( self.branch.suffix, single.weak_edge ) )
    self.assertEqual( self.branch.suffix(), None )
    self.branch.suffix = self.root.weakref()
    self.assertEqual( self.branch.suffix(), self.root )
    
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
    self.assertTrue( 1 in self.branch )
    self.assertTrue( 2 in self.branch )
    self.assertEqual( self.branch[1], self.leaf )
    self.assertEqual( self.branch[2], self.root )
    
    
  def test_const_branch(self):
    
    cb = single.const_edge.from_edge( self.branch )
    
    self.assertEqual( self.branch, cb )
    self.assertEqual( hash( self.branch ), hash( cb ) )
    self.assertNotEqual( single.edge.branch( start = 1, stop = 2 ), cb )
    
    wcb = cb.weakref()
    self.assertTrue( isinstance( wcb, single.const_weak_edge ) )
    self.assertFalse( wcb == cb )
    self.assertTrue( wcb != cb )
    self.assertTrue( isinstance( wcb(), single.const_edge ) )
    self.assertEqual( wcb(), cb )
    
    self.assertEqual( cb.start, 1 )
    self.assertEqual( cb.stop, 2 )
    self.assertRaises( RuntimeError, getattr, cb, "label" )
    
    self.assertTrue( isinstance( cb.parent, single.const_weak_edge ) )
    self.assertEqual( cb.parent(), None )
    self.branch.parent = self.root.weakref()
    self.assertEqual( cb.parent(), self.root )
    
    self.assertTrue( isinstance( cb.suffix, single.const_weak_edge ) )
    self.assertEqual( cb.suffix(), None )
    self.branch.suffix = self.root.weakref()
    self.assertEqual( cb.suffix(), self.root )
    
    self.assertFalse( cb.is_root() )
    self.assertFalse( cb.is_leaf() )
    
    self.assertTrue( cb.is_empty() )
    self.assertEqual( cb.keys(), [] )
    self.assertFalse( 1 in cb )
    self.assertTrue( 1 not in cb )
    self.assertRaises( KeyError, cb.__getitem__, 1 )
    
  
  def test_leaf(self):
    
    wl = self.leaf.weakref()
    self.assertTrue( isinstance( wl, single.weak_edge ) )
    self.assertFalse( wl == self.leaf )
    self.assertTrue( wl != self.leaf )
    self.assertTrue( isinstance( wl(), single.edge ) )
    self.assertEqual( wl(), self.leaf )
    
    self.assertEqual( self.leaf.start, 4 )
    self.leaf.start = 5
    self.assertEqual( self.leaf.start, 5 )
    ld = self.word.length()
    self.assertEqual( self.leaf.stop, ld() )
    self.word.append( 1 )
    self.word.append( 2 )
    self.assertEqual( self.leaf.stop, ld() )
    
    self.assertEqual( self.leaf.label, 6 )
    
    self.assertTrue( isinstance( self.leaf.parent, single.weak_edge ) )
    self.assertEqual( self.leaf.parent(), None )
    self.leaf.parent = self.root.weakref()
    self.assertEqual( self.leaf.parent(), self.root )
    
    self.assertRaises( RuntimeError, getattr, self.leaf, "suffix" )
    self.assertRaises( RuntimeError, setattr, self.leaf, "suffix", self.branch.weakref() )
    
    self.assertFalse( self.leaf.is_root() )
    self.assertTrue( self.leaf.is_leaf() )
    
    self.assertTrue( self.leaf.is_empty() )
    self.assertEqual( self.leaf.keys(), [] )
    self.assertFalse( 1 in self.leaf )
    self.assertTrue( 1 not in self.leaf )
    self.assertRaises( RuntimeError, self.leaf.__getitem__, 1 )
    self.assertRaises( RuntimeError, self.leaf.__setitem__, 1, self.root )
    
    
  def test_const_leaf(self):
    
    cl = single.const_edge.from_edge( self.leaf )
    
    self.assertEqual( self.leaf, cl )
    self.assertEqual( hash( self.leaf ), hash( cl ) )
    self.assertNotEqual(
      single.edge.leaf( start = 4, length = self.word.length(), label = 6 ),
      cl,
      )
    
    wcl = cl.weakref()
    self.assertTrue( isinstance( wcl, single.const_weak_edge ) )
    self.assertFalse( wcl == cl )
    self.assertTrue( wcl != cl )
    self.assertTrue( isinstance( wcl(), single.const_edge ) )
    self.assertEqual( wcl(), cl )
    
    self.assertEqual( cl.start, 4 )
    ld = self.word.length()
    self.assertEqual( cl.stop, ld() )
    self.word.append( 1 )
    self.word.append( 2 )
    self.assertEqual( cl.stop, ld() )
    
    self.assertEqual( cl.label, 6 )
    
    self.assertTrue( isinstance( cl.parent, single.const_weak_edge ) )
    self.assertEqual( cl.parent(), None )
    self.leaf.parent = self.root.weakref()
    self.assertEqual( cl.parent(), self.root )
    
    self.assertRaises( RuntimeError, getattr, cl, "suffix" )
    
    self.assertFalse( cl.is_root() )
    self.assertTrue( cl.is_leaf() )
    
    self.assertTrue( cl.is_empty() )
    self.assertEqual( cl.keys(), [] )
    self.assertFalse( 1 in cl )
    self.assertTrue( 1 not in cl )
    self.assertRaises( KeyError, cl.__getitem__, 1 )
    
    
class TestPreOrder(unittest.TestCase):
  
  def setUp(self):
    
    self.root = single.edge.root()
    self.word = single.word()
    length = self.word.length()
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
    res = list( single.preorder_iteration( root = root ) )
    self.assertEqual( res, [ root ] )
    self.assertTrue( isinstance( res[0], single.edge ) )
    
    cr = single.const_edge.from_edge( root )
    res = list( single.preorder_iteration( root = cr ) )
    self.assertEqual( res, [ cr ] )
    self.assertTrue( isinstance( res[0], single.const_edge ) )
    
    
  def test_single_branch(self):
    
    branch = single.edge.branch( start = 0, stop = 3 )
    res = list( single.preorder_iteration( root = branch ) )
    self.assertEqual( res, [ branch ] )
    self.assertTrue( isinstance( res[0], single.edge ) )
    
    cb = single.const_edge.from_edge( branch )
    res = list( single.preorder_iteration( root = cb ) )
    self.assertEqual( res, [ cb ] )
    self.assertTrue( isinstance( res[0], single.const_edge ) )
    
    
  def test_single_leaf(self):
    
    word = single.word()
    leaf = single.edge.leaf( start = 4, length = word.length(), label = 6 )
    res = list( single.preorder_iteration( root = leaf ) )
    self.assertEqual( res, [ leaf ] )
    self.assertTrue( isinstance( res[0], single.edge ) )
    
    cl = single.const_edge.from_edge( leaf )
    res = list( single.preorder_iteration( root = cl ) )
    self.assertEqual( res, [ cl ] )
    self.assertTrue( isinstance( res[0], single.const_edge ) )
    
    
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
    
    result = list( single.preorder_iteration( root = self.edge_named[ "r_b2" ] ) )
    self.assertEqual( len( result ), 3 )
    self.check_r_b2( result )
    self.assertEqual( result, [] )
  
  
  def test_r_b1_b1(self):
    
    result = list( single.preorder_iteration( root = self.edge_named[ "r_b1_b1" ] ) )
    self.assertEqual( len( result ), 3 )
    self.check_r_b1_b1( result )
    self.assertEqual( result, [] )
  
  
  def test_r_b1(self):
    
    result = list( single.preorder_iteration( root = self.edge_named[ "r_b1" ] ) )
    self.assertEqual( len( result ), 5 )
    self.check_r_b1( result )
    self.assertEqual( result, [] )
    
  
  def test_root(self):
    
    result = list( single.preorder_iteration( root = self.root ) )
    self.assertEqual( len( result ), 10 )
    self.assertTrue( all( isinstance( e, single.edge ) for e in result ) )
    self.check_root( result )
    self.assertEqual( result, [] )
    
    
  def test_const_root(self):
    
    cr = single.const_edge.from_edge( self.root )
    result = list( single.preorder_iteration( root = cr ) )
    self.assertEqual( len( result ), 10 )
    self.assertTrue( all( isinstance( e, single.const_edge ) for e in result ) )
    self.check_root( result )
    self.assertEqual( result, [] )
    

class TestPostOrder(unittest.TestCase):
  
  def setUp(self):
    
    self.root = single.edge.root()
    self.word = single.word()
    length = self.word.length()
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
    res = list( single.postorder_iteration( root = root ) )
    self.assertEqual( res, [ root ] )
    self.assertTrue( isinstance( res[0], single.edge ) )
    
    cr = single.const_edge.from_edge( root )
    res = list( single.postorder_iteration( root = cr ) )
    self.assertEqual( res, [ cr ] )
    self.assertTrue( isinstance( res[0], single.const_edge ) )
    
    
  def test_single_branch(self):
    
    branch = single.edge.branch( start = 0, stop = 3 )
    res = list( single.postorder_iteration( root = branch ) )
    self.assertEqual( res, [ branch ] )
    self.assertTrue( isinstance( res[0], single.edge ) )
    
    cb = single.const_edge.from_edge( branch )
    res = list( single.postorder_iteration( root = cb ) )
    self.assertEqual( res, [ cb ] )
    self.assertTrue( isinstance( res[0], single.const_edge ) )
    
    
  def test_single_leaf(self):
    
    word = single.word()
    leaf = single.edge.leaf( start = 4, length = word.length(), label = 6 )
    res = list( single.postorder_iteration( root = leaf ) )
    self.assertEqual( res, [ leaf ] )
    self.assertTrue( isinstance( res[0], single.edge ) )
    
    cl = single.const_edge.from_edge( leaf )
    res = list( single.postorder_iteration( root = cl ) )
    self.assertEqual( res, [ cl ] )
    self.assertTrue( isinstance( res[0], single.const_edge ) )
    
    
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
    
    result = list( single.postorder_iteration( root = self.edge_named[ "r_b2" ] ) )
    self.assertEqual( len( result ), 3 )
    self.check_r_b2( result )
    self.assertEqual( result, [] )
  
  
  def test_r_b1_b1(self):
    
    result = list( single.postorder_iteration( root = self.edge_named[ "r_b1_b1" ] ) )
    self.assertEqual( len( result ), 3 )
    self.check_r_b1_b1( result )
    self.assertEqual( result, [] )
  
  
  def test_r_b1(self):
    
    result = list( single.postorder_iteration( root = self.edge_named[ "r_b1" ] ) )
    self.assertEqual( len( result ), 5 )
    self.check_r_b1( result )
    self.assertEqual( result, [] )
    
  
  def test_root(self):
    
    result = list( single.postorder_iteration( root = self.root ) )
    self.assertEqual( len( result ), 10 )
    self.assertTrue( all( isinstance( e, single.edge ) for e in result ) )
    self.check_root( result )
    self.assertEqual( result, [] )
    
    
  def test_const_root(self):
    
    cr = single.const_edge.from_edge( self.root )
    result = list( single.postorder_iteration( root = cr ) )
    self.assertEqual( len( result ), 10 )
    self.assertTrue( all( isinstance( e, single.const_edge ) for e in result ) )
    self.check_root( result )
    self.assertEqual( result, [] )
    
    
class TestTree(unittest.TestCase):
  
  def test_tree(self):
    
    tree = single.tree()
    self.assertEqual( tree.in_construction, False )
    w = tree.word
    self.assertTrue( isinstance( w, single.word ) )
    self.assertEqual( w.length()(), 0 )
    r = tree.root
    self.assertTrue( isinstance( r, single.const_edge ) )
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
    
    from scitbx import suffixtree
    tree = single.tree()
    builder = single.ukkonen( tree )
    
    print "a:"
    builder.append( glyph = "a" )
    self.assertTrue( builder.is_valid )
    print suffixtree.tree_string( module = single, tree = tree )
    print
    
    print "an:"
    builder.append( glyph = "n" )
    self.assertTrue( builder.is_valid )
    print suffixtree.tree_string( module = single, tree = tree )
    print
    
    builder.detach()
    self.assertFalse( builder.is_attached )
    self.assertFalse( tree.in_construction )
    
    builder = single.ukkonen( tree )
    self.assertTrue( builder.is_valid )
    self.assertTrue( builder.is_attached )
    self.assertTrue( tree.in_construction )
    
    
    print "ana:"
    builder.append( glyph = "a" )
    self.assertFalse( builder.is_valid )
    self.assertRaises( builder.detach )
    print suffixtree.tree_string( module = single, tree = tree )
    print
    
    print "anan:"
    builder.append( glyph = "n" )
    self.assertFalse( builder.is_valid )
    self.assertRaises( builder.detach )
    print suffixtree.tree_string( module = single, tree = tree )
    print
    
    print "anana:"
    builder.append( glyph = "a" )
    self.assertFalse( builder.is_valid )
    self.assertRaises( builder.detach )
    print suffixtree.tree_string( module = single, tree = tree )
    print
    
    print "ananas:"
    builder.append( glyph = "s" )
    self.assertTrue( builder.is_valid )
    print suffixtree.tree_string( module = single, tree = tree )
    
    builder.detach()
    self.assertFalse( builder.is_attached )
    self.assertFalse( tree.in_construction )
    
    builder = single.ukkonen( tree )
    self.assertTrue( builder.is_valid )
    self.assertTrue( builder.is_attached )
    self.assertTrue( tree.in_construction )
    
    print "ananas$:"
    builder.append( glyph = "$" )
    self.assertTrue( builder.is_valid )
    print suffixtree.tree_string( module = single, tree = tree )
    
    builder.detach()
    
    
    

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


alltests = unittest.TestSuite(
  [
    suite_word,
    suite_edge,
    suite_preorder,
    suite_postorder,
    suite_tree,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )
    
