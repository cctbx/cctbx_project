from __future__ import absolute_import, division, print_function

from mmtbx.geometry import indexing

import unittest

OBJECTS = [
  ( ( 0, 0, 0 ), 1 ),
  ( ( 0.5, 0.5, 0.5 ), 1.5 ),
  ( ( 0.7, 0.7, 0.7 ), 1 ),
  ( ( 0, 0, 3 ), 1.7 ),
  ( ( 1, 1, 1 ), 2 ),
  ( ( 3, 3, 3 ), 1 ),
  ( ( 4, 4, 4 ), "foo" ),
  ( ( 5, 5, 5 ), 3 ),
  ( ( 6, 6, 6 ), 4 ),
  ( ( 7, 7, 7 ), 5 ),
  ]

class TestLinearIndexer(unittest.TestCase):

  def test_1(self):

    indexer = indexing.linear()

    for ( coords, obj ) in OBJECTS:
      indexer.add( object = obj, position = coords )

    self.assertEqual( len( indexer ), len( OBJECTS ) )

    for ( coords, obj ) in OBJECTS:
      closeby = indexer.close_to( centre = coords )
      self.assertEqual( len( closeby ), len( OBJECTS ) )
      self.assertTrue( obj in closeby )
      self.assertFalse( "bar" in closeby )


class TestHashIndexer(unittest.TestCase):

  def split_sorted(self, input_list):
    ints = [x for x in input_list if isinstance(x, (int, float))]
    strs = [x for x in input_list if isinstance(x, str)]
    return sorted(ints) + sorted(strs)

  def test_1(self):

    from mmtbx.geometry.shared_types import voxelizer
    v = voxelizer( base = ( 0, 0, 0 ), step = ( 1, 1, 1 ) )
    indexer = indexing.hash( voxelizer = v, margin = 1 )

    for ( coords, obj ) in OBJECTS:
      indexer.add( object = obj, position = coords )

    self.assertEqual( len( indexer ), len( OBJECTS ) )
    self.assertEqual( indexer.cubes(), 8 )

    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0, 0, 0 ) ) ),
      [ 1, 1, 1.5, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0.5, 0.5, 0.5 ) ) ),
      [ 1, 1, 1.5, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0, 0, 3 ) ) ),
      [ 1.7 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 1, 1, 1 ) ) ),
      [ 1, 1, 1.5, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 2, 2, 2 ) ) ),
      [ 1, 2 ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 3, 3, 3 ) ) ),
      [ 1, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 4, 4, 4 ) ) ),
      [ 1, 3, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 5, 5, 5 ) ) ),
      [ 3, 4, "foo" ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 6, 6, 6 ) ) ),
      [ 3, 4, 5 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 7, 7, 7 ) ) ),
      [ 4, 5 ],
      )


  def test_2(self):

    from mmtbx.geometry.shared_types import voxelizer
    v = voxelizer( base = ( 0, 0, 0 ), step = ( 1, 1, 1 ) )
    indexer = indexing.hash( voxelizer = v, margin = 2 )

    for ( coords, obj ) in OBJECTS:
      indexer.add( object = obj, position = coords )

    self.assertEqual( len( indexer ), len( OBJECTS ) )
    self.assertEqual( indexer.cubes(), 8 )

    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0, 0, 0 ) ) ),
      [ 1, 1, 1.5, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0.5, 0.5, 0.5 ) ) ),
      [ 1, 1, 1.5, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0, 0, 3 ) ) ),
      [ 1.7, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 1, 1, 1 ) ) ),
      [ 1, 1, 1, 1.5, 1.7, 2 ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 2, 2, 2 ) ) ),
      [ 1, 1, 1, 1.5, 1.7, 2, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 3, 3, 3 ) ) ),
      [ 1, 2, 3, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 4, 4, 4 ) ) ),
      [ 1, 3, 4, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 5, 5, 5 ) ) ),
      [ 1, 3, 4, 5, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 6, 6, 6 ) ) ),
      [ 3, 4, 5, "foo" ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 7, 7, 7 ) ) ),
      [ 3, 4, 5 ],
      )


  def test_3(self):

    from mmtbx.geometry.shared_types import voxelizer
    v = voxelizer( base = ( 0, 0, 0 ), step = ( 2, 2, 2 ) )
    indexer = indexing.hash( voxelizer = v, margin = 1 )

    for ( coords, obj ) in OBJECTS:
      indexer.add( object = obj, position = coords )

    self.assertEqual( len( indexer ), len( OBJECTS ) )
    self.assertEqual( indexer.cubes(), 5 )

    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0, 0, 0 ) ) ),
      [ 1, 1, 1, 1.5, 1.7, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0.5, 0.5, 0.5 ) ) ),
      [ 1, 1, 1, 1.5, 1.7, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 0, 0, 3 ) ) ),
      [ 1, 1, 1, 1.5, 1.7, 2 ],
      )
    self.assertEqual(
      sorted( indexer.close_to( centre = ( 1, 1, 1 ) ) ),
      [ 1, 1, 1, 1.5, 1.7, 2 ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 2, 2, 2 ) ) ),
      [ 1, 1, 1, 1.5, 1.7, 2, 3, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 3, 3, 3 ) ) ),
      [ 1, 1, 1, 1.5, 1.7, 2, 3, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 4, 4, 4 ) ) ),
      [ 1, 3, 4, 5, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 5, 5, 5 ) ) ),
      [ 1, 3, 4, 5, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 6, 6, 6 ) ) ),
      [ 3, 4, 5, "foo" ],
      )
    self.assertEqual(
      self.split_sorted( indexer.close_to( centre = ( 7, 7, 7 ) ) ),
      [ 3, 4, 5, "foo" ],
      )


suite_linear_indexer = unittest.TestLoader().loadTestsFromTestCase(
  TestLinearIndexer
  )
suite_hash_indexer = unittest.TestLoader().loadTestsFromTestCase(
  TestHashIndexer
  )


alltests = unittest.TestSuite(
  [
    suite_linear_indexer,
    suite_hash_indexer,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )
