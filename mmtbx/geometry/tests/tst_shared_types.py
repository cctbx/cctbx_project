from __future__ import absolute_import, division, print_function

from mmtbx.geometry import shared_types

import unittest

class TestVoxelizer(unittest.TestCase):

  def setUp(self):

    self.voxelizer = shared_types.voxelizer(
      base = ( 100, 200, 300 ),
      step = ( 2.0, 3.0, 4.0 ),
      )


  def test_1(self):

    self.assertEqual(
      self.voxelizer( vector = ( 100, 200, 300 ) ),
      ( 0, 0, 0 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 101, 201, 301 ) ),
      ( 0, 0, 0 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 101.99, 202, 302 ) ),
      ( 0, 0, 0 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 102, 202.99, 303 ) ),
      ( 1, 0, 0 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 103, 203, 303.99 ) ),
      ( 1, 1, 0 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 104, 204, 304 ) ),
      ( 2, 1, 1 ),
      )

  def test_2(self):

    self.assertEqual(
      self.voxelizer( vector = ( 99.99, 199.99, 299.99 ) ),
      ( -1, -1, -1 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 98.00, 197.00, 296.00 ) ),
      ( -1, -1, -1 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 97.99, 197.00, 296.00 ) ),
      ( -2, -1, -1 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 96.00, 196.99, 296.00 ) ),
      ( -2, -2, -1 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 96.00, 194.00, 295.99 ) ),
      ( -2, -2, -2 ),
      )
    self.assertEqual(
      self.voxelizer( vector = ( 96.00, 194.00, 292.00 ) ),
      ( -2, -2, -2 ),
      )

suite_voxelizer = unittest.TestLoader().loadTestsFromTestCase(
  TestVoxelizer
  )


alltests = unittest.TestSuite(
  [
    suite_voxelizer,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )

