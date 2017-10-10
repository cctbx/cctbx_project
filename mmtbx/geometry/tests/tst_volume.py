from __future__ import division

from builtins import zip
from mmtbx.geometry import volume
from mmtbx.geometry.primitive import sphere, box, bsphere

import math
import unittest

RADIUS = 0.5

FACE_CENTRES = [
  ( -0.5,  0.0,  0.0 ),
  (  0.5,  0.0,  0.0 ),
  (  0.0, -0.5,  0.0 ),
  (  0.0,  0.5,  0.0 ),
  (  0.0,  0.0, -0.5 ),
  (  0.0,  0.0,  0.5 ),
  ]
CORNERS = [
  ( -0.5, -0.5, -0.5 ),
  (  0.5, -0.5, -0.5 ),
  ( -0.5,  0.5, -0.5 ),
  (  0.5,  0.5, -0.5 ),
  ( -0.5, -0.5,  0.5 ),
  (  0.5, -0.5,  0.5 ),
  ( -0.5,  0.5,  0.5 ),
  (  0.5,  0.5,  0.5 ),
  ]

SPHERE_DIAGONAL = RADIUS / math.sqrt( 3 )

SPHERE_CORNERS = [
  ( -SPHERE_DIAGONAL, -SPHERE_DIAGONAL, -SPHERE_DIAGONAL ),
  (  SPHERE_DIAGONAL, -SPHERE_DIAGONAL, -SPHERE_DIAGONAL ),
  ( -SPHERE_DIAGONAL,  SPHERE_DIAGONAL, -SPHERE_DIAGONAL ),
  (  SPHERE_DIAGONAL,  SPHERE_DIAGONAL, -SPHERE_DIAGONAL ),
  ( -SPHERE_DIAGONAL, -SPHERE_DIAGONAL,  SPHERE_DIAGONAL ),
  (  SPHERE_DIAGONAL, -SPHERE_DIAGONAL,  SPHERE_DIAGONAL ),
  ( -SPHERE_DIAGONAL,  SPHERE_DIAGONAL,  SPHERE_DIAGONAL ),
  (  SPHERE_DIAGONAL,  SPHERE_DIAGONAL,  SPHERE_DIAGONAL ),
  ]

class TestPointInOriginSphere(unittest.TestCase):

  def setUp(self):

    self.tester = volume.containment.point_in_origin_sphere()


  def test_origin(self):

    self.assertTrue( self.tester( point = ( 0, 0, 0 ), radius_sq = 0.001 ) )


  def test_face_centres(self):

    for p in FACE_CENTRES:
      self.assertFalse( self.tester( point = p, radius_sq = 0.25 ) )
      self.assertTrue( self.tester( point = p, radius_sq = 0.251 ) )


  def test_corners(self):

    for p in CORNERS:
      self.assertFalse( self.tester( point = p, radius_sq = 0.25 ) )
      self.assertFalse( self.tester( point = p, radius_sq = 0.251 ) )


  def test_sphere_corners(self):

    for p in SPHERE_CORNERS:
      self.assertFalse( self.tester( point = p, radius_sq = 0.25 ) )
      self.assertTrue( self.tester( point = p, radius_sq = 0.251 ) )


class TestPointInOriginDiamond(unittest.TestCase):

  def setUp(self):

    self.tester = volume.containment.point_in_origin_diamond()


  def test_origin(self):

    self.assertTrue( self.tester( point = ( 0, 0, 0 ), radius = 0.001 ) )


  def test_face_centres(self):

    for p in FACE_CENTRES:
      self.assertFalse( self.tester( point = p, radius = 0.5 ) )
      self.assertTrue( self.tester( point = p, radius = 0.51 ) )


  def test_corners(self):

    for p in CORNERS:
      self.assertFalse( self.tester( point = p, radius = 0.5 ) )
      self.assertFalse( self.tester( point = p, radius = 0.51 ) )


  def test_sphere_corners(self):

    for p in SPHERE_CORNERS:
      self.assertFalse( self.tester( point = p, radius = 0.5 ) )
      self.assertFalse( self.tester( point = p, radius = 0.51 ) )


class TestPointInOriginCube(unittest.TestCase):

  def setUp(self):

    self.tester = volume.containment.point_in_origin_cube()


  def test_origin(self):

    self.assertTrue( self.tester( point = ( 0, 0, 0 ), radius = 0.001 ) )


  def test_face_centres(self):

    for p in FACE_CENTRES:
      self.assertFalse( self.tester( point = p, radius = 0.5 ) )
      self.assertTrue( self.tester( point = p, radius = 0.51 ) )


  def test_corners(self):

    for p in CORNERS:
      self.assertFalse( self.tester( point = p, radius = 0.5 ) )
      self.assertTrue( self.tester( point = p, radius = 0.51 ) )


  def test_sphere_corners(self):

    for p in SPHERE_CORNERS:
      self.assertTrue( self.tester( point = p, radius = 0.5 ) )
      self.assertTrue( self.tester( point = p, radius = 0.51 ) )


class TestBetweenSpheres(unittest.TestCase):

  def setUp(self):

    self.tester = volume.overlap.between_spheres()


  def check_with(self, sphere_type):

    b1 = sphere_type( centre = ( 0, 1, 2 ), radius = 0.01 )
    self.assertTrue( self.tester( left = b1, right = b1 ) )
    self.assertTrue(
      self.tester(
        left = sphere_type( centre = ( 1, 2, 3 ), radius = 1 ),
        right = sphere_type( centre = ( 2, 3, 4 ), radius = 1 ),
        )
      )
    self.assertFalse(
      self.tester(
        left = sphere_type( centre = ( 1, 2, 3 ), radius = 0.5 ),
        right = sphere_type( centre = ( 2, 3, 4 ), radius = 0.5 ),
        )
      )
    self.assertFalse(
      self.tester(
        left = sphere_type( centre = ( 1, 2, 3 ), radius = 0.5 ),
        right = sphere_type( centre = ( 2, 2, 3 ), radius = 0.5 ),
        )
      )
    self.assertTrue(
      self.tester(
        left = sphere_type( centre = ( 1, 2, 3 ), radius = 0.51 ),
        right = sphere_type( centre = ( 2, 2, 3 ), radius = 0.5 ),
        )
      )

  def test_with_sphere(self):

    self.check_with( sphere_type = sphere )


  def test_with_bsphere(self):

    self.check_with( sphere_type = bsphere )


class TestBetweenBoxes(unittest.TestCase):

  def setUp(self):

    self.tester = volume.overlap.between_boxes()


  def test_with_bsphere(self):

    b1 = bsphere( centre = ( 0, 1, 2 ), radius = 0.01 )
    self.assertTrue( self.tester( left = b1, right = b1 ) )
    self.assertTrue(
      self.tester(
        left = bsphere( centre = ( 1, 2, 3 ), radius = 1 ),
        right = bsphere( centre = ( 2, 3, 4 ), radius = 1 ),
        )
      )
    self.assertFalse(
      self.tester(
        left = bsphere( centre = ( 1, 2, 3 ), radius = 0.5 ),
        right = bsphere( centre = ( 2, 3, 4 ), radius = 0.5 ),
        )
      )
    self.assertTrue(
      self.tester(
        left = bsphere( centre = ( 1, 2, 3 ), radius = 0.51 ),
        right = bsphere( centre = ( 2, 3, 4 ), radius = 0.5 ),
        )
      )
    self.assertFalse(
      self.tester(
        left = bsphere( centre = ( 1, 2, 3 ), radius = 0.5 ),
        right = bsphere( centre = ( 2, 2, 3 ), radius = 0.5 ),
        )
      )
    self.assertTrue(
      self.tester(
        left = bsphere( centre = ( 1, 2, 3 ), radius = 0.51 ),
        right = bsphere( centre = ( 2, 2, 3 ), radius = 0.5 ),
        )
      )


  def test_with_box(self):

    b = box.from_corners( corner1 = ( -5, -6, 3 ), corner2 = ( 5, -4, 6 ) )
    self.assertTrue( self.tester( left = b, right = b ) )

    directions = [
      ( 1, 0, 0 ), ( 0, 1, 0 ), ( 0, 0, 1 ),
      ( 1, 1, 0 ), ( 1, 0, 1 ), ( 0, 1, 1 ),
      ( 1, -1, 0 ), ( 1, 0, -1 ), ( 0, 1, -1 ),
      ( 1, 1, 1 ), ( 1, 1, -1 ),
      ( 1, -1, 1 ), ( -1, 1, 1 ),
      ]

    for direction in directions:
      self.check_boundary( trial = b, direction = direction )

    self.assertTrue(
      self.tester(
        left = b,
        right = box.from_corners(
          corner1 = ( -4.9, -5.9, 3.1 ),
          corner2 = ( -4.8, -5.8, 3.2 ),
          ),
        )
      )

    self.assertTrue(
      self.tester(
        left = b,
        right = box.from_corners(
          corner1 = ( -0.1, -4.1, 0 ),
          corner2 = ( 0.1, -3.9, 9 ),
          ),
        )
      )


  def check_boundary(self, trial, direction):

    self.assertFalse(
      self.tester(
        left = trial,
        right = box.from_corners(
          corner1 = self.shift( trial.low, self.scale( direction, -2 ) ),
          corner2 = self.shift(
            self.shift( trial.low, ( 1, 1, 1 ) ),
            self.scale( direction, -1.1 ),
            ),
          ),
        )
      )
    self.assertTrue(
      self.tester(
        left = trial,
        right = box.from_corners(
          corner1 = self.shift( trial.low, self.scale( direction, -2 ) ),
          corner2 = self.shift(
            self.shift( trial.low, ( 1, 1, 1 ) ),
            self.scale( direction, -0.9 ),
            ),
          ),
        )
      )
    self.assertFalse(
      self.tester(
        left = trial,
        right = box.from_corners(
          corner1 = self.shift( trial.high, self.scale( direction, 0.1 ) ),
          corner2 = self.shift(
            self.shift( trial.high, ( -1, -1, -1 ) ),
            self.scale( direction, 2 ),
            ),
          ),
        )
      )
    right = box.from_corners(
      corner1 = self.shift( trial.high, self.scale( direction, -0.1 ) ),
      corner2 = self.shift(
        self.shift( trial.high, ( -1, -1, -1 ) ),
        self.scale( direction, 2 ),
        ),
      )
    self.assertTrue(
      self.tester(
        left = trial,
        right = box.from_corners(
          corner1 = self.shift( trial.high, self.scale( direction, -0.1 ) ),
          corner2 = self.shift(
            self.shift( trial.high, ( -1, -1, -1 ) ),
            self.scale( direction, 2 ),
            ),
          ),
        )
      )


  def scale(self, vector, scale):

    return tuple( [ scale * v for v in vector ] )


  def shift(self, vector, shift):

    assert len( vector ) == len( shift )
    return tuple( [ l + r for ( l, r ) in zip( vector, shift ) ] )


suite_point_in_origin_sphere = unittest.TestLoader().loadTestsFromTestCase(
  TestPointInOriginSphere
  )
suite_point_in_origin_diamond = unittest.TestLoader().loadTestsFromTestCase(
  TestPointInOriginDiamond
  )
suite_point_in_origin_cube = unittest.TestLoader().loadTestsFromTestCase(
  TestPointInOriginCube
  )
suite_between_spheres = unittest.TestLoader().loadTestsFromTestCase(
  TestBetweenSpheres
  )
suite_between_boxes = unittest.TestLoader().loadTestsFromTestCase(
  TestBetweenBoxes
  )

alltests = unittest.TestSuite(
  [
    suite_point_in_origin_sphere,
    suite_point_in_origin_diamond,
    suite_point_in_origin_cube,
    suite_between_spheres,
    suite_between_boxes,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

