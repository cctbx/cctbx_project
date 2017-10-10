from __future__ import division

from builtins import zip
from mmtbx.geometry import primitive

import unittest

class GeometryTestCase(unittest.TestCase):

  def assertCoordinatesAlmostEqual(self, left, right, digits):

    self.assertEqual( len( left ), len( right ) )

    for ( e_left, e_right ) in zip( left, right ):
      self.assertAlmostEqual( e_left, e_right, digits )


class TestSphere(GeometryTestCase):

  def setUp(self):

    self.centre = ( 1.1, 2.2, 3.3 )
    self.radius = 4.4
    self.sphere = primitive.sphere( centre = self.centre, radius = self.radius )

  def test_radius(self):

    self.assertAlmostEqual( self.sphere.radius, self.radius, 7 )


  def test_radius_sq(self):

    self.assertAlmostEqual( self.sphere.radius_sq, self.radius ** 2, 7 )


  def test_centre(self):

    self.assertCoordinatesAlmostEqual( self.sphere.centre, self.centre, 7 )


class TestBox(GeometryTestCase):

  def test_from_corners(self):

    corner1 = ( 1.1, 2.4, 3.9 )
    corner2 = ( 1.3, 2.2, 3.3 )
    box = primitive.box.from_corners( corner1 = corner1, corner2 = corner2 )
    self.check( call = box.low, expected = ( 1.1, 2.2, 3.3 ) )
    self.check( call = box.high, expected = ( 1.3, 2.4, 3.9 ) )

    xs = [ 1, 2 ]
    ys = [ 12, 14 ]
    zs = [ 23, 27 ]

    import itertools

    for ( x, y, z ) in itertools.product( xs, ys, zs ):
      box = primitive.box.from_corners(
        corner1 = ( x, y, z ),
        corner2 = (
          self.select_other( array = xs, current = x ),
          self.select_other( array = ys, current = y ),
          self.select_other( array = zs, current = z ),
          )
        )
      self.check( call = box.low, expected = ( min(xs), min(ys), min(zs) ) )
      self.check( call = box.high, expected = ( max(xs), max(ys), max(zs) ) )


  def select_other(self, array, current):

    assert len( array ) == 2
    remaining = [ v for v in array if v != current ]
    assert len( remaining ) == 1
    return remaining[0]


  def test_around_sphere(self):

    centre = ( 1.1, 2.4, 3.9 )
    radius = 4.4
    box = primitive.box.around_sphere( centre = centre, radius = radius )
    self.check( call = box.low, expected = [ v - radius for v in centre ] )
    self.check( call = box.high, expected = [ v + radius for v in centre ] )


  def check(self, call, expected):

    self.assertCoordinatesAlmostEqual( call, expected, 7 )


class TestBSphere(GeometryTestCase):

  def setUp(self):

    self.centre = ( 1.1, 2.2, 3.3 )
    self.radius = 4.4
    self.bsphere = primitive.bsphere(
      centre = self.centre,
      radius = self.radius
      )

  def test_radius(self):

    self.assertAlmostEqual( self.bsphere.radius, self.radius, 7 )


  def test_radius_sq(self):

    self.assertAlmostEqual( self.bsphere.radius_sq, self.radius ** 2, 7 )


  def test_centre(self):

    self.assertCoordinatesAlmostEqual( self.bsphere.centre, self.centre, 7 )


  def test_corners(self):

    self.check(
      call = self.bsphere.low,
      expected = [ v - self.radius for v in self.centre ],
      )
    self.check(
      call = self.bsphere.high,
      expected = [ v + self.radius for v in self.centre ],
      )


  def check(self, call, expected):

    self.assertCoordinatesAlmostEqual( call, expected, 7 )


suite_sphere = unittest.TestLoader().loadTestsFromTestCase(
  TestSphere
  )
suite_box = unittest.TestLoader().loadTestsFromTestCase(
  TestBox
  )
suite_bsphere = unittest.TestLoader().loadTestsFromTestCase(
  TestBSphere
  )

alltests = unittest.TestSuite(
  [
    suite_sphere,
    suite_box,
    suite_bsphere,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

