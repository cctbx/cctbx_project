from __future__ import absolute_import, division, print_function

from mmtbx.geometry import sphere_surface_sampling

import unittest
from six.moves import zip

class TestGoldenSpiral(unittest.TestCase):

  def setUp(self):

    self.sampling = sphere_surface_sampling.golden_spiral( count = 5 )
    self.count = 5
    self.points = [
      ( 0.5999999999999999, -0.8, 0.0 ),
      ( -0.6758097397797129, -0.39999999999999997, 0.6190970809322855 ),
      ( 0.08742572471695988, 5.551115123125783e-17, -0.9961710408648278 ),
      ( 0.5576434272376701, 0.4000000000000002, 0.7273471028736039 ),
      ( -0.590828091189257, 0.8, -0.10450917022758696 )
      ]


  def test_count(self):

    self.assertEqual( self.sampling.count, self.count )


  def test_unit_area(self):

    import math
    self.assertAlmostEqual(
      self.sampling.unit_area,
      4 * math.pi / self.count,
      8,
      )


  def test_points(self):

    got = list( self.sampling.points )

    self.assertEqual( len( got ), len( self.points ) )

    for ( left, right ) in zip( got, self.points ):
      self.assertIterablesAlmostEqual( left, right, 7 )


  """
  def test_transformed(self):

    ( centre_x, centre_y, centre_z ) = ( 0.1, 2.5, -3.6 )
    radius = 3.5

    got = list(
      self.sampling.transformed(
        centre = ( centre_x, centre_y, centre_z ),
        radius = radius,
        )
      )

    self.assertEqual( len( got ), len( self.points ) )

    for ( point, ( x, y, z ) )  in zip( got, self.points ):
      self.assertIterablesAlmostEqual(
        point,
        ( radius * x + centre_x, radius * y + centre_y, radius * z + centre_z ),
        7,
        )
    """

  def assertIterablesAlmostEqual(self, left, right, digits):

    self.assertEqual( len( left ), len( right ) )

    for ( e_left, e_right ) in zip( left, right ):
      self.assertAlmostEqual( e_left, e_right, digits )



suite_golden_spiral = unittest.TestLoader().loadTestsFromTestCase(
    TestGoldenSpiral
    )

alltests = unittest.TestSuite(
  [
      suite_golden_spiral,
      ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

