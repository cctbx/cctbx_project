from scitbx.math import so3_lie_algebra

import scitbx.matrix
import scitbx.math
from scitbx.array_family import flex

import unittest
import math

class TestElement(unittest.TestCase):

  IDENTITY = scitbx.matrix.identity( 3 )
  SAMPLE_SIZE = 10
  SINGULARITY_GUARD = 0.001

  def run_tests_with(self, matrix):

    le = so3_lie_algebra.element.from_rotation_matrix( matrix = matrix )
    self.assertAlmostEqual(
      ( self.IDENTITY - matrix.transpose() * le.exponential() ).norm_sq(),
      0,
      7,
      )

  def test_random_matrices(self):

    for i in range( self.SAMPLE_SIZE ):
      axis = flex.random_double_point_on_sphere()
      angle = flex.random_double() * ( math.pi - self.SINGULARITY_GUARD )
      rotmat = scitbx.matrix.sqr(
        scitbx.math.r3_rotation_axis_and_angle_as_matrix( axis, angle )
        )
      self.run_tests_with( matrix = rotmat )

  def test_identity(self):

    self.run_tests_with( matrix = self.IDENTITY )


class TestAveraging(unittest.TestCase):

  MAX_ERROR = 0.1
  SAMPLE_SIZE = 10

  def test_angle(self):

    axis = scitbx.matrix.col( flex.random_double_point_on_sphere() )
    centre = flex.random_double() * math.pi

    errors = [ ( flex.random_double() - 0.5 )* self.MAX_ERROR
      for i in range( self.SAMPLE_SIZE ) ]

    average = scitbx.math.r3_rotation_average_rotation_via_lie_algebra(
      matrices = [
        scitbx.matrix.sqr(
          scitbx.math.r3_rotation_axis_and_angle_as_matrix( axis, centre + e )
          )
        for e in errors
        ]
      )
    angle = centre + sum( errors ) / len( errors )
    aaf = scitbx.math.r3_rotation_axis_and_angle_from_matrix( average )
    alignment = axis.dot( scitbx.matrix.col( aaf.axis ) )
    self.assertAlmostEqual( abs( alignment ), 1.0, 7 )
    self.assertAlmostEqual( aaf.angle() / alignment, angle, 7 )

  def test_axis(self):

    axis = (0.37394394059075464, 0.49642290523592875, 0.7834093619893614)
    angle = flex.random_double() * math.pi
    rotmat = scitbx.matrix.sqr(
      scitbx.math.r3_rotation_axis_and_angle_as_matrix( axis, angle )
      )

    missets = [
      scitbx.matrix.sqr(
        scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          flex.random_double_point_on_sphere(),
          ( flex.random_double() - 0.5 )* self.MAX_ERROR,
          )
        )
      for i in range( self.SAMPLE_SIZE )
      ]
    matrices = [ misset * rotmat for misset in missets ]
    aver_lie = scitbx.math.r3_rotation_average_rotation_via_lie_algebra(
      matrices = matrices,
      )
    aver_quat = scitbx.matrix.sqr(
      scitbx.math.r3_rotation_average_rotation_matrix_from_matrices(
        *matrices
        )
      )
    diff = aver_quat.transpose() * aver_lie
    self.assertAlmostEqual( ( diff - TestElement.IDENTITY ).norm_sq(), 0, 7 )

suite_element = unittest.TestLoader().loadTestsFromTestCase(
  TestElement
  )
suite_averaging = unittest.TestLoader().loadTestsFromTestCase(
  TestAveraging
  )

alltests = unittest.TestSuite(
  [
    suite_element,
    suite_averaging,
    ]
  )

if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

