import scitbx.matrix

import operator
import math

__doc__ = """
Chirikjian, G. Stochastic models, information theory, and Lie groups.
Volume 2, pp 39-40.
Applied and Numerical Harmonic Analysis, Birkhauser
"""

class element(object):
  """
  An element of the Lie algebra so(3)
  """

  BASES = [
      scitbx.matrix.sqr( ( 0, 0, 0, 0, 0, -1, 0, 1, 0 ) ),
      scitbx.matrix.sqr( ( 0, 0, 1, 0, 0, 0, -1, 0, 0 ) ),
      scitbx.matrix.sqr( ( 0, -1, 0, 1, 0, 0, 0, 0, 0 ) ),
      ]

  def __init__(self, vector):

    self.vector = vector

  @property
  def matrix(self):

    return reduce(
      operator.add,
      [ x * e for ( x, e ) in zip( self.vector, self.BASES ) ],
      )

  def exponential(self):

    x2 = self.vector.length_sq()

    if 0 < x2:
      x = math.sqrt( x2 )
      a1 = math.sin( x ) / x
      a2 = ( 1.0 - math.cos( x ) ) / ( x * x )

    else:
      a1 = 1.0
      a2 = 0.5

    matrix = self.matrix

    return (
      scitbx.matrix.identity( 3 )
      + a1 * matrix
      + a2 * matrix * matrix
      )

  def scalar_multiplied(self, value):

    return self.__class__( vector = self.vector * value )

  def norm_sq(self):

    return 2 * self.vector.norm_sq()

  def __add__(self, other):

    return self.__class__( vector = self.vector + other.vector )

  def __sub__(self, other):

    return self.__class__( vector = self.vector - other.vector )

  def __mul__(self, other):

    if isinstance( other, ( float, int ) ):
      return self.scalar_multiplied( value = other )

    return NotImplemented

  def __rmul__(self, other):

    if isinstance( other, ( float, int ) ):
      return self.scalar_multiplied( value = other )

    return NotImplemented

  @classmethod
  def from_rotation_matrix(cls, matrix):

    import scitbx.math
    return cls.from_lie_matrix(
      matrix = scitbx.math.r3_rotation_matrix_logarithm( matrix ),
      )

  @classmethod
  def from_lie_matrix(cls, matrix):

    return cls(
      scitbx.matrix.col(
        (
          ( matrix[7] - matrix[5] ) / 2.0,
          ( matrix[2] - matrix[6] ) / 2.0,
          ( matrix[3] - matrix[1] ) / 2.0,
          )
        )
      )
