from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from scitbx import matrix
from scitbx.linalg import eigensystem

import operator
import math
from functools import reduce
from six.moves import range
from six.moves import zip

# Exceptions
class MultipleSuperpositionError(Exception):
  """
  Module exception
  """


class NoConvergence(MultipleSuperpositionError):
  """
  Refinement starts to diverge
  """


class BadWeights(MultipleSuperpositionError):
  """
  Zero weights
  """


class BadSiteSets(MultipleSuperpositionError):
  """
  Invalid site sets
  """


# Conversion between rotation representations
def quaternion_to_rotmat(q):

  ( l, m, n, s ) = q

  return matrix.sqr( [
    l*l - m*m - n*n + s*s,
    2.0 * ( l*m - n*s ),
    2.0 * ( l*n + m*s ),
    2.0 * ( l*m + n*s ),
    -l*l + m*m -n*n + s*s,
    2.0 * ( m*n - l*s ),
    2.0 * ( l*n - m*s ),
    2.0 * ( m*n + l*s ),
    -l*l - m*m + n*n + s*s,
    ] )


def quaternion_to_quatmat(q):

  ( l, m, n, s ) = q.elems

  return matrix.sqr( [
     s, -n,  m,  l,
     n,  s, -l,  m,
    -m,  l,  s,  n,
    -l, -m, -n,  s,
    ] )


i_bar = matrix.sqr(
  [ -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1 ]
  )


# Functions to calculate rotational superposition
def best_rotation(moving, reference, weights):

  return top_quaternion(
    p = p_matrix_from( a_sites = moving, b_sites = reference, weights = weights )
    )


def top_quaternion(p):

  eigen = eigensystem.real_symmetric( p.as_flex_double_matrix() )
  return eigen.vectors()[0:4]


def p_matrix_from(a_sites, b_sites, weights):

  # Diamond, Acta Cryst A44, 211-216 (1988)
  assert len( a_sites ) == len( b_sites ) and len( a_sites ) == len( weights )

  ( xa, ya, za ) = component_arrays( a_sites )
  ( xb, yb, zb ) = component_arrays( b_sites )

  m = matrix.sqr(
    [
      flex.sum( weights * xa * xb ),
      flex.sum( weights * xa * yb ),
      flex.sum( weights * xa * zb ),
      flex.sum( weights * ya * xb ),
      flex.sum( weights * ya * yb ),
      flex.sum( weights * ya * zb ),
      flex.sum( weights * za * xb ),
      flex.sum( weights * za * yb ),
      flex.sum( weights * za * zb ),
      ]
    )

  q = m + m.transpose() - 2 * matrix.identity(3) * m.trace()

  # Application of Sum(j, k)[ epsilon( i, j, k ) * m( j, k ) ]
  # where i,j,k = [ 0, 1, 2 ] and epsilon is the Levi-Civitta symbol
  v = ( m[5] - m[7], m[6] - m[2], m[1] - m[3] )

  return matrix.sqr(
    [
      q[0], q[1], q[2], v[0],
      q[3], q[4], q[5], v[1],
      q[6], q[7], q[8], v[2],
      v[0], v[1], v[2], 0,
      ]
    )


# Various vector operations
def component_arrays(vec3_array):

  x, y, z  = zip( *vec3_array )
  return ( flex.double( x ), flex.double( y ), flex.double( z ) )


# Non-member functions
def iterate(superposition, convergence, tolerance=1.0E-8, iterations = 100):

  assert 0 < convergence
  residual = superposition.residual()
  count = 0
  tolerance *= superposition.set_count * superposition.site_count

  while True:
    superposition.full_iteration()
    new_residual = superposition.residual()
    diff = new_residual - residual

    if tolerance < diff:
      raise NoConvergence( "Divergence in refinement", count, new_residual)

    if abs( diff ) < convergence:
      break

    residual = new_residual
    count += 1

    if iterations < count:
      raise NoConvergence( "Iteration limit exceeded", count, new_residual)

  return ( count, new_residual )


def weighted_error_between(left, right, weights):

  differences = left - right

  return flex.mean_weighted( differences.dot( differences ), weights )


def transform_sites(sites, rotation, translation):

  return rotation.elems * sites + translation


class Matrix(object):
  """
  A matrix
  """

  def __init__(self, n, m):

    self.n = n
    self.m = m
    self.store = [ [ None ] * m for i in range( n ) ]


  def set(self, left, right, value):

    assert left < self.n and right < self.m
    self.store[ left ][ right ] = value


  def get(self, left, right):

    assert left < self.n and right < self.m
    return self.store[ left ][ right ]


  def rows(self):

    return self.store


  def off_diagonal_elements_in_rows(self):

    for ( index, row ) in enumerate( self.rows() ):
      yield row[ : index ] + row[ index + 1 : ]


class SymmetricMatrix(object):
  """
  A matrix that is symmetric and only half has to be calculated / stored
  """

  def __init__(self, dimension):

    assert 0 <= dimension
    self.d = dimension
    self.store = [ None ] * ( ( ( self.d + 1 ) * self.d ) // 2 )


  def set(self, left, right, value):

    self.store[ self.transform( left = left, right = right ) ] = value


  def get(self, left, right):

    return self.store[ self.transform( left = left, right = right ) ]


  def column(self, index):

    start1 = self.transform( left = index, right = 0 )
    slice1 = self.store[ start1 : start1 + index + 1 ]
    slice2 = [ self.get( i, index ) for i in range( index + 1, self.d ) ]
    return slice1 + slice2


  def transform(self, left, right):

    if left < right:
      ( left, right ) = ( right, left )

    return ( ( left + 1 ) * left ) // 2 + right


class ErrorOnRotation(object):
  """
  Superposition of two coordinate sets
  """

  def __init__(self, e, p):

    self.e = e
    self.p = p


  def with_rotated_reference(self, q):

    rho = quaternion_to_quatmat( q = q )
    return self.__class__( e = self.e, p = rho * self.p * rho.transpose() )


  def for_quaternion(self, q):

    square = self.e - 2 * ( q.transpose() * self.p * q )[0]

    return ( square if 0 <= square else 0 )


  def best_rotation(self):

    return matrix.col( top_quaternion( p = self.p ) )


  @classmethod
  def from_sites(cls, left, right, weights):

    assert len( left ) == len( right ) == len( weights )
    diffs = left - right
    e = flex.sum( diffs.dot( diffs ) * weights )
    p = p_matrix_from( a_sites = left, b_sites = right, weights = weights )
    return cls( e = e, p = p )


class PairwiseSuperposition(object):
  """
  Simple pairwise superposition
  """

  def __init__(self, reference_sites, moving_sites, weights = None):

    self.nsites = len( reference_sites )
    assert self.nsites == len( moving_sites )

    self.reference_sites = reference_sites
    self.moving_sites = moving_sites

    if weights is None:
      weights = flex.double( [ 1 ] * self.nsites )

    self.set_new_weights( weights = weights )


  def set_new_weights(self, weights):

    if self.nsites != len( weights ):
      raise BadWeights("Weights not equal to number of sites")

    norm = flex.sum( weights )

    if norm <= 0:
      raise BadWeights("Zero weights")

    refmean = self.reference_sites.mean_weighted( weights )
    movmean = self.moving_sites.mean_weighted( weights )
    eor = ErrorOnRotation.from_sites(
      left = self.moving_sites - movmean,
      right = self.reference_sites - refmean,
      weights = weights,
      )
    q = eor.best_rotation()
    self.rmsd = math.sqrt( eor.for_quaternion( q = q ) / norm )
    self.r = quaternion_to_rotmat( q = q )
    self.t = matrix.col( refmean ) - self.r * matrix.col( movmean )


  @property
  def rt(self):

    return matrix.rt( ( self.r, self.t ) )


  @property
  def positional_error_squares(self):

    transformed = transform_sites(
      sites = self.moving_sites,
      rotation = self.r,
      translation = self.t,
      )

    differences = transformed - self.reference_sites
    return differences.dot( differences )


class MultipleSuperposition(object):
  """
  Base class for multiple superpositions
  """

  def __init__(self, site_sets):

    self.set_count = len( site_sets )
    assert 2 <= self.set_count
    self.site_count = len( site_sets[0] )
    assert all( len( s ) == self.site_count for s in site_sets[1:] )
    self.site_sets = site_sets


  def rotations(self):

    return [ matrix.sqr( quaternion_to_rotmat( q ) ) for q in self.qs ]


  def transformations(self):

    rotations = self.rotations()
    return (
      rotations,
      [ -r * matrix.col( t ) for ( r, t ) in zip( rotations, self.centroids ) ],
      )


  def transformed_sites(self):

    ( rotations, translations ) = self.transformations()

    return [
      transform_sites( sites = sites, rotation = rot, translation = tra )
      for ( sites, rot, tra ) in zip( self.site_sets, rotations, translations )
      ]


  def rmsds(self, weights):

    transformed = self.transformed_sites()

    norm = 1.0 / ( self.set_count - 1 )
    return [
      norm * sum(
        [ weighted_error_between( left = l, right = r, weights = weights )
          for r in transformed ]
        )
      for l in transformed
      ]


  def unweighted_rmsd_between(self, left, right):

    ( rotations, translations ) = self.transformations()

    lsites = transform_sites(
      sites = self.site_sets[ left ],
      rotation = rotations[ left ],
      translation = translations[ left ],
      )
    rsites = transform_sites(
      sites = self.site_sets[ right ],
      rotation = rotations[ right ],
      translation = translations[ right ],
      )

    return math.sqrt(
      weighted_error_between(
        left = lsites,
        right = rsites,
        weights = flex.double( [ 1 ] * self.site_count ),
        )
      )


class DiamondAlgorithm(MultipleSuperposition):
  """
  Multiple superposition using Diamond's algorithm
  """

  def __init__(self, site_sets, weights):

    super( DiamondAlgorithm, self ).__init__( site_sets = site_sets )

    # Calculate initial rotations
    self.prepare_new_weights( weights = weights )
    self.calculate_pairwise_error_matrix()
    column = [ r[0] for r in self.pw_matrix.rows() ]
    self.qs = ( [ matrix.col( [ 0, 0, 0, 1 ] ) ]
      + [ eor.best_rotation() for eor in column[1:] ] )
    self.totals = [ None ] * self.set_count
    self.calculate_total_error_vector()


  # Commands
  def set_new_weights(self):

    self.calculate_pairwise_error_matrix()
    self.calculate_total_error_vector()


  def prepare_new_weights(self, weights):

    assert len( weights ) == self.site_count

    if flex.sum( weights ) <= 0:
      raise BadWeights("Zero weights")

    self.new_weight_data = weights


  def full_iteration(self):

    self.qs = [ matrix.col( top_quaternion( p = t.p ) ) for t in self.totals ]
    self.calculate_total_error_vector()


  # Queries
  def rmsd_between(self, left, right):

    if left == right:
      return 0.0

    eor = ( self.pw_matrix.get( left = left, right = right )
      .with_rotated_reference( q = self.qs[ right ] ) )
    square = eor.for_quaternion( q = self.qs[ left ] ) / self.norm

    return math.sqrt( square if 0 <= square else 0 )


  def weighted_rmsds(self):

    norm = 1.0 / ( self.set_count - 1 )
    return [ math.sqrt( r / self.norm * norm ) for r in self.residuals() ]


  def unweighted_rmsds(self):

    transformed = self.transformed_sites()
    weights = flex.double( [ 1 ] * self.site_count )

    norm = 1.0 / ( self.set_count - 1 )
    return [
      math.sqrt(
        norm * sum(
          [ weighted_error_between( left = l, right = r, weights = weights )
              for r in transformed ]
          )
        )
      for l in transformed
      ]


  def rmsd(self):

    return math.sqrt(
       self.residual() / ( self.set_count * ( self.set_count - 1 ) * self.norm )
      )


  def residual(self):

    return sum( self.residuals() )


  def residuals(self):

    return [ t.for_quaternion( q = q ) for ( q, t ) in zip( self.qs, self.totals ) ]


  def average_structure(self):

    return self.average_structure_from_transformed_sites(
      transformeds = self.transformed_sites(),
      )


  def distance_squares_from_average_structure(self):

    transformed = self.transformed_sites()
    average = self.average_structure_from_transformed_sites(
      transformeds = self.transformed_sites(),
      )
    differences = [ t - average for t in transformed ]
    return 1.0 / ( self.set_count - 1 ) * reduce(
      operator.add,
      [ diff.dot( diff ) for diff in differences ]
      )


  # Internal methods
  def average_structure_from_transformed_sites(self, transformeds):

    return 1.0 / self.set_count * reduce( operator.add, transformeds )


  def calculate_pairwise_error_matrix(self):

    self.norm = flex.sum( self.new_weight_data )
    assert 0 < self.norm

    # Calculate centroids
    self.centroids = flex.vec3_double(
      [ s.mean_weighted( self.new_weight_data ) for s in self.site_sets ]
      )
    origin_sets = [ c - o for ( c, o ) in zip( self.site_sets, self.centroids ) ]

    # Calculate starting pairwise squared error matrix
    self.pw_matrix = Matrix( n = self.set_count, m = self.set_count )

    for i in range( self.set_count ):
      for j in range( i + 1, self.set_count ):
        pw = ErrorOnRotation.from_sites(
          left = origin_sets[ i ],
          right = origin_sets[ j ],
          weights = self.new_weight_data,
          )
        self.pw_matrix.set( left = i, right = j, value = pw )
        self.pw_matrix.set(
          left = j,
          right = i,
          value = ErrorOnRotation( e = pw.e, p = i_bar * pw.p * i_bar ),
          )


  def calculate_total_error_vector(self):

    for ( index, row ) in enumerate( self.pw_matrix.rows() ):
      transformed_off_diagonal = [
        eor.with_rotated_reference( q = q )
        for ( i, ( eor, q ) ) in enumerate( zip( row, self.qs ) )
        if i != index
        ]
      sum_e = sum( [ eor.e for eor in transformed_off_diagonal ] )
      sum_p = reduce(
        operator.add,
        [ eor.p for eor in transformed_off_diagonal ],
        )
      self.totals[ index ] = ErrorOnRotation( e = sum_e, p = sum_p )


class WangSnoeyinkAlgorithm(MultipleSuperposition):
  """
  Multiple structure alignment with gaps using the algorithm from
  Wang & Snoeyink
  """

  def __init__(self, site_sets, selections, weights):

    super( WangSnoeyinkAlgorithm, self ).__init__( site_sets = site_sets )

    assert len( selections ) == self.set_count
    assert all( len( s ) == self.site_count for s in selections )

    self.selections = selections
    self.norms = reduce(
      operator.add,
      [ s.as_double() for s in self.selections ],
      )

    if not all( 2 <= v for v in self.norms ):
      raise BadSiteSets("Single positions in supplied site sets")

    self.qs = [ matrix.col( [ 0, 0, 0, 1 ] ) ] * self.set_count
    self.origins = [ ( 0, 0, 0 ) ] * self.set_count
    self.prepare_new_weights( weights = weights )
    self.set_new_weights()


  def set_new_weights(self):

    ( atomic_weights, sums ) = self.next_weight_data

    self.weights = atomic_weights
    self.sums = sums
    self.centroids = [
      s.mean_weighted( w ) for ( s, w ) in zip( self.site_sets, self.weights )
      ]
    self.calculate_superposition_errors()


  def prepare_new_weights(self, weights):

    assert len( weights ) == self.site_count
    atomic_weights = [
      weights.deep_copy().set_selected( ~s, 0 ) for s in self.selections
      ]
    sums = [ flex.sum( w ) for w in atomic_weights ]

    if not all( 0 < s for s in sums ):
      raise BadWeights("Zero weights")

    self.next_weight_data = ( atomic_weights, sums )


  def calculate_superposition_errors(self):

    average = self.average_structure()
    self.origins = [ average.mean_weighted( w ) for w in self.weights ]
    self.errors = [
      ErrorOnRotation.from_sites(
        left = sites - centroid,
        right = average - origin,
        weights = weights,
        )
      for ( sites, centroid, origin, weights ) in zip(
        self.site_sets,
        self.centroids,
        self.origins,
        self.weights,
        )
      ]


  def full_iteration(self):

    self.qs = [ eor.best_rotation() for eor in self.errors ]
    self.calculate_superposition_errors()


  def rmsd_between(self, left, right):

    ( rotations, translations ) = self.transformations()

    lsites = transform_sites(
      sites = self.site_sets[ left ],
      rotation = rotations[ left ],
      translation = translations[ left ],
      )
    rsites = transform_sites(
      sites = self.site_sets[ right ],
      rotation = rotations[ right ],
      translation = translations[ right ],
      )

    return math.sqrt(
      weighted_error_between(
        left = lsites,
        right = rsites,
        weights = self.weights[ left ] * self.weights[ right ],
        )
      )


  def weighted_rmsds(self):

    return self.rmsds( weights = self.weights )


  def unweighted_rmsds(self):

    return self.rmsds( weights = [ s.as_double() for s in self.selections ] )


  def rmsd(self):

    return math.sqrt( 2.0 / ( self.set_count - 1 ) * self.residual() )


  def residual(self):

    return sum( [ eor.for_quaternion( q = q ) / s for ( q, eor, s )
      in zip( self.qs, self.errors, self.sums ) ] )


  def transformations(self):

    ( rotations, translations ) = super( WangSnoeyinkAlgorithm, self ).transformations()

    return (
      rotations,
      [ t + matrix.col( o ) for ( t, o ) in zip( translations, self.origins ) ],
      )


  def average_structure(self):

    return self.average_structure_from_transformed_sites(
      transformeds = self.transformed_sites(),
      )


  def distance_squares_from_average_structure(self):

    transformeds = self.transformed_sites()
    average = self.average_structure_from_transformed_sites(
      transformeds = transformeds,
      )
    differences = [ s - average for s in transformeds ]

    return (
      reduce(
        operator.add,
        [ s.as_double() * d.dot( d ) for ( d, s )
          in zip( differences, self.selections ) ],
        )
      / ( self.norms - 1.0 )
      )


  def rmsds(self, weights):

    transformed = self.transformed_sites()
    norm = 1.0 / ( self.set_count - 1 )

    return [
      math.sqrt(
        norm * sum(
          [ weighted_error_between( left = l, right = r, weights = wl * wr )
            for ( r, wr ) in zip( transformed, weights ) ]
          )
        )
      for ( l, wl ) in zip( transformed, weights )
      ]


  # Internal methods
  def average_structure_from_transformed_sites(self, transformeds):

    return (
      reduce(
        operator.add,
        [ s.set_selected( ~sel, ( 0, 0, 0 ) ) for ( s, sel )
          in zip( transformeds, self.selections ) ]
        )
      / self.norms
      )

# Weighting
class UnitWeightScheme(object):
  """
  Unit weights
  """

  def for_difference_squares(self, squares):

    return flex.double( [ 1.0 ] * len( squares ) )


  def __str__(self):

    return "Unit weighting scheme"


class RobustResistantWeightScheme(object):
  """
  Robust-resistant weights
  """

  def __init__(self, critical_value_square):

    self.critical_value_square = float( critical_value_square )


  def for_difference_squares(self, squares):

    sqr_w = 1.0 - squares / self.critical_value_square
    sqr_w.set_selected( sqr_w < 0.0, 0.0 )
    return sqr_w * sqr_w


  def __str__(self):

    return "Robust-resistant weighting scheme (critical value: %s)" % (
      self.critical_value_square,
      )
