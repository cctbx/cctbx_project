from __future__ import division

from mmtbx.geometry import asa

import iotbx.pdb
import libtbx.load_env

import os.path
import unittest

SOLVENT_RADIUS = 1.4

PDB = os.path.join(
  libtbx.env.under_dist( "mmtbx", "geometry" ),
  "tests",
  "1ahq_trunc.pdb"
  )
ROOT = iotbx.pdb.input( PDB ).construct_hierarchy()
ATOMS = ROOT.atoms()

PARAMS = asa.CalcParams()
SPHERES = asa.convert_to_spheres( atoms = ATOMS, params = PARAMS )

class TestLinearSpheresIndexer(unittest.TestCase):

  def setUp(self):

    self.indexer = asa.indexing.linear_spheres()


  def test_structure(self):

    for sphere in SPHERES:
      self.indexer.add( object = sphere )

    self.assertEqual( len( self.indexer ), len( ATOMS ) )

    for sphere in SPHERES:
      neighbours = list(
        self.indexer.overlapping_with( object = sphere )
        )
      indices = set( n.index for n in neighbours )
      self.assertTrue( sphere.index in indices )

      prefilter = asa.index_filter( index = sphere.index )
      prefiltered_neighbours = list(
        self.indexer.prefiltered_overlapping_with(
          object = sphere,
          prefilter = prefilter,
          )
        )
      self.assertEqual( len( prefiltered_neighbours ), len( neighbours ) - 1 )
      prefiltered_indices = set( n.index for n in prefiltered_neighbours )
      self.assertTrue( sphere.index not in prefiltered_indices )


class TestPythagoreanChecker(unittest.TestCase):

  def setUp(self):

    self.checker = asa.containment.pythagorean_checker()


  def test_add_from_list(self):

    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )
    s = asa.sphere( centre = ( 1, 1, 1 ), radius = 0.1, index = 0 )
    self.checker.add_from_list( neighbours = [ s ] )
    self.assertFalse( self.checker.is_selected( point = ( 1, 1, 1 ) ) )


  def test_add_from_range_1(self):

    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )
    s = asa.sphere( centre = ( 1, 1, 1 ), radius = 0.1, index = 0 )
    indexer = asa.indexing.linear_spheres()
    indexer.add( object = s )

    neighbours = indexer.prefiltered_overlapping_with(
      object = s,
      prefilter = asa.index_filter( index = 0 ),
      )
    self.assertEqual( len( neighbours ), 0 )
    self.checker.add_from_range( neighbours = neighbours )
    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )

    neighbours = indexer.prefiltered_overlapping_with(
      object = s,
      prefilter = asa.index_filter( index = 1 ),
      )
    self.assertEqual( len( neighbours ), 1 )
    self.checker.add_from_range( neighbours = neighbours )
    self.assertFalse( self.checker.is_selected( point = ( 1, 1, 1 ) ) )


  def test_add_from_range_2(self):

    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )
    s = asa.sphere( centre = ( 1, 1, 1 ), radius = 0.1, index = 0 )
    indexer = asa.indexing.linear_spheres()
    indexer.add( object = s )

    neighbours = indexer.overlapping_with(
      object = asa.sphere( centre = ( 2, 2, 2 ), radius = 0.1, index = 1 ),
      )
    self.assertEqual( len( neighbours ), 0 )
    self.checker.add_from_range( neighbours = neighbours )
    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )

    neighbours = indexer.overlapping_with(
      object = asa.sphere( centre = ( 2, 2, 2 ), radius = 1.64, index = 1 ),
      )
    self.assertEqual( len( neighbours ), 1 )
    self.checker.add_from_range( neighbours = neighbours )
    self.assertFalse( self.checker.is_selected( point = ( 1, 1, 1 ) ) )


class TestAccessibleSurfaceArea(unittest.TestCase):

  ACCESSIBLES = [
    431, 322,  43, 228,   0,   4,   1, 213,  41,  70,  10, 302, 140,   0,
    43,   0,   0, 271,   1,   0,   0,   0, 141,  67, 144, 337, 125,   0,
    0,   1,   0,  54,   0,   0,  25,   0,   0,  16,  33, 191,  76, 279,
    409, 119,   7,  21,   0,   0, 119,  61, 164, 108,  49,   0,   0,   9,
    16,  43,   3,   0,   0,  11,   2,   0,  23, 288, 169,   0,   4,  22,
    0,   0,  57,  39,  68, 322,  50,   4, 325, 117,   0,   0,   0,   0,
    0,   0,   0, 120,   1,   0,  67,  16, 323,   0,   0,   0,   0,  26,
    21,  43,  67, 176,  97, 115,   0,   0,   0,   6,   0,  60,   5, 247,
    0,   0,  75,   0,   0,   0,   0,   0,   0,  54,  11,  35, 117,   0,
    0,   0,   0,  23,   0,  22,  46, 222,   0,   0,  38,   0, 134,   5,
    41,   0,   0,   4,   0, 121, 101, 173,   0,   0,   0, 126,   6,   1,
    0,  84,   0,   0,  70,  17, 207,   0,   0,  12,   2, 115, 115,  48,
    25,  90, 366,  30,   0, 307, 315,   0,  56,   0,  59,  11,  88,  55,
    108,  46,   0, 207, 152,   0,   0,   0,  39,   0,   0,   0,   0, 128,
    0,   0,  43,   0,   0,   0,   0, 156,   0,  45,  59,  11,  56,  23,
    55,   0, 168,   0, 109, 291, 264,   0,  12,   0,   0, 134,   0,   0,
    11,   0,   0,   0, 118,   0,  51,   0,   0,   0, 122,   0,   0,   0,
    75,   0,  10,   0,   0,   0,   0, 186,   0,  20,   0,   0,   5, 222,
    44,  14,   3, 144,   4, 167,  47,  48,   0,  47,   5,   0, 118,   0,
    87,  37,  12,   0, 247, 272, 333,   0,   0,   0,  47,  40, 117,   0,
    391, 117,   0,  23,  20,  60,  53,   0,   0,  14,   0,  20, 259,  11,
    70,  16, 253, 463, 119,   0,  77,   3, 152, 252,  29,  12,  59,   0,
    0,   0, 114,  91,  51, 141,  55,   0, 334, 244,   0,  26,   0,  58,
    79,  44, 325,   0, 119,   0,   0,   0,   0,   4,   0,  37, 223, 299,
    0,   0,   0,   0,   0,   9, 115,  93,   0,   0,   0,   0,   0,  50,
    134,  60,   0,   0,   0,   1, 133,  26,  42,   0,  41,   0,   0,   0,
    81,   0,  31,   7,  38, 308,   0,   0,   0,   0,   0,  26,  21,   7,
    145, 256,  16,   0, 117, 306,   1,   0,   9, 218,  15,  72, 133, 127,
    0,  23,   0,   0,   0,   2, 153,   0,   0, 137,   0,   0,   1,  47,
    125, 237, 101,   0,   0,  23, 153, 181,  63, 337,  65,  68, 338, 206,
    11,  28,   0,   0,  75,   0,   0,   1,   1,   1,  53,  26, 346,   0,
    49,   3,  35,   0,   0, 183,  23, 172, 126,  31,  99,  35, 147, 130,
    226,   5,   5,   0,   0, 129, 201,  63, 209, 262,  46,   0,  19,   0,
    137,  52,  56, 292, 131,   0,   0,   0,   0, 145,  11,   2,   1,  21,
    88, 137, 276,   0,   0,  20,   0, 228,   0,  73,  82, 183,   9,   0,
    285, 299, 326,   0,  62, 149, 407, 174,  55,   0, 178,  27, 431, 427,
    383, 363,   0,   6,  36,   0,   0,   0,   0,   2,   0,   0, 210, 187,
    231,  73, 213, 246,   0,   8,   3,   0,   0, 135,   0,   0,   0,   0,
    0,   0,   0,   0,   0, 116,   0,   0,   0,   0,   0, 124,  53,  90,
    84, 219, 410, 960, 957, 960, 960, 960, 948, 586,   0, 960, 958, 960,
    635, 960, 542, 726, 960, 511, 828, 960, 875,
    ]

  def check_indexer(self, indexer, type):

    self.assertEqual( indexer.factory, type )
    self.assertTrue( isinstance( indexer.regular, type ) )
    self.assertEqual( len( indexer.regular ), 0 )


  def test_get_linear_indexer(self):

    indexer = asa.get_linear_indexer_for( atoms = ATOMS )
    self.assertTrue( isinstance( indexer, asa.Indexer ) )
    self.check_indexer( indexer = indexer, type = asa.indexing.linear_spheres )


  def test_get_optimal_indexer(self):

    indexer = asa.get_optimal_indexer_for( atoms = ATOMS )
    self.assertTrue( isinstance( indexer, asa.Indexer ) )
    self.check_indexer( indexer = indexer, type = asa.indexing.linear_spheres )


  def test_single(self):

    indexer = asa.indexing.linear_spheres()

    for sphere in SPHERES:
      indexer.add( object = sphere )

    results = [
      asa.single( indexer = indexer, sphere = sphere, sampling = PARAMS.sampling )
      for sphere in SPHERES
      ]

    self.assertEqual( results, self.ACCESSIBLES )


  def test_asa_calculation(self):

    results = asa.calculate( atoms = ATOMS, params = PARAMS )

    self.assertEqual( len( results ), len( self.ACCESSIBLES ) )
    self.assertEqual( len( results ), len( SPHERES ) )

    import math
    unit = 4 * math.pi / PARAMS.sampling.count

    for ( index, ( accessible, sphere, res ) ) in enumerate( zip(
    self.ACCESSIBLES, SPHERES, results ) ):
      self.assertAlmostEqual( res, unit * accessible * sphere.radius_sq, 7, "%s does not match" % index )



suite_linear_spheres_indexer = unittest.TestLoader().loadTestsFromTestCase(
  TestLinearSpheresIndexer
  )
suite_pythagorean_checker = unittest.TestLoader().loadTestsFromTestCase(
  TestPythagoreanChecker
  )
suite_accessible_surface_area = unittest.TestLoader().loadTestsFromTestCase(
  TestAccessibleSurfaceArea
  )


alltests = unittest.TestSuite(
  [
    suite_linear_spheres_indexer,
    suite_pythagorean_checker,
    suite_accessible_surface_area,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )

