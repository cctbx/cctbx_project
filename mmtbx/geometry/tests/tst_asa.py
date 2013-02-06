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

print ", ".join( [ a.name for a in ATOMS[:20 ] ] )

PARAMS = asa.CalcParams()
TABLE = PARAMS.van_der_waals_radii
RADII = [
  TABLE[ atom.determine_chemical_element_simple().strip().capitalize() ]
  for atom in ATOMS
  ]
SPHERES = [
  asa.sphere.create( centre = atom.xyz, radius = radius + PARAMS.probe )
  for ( atom, radius ) in zip( ATOMS, RADII )
  ]

class TestSphere(unittest.TestCase):

  def test_creation(self):
    
    s1 = asa.sphere.create( centre = ( 0, 0, 0 ), radius = 3 )
    s2 = asa.sphere.create( centre = ( 0, 0, 0 ), radius = 3 )
    
    self.assertTrue( s1.index != s2.index )
    
    s3 = asa.sphere( centre = ( 0, 0, 0 ), radius = 3, index = 0 )
    self.assertEqual( s3.index, 0 )
    
  def test_equality(self):
    
    s1 = asa.sphere.create( centre = ( 0, 0, 0 ), radius = 3 )
    s2 = asa.sphere.create( centre = ( 0, 0, 0 ), radius = 3 )
    self.assertTrue( s1 != s2 )
    
    s3 = asa.sphere( centre = ( 1, 1, 1 ), radius = 2, index = s1.index )
    self.assertTrue( s1 == s3 )
    self.assertTrue( s2 != s3 )
    

class TestLinearSpheresIndexer(unittest.TestCase):

  def setUp(self):

    self.indexer = asa.indexing.linear_spheres()


  def test_structure(self):

    for sphere in SPHERES:
      self.indexer.add( object = sphere )

    self.assertEqual( len( self.indexer ), len( SPHERES ) )

    for sphere in SPHERES:
      neighbours = list(
        self.indexer.overlapping_with( object = sphere )
        )
      indices = set( n.index for n in neighbours )
      self.assertTrue( sphere.index not in indices )


class TestPythagoreanChecker(unittest.TestCase):

  def setUp(self):

    self.checker = asa.accessibility.pythagorean_checker()


  def test_add_from_list(self):

    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )
    s = asa.sphere.create( centre = ( 1, 1, 1 ), radius = 0.1 )
    self.checker.add_from_list( neighbours = [ s ] )
    self.assertFalse( self.checker.is_selected( point = ( 1, 1, 1 ) ) )
    self.assertTrue( self.checker.is_selected( point = ( 1.2, 1, 1 ) ) )


  def test_add_from_range_1(self):

    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )
    
    s = asa.sphere.create( centre = ( 1, 1, 1 ), radius = 0.1 )
    indexer = asa.indexing.linear_spheres()
    indexer.add( object = s )

    neighbours = indexer.overlapping_with( object = s )
    self.assertEqual( len( neighbours ), 0 )
    self.checker.add_from_range( neighbours = neighbours )
    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )

    neighbours = indexer.overlapping_with(
      object = asa.sphere.create( centre = ( 1, 1, 1 ), radius = 0.1 ),
      )
    self.assertEqual( len( neighbours ), 1 )
    self.checker.add_from_range( neighbours = neighbours )
    self.assertFalse( self.checker.is_selected( point = ( 1, 1, 1 ) ) )


  def test_add_from_range_2(self):

    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )
    
    s = asa.sphere.create( centre = ( 1, 1, 1 ), radius = 0.1 )
    indexer = asa.indexing.linear_spheres()
    indexer.add( object = s )

    neighbours = indexer.overlapping_with(
      object = asa.sphere.create( centre = ( 2, 2, 2 ), radius = 0.1 ),
      )
    self.assertEqual( len( neighbours ), 0 )
    self.checker.add_from_range( neighbours = neighbours )
    self.assertTrue( self.checker.is_selected( point = ( 1, 1, 1 ) ) )

    neighbours = indexer.overlapping_with(
      object = asa.sphere.create( centre = ( 2, 2, 2 ), radius = 1.64 ),
      )
    self.assertEqual( len( neighbours ), 1 )
    self.checker.add_from_range( neighbours = neighbours )
    self.assertFalse( self.checker.is_selected( point = ( 1, 1, 1 ) ) )


class TestAccessibleSurfaceArea(unittest.TestCase):

  ACCESSIBLES = [
    431, 322,  43, 228,   0,   4,   1, 213,  41,  70,  10, 302, 140, 140, 0, 43,
    0,  0, 271,   1,  0,   0,   0, 141,  67, 144, 337, 125,   0,
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


  def test_asa_calculation_simple_linear(self):

    result = asa.calculate(
      atoms = ATOMS,
      params = PARAMS,
      indexer_selector = asa.get_linear_indexer_for,
      calculator = asa.SimpleSurfaceCalculator,
      )
    self.asa_result_check( result = result )
    
    
  def test_asa_calculation_aa_linear(self):

    result = asa.calculate(
      atoms = ATOMS,
      params = PARAMS,
      indexer_selector = asa.get_linear_indexer_for,
      calculator = asa.AltlocAveragedCalculator,
      )
    self.asa_result_check( result = result )
    
    
  def asa_result_check(self, result):

    self.assertEqual( len( result.values ), len( self.ACCESSIBLES ) )
    self.assertEqual( len( result.values ), len( SPHERES ) )
    
    self.assertEqual( len( result.points ), len( self.ACCESSIBLES ) )
    self.assertEqual( result.points, self.ACCESSIBLES )
    
    self.assertEqual( len( result.areas ), len( self.ACCESSIBLES ) )
    
    for ( count, sphere, area ) in zip( self.ACCESSIBLES, SPHERES, result.areas ):
      self.assertAlmostEqual(
        area,
        PARAMS.sampling.unit_area * count * sphere.radius_sq,
        7,
        )


class TestVoxelizer(unittest.TestCase):

  def setUp(self):

    self.voxelizer = asa.indexing.voxelizer(
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


suite_sphere = unittest.TestLoader().loadTestsFromTestCase(
  TestSphere
  )
suite_linear_spheres_indexer = unittest.TestLoader().loadTestsFromTestCase(
  TestLinearSpheresIndexer
  )
suite_pythagorean_checker = unittest.TestLoader().loadTestsFromTestCase(
  TestPythagoreanChecker
  )
suite_accessible_surface_area = unittest.TestLoader().loadTestsFromTestCase(
  TestAccessibleSurfaceArea
  )
suite_voxelizer = unittest.TestLoader().loadTestsFromTestCase(
  TestVoxelizer
  )


alltests = unittest.TestSuite(
  [
    suite_sphere,
    suite_linear_spheres_indexer,
    suite_pythagorean_checker,
    suite_accessible_surface_area,
    suite_voxelizer,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )

