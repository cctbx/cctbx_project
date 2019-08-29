from __future__ import absolute_import, division, print_function

from mmtbx.geometry import asa
from mmtbx.geometry import altloc

import iotbx.pdb
import libtbx.load_env

import os.path
import unittest
from six.moves import zip
from six.moves import range

SOLVENT_RADIUS = 1.4

PDB = os.path.join(
  libtbx.env.under_dist( "mmtbx", "geometry" ),
  "tests",
  "1ahq_trunc.pdb"
  )
ROOT = iotbx.pdb.input( PDB ).construct_hierarchy( sort_atoms = False )
ATOMS = ROOT.atoms()

from mmtbx.geometry import sphere_surface_sampling
SAMPLING_POINTS = sphere_surface_sampling.golden_spiral( count = 960 )

from cctbx.eltbx import van_der_waals_radii
TABLE = van_der_waals_radii.vdw.table

PROBE = 1.4

RADII = [
  TABLE[ atom.determine_chemical_element_simple().strip().capitalize() ]
  for atom in ATOMS
  ]
SPHERES = [
  asa.sphere( centre = atom.xyz, radius = radius + PROBE, index = index )
  for ( index, ( atom, radius ) ) in enumerate( zip( ATOMS, RADII ) )
  ]

DESCRIPTIONS = [
  altloc.Description(
    data = s,
    coordinates = s.centre,
    altid = altloc.altid_for( atom = atom ),
    )
    for ( s, atom ) in zip( SPHERES, ATOMS )
  ]

class TestSphere(unittest.TestCase):

  def test_creation(self):

    s1 = asa.sphere( centre = ( 0, 0, 0 ), radius = 3, index = 1 )
    s2 = asa.sphere( centre = ( 0, 0, 0 ), radius = 3, index = 2 )

    self.assertTrue( s1.index != s2.index )


  def test_low(self):

    self.assertIterablesAlmostEqual(
      asa.sphere( centre = ( 1, 2, 3 ), radius = 4, index = 0 ).low,
      ( -3, -2, -1 ),
      7,
      )


  def test_high(self):

    self.assertIterablesAlmostEqual(
      asa.sphere( centre = ( 1, 2, 3 ), radius = 4, index = 0 ).high,
      ( 5, 6, 7 ),
      7,
      )


  def assertIterablesAlmostEqual(self, left, right, digits):

    self.assertEqual( len( left ), len( right ) )

    for ( e_left, e_right ) in zip( left, right ):
      self.assertAlmostEqual( e_left, e_right, digits )


class TestTransformation(unittest.TestCase):

  def setUp(self):

    ( self.x, self.y, self.z ) = ( 0.1, 2.5, -3.6 )
    self.radius = 3.5

    self.transformation = asa.accessibility.transformation(
      centre = ( self.x, self.y, self.z ),
      radius = self.radius,
      )


  def transform(self, site):

    ( x, y, z ) = site
    return ( self.radius * x + self.x, self.radius * y + self.y, self.radius * z + self.z )


  def test_single_site(self):

    site = ( 1, 2, 3 )

    self.assertIterablesAlmostEqual(
      self.transformation( point = site ),
      self.transform( site = site ),
      7,
      )


  def test_range(self):

    from mmtbx.geometry import sphere_surface_sampling
    sampling = sphere_surface_sampling.golden_spiral( count = 960 )
    transformed = asa.accessibility.transform(
      range = sampling.points,
      transformation = self.transformation,
      )
    self.assertEqual( len( sampling.points ), len( transformed ) )

    for ( o, t ) in zip( sampling.points, transformed ):
      self.assertIterablesAlmostEqual( t, self.transform( site = o ), 7 )


  def get_transformed_points(self):

    from mmtbx.geometry import sphere_surface_sampling
    sampling = sphere_surface_sampling.golden_spiral( count = 5 )
    return (
      list( sampling.points ),
      asa.accessibility.transform(
        range = sampling.points,
        transformation = self.transformation,
        )
      )

  def test_transformed_range_lifetime(self):

    ( points, transformed ) = self.get_transformed_points()
    for ( o, t ) in zip( points, transformed ):
      self.assertIterablesAlmostEqual( t, self.transform( site = o ), 7 )


  def assertIterablesAlmostEqual(self, left, right, digits):

    self.assertEqual( len( left ), len( right ) )

    for ( e_left, e_right ) in zip( left, right ):
      self.assertAlmostEqual( e_left, e_right, digits )


class TestLinearSpheresIndexer(unittest.TestCase):

  def setUp(self):

    self.indexer = asa.indexing.linear_spheres()


  def test_structure(self):

    for sphere in SPHERES:
      self.indexer.add( object = sphere, position = sphere.centre )

    self.assertEqual( len( self.indexer ), len( SPHERES ) )

    for sphere in SPHERES:
      closeby = self.indexer.close_to( centre = sphere.centre )
      self.assertEqual( len( closeby ), len( SPHERES ) )


class TestHashSpheresIndexer(unittest.TestCase):

  def setUp(self):

    self.indexer = asa.indexing.hash_spheres(
      voxelizer = asa.get_voxelizer_for( descriptions = DESCRIPTIONS ),
      margin = 1,
      )


  def test_structure(self):

    for sphere in SPHERES:
      self.indexer.add( object = sphere, position = sphere.centre )

    self.assertEqual( len( self.indexer ), len( SPHERES ) )


class TestPredicate(unittest.TestCase):

  def setUp(self):

    self.sphere = asa.sphere( centre = ( 1, 1, 1 ), radius = 2, index = 0 )
    self.predicate = asa.accessibility.overlap_equality_predicate( object = self.sphere )
    self.overlap1 = asa.sphere( centre = ( 1, 1, 1 ), radius = 0.1, index = 1 )
    self.overlap2 = asa.sphere( centre = ( 2, 2, 2 ), radius = 0.1, index = 1 )
    self.no_overlap = asa.sphere( centre = ( -1, -1, -1 ), radius = 0.1, index = 1 )


  def test_filtering(self):

    self.assertFalse( self.predicate( other = self.sphere ) )
    self.assertTrue( self.predicate( other =  self.overlap1 ) )
    self.assertTrue( self.predicate( other =  self.overlap2 ) )
    self.assertFalse( self.predicate( other = self.no_overlap ) )


  def make_indexer_linear(self):

    indexer = asa.indexing.linear_spheres()
    indexer.add( object = self.sphere, position = self.sphere.centre )
    indexer.add( object = self.overlap1, position = self.overlap1.centre )
    indexer.add( object = self.overlap2, position = self.overlap2.centre )
    indexer.add( object = self.no_overlap, position = self.no_overlap.centre )

    return indexer


  def test_filtering_linear(self):

    indexer = self.make_indexer_linear()

    range = asa.accessibility.filter(
      range = indexer.close_to( centre = self.sphere.centre ),
      predicate = self.predicate,
      )
    self.assertFalse( range.empty() )
    self.assertEqual( len( range ), 2 )
    self.assertTrue( self.overlap1.index in set( s.index for s in range ) )
    self.assertTrue( self.overlap2.index in set( s.index for s in range ) )

    range = asa.accessibility.filter(
      range = indexer.close_to( centre = self.no_overlap.centre ),
      predicate = asa.accessibility.overlap_equality_predicate( object = self.no_overlap ),
      )
    self.assertTrue( range.empty() )
    self.assertEqual( len( range ), 0 )


  def get_overlapping_spheres(self, sphere):

    return asa.accessibility.filter(
      range = self.make_indexer_linear().close_to( centre = sphere.centre ),
      predicate = asa.accessibility.overlap_equality_predicate( object = sphere ),
      )


  def test_filtering_lifetime(self):

    range = self.get_overlapping_spheres( sphere = self.sphere )
    self.assertEqual( len( range ), 2 )
    self.assertTrue( self.overlap1.index in set( s.index for s in range ) )
    self.assertTrue( self.overlap2.index in set( s.index for s in range ) )


class TestPythagoreanChecker(unittest.TestCase):

  def setUp(self):

    self.checker = asa.accessibility.pythagorean_checker()


  def test_add_from_list(self):

    self.assertTrue( self.checker( point = ( 1, 1, 1 ) ) )
    s = asa.sphere( centre = ( 1, 1, 1 ), radius = 0.1, index = 0 )
    self.checker.add( neighbours = [ s ] )
    self.assertFalse( self.checker( point = ( 1, 1, 1 ) ) )
    self.assertTrue( self.checker( point = ( 1.2, 1, 1 ) ) )


  def test_add_from_range_1(self):

    self.assertTrue( self.checker( point = ( 1, 1, 1 ) ) )

    s = asa.sphere( centre = ( 1, 1, 1 ), radius = 0.1, index = 0 )
    indexer = asa.indexing.linear_spheres()
    indexer.add( object = s, position = s.centre )

    neighbours = asa.accessibility.filter(
      range = indexer.close_to( centre = s.centre ),
      predicate = asa.accessibility.overlap_equality_predicate( object = s ),
      )

    self.assertEqual( len( neighbours ), 0 )
    self.checker.add( neighbours = neighbours )
    self.assertTrue( self.checker( point = ( 1, 1, 1 ) ) )

    neighbours = asa.accessibility.filter(
      range = indexer.close_to( centre = s.centre ),
      predicate = asa.accessibility.overlap_equality_predicate(
        object = asa.sphere( centre = ( 1, 1, 1 ), radius = 0.1, index = 1 )
        ),
      )
    self.assertEqual( len( neighbours ), 1 )
    self.checker.add( neighbours = neighbours )
    self.assertFalse( self.checker( point = ( 1, 1, 1 ) ) )


  def test_add_from_range_2(self):

    self.assertTrue( self.checker( point = ( 1, 1, 1 ) ) )

    s = asa.sphere( centre = ( 1, 1, 1 ), radius = 0.1, index = 0 )
    indexer = asa.indexing.linear_spheres()
    indexer.add( object = s, position = s.centre )

    neighbours = asa.accessibility.filter(
      range = indexer.close_to( centre = s.centre ),
      predicate = asa.accessibility.overlap_equality_predicate(
        object = asa.sphere( centre = ( 2, 2, 2 ), radius = 0.1, index = 1 )
        ),
      )
    self.assertEqual( len( neighbours ), 0 )
    self.checker.add( neighbours = neighbours )
    self.assertTrue( self.checker( point = ( 1, 1, 1 ) ) )

    neighbours = asa.accessibility.filter(
      range = indexer.close_to( centre = s.centre ),
      predicate = asa.accessibility.overlap_equality_predicate(
        object = asa.sphere( centre = ( 2, 2, 2 ), radius = 1.64, index = 2 )
        ),
      )
    self.assertEqual( len( neighbours ), 1 )
    self.checker.add( neighbours = neighbours )
    self.assertFalse( self.checker( point = ( 1, 1, 1 ) ) )


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

    self.assertTrue( isinstance( indexer.regular, type ) )
    self.assertEqual(
      len( indexer.regular ),
      len( [ a for a in ATOMS if not a.parent().altloc ] ),
      )


  def test_get_linear_indexer(self):

    indexer = asa.get_linear_indexer_for( descriptions = DESCRIPTIONS )
    self.assertTrue( isinstance( indexer, altloc.Indexer ) )
    self.check_indexer( indexer = indexer, type = asa.indexing.linear_spheres )


  def test_get_hash_indexer(self):

    indexer = asa.get_hash_indexer_for( descriptions = DESCRIPTIONS )
    self.assertTrue( isinstance( indexer, altloc.Indexer ) )
    self.check_indexer( indexer = indexer, type = asa.indexing.hash_spheres )


  def test_get_optimal_indexer(self):

    indexer = asa.get_optimal_indexer_for( descriptions = DESCRIPTIONS )
    self.assertTrue( isinstance( indexer, altloc.Indexer ) )
    self.check_indexer( indexer = indexer, type = asa.indexing.linear_spheres )


  def test_get_voxelizer_for(self):

    voxelizer = asa.get_voxelizer_for( descriptions = DESCRIPTIONS )
    from mmtbx.geometry import shared_types
    self.assertTrue( isinstance( voxelizer, shared_types.voxelizer ) )


  def test_asa_calculation_simple_linear(self):

    result = asa.calculate(
      atoms = ATOMS,
      indexer_selector = asa.get_linear_indexer_for,
      calculation = asa.simple_surface_calculation,
      )
    self.asa_result_check( result = result )


  def test_asa_calculation_aa_linear(self):

    result = asa.calculate(
      atoms = ATOMS,
      indexer_selector = asa.get_linear_indexer_for,
      calculation = asa.altloc_averaged_calculation,
      )
    self.asa_result_check( result = result )


  def test_asa_calculation_aa_hash(self):

    result = asa.calculate(
      atoms = ATOMS,
      indexer_selector = asa.get_hash_indexer_for,
      calculation = asa.altloc_averaged_calculation,
      )
    self.asa_result_check( result = result )


  def test_asa_calculation_simple_hash(self):

    result = asa.calculate(
      atoms = ATOMS,
      indexer_selector = asa.get_hash_indexer_for,
      calculation = asa.simple_surface_calculation,
      )
    self.asa_result_check( result = result )


  def test_asa_calculator0(self):

    myatoms = ATOMS[:2]
    myradii = ( 2, -2 )

    calc = asa.calculator(
      coordinate_adaptor = asa.coordinate_adaptor( array = myatoms.extract_xyz() ),
      radius_adaptor = asa.radius_adaptor( array = myradii ),
      probe = PROBE,
      )
    self.assertEqual( calc.accessible_surface_points( index = 0 ), 960 )
    self.assertRaises( RuntimeError, calc.accessible_surface_points, index = 1 )


  def test_asa_calculator1(self):

    myatoms = ATOMS[:12]
    myatoms.extend( ATOMS[13:] )
    mycoords = myatoms.extract_xyz()
    myradii = RADII[:12] + RADII[13:]

    calc = asa.calculator(
      coordinate_adaptor = asa.coordinate_adaptor( array = mycoords ),
      radius_adaptor = asa.radius_adaptor( array = myradii ),
      probe = PROBE,
      )
    myvalues = self.ACCESSIBLES[:12] + self.ACCESSIBLES[13:]

    for ( index, count, radius ) in zip( range( len( myatoms )), myvalues, myradii ):
      self.assertEqual( calc.accessible_surface_points( index = index ), count )
      self.assertAlmostEqual(
        calc.accessible_surface_area( index = index ),
        SAMPLING_POINTS.unit_area * count * ( radius + PROBE ) ** 2,
        7,
        )

    self.assertTrue(
      calc.is_overlapping_sphere( centre = mycoords[0], radius = myradii[0] )
      )
    self.assertFalse(
      calc.is_overlapping_sphere( centre = ( 0, 0, 0 ), radius = 10 )
      )
    self.assertTrue(
      calc.is_overlapping_sphere( centre = ( 0, 0, 0 ), radius = 15 )
      )


  def test_asa_calculator2(self):

    myatoms = ATOMS[:12]
    myatoms.extend( ATOMS[13:] )
    myradii = RADII[:12] + RADII[13:]

    calc = asa.calculator(
      coordinate_adaptor = asa.coordinate_adaptor(
        array = myatoms,
        transformation = lambda a: a.xyz,
        ),
      radius_adaptor = asa.radius_adaptor( array = myradii ),
      probe = PROBE,
      )
    myvalues = self.ACCESSIBLES[:12] + self.ACCESSIBLES[13:]

    for ( index, count, radius ) in zip( range( len( myatoms )), myvalues, myradii ):
      self.assertEqual( calc.accessible_surface_points( index = index ), count )
      self.assertAlmostEqual(
        calc.accessible_surface_area( index = index ),
        SAMPLING_POINTS.unit_area * count * ( radius + PROBE ) ** 2,
        7,
        )


  def test_is_overlapping_sphere(self):

    mycoord = ATOMS.extract_xyz()[0]
    myradius = 2

    calc = asa.calculator(
      coordinate_adaptor = asa.coordinate_adaptor( array = [ mycoord ] ),
      radius_adaptor = asa.radius_adaptor( array = [ myradius ] ),
      probe = PROBE,
      )
    self.assertTrue(
      calc.is_overlapping_sphere( centre = mycoord, radius = myradius )
      )

    (x, y, z ) = mycoord
    diff = 1.9 * myradius
    self.assertTrue(
      calc.is_overlapping_sphere( centre = ( x + diff, y, z ), radius = myradius )
      )
    self.assertTrue(
      calc.is_overlapping_sphere( centre = ( x, y + diff, z ), radius = myradius )
      )
    self.assertTrue(
      calc.is_overlapping_sphere( centre = ( x, y, z + diff ), radius = myradius )
      )

    diff = 2.1 * myradius
    self.assertFalse(
      calc.is_overlapping_sphere( centre = ( x + diff, y, z ), radius = myradius )
      )
    self.assertFalse(
      calc.is_overlapping_sphere( centre = ( x, y + diff, z ), radius = myradius )
      )
    self.assertFalse(
      calc.is_overlapping_sphere( centre = ( x, y, z + diff ), radius = myradius )
      )

    import math
    diff = 1.9 * myradius / math.sqrt( 3 )
    self.assertTrue(
      calc.is_overlapping_sphere( centre = ( x + diff, y + diff, z + diff ), radius = myradius )
      )

    diff = 2.1  * myradius / math.sqrt( 3 )
    self.assertFalse(
      calc.is_overlapping_sphere( centre = ( x + diff, y + diff, z + diff ), radius = myradius )
      )


  def test_asa_const_radius_calculator(self):

    myatoms = ATOMS[:2]

    calc = asa.const_radius_calculator(
      coordinate_adaptor = asa.coordinate_adaptor( array = myatoms.extract_xyz() ),
      radius = 2,
      probe = PROBE,
      )
    self.assertEqual( calc.accessible_surface_points( index = 0 ), 583 )
    self.assertEqual( calc.accessible_surface_points( index = 1 ), 586 )


  def asa_result_check(self, result):

    self.assertEqual( len( result.values ), len( self.ACCESSIBLES ) )
    self.assertEqual( len( result.values ), len( SPHERES ) )

    self.assertEqual( len( result.points ), len( self.ACCESSIBLES ) )
    self.assertEqual( result.points, self.ACCESSIBLES )

    self.assertEqual( len( result.areas ), len( self.ACCESSIBLES ) )

    for ( count, sphere, area ) in zip( self.ACCESSIBLES, SPHERES, result.areas ):
      self.assertAlmostEqual(
        area,
        SAMPLING_POINTS.unit_area * count * sphere.radius_sq,
        7,
        )


suite_sphere = unittest.TestLoader().loadTestsFromTestCase(
  TestSphere
  )
suite_transformation = unittest.TestLoader().loadTestsFromTestCase(
  TestTransformation
  )
suite_linear_spheres_indexer = unittest.TestLoader().loadTestsFromTestCase(
  TestLinearSpheresIndexer
  )
suite_hash_spheres_indexer = unittest.TestLoader().loadTestsFromTestCase(
  TestHashSpheresIndexer
  )
suite_predicate = unittest.TestLoader().loadTestsFromTestCase(
  TestPredicate
  )
suite_pythagorean_checker = unittest.TestLoader().loadTestsFromTestCase(
  TestPythagoreanChecker
  )
suite_accessible_surface_area = unittest.TestLoader().loadTestsFromTestCase(
  TestAccessibleSurfaceArea
  )


alltests = unittest.TestSuite(
  [
    suite_sphere,
    suite_transformation,
    suite_linear_spheres_indexer,
    suite_hash_spheres_indexer,
    suite_predicate,
    suite_pythagorean_checker,
    suite_accessible_surface_area,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )
