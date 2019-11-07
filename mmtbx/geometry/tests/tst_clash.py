from __future__ import absolute_import, division, print_function

from mmtbx.geometry import clash

import unittest
from six.moves import zip

class TestAltlocStrategy(unittest.TestCase):

  def setUp(self):

    self.regular1 = clash.altloc_strategy.regular()
    self.regular2 = clash.altloc_strategy.regular()
    self.altloc1 = clash.altloc_strategy.alternate( identifier = "A" )
    self.altloc2 = clash.altloc_strategy.alternate( identifier = "A" )
    self.altloc3 = clash.altloc_strategy.alternate( identifier = "B" )


  def test_equality(self):

    self.assertFalse( self.regular1 == self.regular2 )
    self.assertFalse( self.regular1 == self.altloc1 )
    self.assertFalse( self.regular1 == self.altloc2 )
    self.assertFalse( self.regular1 == self.altloc3 )

    self.assertFalse( self.altloc1 == self.altloc2 )
    self.assertFalse( self.altloc1 == self.altloc3 )


  def run_interaction_tests(self, left, right, outcome):

    self.assertEqual( left.is_interacting_with( other = right ), outcome )
    self.assertEqual( right.is_interacting_with( other = left ), outcome )


  def test_interaction_regular(self):

    self.run_interaction_tests( self.regular1, self.regular1, True )
    self.run_interaction_tests( self.regular1, self.regular2, True )
    self.run_interaction_tests( self.regular1, self.altloc1, True )
    self.run_interaction_tests( self.regular1, self.altloc2, True )
    self.run_interaction_tests( self.regular1, self.altloc3, True )


  def test_interaction_altloc(self):

    self.run_interaction_tests( self.altloc1, self.altloc1, True )
    self.run_interaction_tests( self.altloc1, self.altloc2, True )
    self.run_interaction_tests( self.altloc1, self.altloc3, False )


class TestSphere(unittest.TestCase):

  def setUp(self):

    self.centre = ( 1, 2, 3 )
    self.radius = 4
    self.molecule = 5
    self.atom = 6
    self.symop = clash.sgtbx.rt_mx( ( 1, 0, 0, 0, 1, 0, 0, 0, 1 ), ( 0, 0, 1 ) )

    self.regular = clash.altloc_strategy.regular()
    self.altloc_a = clash.altloc_strategy.alternate( identifier = "A" )
    self.altloc_b = clash.altloc_strategy.alternate( identifier = "B" )


  def assertVectorsAlmostEqual(self, left, right, precision = 7):

    self.assertEqual( len( left ), len( right ) )

    for ( l, r ) in zip( left, right ):
      self.assertAlmostEqual( l, r, precision )


  def check_sphere(self, sphere):

    self.assertVectorsAlmostEqual( self.centre, sphere.centre )
    self.assertAlmostEqual( self.radius, sphere.radius )
    self.assertEqual( self.molecule, sphere.molecule )
    self.assertEqual( self.atom, sphere.atom )
    self.assertEqual( self.symop, sphere.symop )


  def check_altloc(self, sphere, interacting = [], noninteracting = []):

    altloc = sphere.altloc

    for a in interacting:
      self.assertTrue( altloc.is_interacting_with( other = a ) )

    for a in noninteracting:
      self.assertFalse( altloc.is_interacting_with( other = a ) )


  def test_regular(self):

    sphere = clash.sphere(
      centre = self.centre,
      radius = self.radius,
      molecule = self.molecule,
      atom = self.atom,
      altloc = self.regular,
      symop = self.symop,
      )
    self.check_sphere( sphere = sphere )
    self.check_altloc(
      sphere = sphere,
      interacting = [ self.regular, self.altloc_a, self.altloc_b],
      )


  def test_altloc_a(self):

    sphere = clash.sphere(
      centre = self.centre,
      radius = self.radius,
      molecule = self.molecule,
      atom = self.atom,
      altloc = self.altloc_a,
      symop = self.symop,
      )
    self.check_sphere( sphere = sphere )
    self.check_altloc(
      sphere = sphere,
      interacting = [ self.regular, self.altloc_a ],
      noninteracting = [ self.altloc_b ],
      )


  def test_altloc_b(self):

    sphere = clash.sphere(
      centre = self.centre,
      radius = self.radius,
      molecule = self.molecule,
      atom = self.atom,
      altloc = self.altloc_b,
      symop = self.symop,
      )
    self.check_sphere( sphere = sphere )
    self.check_altloc(
      sphere = sphere,
      interacting = [ self.regular, self.altloc_b],
      noninteracting = [ self.altloc_a ],
      )


  def test_altloc_c(self):

    sphere = clash.sphere(
      centre = self.centre,
      radius = self.radius,
      molecule = self.molecule,
      atom = self.atom,
      altloc = clash.altloc_strategy.alternate( identifier = "C" ),
      symop = self.symop,
      )
    self.check_sphere( sphere = sphere )
    self.check_altloc(
      sphere = sphere,
      interacting = [ self.regular ],
      noninteracting = [ self.altloc_a, self.altloc_b ],
      )


class TestOverlapInteractionPredicate(unittest.TestCase):

  def setUp(self):

    self.centre = ( 0, 0, 0 )
    self.radius = 1.7
    self.molecule = 0
    self.atom = 0
    self.symop = clash.sgtbx.rt_mx( ( 1, 0, 0, 0, 1, 0, 0, 0, 1 ), ( 0, 0, 0 ) )
    self.symop2 = clash.sgtbx.rt_mx( ( 1, 0, 0, 0, 1, 0, 0, 0, 1 ), ( 0, 0, 1 ) )
    self.altloc = clash.altloc_strategy.alternate( identifier = "A" )

    self.predicate = clash.overlap_interaction_predicate(
      object = clash.sphere(
        centre = self.centre,
        radius = self.radius,
        molecule = self.molecule,
        atom = self.atom,
        altloc = self.altloc,
        symop = self.symop,
        ),
      tolerance = 0,
      )


  def test_same_atom(self):

    # Exact same
    self.assertFalse(
      self.predicate(
        other = clash.sphere(
          centre = self.centre,
          radius = self.radius,
          molecule = self.molecule,
          atom = self.atom,
          altloc = self.altloc,
          symop = self.symop,
          )
        )
      )

    # Symmetry equivalent
    self.assertTrue(
      self.predicate(
        other = clash.sphere(
          centre = self.centre,
          radius = self.radius,
          molecule = self.molecule,
          atom = self.atom,
          altloc = self.altloc,
          symop = self.symop2,
          )
        )
      )
    self.assertFalse(
      self.predicate(
        other = clash.sphere(
          centre = ( 0, 0, 3 * self.radius ),
          radius = self.radius,
          molecule = self.molecule,
          atom = self.atom,
          altloc = self.altloc,
          symop = self.symop2,
          )
        )
      )


  def test_same_molecule(self):

    # Exact same
    self.assertFalse(
      self.predicate(
        other = clash.sphere(
          centre = self.centre,
          radius = self.radius,
          molecule = self.molecule,
          atom = self.atom + 1,
          altloc = self.altloc,
          symop = self.symop,
          )
        )
      )

    # Symmetry equivalent
    self.assertTrue(
      self.predicate(
        other = clash.sphere(
          centre = self.centre,
          radius = self.radius,
          molecule = self.molecule,
          atom = self.atom + 1,
          altloc = self.altloc,
          symop = self.symop2,
          )
        )
      )
    self.assertFalse(
      self.predicate(
        other = clash.sphere(
          centre = ( 0, 0, 3 * self.radius ),
          radius = self.radius,
          molecule = self.molecule,
          atom = self.atom + 1,
          altloc = self.altloc,
          symop = self.symop2,
          )
        )
      )
    self.assertFalse(
      self.predicate(
        other = clash.sphere(
          centre = self.centre,
          radius = self.radius,
          molecule = self.molecule,
          atom = self.atom + 1,
          altloc = clash.altloc_strategy.alternate( identifier = "B" ),
          symop = self.symop2,
          )
        )
      )


  def test_other_molecules(self):

    self.assertTrue(
      self.predicate(
        other = clash.sphere(
          centre = self.centre,
          radius = self.radius,
          molecule = self.molecule + 1,
          atom = self.atom,
          altloc = self.altloc,
          symop = self.symop,
          )
        )
      )
    self.assertTrue(
      self.predicate(
        other = clash.sphere(
          centre = self.centre,
          radius = self.radius,
          molecule = self.molecule + 1,
          atom = self.atom,
          altloc = self.altloc,
          symop = self.symop2,
          )
        )
      )
    self.assertFalse(
      self.predicate(
        other = clash.sphere(
          centre = self.centre,
          radius = self.radius,
          molecule = self.molecule + 1,
          atom = self.atom,
          altloc = clash.altloc_strategy.alternate( identifier = "B" ),
          symop = self.symop,
          )
        )
      )


class TestLinearSpheresIndexer(unittest.TestCase):

  def setUp(self):

    self.indexer = clash.indexing.linear_spheres()

    self.centre = ( 1, 2, 3 )
    self.radius = 4
    self.molecule = 5
    self.atom = 6
    self.altloc = clash.altloc_strategy.regular()
    self.symop = clash.sgtbx.rt_mx( ( 1, 0, 0, 0, 1, 0, 0, 0, 1 ), ( 0, 0, 1 ) )


  def test_structure(self):

    self.assertEqual( len( self.indexer ), 0 )

    sphere = clash.sphere(
      centre = self.centre,
      radius = self.radius,
      molecule = self.molecule,
      atom = self.atom,
      altloc = self.altloc,
      symop = self.symop,
      )

    self.indexer.add( object = sphere, position = sphere.centre )

    self.assertEqual( len( self.indexer ), 1 )

    closeby = self.indexer.close_to( centre = sphere.centre )
    self.assertEqual( len( closeby ), 1 )

    spheres = list( closeby )

    self.assertEqual( sphere.atom, spheres[0].atom )
    self.assertEqual( sphere.molecule, spheres[0].molecule )
    self.assertEqual( sphere.symop, spheres[0].symop )


class TestHashSpheresIndexer(unittest.TestCase):

  def setUp(self):

    import mmtbx.geometry.shared_types
    self.indexer = clash.indexing.hash_spheres(
      voxelizer = mmtbx.geometry.shared_types.voxelizer(
        base = ( 0, 0, 0 ),
        step = ( 1, 1, 1 ),
        ),
      margin = 1,
      )

    self.centre = ( 1, 2, 3 )
    self.radius = 4
    self.molecule = 5
    self.atom = 6
    self.altloc = clash.altloc_strategy.regular()
    self.symop = clash.sgtbx.rt_mx( ( 1, 0, 0, 0, 1, 0, 0, 0, 1 ), ( 0, 0, 1 ) )


  def test_structure(self):

    self.assertEqual( len( self.indexer ), 0 )

    sphere = clash.sphere(
      centre = self.centre,
      radius = self.radius,
      molecule = self.molecule,
      atom = self.atom,
      altloc = self.altloc,
      symop = self.symop,
      )

    self.indexer.add( object = sphere, position = sphere.centre )

    self.assertEqual( len( self.indexer ), 1 )

    closeby = self.indexer.close_to( centre = sphere.centre )
    self.assertEqual( len( closeby ), 1 )

    spheres = list( closeby )

    self.assertEqual( sphere.atom, spheres[0].atom )
    self.assertEqual( sphere.molecule, spheres[0].molecule )
    self.assertEqual( sphere.symop, spheres[0].symop )

    closeby = self.indexer.close_to( centre = [ c + 3 for c in self.centre ] )
    self.assertEqual( len( closeby ), 0 )


suite_altloc_strategy = unittest.TestLoader().loadTestsFromTestCase(
  TestAltlocStrategy
  )
suite_sphere = unittest.TestLoader().loadTestsFromTestCase(
  TestSphere
  )
suite_overlap_interaction_predicate = unittest.TestLoader().loadTestsFromTestCase(
  TestOverlapInteractionPredicate
  )
suite_linear_spheres_indexer = unittest.TestLoader().loadTestsFromTestCase(
  TestLinearSpheresIndexer
  )
suite_hash_spheres_indexer = unittest.TestLoader().loadTestsFromTestCase(
  TestHashSpheresIndexer
  )


alltests = unittest.TestSuite(
  [
    suite_altloc_strategy,
    suite_sphere,
    suite_overlap_interaction_predicate,
    suite_linear_spheres_indexer,
    suite_hash_spheres_indexer,
    ]
  )


def load_tests(loader, tests, pattern):

    return alltests


if __name__ == "__main__":
    unittest.TextTestRunner( verbosity = 2 ).run( alltests )

