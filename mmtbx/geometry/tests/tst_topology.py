from __future__ import absolute_import, division, print_function

from mmtbx.geometry import topology

import unittest

class TestAtom(unittest.TestCase):

  def test_1(self):

    foo = object()
    bar = object()

    a = topology.Atom( foo = foo, bar = bar )
    self.assertEqual( a.foo, foo )
    self.assertEqual( a.bar, bar )


class TestMolecule(unittest.TestCase):

  def test_0(self):

    m = topology.Molecule()
    self.assertEqual( m.size(), 0 )
    self.assertEqual( m.atoms, [] )
    self.assertEqual( m.atom_for, {} )
    self.assertEqual( m.descriptor_for, {} )
    self.assertEqual( list( m.graph.vertices() ), [] )
    self.assertEqual( list( m.graph.edges() ), [] )


  def test_1(self):

    m = topology.Molecule()
    a = topology.Atom()
    m.add( atom = a,  xyz = ( 0, 0, 0 ) )
    self.assertEqual( m.size(), 1 )
    self.assertEqual( m.atoms, [ a ] )
    self.assertEqual( len( m.atom_for ), 1 )
    self.assertTrue( a in m.atom_for.values() )
    self.assertEqual( len( m.descriptor_for ), 1 )
    self.assertTrue( a in m.descriptor_for )
    self.assertEqual( len( list( m.graph.vertices() ) ), 1 )
    self.assertEqual( list( m.graph.edges() ), [] )


  def test_2(self):

    m = topology.Molecule()
    a1 = topology.Atom()
    a2 = topology.Atom()
    m.add( atom = a1, xyz = ( 0, 0, 0 ) )
    m.add( atom = a2, xyz = ( 1, 1, 1 ) )
    self.assertEqual( m.size(), 2 )
    self.assertEqual( set( m.atoms ), set( [ a1, a2 ] ) )
    self.assertEqual( len( m.atom_for ), 2 )
    self.assertTrue( a1 in m.atom_for.values() )
    self.assertTrue( a2 in m.atom_for.values() )
    self.assertEqual( len( m.descriptor_for ), 2 )
    self.assertTrue( a1 in m.descriptor_for )
    self.assertTrue( a2 in m.descriptor_for )
    self.assertEqual( len( list( m.graph.vertices() ) ), 2 )
    self.assertEqual( len( list( m.graph.edges() ) ), 1 )
    edge = next(m.graph.edges())
    self.assertAlmostEqual( m.graph.edge_weight( edge = edge ), 1.73205, 5 )


class TestCompound(unittest.TestCase):

  def test_0(self):

    m = topology.Compound.create()
    self.assertEqual( m.atoms, [] )
    self.assertEqual( m.atom_for, {} )
    self.assertEqual( m.descriptor_for, {} )
    self.assertEqual( list( m.graph.vertices() ), [] )
    self.assertEqual( list( m.graph.edges() ), [] )

  def test_1(self):

    m = topology.Compound.create()
    a = topology.Atom()
    m.add_atom( atom = a )
    self.assertEqual( m.atoms, [ a ] )
    self.assertEqual( len( m.atom_for ), 1 )
    self.assertTrue( a in m.atom_for.values() )
    self.assertEqual( len( m.descriptor_for ), 1 )
    self.assertTrue( a in m.descriptor_for )
    self.assertEqual( len( list( m.graph.vertices() ) ), 1 )
    self.assertEqual( list( m.graph.edges() ), [] )
    self.assertEqual( m.distances_from( atom = a ), { a: 0 } )
    self.assertEqual( m.connected_segments(), [ [ a ] ] )


  def test_2(self):

    m = topology.Compound.create()
    a1 = topology.Atom()
    a2 = topology.Atom()
    m.add_atom( atom = a1 )
    m.add_atom( atom = a2 )
    self.assertEqual( set( m.atoms ), set( [ a1, a2 ] ) )
    self.assertEqual( len( m.atom_for ), 2 )
    self.assertTrue( a1 in m.atom_for.values() )
    self.assertTrue( a2 in m.atom_for.values() )
    self.assertEqual( len( m.descriptor_for ), 2 )
    self.assertTrue( a1 in m.descriptor_for )
    self.assertTrue( a2 in m.descriptor_for )
    self.assertEqual( len( list( m.graph.vertices() ) ), 2 )
    self.assertEqual( len( list( m.graph.edges() ) ), 0 )
    self.assertEqual( m.distances_from( atom = a1 ), { a1: 0, a2: None } )
    self.assertEqual( m.distances_from( atom = a2 ), { a2: 0, a1: None } )
    self.assertEqual(
      set( frozenset( s ) for s in m.connected_segments() ),
      set( [ frozenset( [ a1 ] ), frozenset( [ a2 ] ) ] ),
      )

    m.add_bond( left = a1, right = a2 )
    self.assertEqual( len( list( m.graph.vertices() ) ), 2 )
    self.assertEqual( len( list( m.graph.edges() ) ), 1 )
    self.assertEqual( m.distances_from( atom = a1 ), { a1: 0, a2: 1 } )
    self.assertEqual( m.distances_from( atom = a2 ), { a2: 0, a1: 1 } )
    self.assertEqual(
      set( frozenset( s ) for s in m.connected_segments() ),
      set( [ frozenset( [ a1, a2 ] ) ] ),
      )

    ss1 = m.subset( atoms = [ a1 ] )
    self.assertEqual( len( ss1.atom_for ), 1 )
    self.assertTrue( a1 in ss1.atom_for.values() )
    self.assertEqual( len( ss1.descriptor_for ), 1 )
    self.assertTrue( a1 in ss1.descriptor_for )
    self.assertEqual( len( list( ss1.graph.vertices() ) ), 1 )
    self.assertEqual( len( list( ss1.graph.edges() ) ), 0 )

    ss2 = m.subset( atoms = [ a2 ] )
    self.assertEqual( len( ss2.atom_for ), 1 )
    self.assertTrue( a2 in ss2.atom_for.values() )
    self.assertEqual( len( ss2.descriptor_for ), 1 )
    self.assertTrue( a2 in ss2.descriptor_for )
    self.assertEqual( len( list( ss2.graph.vertices() ) ), 1 )
    self.assertEqual( len( list( ss2.graph.edges() ) ), 0 )


  def test_3(self):

    atoms = [
      topology.Atom( name = "N", element =  "N", xyz = ( 11.498, 10.510, 10.231 ) ),
      topology.Atom( name = "CA", element =  "C", xyz = ( 12.730, 11.073, 10.769 ) ),
      topology.Atom( name = "C", element =  "C", xyz = ( 13.674, 9.966, 11.221 ) ),
      topology.Atom( name = "O", element =  "O", xyz = ( 13.739, 8.902, 10.605 ) ),
      topology.Atom( name = "CB", element =  "C", xyz = ( 12.421, 12.004, 11.944 ) ),
      topology.Atom( name = "CG", element =  "C", xyz = ( 11.478, 13.179, 11.661 ) ),
      topology.Atom( name = "CD1", element =  "C", xyz = ( 11.043, 13.834, 12.963 ) ),
      topology.Atom( name = "CD2", element =  "C", xyz = ( 12.126, 14.201, 10.736 ) ),
      ]
    compound = topology.Compound.from_structure( atoms = atoms, tolerance = 0.1 )

    self.assertEqual(
      set( frozenset( [ l.name, r.name ] ) for ( l, r ) in compound.bonds ),
      set(
        [ frozenset( [ "N", "CA" ] ), frozenset( [ "CA", "C" ] ),
          frozenset( [ "C", "O" ] ), frozenset( [ "CA", "CB" ] ),
          frozenset( [ "CB", "CG" ] ), frozenset( [ "CG", "CD1" ] ),
          frozenset( [ "CG", "CD2" ] ),
          ]
        )
      )


class TestMcGregorMatch(unittest.TestCase):

  def test_asn_leu(self):

    l_ca = topology.Atom( label = "CA" )
    l_cb = topology.Atom( label = "C" )
    l_cg = topology.Atom( label = "C" )
    l_cd1 = topology.Atom( label = "C" )
    l_cd2 = topology.Atom( label = "C" )
    leu = topology.Molecule()
    leu.add( atom = l_ca, xyz = ( -1.0085, -0.590773,  0.814318 ) )
    leu.add( atom = l_cb,  xyz = (  0.0275, -0.557773, -0.314682 ) )
    leu.add( atom = l_cg,  xyz = (  1.2335,  0.374227, -0.138682 ) )
    leu.add( atom = l_cd1, xyz = (  2.3065,  0.046227, -1.16768  ) )
    leu.add( atom = l_cd2, xyz = (  0.8395,  1.84323,  -0.230682 ) )

    a_ca = topology.Atom( label = "CA" )
    a_cb = topology.Atom( label = "C" )
    a_cg = topology.Atom( label = "C" )
    a_od1 = topology.Atom( label = "C" )
    a_nd2 = topology.Atom( label = "C" )
    asn = topology.Molecule()
    asn.add( atom = a_ca, xyz = ( -1.03327, -0.544348,  0.860946 ) )
    asn.add( atom = a_cb,  xyz = (  0.10486, -0.548357, -0.164901 ) )
    asn.add( atom = a_cg,  xyz = (  0.990984, 0.682823, -0.070521 ) )
    asn.add( atom = a_od1, xyz = (  1.39496,  1.24684,  -1.08724  ) )
    asn.add( atom = a_nd2, xyz = (  1.29745,  1.10599,   1.15228  ) )

    res = topology.McGregorMatch(
      molecule1 = leu,
      molecule2 = asn,
      is_valid = lambda match: any( m.label == "CA" for m in match ),
      vertex_equality = lambda l, r: l.label == r.label,
      edge_equality = lambda l, r: abs( l - r ) < 0.1
      )
    self.assertEqual( res.length(), 3 )
    mapping = res.remapped()
    self.assertTrue( ( l_ca, a_ca ) in mapping )
    self.assertTrue( ( l_cb, a_cb ) in mapping )
    self.assertTrue( ( l_cg, a_cg ) in mapping )
    self.assertTrue( ( l_cd1, a_od1 ) not in mapping )


class TestRascalMatch(unittest.TestCase):

  def test_asn_leu(self):

    l_ca = topology.Atom( label = "CA" )
    l_cb = topology.Atom( label = "C" )
    l_cg = topology.Atom( label = "C" )
    l_cd1 = topology.Atom( label = "C" )
    l_cd2 = topology.Atom( label = "C" )
    leu = topology.Molecule()
    leu.add( atom = l_ca, xyz = ( -1.0085, -0.590773,  0.814318 ) )
    leu.add( atom = l_cb, xyz = (  0.0275, -0.557773, -0.314682 ) )
    leu.add( atom = l_cg, xyz = (  1.2335,  0.374227, -0.138682 ) )
    leu.add( atom = l_cd1, xyz = (  2.3065,  0.046227, -1.16768  ) )
    leu.add( atom = l_cd2, xyz = (  0.8395,  1.84323,  -0.230682 ) )

    a_ca = topology.Atom( label = "CA" )
    a_cb = topology.Atom( label = "C" )
    a_cg = topology.Atom( label = "C" )
    a_od1 = topology.Atom( label = "C" )
    a_nd2 = topology.Atom( label = "C" )
    asn = topology.Molecule()
    asn.add( atom = a_ca, xyz = ( -1.03327, -0.544348,  0.860946 ) )
    asn.add( atom = a_cb,  xyz = (  0.10486, -0.548357, -0.164901 ) )
    asn.add( atom = a_cg,  xyz = (  0.990984, 0.682823, -0.070521 ) )
    asn.add( atom = a_od1, xyz = (  1.39496,  1.24684,  -1.08724  ) )
    asn.add( atom = a_nd2, xyz = (  1.29745,  1.10599,   1.15228  ) )

    m = topology.RascalMatch(
      molecule1 = leu,
      molecule2 = asn,
      vertex_equality = lambda l, r: l.label == r.label,
      edge_equality = lambda l, r: abs( l - r ) <= 0.1,
      )
    self.assertEqual( m.count(), 1 )
    self.assertEqual( m.length(), 3 )
    mapping = m.remapped()[0]
    self.assertEqual( len( mapping ), 3 )
    self.assertTrue( ( l_ca, a_ca ) in mapping )
    self.assertTrue( ( l_cb, a_cb ) in mapping )
    self.assertTrue( ( l_cg, a_cg ) in mapping )
    self.assertTrue( ( l_cd1, a_od1 ) not in mapping )


class TestGreedyMatch(unittest.TestCase):

  def test_asn_leu(self):

    l_ca = topology.Atom( label = "CA" )
    l_cb = topology.Atom( label = "C" )
    l_cg = topology.Atom( label = "C" )
    l_cd1 = topology.Atom( label = "C" )
    l_cd2 = topology.Atom( label = "C" )
    leu = topology.Molecule()
    leu.add( atom = l_ca, xyz = ( -1.0085, -0.590773,  0.814318 ) )
    leu.add( atom = l_cb, xyz = (  0.0275, -0.557773, -0.314682 ) )
    leu.add( atom = l_cg, xyz = (  1.2335,  0.374227, -0.138682 ) )
    leu.add( atom = l_cd1, xyz = (  2.3065,  0.046227, -1.16768  ) )
    leu.add( atom = l_cd2, xyz = (  0.8395,  1.84323,  -0.230682 ) )

    a_ca = topology.Atom( label = "CA" )
    a_cb = topology.Atom( label = "C" )
    a_cg = topology.Atom( label = "C" )
    a_od1 = topology.Atom( label = "C" )
    a_nd2 = topology.Atom( label = "C" )
    asn = topology.Molecule()
    asn.add( atom = a_ca, xyz = ( -1.03327, -0.544348,  0.860946 ) )
    asn.add( atom = a_cb,  xyz = (  0.10486, -0.548357, -0.164901 ) )
    asn.add( atom = a_cg,  xyz = (  0.990984, 0.682823, -0.070521 ) )
    asn.add( atom = a_od1, xyz = (  1.39496,  1.24684,  -1.08724  ) )
    asn.add( atom = a_nd2, xyz = (  1.29745,  1.10599,   1.15228  ) )

    m = topology.GreedyMatch(
      molecule1 = leu,
      molecule2 = asn,
      vertex_equality = lambda l, r: l.label == r.label,
      edge_equality = lambda l, r: abs( l - r ) <= 0.1,
      )
    self.assertEqual( m.count(), 1 )
    self.assertEqual( m.length(), 3 )
    mapping = m.remapped()[0]
    self.assertEqual( len( mapping ), 3 )
    self.assertTrue( ( l_ca, a_ca ) in mapping )
    self.assertTrue( ( l_cb, a_cb ) in mapping )
    self.assertTrue( ( l_cg, a_cg ) in mapping )
    self.assertTrue( ( l_cd1, a_od1 ) not in mapping )


suite_atom = unittest.TestLoader().loadTestsFromTestCase(
  TestAtom
  )
suite_molecule = unittest.TestLoader().loadTestsFromTestCase(
  TestMolecule
  )
suite_compound = unittest.TestLoader().loadTestsFromTestCase(
  TestCompound
  )
suite_mcgregor_match = unittest.TestLoader().loadTestsFromTestCase(
  TestMcGregorMatch
  )
suite_rascal_match= unittest.TestLoader().loadTestsFromTestCase(
  TestRascalMatch
  )
suite_greedy_match= unittest.TestLoader().loadTestsFromTestCase(
  TestGreedyMatch
  )


alltests = unittest.TestSuite(
  [
    suite_atom,
    suite_molecule,
    suite_compound,
    suite_mcgregor_match,
    suite_rascal_match,
    suite_greedy_match,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

