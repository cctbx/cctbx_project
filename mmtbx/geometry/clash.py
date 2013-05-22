from __future__ import division

from cctbx import sgtbx # import dependency
import mmtbx.geometry.primitive # import dependency
import mmtbx.geometry.shared_types # import dependency

import boost.python
ext = boost.python.import_ext( "mmtbx_geometry_clash_ext" )
from mmtbx_geometry_clash_ext import *

def to_xray_structure(atoms, symmetry):
  """
  A utility function to convert atoms into cctbx.xray.structure
  """

  from cctbx import xray
  from cctbx.array_family import flex
  cell = symmetry.unit_cell()

  return xray.structure(
      crystal_symmetry = symmetry,
      scatterers = flex.xray_scatterer(
        [
          xray.scatterer(
            label = a.pdb_label_columns(),
            site = cell.fractionalize( a.xyz ),
            )
          for a in atoms
          ]
        )
      )


def altloc_strategy_from_atom(atom):

  altloc = atom.parent().altloc

  if altloc:
    return altloc_strategy.alternate( identifier = altloc )

  else:
    return altloc_strategy.regular()


def radius_from_atom(atom, table):

  return table[ atom.determine_chemical_element_simple().strip().capitalize() ]


def largest_atom_in(root, radii_table):

  radii = [
    radius_from_atom( atom = a, table = radii_table ) for a in root.atoms()
    ]

  if radii:
    return max( radii )

  else:
    return 0.0


def max_radius_in(roots, vdw_table):

  max_radii = [
      largest_atom_in( root = r, radii_table = vdw_table ) for r in roots
      ]
  return max( max_radii ) if max_radii else 0.0


def voxelizer_for(symmetry, max_radius, margin = 0.5):

  ( xs, ys, zs ) = zip(
    *symmetry.space_group_info().direct_space_asu().shape_vertices()
    )

  centre = (
    ( min( xs ) + max( xs ) ) / 2.0,
    ( min( ys ) + max( ys ) ) / 2.0,
    ( min( zs ) + max( zs ) ) / 2.0,
    )

  step = 2 * max_radius + margin

  return mmtbx.geometry.shared_types.voxelizer(
    base = symmetry.unit_cell().orthogonalize( centre ),
    step = ( step, step, step ),
    )


def linear_indexer_for(symmetry, max_radius, margin):

  return indexing.linear_spheres()


def hash_indexer_for(symmetry, max_radius, margin):

  return indexing.hash_spheres(
    voxelizer = voxelizer_for(
      symmetry = symmetry,
      max_radius = max_radius,
      margin = margin,
      )
    )


class SymmetryEquivalent(object):
  """
  A potential symmetry equivalent of an atom
  """

  def __init__(self, atom, symop):

    self.atom = atom
    self.symop = symop


  @property
  def atom_group(self):

    return self.atom.parent()


  @property
  def residue_group(self):

    return self.atom_group.parent()


  @property
  def chain(self):

    return self.residue_group.parent()


  @property
  def molecule(self):

    return self.chain.parent()

  @property
  def root(self):

    return self.molecule.parent()


  def __eq__(self, other):

    return ( self.atom == other.atom ) and ( self.symop == other.symop )


  def __ne__(self, other):

    return not ( self == other )


  def __hash__(self):

    return hash( ( self.atom, self.symop ) )


class ClashingPair(object):
  """
  Clashing atom pair
  """

  def __init__(self, left, right):

    self.left = left
    self.right = right


  def is_self_clash(self):

    return self.left.root == self.right.root


  def frozenset(self):

    return frozenset( ( self.left, self.right ) )


  def atom_pair(self):

    return frozenset( ( self.left.atom, self.right.atom ) )


  def __hash__(self):

    return hash( self.frozenset() )


  def __eq__(self, other):

    return self.frozenset() == other.frozenset()


  def __ne__(self, other):

    return not ( self == other )


class ClashPoint(object):
  """
  Clashes with a central atom
  """

  def __init__(self, centre, others):

    self.centre = centre
    self.others = others


  def pairs(self):

    return ( ClashingPair( left = self.centre, right = o ) for o in self.others )


  def count(self):

    return len( self.others )


class ClashCollection(object):
  """
  All clashes in a model
  """

  def __init__(self, points):

    self.points = points


  def pairs(self):

    for point in self.points:
      for pair in point.pairs():
        yield pair


class DistanceCalculation(object):
  """
  Calculates clash distance for atom pair
  """

  def __init__(self, cell):

    self.cell = cell


  def __call__(self, pair):

    ( l_x, l_y, l_z ) = pair.left.symop * self.cell.fractionalize( pair.left.atom.xyz )
    ( r_x, r_y, r_z ) = pair.right.symop * self.cell.fractionalize( pair.right.atom.xyz )

    diff = self.cell.orthogonalize( ( l_x - r_x, l_y - r_y, l_z - r_z ) )

    import math
    return math.sqrt( sum( [ d * d for d in diff ] ) )


def calculate(roots, symmetry, vdw_table, ifactory = hash_indexer_for, margin = 1.0, tolerance = 0.5):

  max_radius = max_radius_in( roots = roots, vdw_table = vdw_table )
  asu_buffer = 2 * max_radius + margin

  asu_spheres = []
  indexer = ifactory(
    symmetry = symmetry,
    max_radius = max_radius,
    margin = margin,
    )

  for ( mid, root ) in enumerate( roots ):
    structure = to_xray_structure( atoms = root.atoms(), symmetry = symmetry )
    asu_mapping = structure.asu_mappings( buffer_thickness = asu_buffer )
    mappings = asu_mapping.mappings()
    assert len( root.atoms() ) == len( mappings )
    root_asu_spheres = []

    for ( aid, ( atom, mapping ) ) in enumerate( zip( root.atoms(), mappings ) ):
      radius = radius_from_atom( atom = atom, table = vdw_table )
      altloc = altloc_strategy_from_atom( atom = atom )
      considered = [
        sphere(
          centre = m.mapped_site(),
          radius = radius,
          molecule = mid,
          atom = aid,
          altloc = altloc,
          symop = asu_mapping.get_rt_mx( m ),
          )
        for m in mapping
        ]
      assert 1 <= len( considered )
      root_asu_spheres.append( considered[0] )
      indexer.add( object = considered[0] )

      for sph in considered[1:]:
        indexer.add( object = sph )

    asu_spheres.append( root_asu_spheres )

  clashes = []
  root_atoms = [ r.atoms() for r in roots ]

  for root_asu_spheres in asu_spheres:
    clash_points = []

    for sph in root_asu_spheres:
      interacting = filter(
        range = indexer.close_to( object = sph ),
        predicate = overlap_interaction_predicate( object = sph, tolerance = tolerance ),
        )

      if interacting:
        clash_points.append(
          ClashPoint(
            centre = SymmetryEquivalent(
              atom = root_atoms[ sph.molecule ][ sph.atom ],
              symop = sph.symop,
              ),
            others = [
              SymmetryEquivalent(
                atom = root_atoms[ s.molecule ][ s.atom ],
                symop = s.symop,
                )
              for s in interacting
              ]
            )
          )

    clashes.append( ClashCollection( points = clash_points ) )

  return clashes

