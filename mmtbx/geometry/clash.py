from __future__ import absolute_import, division, print_function

from cctbx import sgtbx # import dependency
import mmtbx.geometry.shared_types # import dependency

import boost_adaptbx.boost.python as bp
from functools import reduce
from six.moves import zip
ext = bp.import_ext( "mmtbx_geometry_clash_ext" )
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


def element_of(atom):

  return atom.element.strip().capitalize()


def linear_indexer_for(params):

  return indexing.linear_spheres()


def hash_indexer_for(params):

  return indexing.hash_spheres( voxelizer = params.voxelizer(), margin = 1 )


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


class Model(object):
  """
  A model, providing fast access to atoms by memory_id
  """

  def __init__(self, root):

    self.identifier = root.memory_id()
    self.atom_for = dict( ( a.memory_id(), a ) for a in root.atoms() )


  def atoms(self):

    return list(self.atom_for.values())


  def elements(self):

    return set( element_of( atom = a ) for a in self.atoms() )


  def __hash__(self):

    return hash( self.identifier )


  def __eq__(self, other):

    return self.identifier == other.identifier


  def __ne__(self, other):

    return not( self == other )


  def __getitem__(self, key):

    return self.atom_for[ key ]


class Case(object):
  """
  Calculation parameters that are transferable
  """

  def __init__(self, symmetry, vdw_radius_for, max_radius, margin):

    self.symmetry = symmetry
    self.vdw_radius_for = vdw_radius_for
    self.max_radius = max_radius
    self.margin = margin


  def asu_mapping_for(self, atoms):

    structure = to_xray_structure( atoms = atoms, symmetry = self.symmetry )
    return structure.asu_mappings(
      buffer_thickness = 2 * self.max_radius + self.margin
      )


  def voxelizer(self):

    vertices = self.symmetry.space_group_info().direct_space_asu().shape_vertices()
    xs, ys, zs = zip( *vertices )

    centre = (
      ( min( xs ) + max( xs ) ) / 2.0,
      ( min( ys ) + max( ys ) ) / 2.0,
      ( min( zs ) + max( zs ) ) / 2.0,
      )

    step = 2 * self.max_radius

    return mmtbx.geometry.shared_types.voxelizer(
      base = self.symmetry.unit_cell().orthogonalize( centre ),
      step = ( step, step, step ),
      )

  def radius_for(self, atom):

    return self.vdw_radius_for[ element_of( atom = atom ) ]


  def sphere_model_group_for(self, model):

    atoms = model.atoms()
    asu_mapping = self.asu_mapping_for( atoms = atoms )
    mappings = asu_mapping.mappings()
    assert len( atoms ) == len( mappings )
    asu_spheres = []
    other_spheres = []

    for ( atom, mapping ) in zip( atoms, mappings ):
      radius = self.radius_for( atom = atom )
      altloc = altloc_strategy_from_atom( atom = atom )
      aid = atom.memory_id()
      considered = [
        sphere(
          centre = m.mapped_site(),
          radius = radius,
          molecule = model.identifier,
          atom = aid,
          altloc = altloc,
          symop = asu_mapping.get_rt_mx( m ),
          )
        for m in mapping
        ]
      assert 1 <= len( considered )
      asu_spheres.append( considered[0] )
      other_spheres.extend( considered[1:] )

    return SphereModelGroup(
      model = model,
      group = [ asu_spheres, other_spheres ],
      )


  @classmethod
  def create(cls, symmetry, models, vdw_radius_for, margin = 1.0):

    import operator
    elements = reduce(
      operator.or_,
      [ m.elements() for m in models ],
      set()
      )

    return cls(
      symmetry = symmetry,
      vdw_radius_for = vdw_radius_for,
      max_radius = (
        max( [ vdw_radius_for[ e ] for e in elements ] ) if elements else 0.0
        ),
      margin = margin,
      )


class AtomDescriptor(object):
  """
  A descriptor that allows looking up the original object is necessary
  """

  def __init__(self, model, sphere):

    self.model = model
    self.sphere = sphere


  @property
  def atom(self):

    return self.model[ self.sphere.atom ]


class SphereModel(object):
  """
  Sphere representation of a model
  """

  def __init__(self, model, spheres):

    self.model = model
    self.spheres = spheres


  def descriptors(self):

    return (
      AtomDescriptor( model = self.model, sphere = s ) for s in self.spheres
      )


class SphereModelGroup(object):
  """
  Partitioned sphere model
  """

  def __init__(self, model, group):

    self.model = model
    self._group = group


  def group(self):

    return [
      SphereModel( model = self.model, spheres = s ) for s in self._group
      ]


class ASUContent(object):
  """
  Contains the asymmetric unit
  """

  def __init__(self, case, ifactory = hash_indexer_for):

    self.case = case
    self.indexer = ifactory( params = self.case )
    self.model_with = {}
    self.asu_transforms = []


  def add(self, sequence):

    self.model_with[ sequence.model.identifier ] = sequence.model
    ( asu, other ) = sequence.group()
    self.asu_transforms.append( asu )

    self._insert( sphere_model = asu )
    self._insert( sphere_model = other )


  def clashing_with(self, sphere, tolerance):

    return list(filter(
      range = self.indexer.close_to( object = sphere ),
      predicate = overlap_interaction_predicate(
        object = sphere,
        tolerance = tolerance,
        ),
      ))


  def descriptor_for(self, sphere):

    return AtomDescriptor(
      model = self.model_with[ sphere.molecule ],
      sphere = sphere,
      )


  def _insert(self, sphere_model):

    for d in sphere_model.descriptors():
      self.indexer.add( object = d.sphere, position = d.sphere.centre )


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


class ClashCalculation(object):
  """
  Calculates clashpoints
  """

  def __init__(self, tolerance = 0.5):

    self.tolerance = tolerance


  def __call__(self, asu, descriptor, accumulator):

    clashing = asu.clashing_with(
      sphere = descriptor.sphere,
      tolerance = self.tolerance,
      )

    if clashing:
      accumulator(
        ClashPoint(
          centre = SymmetryEquivalent(
            atom = descriptor.atom,
            symop = descriptor.sphere.symop,
            ),
          others = [
            SymmetryEquivalent(
              atom = asu.descriptor_for( sphere = s ).atom,
              symop = s.symop,
              )
            for s in clashing
            ]
          )
        )


class ClashCountCalculation(object):
  """
  Calculates number of clashes
  """

  def __init__(self, tolerance = 0.5):

    self.tolerance = tolerance


  def __call__(self, asu, descriptor, accumulator):

    clashing = asu.clashing_with(
      sphere = descriptor.sphere,
      tolerance = self.tolerance,
      )
    accumulator( len( clashing ) )


class ElementAccumulator(object):
  """
  Collects all elements in a list
  """

  def __init__(self):

    self.result = []


  def __call__(self, element):

    self.result.append( element )


class ReducingAccumulator(object):
  """
  Sums up the elements as they come in
  """

  def __init__(self, initial, operation):

    self.result = initial
    self.operation = operation


  def __call__(self, element):

    self.result = self.operation( self.result, element )


def calculate_contents(roots, symmetry, vdw_radius_for, ifactory, margin):

  models = [ Model( root = r ) for r in roots ]
  asu = ASUContent(
    case = Case.create(
      symmetry = symmetry,
      models = models,
      vdw_radius_for = vdw_radius_for,
      margin = margin,
      ),
    ifactory = ifactory,
    )

  for m in models:
    asu.add( sequence = asu.case.sphere_model_group_for( model = m ) )

  return asu


def calculate(
  roots,
  symmetry,
  vdw_table,
  ifactory = hash_indexer_for,
  margin = 1.0,
  calculation = None,
  accumulator = None,
  ):

  asu = calculate_contents(
    roots = roots,
    symmetry = symmetry,
    vdw_radius_for = vdw_table,
    ifactory = ifactory,
    margin = margin,
    )

  if calculation is None:
    calculation = ClashCalculation()

  if accumulator is None:
    accumulator = ElementAccumulator() # avoid changing mutable default

  for sphere_model in asu.asu_transforms:
    for d in sphere_model.descriptors():
      calculation( asu = asu, descriptor = d, accumulator = accumulator )

  return accumulator.result
