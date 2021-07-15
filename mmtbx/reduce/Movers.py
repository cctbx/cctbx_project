##################################################################################
#                Copyright 2021  Richardson Lab at Duke University
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import math
import scitbx.matrix, scitbx.math
from scitbx.array_family import flex
from iotbx import pdb
import traceback

##################################################################################
# This is a set of classes that implement Reduce's "Movers".  These are sets of
# atoms that have more than one potential set of locations.
#  - They may have a-priori preferences for some locations over others.
#  - They may have both coarse locations that they test and a set of fine
#    locations around each coarse location.
#  - They may have final tune-up behaviors once their final locations are selected.
#
# There are two basic types of Movers:
#  - Rotator: One or more Hydrogen atoms that have a set of orientations spinning
#    around a common axis.
#  - Flipper: A structure that has pairs of atom groups that have similar densities, such
#    that flipping them across a center axis produces similar density fits.
#    They sometimes have optional hydrogen placements in the atom groups.
#
# All Movers have the following methods, parameters, and return types:
#  - type PositionReturn: ( flex<atom> atoms, flex<flex<vec3>> positions, flex<float> preferenceEnergies)
#  - PositionReturn CoarsePositions()
#     The first return lists all of the hierarchy atoms that move.
#     The second return has a vector that has some number of positions, each of
#       which is a vector the same length as the first return
#       with a vec3 position for each atom.
#     The third return has a vector of the same length as the number of positions in
#       each element of the third return; it indicates the relative favorability of
#       each possible location and should be added to the total energy by an
#       optimizer to break ties in otherwise-equivalent orientations.
#  - PositionReturn FinePositions(coarseIndex)
#     The coarseIndex indicates the index (0 for the first) of the relevant coarse
#       orientation.
#     The return values are the same as for CoarsePositions and they list potential
#       fine motions around the particular coarse position (not including the position
#       itself).  This function can be used by optimizers that wish to do heavy-weight
#       operations at a coarse resolution and then lightweight operations at a finer
#       scale; other optimizers may choose to ask for all of the fine positions and
#       append them to the coarse positions and globally optimize.
#     Note: Some Movers will return empty arrays.
#  - type FixUpReturn: ( flex<atom> atoms, flex<vec3> newPositions )
#  - FixUpReturn FixUp(coarseIndex)
#     The coarseIndex indicates the index (0 for the first) of the relevant coarse
#       orientation that was finally chosen.
#     The first return lists the atoms that should be repositioned.
#     The second return lists the new locations for each atom.
#     Note: This function may ask to modify atoms other than the ones it reported in its
#       coarse and fine position functions.  This is to attempt to optimizing bonds
#       for Flippers.  None of these depend on fine index.
#     Note: Some Movers will return empty arrays.  This indicates that no fix-up is
#       needed and the atoms are in their final positions.
#
# The caller is responsible for moving the specified atoms to their final positions,
# whether based on coarse or fine positions.  After applying the coarse or fine
# adjustment, they must call FixUp() and apply any final changes (which may require
# removing Hydrogens).
#
# The InteractionGraph.py script provides functions for determining which pairs of
# Movers have overlaps between movable atoms.
#

##################################################################################
class PositionReturn:
  # Return type from CoarsePosition() and FinePosition() calls.
  def __init__(self, atoms, positions, preferenceEnergies):
    self.atoms = atoms
    self.positions = positions
    self.preferenceEnergies = preferenceEnergies

##################################################################################
class FixUpReturn:
  # Return type from FixUp() calls.
  def __init__(self, atoms, newPositions):
    self.atoms = atoms
    self.newPositions = newPositions

##################################################################################
class MoverNull:
  '''A trivial Mover that returns a single result atom at a single location.
     Useful as a simple and fast test case for programs that use Movers.
     It also serves as a basic example of how to develop a new Mover.
  '''
  def __init__(self, atom):
    self._atom = atom
  def CoarsePositions(self):
    # returns: The original atom at its original position.
    return PositionReturn([ self._atom ],
        [ [ [ self._atom.xyz[0], self._atom.xyz[1], self._atom.xyz[2] ] ] ],
        [ 0.0 ])
  def FinePositions(self, coarseIndex):
    # returns: No fine positions for any coarse position.
    return PositionReturn([], [], [])
  def FixUp(self, coarseIndex):
    # No fixups for any coarse index.
    return FixUpReturn([], [])

##################################################################################
# @todo Consider having another constructor parameter that gives the initial rotation
# offset to place the hydrogens in a good starting location, then derived classes can
# indicate this rotation without actually having to change the positions of the Hydrogens
# that are rotating.  This offset will be added to all resulting rotations, and subtracted
# before applying the preference function.
class _MoverRotator:
  def __init__(self, atoms, axis, coarseRange, coarseStepDegrees = 30.0,
                doFineRotations = True, fineStepDegrees = 1.0,
                preferenceFunction = None, preferredOrientationScale = 1.0):
    """Constructs a Rotator, to be called by a derived class or test program but
       not usually user code.
       A base class for types of Movers that rotate one or more atoms around a single
       axis.  The derived class must determine the initial position for each atom, and
       a preferred-energies function that maps from an angle that is 0 at the initial
       location to an energy.  The range of coarse rotations is also specified, along
       with a flag telling whether fine rotations are allowed.  The coarse step size
       can also be specified, overriding the value from coarseStepDegrees.
       :param atoms: flex array of atoms that will be rotated.  These should be
       placed in one of the preferred orientations if there are any.
       :param axis: flex array with two scitbx::vec3<double> points, the first
       of which is the origin and the second is a vector pointing in the direction
       of the axis of rotation.  Positive rotations will be right-handed rotations
       around this axis.
       :param coarseRange: Range in degrees that the atoms can rotate around the axis
       from the starting position.  The range is considered closed on the negative
       end and open on the positive: [-coarseRange..coarseRange).  For example a
       range of 180 will rotate all the way around, from -180 (inclusive) to just
       slightly less than +180.
       :param coarseStepDegrees: The coarse step to take.
       To test only two fixed orientations, doFineRotations should be set
       to False, coarseRange to 180, and coarseStepDegrees to 180.
       :param doFineRotations: Specifies whether fine rotations are allowed for
       this instance.  To test only two fixed orientations, this should be set
       to False, coarseRange to 180, and coarseStepDegrees to 180.
       :param fineStepDegrees: The fine step to take.
       :param preferenceFunction: A function that takes a floating-point
       argument in degrees of rotation about the specified axis and produces
       a floating-point energy value that is multiplied by the value of
       preferredOrientationScale and then added to the energy
       returned for this orientation by CoarsePositions() and FinePositions().
       For the default value of None, no addition is performed.
       :param preferredOrientationScale: How much to scale the preferred-orientation
       energy by before adding it to the total.
    """
    self._atoms = atoms
    self._axis = axis
    self._coarseRange = coarseRange
    self._coarseStepDegrees = coarseStepDegrees
    self._doFineRotations = doFineRotations
    self._preferenceFunction = preferenceFunction
    self._fineStepDegrees = fineStepDegrees
    self._preferredOrientationScale = preferredOrientationScale

    # Make a list of coarse angles to try based on the coarse range (inclusive
    # for negative and exclusive for positive) and the step size.  We always
    # try 0 degrees.
    # Store the angles for use by CoarsePositions() and FinePositions()
    self._coarseAngles = [0.0]
    curStep = self._coarseStepDegrees
    while curStep <= self._coarseRange:
      # The interval is closed on the left, so we place
      self._coarseAngles.append(-curStep)
      # We only place on the right if we're strictly less than because
      # the interval is open on the right.
      if curStep < self._coarseRange:
        self._coarseAngles.append(curStep)
      curStep += self._coarseStepDegrees

    # Make a list of fine angles to try to cover half the coarse step size (inclusive
    # for negative and exclusive for positive) using the fine step size.  We never
    # try +0 degrees because that is the location of the coarse rotation.
    self._fineAngles = []
    curStep = self._fineStepDegrees
    fineRange = self._coarseStepDegrees / 2
    while curStep <= fineRange:
      # The interval is closed on the left, so we place
      self._fineAngles.append(-curStep)
      # We only place on the right if we're strictly less than because
      # the interval is open on the right.
      if curStep < fineRange:
        self._fineAngles.append(curStep)
      curStep += self._fineStepDegrees

  def _preferencesFor(self, angles, scale):
    """Return the preference energies for the specified angles, scaled by the scale.
       :param angles: List of angles to compute the preferences at.
       :param scale: Scale factor to apply to the resulting energy function that is
       evaluated at the angles.
       :returns Energy function evaluated at the angles multiplied by the scale factor.
       If there is no energy function all entries are 0.0.
    """
    preferences = [0.0] * len(angles)
    if self._preferenceFunction is not None:
      for i, a in enumerate(angles):
        preferences[i] = self._preferenceFunction(a) * scale
    return preferences

  def _posesFor(self, angles):
    """Return the atom poses for the specified angles.
       :param angles: List of angles to compute the poses at.
       :returns List of flex.vec3_double entries, where each entry records the
       positions of all atoms when rotated about the axis by the associated angle.
       There is one entry in the list for each angle.
    """
    # @todo Turn this into a flex array of flex arrays rather than a list of flex arrays.
    poses = []
    for agl in angles:
      atoms = flex.vec3_double()
      for atm in self._atoms:
        atoms.append(_rotateAroundAxis(atm, self._axis, agl))
      poses.append(atoms)
    return poses;

  def CoarsePositions(self):

    # Return the atoms, coarse-angle poses, and coarse-angle preferences
    return PositionReturn(self._atoms,
      self._posesFor(self._coarseAngles),
      self._preferencesFor(self._coarseAngles, self._preferredOrientationScale))

  def FinePositions(self, coarseIndex):
    if not self._doFineRotations:
      # No fine positions for any coarse position.
      return PositionReturn([], [], [])

    # We add the range of values to the coarse angle we're starting with to provide
    # the list of fine angles to try.
    angles = []
    try:
      ca = self._coarseAngles[coarseIndex]
    except Exception as e:
      raise ValueError("MoverRotator.FinePositions(): Bad coarseIndex: "+str(e))
    for fa in self._fineAngles:
      angle = fa + ca
      angles.append(angle)

    # Return the atoms and poses along with the preferences.
    return PositionReturn(self._atoms,
      self._posesFor(angles),
      self._preferencesFor(angles, self._preferredOrientationScale))

  def FixUp(self, coarseIndex):
    # No fixups for any coarse index.
    return FixUpReturn([], [])

##################################################################################
class MoverSingleHydrogenRotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, potentialAcceptors = [],
                  coarseStepDegrees = 10.0,
                  fineStepDegrees = 1.0, preferredOrientationScale = 1.0):
    """ A Mover that rotates a single Hydrogen around an axis from its bonded partner
       to the single bonded partner of its partner.  This is designed for use with OH,
       SH, and SeH configurations.  For partner-partner atoms that are bonded to a
       tetrahedron, the starting orientation is aligned between two of the edges of
       the tetrahedron.  For partner-partner atoms that are bonded to two atoms, the
       starting orientation is in the plane of these atoms and the partner-partner atom.
       :param atom: Hydrogen atom that will be rotated.  It must be bonded to a single
       atom, and that atom must also be bonded to a single other atom.  NOTE: As a side
       effect, this Hydrogen is immediately rotated to lie either in the plane of the two
       friends of the partner atom or to be between two of friends if there are three.
       :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
       structure that the atom from the first parameter interacts with that lists all of the
       bonded atoms.  Can be obtained by calling getBondedNeighborLists().
       :param potentialAcceptors: A flex array of atoms that are nearby potential acceptors.
       Coarse orientations are added that aim the hydrogen in the direction of these potential
       partners.
       :param coarseStepDegrees: The coarse step to take.
       :param fineStepDegrees: The fine step to take.
       :param preferredOrientationScale: How much to scale the preferred-orientation
       energy by before adding it to the total.
    """

    # Check the conditions to make sure we've been called with a valid atom.  This is a hydrogen with
    # a single bonded neighbor that has a single other bonded partner that has 2-3 other bonded friends.
    # Find the friends bonded to the partner besides the neighbor, which will be used to
    # determine the initial orientation for the hydrogen.
    if atom.element != "H":
      raise ValueError("MoverSingleHydrogenRotator(): Atom is not a hydrogen")
    neighbors = bondedNeighborLists[atom]
    if len(neighbors) != 1:
      raise ValueError("MoverSingleHydrogenRotator(): Atom does not have a single bonded neighbor")
    neighbor = neighbors[0]
    partners = bondedNeighborLists[neighbor]
    if len(partners) != 2:  # Me and one other
      raise ValueError("MoverSingleHydrogenRotator(): Atom does not have a single bonded neighbor's neighbor")
    partner = partners[0]
    if partner.i_seq == atom.i_seq:
      partner = partners[1]
    bonded = bondedNeighborLists[partner]
    friends = []
    for b in bonded:
      if b.i_seq != neighbor.i_seq:
        friends.append(b)

    # Determine the preference function (180 or 120) based on friends bonding structure
    # @todo Consider parameterizing the magic constant of 0.1 for the preference magnitude
    if len(friends) == 2:
      # There are two neighbors, so our function repeats every 180 degrees
      def preferenceFunction(degrees): return 0.1 + 0.1 * math.cos(degrees * (math.pi/180) * (360/180))
    elif len(friends) == 3:
      # There are three neighbors, so our function repeats every 120 degrees
      def preferenceFunction(degrees): return 0.1 + 0.1 * math.cos(degrees * (math.pi/180) * (360/120))
    else:
      raise ValueError("MoverSingleHydrogenRotator(): Atom's bonded neighbor's neighbor does not have 2-3 other bonds "+
      "it has "+str(len(friends)))

    # Determine the axis to rotate around, which starts at the partner atom and points at the neighbor.
    normal = (_rvec3(neighbor.xyz) - _rvec3(partner.xyz)).normalize()
    axis = flex.vec3_double([partner.xyz, normal])

    # Move the atom so that it is in one of the two preferred locations.  The preferred location depends on
    # whether there are two or three friends, it wants to be in the plane if there are two of them and it
    # wants to be between two edges if there are three.  It turns out that rotating it to point away from
    # one of the friends works in both of those cases.
    atom.xyz = _rotateOppositeFriend(atom, axis, partner, friends)

    # Make a list that contains just the single atom.
    atoms = [ atom ]

    # Construct our parent class, which will do all of the actual work based on our inputs.
    _MoverRotator.__init__(self, atoms, axis, 180, coarseStepDegrees = coarseStepDegrees,
      fineStepDegrees = fineStepDegrees, preferenceFunction = preferenceFunction, 
      preferredOrientationScale = preferredOrientationScale)

    # Now add orientations that point in the direction of the potential acceptors.
    # @todo The original C++ code aimed only at these (or near them for clashes) and in a
    # direction far from clashes.  We may need to do this for speed reasons and to reduce the
    # number of elements in each clique, but for now we try all coarse orientations.

    # Compute the dihedral angle from the Hydrogen to the potential acceptor through
    # the partner and neighbor.  This is the amount to rotate the hydrogen by.
    for a in potentialAcceptors:
      sites = [ atom.xyz, partner.xyz, neighbor.xyz, a.xyz ]
      degrees = scitbx.math.dihedral_angle(sites=sites, deg=True)
      self._coarseAngles.append(degrees)

##################################################################################
class MoverNH3Rotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, coarseStepDegrees = 30.0, fineStepDegrees = 1.0,
                preferredOrientationScale = 1.0):
    """ A Mover that rotates three Hydrogens around an axis from their bonded Nitrogen neighbor
       to the single bonded partner of its partner.  This is designed for use with NH3+,
       whose partner-partner atoms are bonded to a tetrahedron (two Carbons and a Hydrogen)
       of friends.
       The starting orientation has one of the Hydrogens aligned between two of the edges of
       the tetrahedron.
       :param atom: Nitrogen atom bonded to the three Hydrogens that will be rotated.
       It must be bonded to three Hydrogens and a single other
       atom, and the other atom must be bonded to three other atoms.  NOTE: As a side
       effect, the Hydrogens are immediately rotated to lie between two of the friends.
       :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
       structure that the atom from the first parameter interacts with that lists all of the
       bonded atoms.  Can be obtained by calling getBondedNeighborLists().
       the value.
       :param fineStepDegrees: The fine step to take.
       :param preferredOrientationScale: How much to scale the preferred-orientation
       energy by before adding it to the total.
    """

    # The Nitrogen is the neighbor in these calculations, making this code symmetric with the other
    # class code.
    neighbor = atom

    # Check the conditions to make sure we've been called with a valid atom.  This is a Nitrogen with
    # three hydrogens bonded and a single bonded neighbor that has 3 other bonded friends.
    # Find the friends bonded to the partner besides the neighbor, which will be used to
    # determine the initial orientation for the hydrogens.
    if neighbor.element != "N":
      raise ValueError("MoverNH3Rotator(): atom is not a Nitrogen")
    partners = bondedNeighborLists[neighbor]
    if len(partners) != 4:
      raise ValueError("MoverNH3Rotator(): atom does not have four bonded neighbors")
    hydrogens = []
    for a in partners:
      if a.element == "H":
        hydrogens.append(a)
      else:
        partner = a
    if len(hydrogens) != 3:
      raise ValueError("MoverNH3Rotator(): atom does not have three bonded hydrogens")
    bonded = bondedNeighborLists[partner]
    friends = []
    for b in bonded:
      if b.i_seq != neighbor.i_seq:
        friends.append(b)
    if len(friends) != 3:
      raise ValueError("MoverNH3Rotator(): Partner does not have three bonded friends")

    # Set the preference function to like 120-degree rotations away from the starting location.
    # @todo Consider parameterizing the magic constant of 0.1 for the preference magnitude
    # (for example, we might use preferredOrientationScale for this and not set a separate one?)
    def preferenceFunction(degrees): return 0.1 + 0.1 * math.cos(degrees * (math.pi/180) * (360/120))

    # Determine the axis to rotate around, which starts at the partner and points at the neighbor.
    normal = (_rvec3(neighbor.xyz) - _rvec3(partner.xyz)).normalize()
    axis = flex.vec3_double([partner.xyz, normal])

    # Move the Hydrogens so that they are in one of the preferred locations by rotating one of them to
    # point away from one of the friends.  The other two are located at +120 and -120 degrees rotated
    # around the axis from the first.
    hydrogens[0].xyz = _rotateOppositeFriend(hydrogens[0], axis, partner, friends)
    hydrogens[1].xyz = _rotateAroundAxis(hydrogens[0], axis, 120)
    hydrogens[2].xyz = _rotateAroundAxis(hydrogens[0], axis, -120)

    # Construct our parent class, which will do all of the actual work based on our inputs.
    _MoverRotator.__init__(self, hydrogens, axis, 180, coarseStepDegrees,
      fineStepDegrees = fineStepDegrees, preferenceFunction = preferenceFunction,
      preferredOrientationScale = preferredOrientationScale)

##################################################################################
class MoverAromaticMethylRotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists):
    """ A Mover that rotates three Hydrogens around an axis from their bonded Carbon neighbor
       to the single bonded partner of its partner.  This is designed for use with Aromatic
       CH3 (Methly) groups, whose partner-partner atoms are bonded to an aromatic ring, having
       two friends.
       The starting orientation has one of the Hydrogens pointing away from the plane of the ring.
       The only other preferred orientation is having that hydrogen point out the other side
       of the ring, so we only do 180-degree coarse and no fine rotation.
       :param atom: Carbon atom bonded to the three Hydrogens that will be rotated.
       It must be bonded to three Hydrogens and a single other
       atom, and the other atom must be bonded to two other atoms.  NOTE: As a side
       effect, the Hydrogens are immediately rotated to lie perpendicular to the friends.
       :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
       structure that the atom from the first parameter interacts with that lists all of the
       bonded atoms.  Can be obtained by calling getBondedNeighborLists().
    """

    # The Carbon is the neighbor in these calculations, making this code symmetric with the other
    # class code.
    neighbor = atom

    # Check the conditions to make sure we've been called with a valid atom.  This is a Carbon with
    # three hydrogens bonded and a single bonded neighbor that has 2 other bonded friends.
    # Find the friends bonded to the partner besides the neighbor, which will be used to
    # determine the initial orientation for the hydrogens.
    if neighbor.element != "C":
      raise ValueError("MoverAromaticMethylRotator(): atom is not a Carbon")
    partners = bondedNeighborLists[neighbor]
    if len(partners) != 4:
      raise ValueError("MoverAromaticMethylRotator(): atom does not have four bonded neighbors")
    hydrogens = []
    for a in partners:
      if a.element == "H":
        hydrogens.append(a)
      else:
        partner = a
    if len(hydrogens) != 3:
      raise ValueError("MoverAromaticMethylRotator(): atom does not have three bonded hydrogens")
    bonded = bondedNeighborLists[partner]
    friends = []
    for b in bonded:
      if b.i_seq != neighbor.i_seq:
        friends.append(b)
    if len(friends) != 2:
      raise ValueError("MoverAromaticMethylRotator(): Partner does not have two bonded friends")

    # Determine the axis to rotate around, which starts at the partner and points at the neighbor.
    normal = (_rvec3(neighbor.xyz) - _rvec3(partner.xyz)).normalize()
    axis = flex.vec3_double([partner.xyz, normal])

    # Move the Hydrogens so that they are in one of the preferred locations by rotating one of them to
    # point away from one of the friends and then rotating it 90 degrees.  The other two are located
    # at +120 and -120 degrees rotated around the axis from the first.
    hydrogens[0].xyz = _rotateOppositeFriend(hydrogens[0], axis, partner, friends)
    hydrogens[0].xyz = _rotateAroundAxis(hydrogens[0], axis, 90)
    hydrogens[1].xyz = _rotateAroundAxis(hydrogens[0], axis, 120)
    hydrogens[2].xyz = _rotateAroundAxis(hydrogens[0], axis, -120)

    # Construct our parent class, which will do all of the actual work based on our inputs.
    # We have a coarse step size of 180 degrees and a range of 180 degrees and do not
    # allow fine rotations.
    _MoverRotator.__init__(self, hydrogens, axis, 180, 180, doFineRotations = False)

##################################################################################
class MoverTetrahedralMethylRotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, coarseStepDegrees = 30.0,
                  fineStepDegrees = 1.0, preferredOrientationScale = 1.0):
    """ A Mover that rotates three Hydrogens around an axis from their bonded Carbon neighbor
       to the single bonded partner of its partner.  This is designed for use with tetrahedral
       partners whose partner-partner atoms are bonded to three friends.
       The starting orientation has the Hydrogens pointing between the side of the tetrahedron.
       It can rotate to any angle to optimize for hydrogen bonds.
       Note: Reduce does not normally rotate these groups because it is not worth the computational
       cost to do so (and because this can mask mis-placed atoms by forming spurious hydrogen
       bonds), but it will construct them so that they will be aligned staggered to the
       tetrahedron.
       @todo Construct these in Reduce so that they will be staggered but do not optimize them.
       :param atom: Carbon atom bonded to the three Hydrogens that will be rotated.
       It must be bonded to three Hydrogens and a single other
       atom, and the other atom must be bonded to three other atoms.  NOTE: As a side
       effect, the Hydrogens are immediately rotated to lie staggered.
       :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
       structure that the atom from the first parameter interacts with that lists all of the
       bonded atoms.  Can be obtained by calling getBondedNeighborLists().
       :param coarseStepDegrees: The coarse step to take.
       :param fineStepDegrees: The fine step to take.
       :param preferredOrientationScale: How much to scale the preferred-orientation
       energy by before adding it to the total.
    """

    # The Carbon is the neighbor in these calculations, making this code symmetric with the other
    # class code.
    neighbor = atom

    # Check the conditions to make sure we've been called with a valid atom.  This is a Carbon with
    # three hydrogens bonded and a single bonded neighbor that has 2 other bonded friends.
    # Find the friends bonded to the partner besides the neighbor, which will be used to
    # determine the initial orientation for the hydrogens.
    if neighbor.element != "C":
      raise ValueError("MoverTetrahedralMethylRotator(): atom is not a Carbon")
    partners = bondedNeighborLists[neighbor]
    if len(partners) != 4:
      raise ValueError("MoverTetrahedralMethylRotator(): atom does not have four bonded neighbors")
    hydrogens = []
    for a in partners:
      if a.element == "H":
        hydrogens.append(a)
      else:
        partner = a
    if len(hydrogens) != 3:
      raise ValueError("MoverTetrahedralMethylRotator(): atom does not have three bonded hydrogens")
    bonded = bondedNeighborLists[partner]
    friends = []
    for b in bonded:
      if b.i_seq != neighbor.i_seq:
        friends.append(b)
    if len(friends) != 3:
      raise ValueError("MoverTetrahedralMethylRotator(): Partner does not have two bonded friends")

    # Determine the axis to rotate around, which starts at the partner and points at the neighbor.
    normal = (_rvec3(neighbor.xyz) - _rvec3(partner.xyz)).normalize()
    axis = flex.vec3_double([partner.xyz, normal])

    # Move the Hydrogens so that they are in one of the preferred locations by rotating one of them to
    # point away from one of the friends.  The other two are located at +120 and -120 degrees rotated
    # around the axis from the first.
    hydrogens[0].xyz = _rotateOppositeFriend(hydrogens[0], axis, partner, friends)
    hydrogens[1].xyz = _rotateAroundAxis(hydrogens[0], axis, 120)
    hydrogens[2].xyz = _rotateAroundAxis(hydrogens[0], axis, -120)

    # Set the preference function to like 120-degree rotations away from the starting location.
    # @todo Consider parameterizing the magic constant of 0.1 for the preference magnitude
    def preferenceFunction(degrees): return 0.1 + 0.1 * math.cos(degrees * (math.pi/180) * (360/120))

    # Construct our parent class, which will do all of the actual work based on our inputs.
    # We have a coarse step size of 180 degrees and a range of 180 degrees and do not
    # allow fine rotations.
    _MoverRotator.__init__(self, hydrogens, axis, 180, fineStepDegrees = fineStepDegrees,
      preferenceFunction = preferenceFunction,
      preferredOrientationScale = preferredOrientationScale)

##################################################################################
class MoverNH2Flip:
  def __init__(self, nh2Atom, caAtomName, bondedNeighborLists):
    """Constructs a Mover that will handle flipping an NH2 with an O, both of which
       are attached to the same Carbon atom (and each of which has no other bonds).
       This Mover uses a simple swap of the center positions of the heavy atoms (with
       repositioning of the Hydrogens to lie in the plane with the other three atoms)
       for its testing, but during FixUp it adjusts the bond lengths of the Oxygen
       and the Nitrogen; per Protein Science Vol 27:293-315.
       This handles flips for Asn and Gln.
       :param nh2Atom: Nitrogen atom that is attached to two Hydrogens.
       :param ca2AtomName: Name of the Alpha Carbon for the amino acid -- the atom around which the
       rigid body is rotated for the final docking motion.  All atoms linked to before
       this one is reached will be included in the list of movable atoms -- there is one
       more of these for Gln than for Asn.
       :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
       structure that the atom from the first parameter interacts with that lists all of the
       bonded atoms.  Can be obtained by calling getBondedNeighborLists().
    """

    # Verify that we've been run on a valid structure and get a list of all of the
    # atoms up to and including the pivot atom.
    if nh2Atom.element != "N":
      raise ValueError("MoverNH2Flip(): nh2Atom is not a Nitrogen")
    partners = bondedNeighborLists[nh2Atom]
    if len(partners) != 3:
      raise ValueError("MoverNH2Flip(): nh2Atom does not have three bonded neighbors")
    nh2Hydrogens = []
    hinge = None
    for a in partners:
      if a.element == "H":
        nh2Hydrogens.append(a)
      else:
        hinge = a
    if len(nh2Hydrogens) != 2:
      raise ValueError("MoverNH2Flip(): nh2Atom does not have two bonded hydrogens")
    if hinge is None:
      raise ValueError("MoverNH2Flip(): nh2Atom does not have bonded (hinge) Carbon friend")

    bonded = bondedNeighborLists[hinge]
    oyxgen = None
    pivot = None
    for b in bonded:
      if b.element == "O":
        oxygen = b
      elif b.element == "C":
        pivot = b
    if pivot is None:
      raise ValueError("MoverNH2Flip(): Hinge does not have bonded (pivot) Carbon friend")
    if oxygen is None:
      raise ValueError("MoverNH2Flip(): Hinge does not have bonded oxygen friend")
    if len(bondedNeighborLists[oxygen]) != 1:
      raise ValueError("MoverNH2Flip(): Oxygen has more than one bonded neighbor")
    self._atoms = [ nh2Hydrogens[0], nh2Hydrogens[1], nh2Atom, oxygen, hinge, pivot ]

    # Find any linking atoms between the pivot atom and the Alpha Carbon and add them to the list.
    # Also add hydrogens that are bound to the pivot and linker atoms.
    linkers = []
    prevID = hinge.i_seq
    foundAlpha = False
    cur = pivot
    linkerHydrogens = []
    while not foundAlpha:
      bonded = bondedNeighborLists[cur]
      link = None
      if len(bonded) != 4:
        raise ValueError("MoverNH2Flip(): Linker chain has an element with other than four bonds")
      for b in bonded:
        if b.i_seq == prevID:
          continue
        elif b.element == "C":
          link = b
        elif b.element == "H":
          linkerHydrogens.append(b)
      if link is None:
        raise ValueError("MoverNH2Flip(): Did not find Carbon in linker chain step")
      if link.name.strip().upper() == caAtomName.strip().upper():
        caAtom = link
        foundAlpha = True
      else:
        linkers.append(link)
        prevID = cur.i_seq
        cur = link
    if len(linkerHydrogens) != 2*(len(linkers)+1):
      raise ValueError("MoverNH2Flip(): Linker carbons do not have two Hydrogens each: "+
        str(len(linkers))+","+str(len(linkerHydrogens)))
    self._atoms.extend(linkers)
    self._atoms.extend(linkerHydrogens)

    #########################
    # Compute the new positions for the Hydrogens such that they are at the same distance from
    # the Oxygen as one of them is from the Nitrogen and located at +/-120 degrees from the
    # Carbon-Oxygen bond in the plane of the Nitrogen, Carbon, and Oxygen.
    cToO = _lvec3(oxygen.xyz) - _lvec3(hinge.xyz)
    nToO = _rvec3(nh2Atom.xyz) - _rvec3(hinge.xyz)

    # Normal to the plane containing Nitrogen, Carbon, and Oxygen
    normal = _lvec3(scitbx.matrix.cross_product_matrix(cToO) * nToO).normalize()

    hBond0Len = (_rvec3(nh2Hydrogens[0].xyz) - _rvec3(nh2Atom.xyz)).length()
    hBond1Len = (_rvec3(nh2Hydrogens[1].xyz) - _rvec3(nh2Atom.xyz)).length()
    newH0 = _lvec3(oxygen.xyz) + ((-cToO.normalize()) * hBond0Len).rotate_around_origin(normal, 120 * math.pi/180)
    newH1 = _lvec3(oxygen.xyz) + ((-cToO.normalize()) * hBond1Len).rotate_around_origin(normal,-120 * math.pi/180)

    #########################
    # Compute the list of positions for all of the atoms. This consists of the original
    # location and the flipped location where we swap the locations of the two heavy atoms
    # and bring the Hydrogens along for the ride.

    startPos = []
    for a in self._atoms:
      startPos.append(a.xyz)

    newPos = startPos.copy()
    newPos[0] = newH0
    newPos[1] = newH1
    newPos[2] = oxygen.xyz
    newPos[3] = nh2Atom.xyz

    self._coarsePositions = [ startPos, newPos ]

    #########################
    # Compute the list of Fixup returns.
    hingeIndex = 4
    firstDockIndex = 3
    secondDockIndex = 2
    movable = _rotateHingeDock(self._atoms, hingeIndex, firstDockIndex, secondDockIndex, caAtom)

    # No fix-up for coarse position 0, do the above adjustment for position 1
    self._fixUpPositions = [ [], movable ]
    
  def CoarsePositions(self):
    # returns: The two possible coarse positions with 0 energy offset for either.
    return PositionReturn(self._atoms, self._coarsePositions, [0.0, 0.0])

  def FinePositions(self, coarseIndex):
    # returns: No fine positions for any coarse position.
    return PositionReturn([], [], [])

  def FixUp(self, coarseIndex):
    # Return the appropriate fixup
    return FixUpReturn(self._atoms, self._fixUpPositions[coarseIndex])

##################################################################################
class MoverHistidineFlip:
  def __init__(self, ne2Atom, bondedNeighborLists):
    """Constructs a Mover that will handle flipping a Histidine ring.
       This Mover uses a simple swap of the center positions of the heavy atoms (with
       repositioning of the Hydrogens to lie in the same directions)
       for its testing, but during FixUp it adjusts the bond lengths per
       Protein Science Vol 27:293-315.
       :param ne2Atom: NE2 atom within the Histidine ring.
       :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
       structure that the atom from the first parameter interacts with that lists all of the
       bonded atoms.  Can be obtained by calling getBondedNeighborLists().
    """

    # Verify that we've been run on a valid structure and get a list of all of the
    # atoms up to and including the pivot atom.
    if ne2Atom.element != "N":
      raise ValueError("MoverHistidineFlip(): ne2Atom is not a Nitrogen")
    partners = bondedNeighborLists[ne2Atom]
    if len(partners) != 3:
      raise ValueError("MoverHistidineFlip(): ne2Atom does not have three bonded neighbors")
    hydrogens = []
    carbons = []
    for a in partners:
      if a.element == "H":
        hydrogens.append(a)
      elif a.element == "C":
        carbons.append(a)
    if len(hydrogens) != 1:
      raise ValueError("MoverHistidineFlip(): ne2Atom does not have one bonded hydrogen")
    if len(carbons) != 2:
      raise ValueError("MoverHistidineFlip(): ne2Atom does not have two bonded carbons")
    ne2HAtom = hydrogens[0]

    # Determine if the first Carbon is CE1, which is bonded to two Carbons.  If not,
    # then the second one must be.  Fill in CD2 and CE1 based on this information, then we
    # can check them and continue parsing.
    cTest = carbons[0]
    bonded = bondedNeighborLists[cTest]
    if len(bonded) != 3:
      raise ValueError("MoverHistidineFlip(): ne2Atom neighbor does not have three bonded neighbors")
    cs = 0
    for b in bonded:
      if b.element == "C":
        cs += 1
    if cs == 0:
      ce1Atom = carbons[0]
      cd2Atom = carbons[1]
    else:
      ce1Atom = carbons[1]
      cd2Atom = carbons[0]

    # Find the Hydrogen on CE1.
    bonded = bondedNeighborLists[ce1Atom]
    if len(bonded) != 3:
      raise ValueError("MoverHistidineFlip(): CE1 does not have three bonded neighbors")
    ce1HAtom = None
    for b in bonded:
      if b.element == "H":
        ce1HAtom = b
    if ce1HAtom is None:
      raise ValueError("MoverHistidineFlip(): Could not find Hydrogen attached to CE1")

    # Find the Hydrogen on CD2.
    bonded = bondedNeighborLists[cd2Atom]
    if len(bonded) != 3:
      raise ValueError("MoverHistidineFlip(): CD2 does not have three bonded neighbors")
    cd2HAtom = None
    for b in bonded:
      if b.element == "H":
        cd2HAtom = b
    if cd2HAtom is None:
      raise ValueError("MoverHistidineFlip(): Could not find Hydrogen attached to CD2")

    # Find CG on the other side of CD2.
    cgAtom = None
    bonded = bondedNeighborLists[cd2Atom]
    if len(bonded) < 2:
      raise ValueError("MoverHistidineFlip(): CD2 does not have at least two bonded neighbors")
    for b in bonded:
      if b.i_seq != cd2Atom.i_seq and b.element == "C":
        cgAtom = b
    if cgAtom is None:
      raise ValueError("MoverHistidineFlip(): Could not find CG")

    # Find CB and ND1 on the other side of CG
    cbAtom = None
    nd1Atom = None
    bonded = bondedNeighborLists[cgAtom]
    if len(bonded) != 3:
      raise ValueError("MoverHistidineFlip(): CG does not have three bonded neighbors")
    for b in bonded:
      if b.i_seq == cgAtom.i_seq:
        continue
      elif b.element == "N":
        nd1Atom = b
      elif b.element == "C":
        cbAtom = b
    if nd1Atom is None:
      raise ValueError("MoverHistidineFlip(): Could not find ND1")
    if cbAtom is None:
      raise ValueError("MoverHistidineFlip(): Could not find CB")

    # Find the Hydrogen attached to ND1
    bonded = bondedNeighborLists[nd1Atom]
    if len(bonded) != 3:
      raise ValueError("MoverHistidineFlip(): ND1 does not have three bonded neighbors")
    nd1HAtom = None
    for b in bonded:
      if b.element == "H":
        nd1HAtom = b
    if nd1HAtom is None:
      raise ValueError("MoverHistidineFlip(): ND1 does not have a bonded hydrogen")

    # Find CA on the other side of CB and find CBs Hydrogens
    caAtom = None
    cbHydrogens = []
    bonded = bondedNeighborLists[cbAtom]
    if len(bonded) != 4:
      raise ValueError("MoverHistidineFlip(): CB does not have four bonded neighbors")
    for b in bonded:
      if b.i_seq == cgAtom.i_seq:
        continue
      elif b.element == "C":
        caAtom = b
      elif b.element == "H":
        cbHydrogens.append(b)
    if caAtom is None:
      raise ValueError("MoverHistidineFlip(): Could not find CA")
    if len(cbHydrogens) != 2:
      raise ValueError("MoverHistidineFlip(): Could not find Hydrogens on CB")

    self._atoms = [ ne2Atom, ne2HAtom, ce1Atom, ce1HAtom, nd1Atom, nd1HAtom, cd2Atom, cd2HAtom,
      cgAtom, cbAtom, cbHydrogens[0], cbHydrogens[1] ]

    #########################
    # Compute the new positions for the Hydrogens such that they are at the same distance from
    # their swapped parent atoms and in the direction of the Hydrogens from the original atoms at
    # each location.  We swap ND1 with CD2 and NE2 with CE1.
    nd1HVec = _lvec3(nd1HAtom.xyz) - _lvec3(nd1Atom.xyz)
    ne2HVec = _lvec3(ne2HAtom.xyz) - _lvec3(ne2Atom.xyz)
    ce1HVec = _lvec3(ce1HAtom.xyz) - _lvec3(ce1Atom.xyz)
    cd2HVec = _lvec3(cd2HAtom.xyz) - _lvec3(cd2Atom.xyz)

    nd1HNew = _lvec3(cd2Atom.xyz) + nd1HVec.length() * cd2HVec.normalize()
    cd2HNew = _lvec3(nd1Atom.xyz) + cd2HVec.length() * nd1HVec.normalize()
    ce1HNew = _lvec3(ne2Atom.xyz) + ce1HVec.length() * ne2HVec.normalize()
    ne2HNew = _lvec3(ce1Atom.xyz) + ne2HVec.length() * ce1HVec.normalize()

    #########################
    # Compute the list of positions for all of the atoms. This consists of the original
    # location and the flipped location where we swap the locations of the two pairs of heavy atoms
    # and bring the Hydrogens along for the ride.

    startPos = []
    for a in self._atoms:
      startPos.append(a.xyz)

    newPos = startPos.copy()
    newPos[0] = ce1Atom.xyz   # ne2 swapped to this location
    newPos[1] = ne2HNew
    newPos[2] = ne2Atom.xyz   # ce1 swapped to this location
    newPos[3] = ce1HNew
    newPos[4] = cd2Atom.xyz   # nd1 swapped to this location
    newPos[5] = nd1HNew
    newPos[6] = nd1Atom.xyz   # cd2 swapped to this location
    newPos[7] = cd2HNew

    self._coarsePositions = [ startPos, newPos ]

    #########################
    # Compute the list of Fixup returns.
    hingeIndex = 8
    firstDockIndex = 0
    secondDockIndex = 2
    movable = _rotateHingeDock(self._atoms, hingeIndex, firstDockIndex, secondDockIndex, caAtom)

    # No fix-up for coarse position 0, do the above adjustment for position 1
    self._fixUpPositions = [ [], movable ]

  def CoarsePositions(self):
    # returns: The two possible coarse positions with 0 energy offset for either.
    return PositionReturn(self._atoms, self._coarsePositions, [0.0, 0.0])

  def FinePositions(self, coarseIndex):
    # returns: No fine positions for any coarse position.
    return PositionReturn([], [], [])

  def FixUp(self, coarseIndex):
    # Return the appropriate fixup
    return FixUpReturn(self._atoms, self._fixUpPositions[coarseIndex])


##################################################################################
# Test function to verify that all Movers behave properly.

def Test():
  """Test function for all classes provided above.
  :returns Empty string on success, string describing the problem on failure.
  :returns Empty string on success, string describing the problem on failure.
  """

  # Test the _MoverRotator class.
  try:
    # Construct a _MoverRotator with three atoms, each at +1 in Z with one at +1 in X and 0 in Y
    # and the other two at +/-1 in Y and 0 in X.  It will rotate around the Z axis with an offset
    # of -2 for the axis start location by 180 degrees with a coarse step size of 90 and a
    # preference function that always returns 1.
    atoms = []
    coords = [[ 1.0, 0.0, 1.0],
              [ 0.0, 1.0, 1.0],
              [ 0.0,-1.0, 1.0]]
    for i in range(3):
      a = pdb.hierarchy.atom()
      a.xyz = coords[i]
      atoms.append(a)
    axis = flex.vec3_double([ [ 0.0, 0.0,-2.0], [ 0.0, 0.0, 1.0] ])
    def prefFunc(x):
      return 1.0
    rot = _MoverRotator(atoms,axis, 180, 90.0, True, preferenceFunction = prefFunc)

    # See if the results of each of the functions are what we expect in terms of sizes and locations
    # of atoms and preferences.  We'll use the default fine step size.
    # The first coarse rotation should be by -90 degrees, moving the first atom to (0, -1, 1)
    coarse = rot.CoarsePositions()
    if len(coarse.atoms) != 3:
      return "Movers.Test() _MoverRotator basic: Expected 3 atoms for CoarsePositions, got "+str(len(coarse.atoms))
    atom0pos1 = coarse.positions[1][0]
    if (_lvec3(atom0pos1) - _lvec3([0,-1,1])).length() > 1e-5:
      return "Movers.Test() _MoverRotator basic: Expected location = (0,-1,1), got "+str(atom0pos1)

    # The first fine rotation (index 0) around the second coarse index (index 1) should be to -91 degrees,
    # moving the first atom to the appropriately rotated location around the Z axis
    rad = -91 / 180 * math.pi
    x = math.cos(rad)
    y = math.sin(rad)
    z = 1
    fine = rot.FinePositions(1)
    atom0pos1 = fine.positions[0][0]
    if (_lvec3(atom0pos1) - _lvec3([x,y,z])).length() > 1e-5:
      return "Movers.Test() _MoverRotator basic: Expected fine location = "+str([x,y,z])+", got "+str(atom0pos1)

    # The preference function should always return 1.
    for p in fine.preferenceEnergies:
      if p != 1:
        return "Movers.Test() _MoverRotator basic: Expected preference energy = 1, got "+str(p)

    # Test different preference scale value.
    rot = _MoverRotator(atoms,axis, 180, 90.0, True, preferenceFunction = prefFunc, preferredOrientationScale = 2)
    coarse = rot.CoarsePositions()
    for p in coarse.preferenceEnergies:
      if p != 2:
        return "Movers.Test() _MoverRotator Scaled preference: Expected preference energy = 2, got "+str(p)

    # Test None preference function and a sinusoidal one.
    rot = _MoverRotator(atoms,axis, 180, 90.0, True, preferenceFunction = None)
    coarse = rot.CoarsePositions()
    for p in coarse.preferenceEnergies:
      if p != 0:
        return "Movers.Test() _MoverRotator None preference function: Expected preference energy = 0, got "+str(p)
    def prefFunc2(x):
      return math.cos(x * math.pi / 180)
    rot = _MoverRotator(atoms,axis, 180, 180, True, preferenceFunction = prefFunc2)
    coarse = rot.CoarsePositions()
    expected = [1, -1]
    for i,p  in enumerate(coarse.preferenceEnergies):
      val = expected[i]
      if p != val:
        return "Movers.Test() _MoverRotator Sinusoidal preference function: Expected preference energy = "+str(val)+", got "+str(p)

    # Test coarseStepDegrees default behavior.
    rot = _MoverRotator(atoms,axis, 180)
    coarse = rot.CoarsePositions()
    if len(coarse.positions) != 12:
      return "Movers.Test() _MoverRotator Default coarse step: Expected 12, got "+str(len(coarse.positions))

    # Test fineStepDegrees setting.
    rot = _MoverRotator(atoms,axis, 180, fineStepDegrees = 2)
    fine = rot.FinePositions(0)
    # +/- 15 degrees in 1-degree steps, but we don't do the +15 because it will be handled by the next
    # rotation up.
    if len(fine.positions) != 14:
      return "Movers.Test() _MoverRotator setting fine step: Expected 14, got "+str(len(fine.positions))

    # Test doFineRotations = False and 180 degree coarseStepDegrees.
    rot = _MoverRotator(atoms,axis, 180, 180, False)
    coarse = rot.CoarsePositions()
    if len(coarse.positions) != 2:
      return "Movers.Test() _MoverRotator 180 coarse steps: Expected 2, got "+str(len(coarse.positions))
    fine = rot.FinePositions(0)
    if len(fine.positions) != 0:
      return "Movers.Test() _MoverRotator 180 coarse steps: Expected 0, got "+str(len(fine.positions))

  except Exception as e:
    return "Movers.Test() _MoverRotator basic: Exception during test of _MoverRotator: "+str(e)+"\n"+traceback.format_exc()

  # Test the MoverSingleHydrogenRotator class.
  try:
    # Construct a MoverSingleHydrogenRotator that has an atom that starts out at 45 degrees around Z
    # that is bonded to a neighbor and partner that are vertical and then partner is bonded to two
    # friends that are in the Y=0 plane.  This should cause us to get the atom rotated to lie in
    # the Y=0 plane at a distance of sqrt(2) and a fitness function that prefers the orientations
    # that are in this plane.
    h = pdb.hierarchy.atom()
    h.element = "H"
    h.xyz = [ 1.0, 1.0, 1.0 ]

    n = pdb.hierarchy.atom()
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.xyz = [-1.0, 0.0,-2.0 ]

    # Build the hierarchy so we can reset the i_seq values.
    ag = pdb.hierarchy.atom_group()
    ag.append_atom(h)
    ag.append_atom(n)
    ag.append_atom(p)
    ag.append_atom(f1)
    ag.append_atom(f2)
    rg = pdb.hierarchy.residue_group()
    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    bondedNeighborLists = {}
    bondedNeighborLists[h] = [ n ]
    bondedNeighborLists[n] = [ h, p ]
    bondedNeighborLists[p] = [ n, f1, f2 ]
    bondedNeighborLists[f1] = [ p ]
    bondedNeighborLists[f2] = [ p ]

    # Add a non-bonded potential acceptor atom at 13 degrees rotation towards the Y axis from
    # the X axis.
    acc = pdb.hierarchy.atom()
    acc.xyz = [ math.cos(13*math.pi/180), math.sin(13*math.pi/180), 1.0 ]

    mover = MoverSingleHydrogenRotator(h, bondedNeighborLists, [acc])

    # Check for hydrogen rotated into Y=0 plane at a distance of sqrt(2) from Z axis
    if h.xyz[2] != 1 or abs(abs(h.xyz[0])-math.sqrt(2)) > 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator pair: bad H placement"

    # Check fitness function preferring 0 and 180 rotations
    zero = mover._preferenceFunction(0)
    ninety = mover._preferenceFunction(90)
    oneEighty = mover._preferenceFunction(180)
    if abs(zero - oneEighty) > 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator pair: bad preference function"
    if zero - ninety < 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator pair: bad preference function"

    # Check that one of the orientations has a dihedral angle of 0 compared to the original
    # 13-degree-off potential acceptor atom.
    found = False
    for pos in mover.CoarsePositions().positions:
      sites = [ pos[0], p.xyz, n.xyz, acc.xyz ]
      dihedral = scitbx.math.dihedral_angle(sites=sites, deg=True)
      if abs(dihedral) < 1e-5:
        found = True
        break
    if not found:
      return "Movers.Test() MoverSingleHydrogenRotator pair: no orientation towards acceptor"

  except Exception as e:
    return "Movers.Test() MoverSingleHydrogenRotator pair: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  try:
    # Construct a MoverSingleHydrogenRotator that has an atom that starts out at 45 degrees around Z
    # that is bonded to a neighbor and partner that are vertical and then partner is bonded to three
    # friends one of which is on the +X axis and the other two are +/-120 away from it.
    # This should cause us to get the atom rotated to be between two of the friends and a fitness
    # function that prefers its location and +/-120 degrees from it.
    h = pdb.hierarchy.atom()
    h.element = "H"
    h.xyz = [ 1.0, 1.0, 1.0 ]

    n = pdb.hierarchy.atom()
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.xyz = _rotateAroundAxis(f1, axis, -120)

    f3 = pdb.hierarchy.atom()
    f3.xyz = _rotateAroundAxis(f1, axis, 120)

    # Build the hierarchy so we can reset the i_seq values.
    ag = pdb.hierarchy.atom_group()
    ag.append_atom(h)
    ag.append_atom(n)
    ag.append_atom(p)
    ag.append_atom(f1)
    ag.append_atom(f2)
    ag.append_atom(f3)
    rg = pdb.hierarchy.residue_group()
    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    bondedNeighborLists = {}
    bondedNeighborLists[h] = [ n ]
    bondedNeighborLists[n] = [ h, p ]
    bondedNeighborLists[p] = [ n, f1, f2, f3 ]
    bondedNeighborLists[f1] = [ p ]
    bondedNeighborLists[f2] = [ p ]
    bondedNeighborLists[f3] = [ p ]

    mover = MoverSingleHydrogenRotator(h, bondedNeighborLists)

    # Check for a hydrogen on the -X axis at a distance of sqrt(2) from Z axis,
    # or +/-120 from there.
    angles = [0, 120, -120]
    angle = None
    for a in angles:
      rotated = _rotateAroundAxis(h, axis, a)
      if rotated[2] == 1 and rotated[0]+math.sqrt(2) < 1e-5:
        angle = a
    if angle == None:
      return "Movers.Test() MoverSingleHydrogenRotator triple: bad H placement"

    # Check fitness function preferring 180 and +/- 120 from there rotations away from
    # the angle away from 180 degrees.
    zero = mover._preferenceFunction(angle+0)
    oneEighty = mover._preferenceFunction(angle+180)
    off1 = mover._preferenceFunction(angle+180+120)
    off2 = mover._preferenceFunction(angle+180-120)
    if abs(off1 - oneEighty) > 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator triple: bad preference function"
    if abs(off2 - oneEighty) > 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator triple: bad preference function"
    if zero - oneEighty < 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator triple: bad preference function"

  except Exception as e:
    return "Movers.Test() MoverSingleHydrogenRotator triple: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  # Test the MoverNH3Rotator class.
  try:
    # Construct a MoverNH3Rotator that has one hydrogen start out at 45 degrees around Z and the
    # other two at +/-120 degrees from that one.
    # They are bonded to a Nitrogen and partner that are vertical and then partner is bonded to three
    # friends with one on the +X axis and the others +/-120.  This should cause us to get the hydrogens at
    # 180 and 120 away from that and a fitness function that prefers the orientations 120 degrees
    # apart.
    axis = flex.vec3_double([ [0,0,0], [0,0,1] ])
    h1 = pdb.hierarchy.atom()
    h1.element = "H"
    h1.xyz = [ 1.0, 1.0, 1.0 ]

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.xyz = _rotateAroundAxis(h1, axis, -120)

    h3 = pdb.hierarchy.atom()
    h3.element = "H"
    h3.xyz = _rotateAroundAxis(h1, axis, 120)

    n = pdb.hierarchy.atom()
    n.element = "N"
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.xyz = _rotateAroundAxis(f1, axis, -120)

    f3 = pdb.hierarchy.atom()
    f3.xyz = _rotateAroundAxis(f1, axis, 120)

    # Build the hierarchy so we can reset the i_seq values.
    ag = pdb.hierarchy.atom_group()
    ag.append_atom(h1)
    ag.append_atom(h2)
    ag.append_atom(h3)
    ag.append_atom(n)
    ag.append_atom(p)
    ag.append_atom(f1)
    ag.append_atom(f2)
    ag.append_atom(f3)
    rg = pdb.hierarchy.residue_group()
    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    bondedNeighborLists = {}
    bondedNeighborLists[h1] = [ n ]
    bondedNeighborLists[h2] = [ n ]
    bondedNeighborLists[h3] = [ n ]
    bondedNeighborLists[n] = [ h1, h2, h3, p ]
    bondedNeighborLists[p] = [ n, f1, f2, f3 ]
    bondedNeighborLists[f1] = [ p ]
    bondedNeighborLists[f2] = [ p ]
    bondedNeighborLists[f3] = [ p ]

    mover = MoverNH3Rotator(n, bondedNeighborLists)

    # Check for a hydrogen on the -X axis at a distance of sqrt(2) from Z axis
    found = False
    for h in [h1, h2, h3]:
      if h.xyz[2] == 1 and h.xyz[0]+math.sqrt(2) < 1e-5:
        found = True
    if not found:
      return "Movers.Test() MoverNH3Rotator basic: bad H placement"

    # Check fitness function preferring 180 and +/- 120 from there rotations
    zero = mover._preferenceFunction(0)
    oneEighty = mover._preferenceFunction(180)
    off1 = mover._preferenceFunction(180+120)
    off2 = mover._preferenceFunction(180-120)
    if abs(off1 - oneEighty) > 1e-5:
      return "Movers.Test() MoverNH3Rotator basic: bad preference function"
    if abs(off2 - oneEighty) > 1e-5:
      return "Movers.Test() MoverNH3Rotator basic: bad preference function"
    if zero - oneEighty < 1e-5:
      return "Movers.Test() MoverNH3Rotator basic: bad preference function"

  except Exception as e:
    return "Movers.Test() MoverNH3Rotator basic: Exception during test: "+str(e)+"\n"+traceback.format_exc()
    
  # Test the MoverAromaticMethylRotator class.
  try:
    # Construct a MoverAromaticMethylRotator that has one hydrogen start out at 45 degrees around Z and the
    # other two at +/-120 degrees from that one.
    # They are bonded to a Carbon and partner that are vertical and then partner is bonded to two
    # friends that are in the Y=0 plane.  This should cause us to get the one of the hydrogens at
    # +90 or -90 and the others 120 away from the first with only two coarse choices and no fine choices.
    axis = flex.vec3_double([ [0,0,0], [0,0,1] ])
    h1 = pdb.hierarchy.atom()
    h1.element = "H"
    h1.xyz = [ 1.0, 1.0, 1.0 ]

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.xyz = _rotateAroundAxis(h1, axis, -120)

    h3 = pdb.hierarchy.atom()
    h3.element = "H"
    h3.xyz = _rotateAroundAxis(h1, axis, 120)

    n = pdb.hierarchy.atom()
    n.element = "C"
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f1.xyz = [-1.0, 0.0,-2.0 ]

    # Build the hierarchy so we can reset the i_seq values.
    ag = pdb.hierarchy.atom_group()
    ag.append_atom(h1)
    ag.append_atom(h2)
    ag.append_atom(h3)
    ag.append_atom(n)
    ag.append_atom(p)
    ag.append_atom(f1)
    ag.append_atom(f2)
    rg = pdb.hierarchy.residue_group()
    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    bondedNeighborLists = {}
    bondedNeighborLists[h1] = [ n ]
    bondedNeighborLists[h2] = [ n ]
    bondedNeighborLists[h3] = [ n ]
    bondedNeighborLists[n] = [ h1, h2, h3, p ]
    bondedNeighborLists[p] = [ n, f1, f2 ]
    bondedNeighborLists[f1] = [ p ]
    bondedNeighborLists[f2] = [ p ]

    mover = MoverAromaticMethylRotator(n, bondedNeighborLists)

    # Check for a hydrogen on the +/-Y axis at a distance of sqrt(2) from the Z axis
    found = False
    for h in [h1, h2, h3]:
      if h.xyz[2] == 1 and abs(abs(h.xyz[1])-math.sqrt(2)) < 1e-5:
        found = True
    if not found:
      return "Movers.Test() MoverAromaticMethylRotator basic: bad H placement"

    # Check that we get two coarse and no fine orientations
    coarse = mover.CoarsePositions().positions
    if len(coarse) != 2:
      return "Movers.Test() MoverAromaticMethylRotator basic: bad coarse count: "+str(len(coarse))
    fine = mover.FinePositions(0).positions
    if len(fine) != 0:
      return "Movers.Test() MoverAromaticMethylRotator basic: bad fine count: "+str(len(fine))

  except Exception as e:
    return "Movers.Test() MoverAromaticMethylRotator basic: Exception during test: "+str(e)+"\n"+traceback.format_exc()
    
  # Test the MoverTetrahedralMethylRotator class.
  try:
    # Construct a MoverTetrahedralMethylRotator that has one hydrogen start out at 45 degrees around Z and the
    # other two at +/-120 degrees from that one.
    # They are bonded to a Carbon and partner that are vertical and then partner is bonded to three
    # friends with one in the +X direction.  This should cause us to get the one of the hydrogens at
    # 180 the others 120 away from the first.
    axis = flex.vec3_double([ [0,0,0], [0,0,1] ])
    h1 = pdb.hierarchy.atom()
    h1.element = "H"
    h1.xyz = [ 1.0, 1.0, 1.0 ]

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.xyz = _rotateAroundAxis(h1, axis, -120)

    h3 = pdb.hierarchy.atom()
    h3.element = "H"
    h3.xyz = _rotateAroundAxis(h1, axis, 120)

    n = pdb.hierarchy.atom()
    n.element = "C"
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.xyz = _rotateAroundAxis(f1, axis, -120)

    f3 = pdb.hierarchy.atom()
    f3.xyz = _rotateAroundAxis(f1, axis,  120)

    # Build the hierarchy so we can reset the i_seq values.
    ag = pdb.hierarchy.atom_group()
    ag.append_atom(h1)
    ag.append_atom(h2)
    ag.append_atom(h3)
    ag.append_atom(n)
    ag.append_atom(p)
    ag.append_atom(f1)
    ag.append_atom(f2)
    ag.append_atom(f3)
    rg = pdb.hierarchy.residue_group()
    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    bondedNeighborLists = {}
    bondedNeighborLists[h1] = [ n ]
    bondedNeighborLists[h2] = [ n ]
    bondedNeighborLists[h3] = [ n ]
    bondedNeighborLists[n] = [ h1, h2, h3, p ]
    bondedNeighborLists[p] = [ n, f1, f2, f3 ]
    bondedNeighborLists[f1] = [ p ]
    bondedNeighborLists[f2] = [ p ]
    bondedNeighborLists[f3] = [ p ]

    mover = MoverTetrahedralMethylRotator(n, bondedNeighborLists)

    # Check for a hydrogen on the +/-Y axis at a distance of sqrt(2) from the Z axis
    found = False
    for h in [h1, h2, h3]:
      if h.xyz[2] == 1 and h.xyz[0]+math.sqrt(2) < 1e-5:
        found = True
    if not found:
      return "Movers.Test() MoverTetrahedralMethylRotator basic: bad H placement"

    # Check fitness function preferring 180 and +/- 120 from there rotations.
    zero = mover._preferenceFunction(0)
    oneEighty = mover._preferenceFunction(180)
    off1 = mover._preferenceFunction(180+120)
    off2 = mover._preferenceFunction(180-120)
    if abs(off1 - oneEighty) > 1e-5:
      return "Movers.Test() MoverTetrahedralMethylRotator: bad preference function"
    if abs(off2 - oneEighty) > 1e-5:
      return "Movers.Test() MoverTetrahedralMethylRotator: bad preference function"
    if zero - oneEighty < 1e-5:
      return "Movers.Test() MoverTetrahedralMethylRotator: bad preference function"

  except Exception as e:
    return "Movers.Test() MoverTetrahedralMethylRotator basic: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  # Test the MoverNH2Flip class with no linker (similar to Asn).
  try:
    # Test behavior with offsets for each atom so that we test the generic case.
    # Construct a MoverNH2Flip that has the N, H's and Oxygen located (non-physically)
    # slightly in the +Y direction out of the X-Z plane, with the Hydrogens 120 around the
    # same offset axis.
    # They are bonded to a Carbon and friend (pivot) and alpha carbon that are on the Z axis.
    # Non-physically, atoms are 1 unit apart.
    p = pdb.hierarchy.atom()
    p.xyz = [ 0.0, 0.0, 0.0 ]

    f = pdb.hierarchy.atom()
    f.element = "C"
    f.xyz = [ 0.0, 0.0,-1.0 ]

    fh1 = pdb.hierarchy.atom()
    fh1.element = "H"
    fh1.xyz = [ 1.0, 0.0,-1.0 ]

    fh2 = pdb.hierarchy.atom()
    fh2.element = "H"
    fh2.xyz = [-1.0, 0.0,-1.0 ]

    ca = pdb.hierarchy.atom()
    ca.element = "C"
    ca.name = "CA"
    ca.xyz = [ 0.0, 0.0,-2.0 ]

    # Nitrogen and Oxygen are +/-120 degrees from carbon-carbon bond
    axis = flex.vec3_double([ [0,0,0], [0,1,0] ])
    n = pdb.hierarchy.atom()
    n.element = "N"
    n.xyz = _rotateAroundAxis(f, axis,-120) + _lvec3([0,0.01,0]) + _lvec3([ 0.002, 0.003,-0.004])

    o = pdb.hierarchy.atom()
    o.element = "O"
    o.xyz = _rotateAroundAxis(f, axis, 120) + _lvec3([0,0.01,0]) + _lvec3([-0.003, 0.002, 0.003])

    # Hydrogens are +/-120 degrees from nitrogen-carbon bond
    axis = flex.vec3_double([ n.xyz, [0,1,0] ])
    h1 = pdb.hierarchy.atom()
    h1.element = "H"
    h1.xyz = _rotateAroundAxis(p, axis,-120) + _lvec3([0,0.01,0]) + _lvec3([-0.008, 0.001, 0.008])

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.xyz = _rotateAroundAxis(p, axis, 120) + _lvec3([0,0.01,0]) + _lvec3([ 0.007,-0.001, 0.007])

    # Build the hierarchy so we can reset the i_seq values.
    ag = pdb.hierarchy.atom_group()
    ag.append_atom(h1)
    ag.append_atom(h2)
    ag.append_atom(n)
    ag.append_atom(o)
    ag.append_atom(p)
    ag.append_atom(f)
    ag.append_atom(fh1)
    ag.append_atom(fh2)
    ag.append_atom(ca)
    rg = pdb.hierarchy.residue_group()
    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    bondedNeighborLists = {}
    bondedNeighborLists[h1] = [ n ]
    bondedNeighborLists[h2] = [ n ]
    bondedNeighborLists[n] = [ h1, h2, p ]
    bondedNeighborLists[o] = [ p ]
    bondedNeighborLists[p] = [ n, o, f ]
    bondedNeighborLists[f] = [ p, ca, fh1, fh2 ]
    bondedNeighborLists[fh1] = [ f ]
    bondedNeighborLists[fh2] = [ f ]
    bondedNeighborLists[ca] = [ f ]

    mover = MoverNH2Flip(n, ca.name, bondedNeighborLists)

    # Ensure that the coarse-flip results meet the expections:
    # 1) N and O are flipped in position
    # 2) H remain at the same distance from the new N.

    coarse = mover.CoarsePositions()
    if len(coarse.positions) != 2:
      return "Movers.Test() MoverNH2Flip basic: Did not find two locations: "+str(len(coarse.positions))
    newPos = coarse.positions[1]
    dist = (_lvec3(newPos[2]) - _lvec3(o.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverNH2Flip basic: Nitrogen moved incorrectly: "+str(dist)
    dist = (_lvec3(newPos[3]) - _lvec3(n.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverNH2Flip basic: Oxygen moved incorrectly: "+str(dist)

    dHydrogen = (_lvec3(newPos[0]) - _lvec3(newPos[2])).length()
    oldDHydrogen = (_lvec3(h1.xyz)-_lvec3(n.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverNH2Flip basic: Bad coarse hydrogen1 motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (_lvec3(newPos[1]) - _lvec3(newPos[2])).length()
    oldDHydrogen = (_lvec3(h2.xyz)-_lvec3(n.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverNH2Flip basic: Bad coarse hydrogen2 motion: "+str(dHydrogen-oldDHydrogen)

    # Ensure that the Fixup results meet the specifications:
    # 1) New Oxygen on the line from the alpha carbon to the old Nitrogen
    # 2) New plane of Oxygen, Nitrogen, Alpha Carbon matches old plane, but flipped
    # 3) Carbons and pivot Hydrogens move slightly due to rigid-body motion

    fixed = mover.FixUp(1).newPositions
    newODir = (fixed[3] - _lvec3(f.xyz)).normalize()
    oldNDir = (_rvec3(n.xyz) - _rvec3(f.xyz)).normalize()
    if (newODir * oldNDir)[0] < 0.9999:
      return "Movers.Test() MoverNH2Flip basic: Bad oxygen alignment: "+str((newODir * oldNDir)[0])

    newNDir = (fixed[2] - _lvec3(f.xyz)).normalize()
    oldODir = (_rvec3(o.xyz) - _rvec3(f.xyz)).normalize()
    newNormal = (scitbx.matrix.cross_product_matrix(_lvec3(newNDir)) * _rvec3(newODir)).normalize()
    oldNormal = (scitbx.matrix.cross_product_matrix(_lvec3(oldNDir)) * _rvec3(oldODir)).normalize()
    dot = (_lvec3(newNormal) * _rvec3(oldNormal))[0]
    if dot > -0.99999:
      return "Movers.Test() MoverNH2Flip basic: Bad plane alignment: "+str(dot)

    dCarbon = (fixed[4] - _lvec3(p.xyz)).length()
    if dCarbon < 0.001 or dCarbon > 0.1:
      return "Movers.Test() MoverNH2Flip basic: Bad hinge motion: "+str(dCarbon)

    dCarbon = (fixed[5] - _lvec3(f.xyz)).length()
    if dCarbon < 0.0005 or dCarbon > 0.1:
      return "Movers.Test() MoverNH2Flip basic: Bad pivot motion: "+str(dCarbon)

    dHydrogen = (fixed[6] - _lvec3(fh1.xyz)).length()
    if dHydrogen < 0.0005 or dHydrogen > 0.1:
      return "Movers.Test() MoverNH2Flip basic: Bad pivot hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixed[7] - _lvec3(fh2.xyz)).length()
    if dHydrogen < 0.0005 or dHydrogen > 0.1:
      return "Movers.Test() MoverNH2Flip basic: Bad pivot hydrogen motion: "+str(dHydrogen)

  except Exception as e:
    return "Movers.Test() MoverNH2Flip basic: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  # Test the MoverNH2Flip class with a linker (similar to Gln).
  try:
    # Test behavior with offsets for each atom so that we test the generic case.
    # Construct a MoverNH2Flip that has the N, H's and Oxygen located (non-physically)
    # slightly in the +Y direction out of the X-Z plane, with the Hydrogens 120 around the
    # same offset axis.
    # They are bonded to a Carbon and friend (pivot) and linker and alpha carbon that are
    # on the Z axis.
    # Non-physically, atoms are 1 unit apart.
    p = pdb.hierarchy.atom()
    p.xyz = [ 0.0, 0.0, 0.0 ]

    f = pdb.hierarchy.atom()
    f.element = "C"
    f.xyz = [ 0.0, 0.0,-1.0 ]

    fh1 = pdb.hierarchy.atom()
    fh1.element = "H"
    fh1.xyz = [ 1.0, 0.0,-1.0 ]

    fh2 = pdb.hierarchy.atom()
    fh2.element = "H"
    fh2.xyz = [-1.0, 0.0,-1.0 ]

    ln = pdb.hierarchy.atom()
    ln.element = "C"
    ln.xyz = [ 0.0, 0.0,-2.0 ]

    lnh1 = pdb.hierarchy.atom()
    lnh1.element = "H"
    lnh1.xyz = [ 1.0, 0.0,-2.0 ]

    lnh2 = pdb.hierarchy.atom()
    lnh2.element = "H"
    lnh2.xyz = [-1.0, 0.0,-2.0 ]

    ca = pdb.hierarchy.atom()
    ca.element = "C"
    ca.name = "CA"
    ca.xyz = [ 0.0, 0.0,-3.0 ]

    # Nitrogen and Oxygen are +/-120 degrees from carbon-carbon bond
    axis = flex.vec3_double([ [0,0,0], [0,1,0] ])
    n = pdb.hierarchy.atom()
    n.element = "N"
    n.xyz = _rotateAroundAxis(f, axis,-120) + _lvec3([0,0.01,0]) + _lvec3([ 0.002, 0.003,-0.004])

    o = pdb.hierarchy.atom()
    o.element = "O"
    o.xyz = _rotateAroundAxis(f, axis, 120) + _lvec3([0,0.01,0]) + _lvec3([-0.003, 0.002, 0.003])

    # Hydrogens are +/-120 degrees from nitrogen-carbon bond
    axis = flex.vec3_double([ n.xyz, [0,1,0] ])
    h1 = pdb.hierarchy.atom()
    h1.element = "H"
    h1.xyz = _rotateAroundAxis(p, axis,-120) + _lvec3([0,0.01,0]) + _lvec3([-0.008, 0.001, 0.008])

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.xyz = _rotateAroundAxis(p, axis, 120) + _lvec3([0,0.01,0]) + _lvec3([ 0.007,-0.001, 0.007])

    # Build the hierarchy so we can reset the i_seq values.
    ag = pdb.hierarchy.atom_group()
    ag.append_atom(h1)
    ag.append_atom(h2)
    ag.append_atom(n)
    ag.append_atom(o)
    ag.append_atom(p)
    ag.append_atom(f)
    ag.append_atom(fh1)
    ag.append_atom(fh2)
    ag.append_atom(ln)
    ag.append_atom(lnh1)
    ag.append_atom(lnh2)
    ag.append_atom(ca)
    rg = pdb.hierarchy.residue_group()
    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    bondedNeighborLists = {}
    bondedNeighborLists[h1] = [ n ]
    bondedNeighborLists[h2] = [ n ]
    bondedNeighborLists[n] = [ h1, h2, p ]
    bondedNeighborLists[o] = [ p ]
    bondedNeighborLists[p] = [ n, o, f ]
    bondedNeighborLists[f] = [ p, ln, fh1, fh2 ]
    bondedNeighborLists[fh1] = [ f ]
    bondedNeighborLists[fh2] = [ f ]
    bondedNeighborLists[ln] = [ f, ca, lnh1, lnh2 ]
    bondedNeighborLists[lnh1] = [ ln ]
    bondedNeighborLists[lnh2] = [ ln ]
    bondedNeighborLists[ca] = [ ln ]

    mover = MoverNH2Flip(n, ca.name, bondedNeighborLists)
    fixed = mover.FixUp(1).newPositions

    # Ensure that the results meet the specifications:
    # 1) New Oxygen on the line from the alpha carbon to the old Nitrogen
    # 2) New plane of Oxygen, Nitrogen, Alpha Carbon matches old plane, but flipped
    # 3) Pivot and linker Carbons and Hydrogens move slightly due to rigid-body motion

    newODir = (fixed[3] - _lvec3(ca.xyz)).normalize()
    oldNDir = (_rvec3(n.xyz) - _rvec3(ca.xyz)).normalize()
    if (newODir * oldNDir)[0] < 0.9999:
      return "Movers.Test() MoverNH2Flip linked: Bad oxygen alignment: "+str((newODir * oldNDir)[0])

    newNDir = (fixed[2] - _lvec3(ca.xyz)).normalize()
    oldODir = (_rvec3(o.xyz) - _rvec3(ca.xyz)).normalize()
    newNormal = (scitbx.matrix.cross_product_matrix(_lvec3(newNDir)) * _rvec3(newODir)).normalize()
    oldNormal = (scitbx.matrix.cross_product_matrix(_lvec3(oldNDir)) * _rvec3(oldODir)).normalize()
    dot = (_lvec3(newNormal) * _rvec3(oldNormal))[0]
    if dot > -0.99999:
      return "Movers.Test() MoverNH2Flip linked: Bad plane alignment: "+str(dot)

    dCarbon = (fixed[4] - _lvec3(p.xyz)).length()
    if dCarbon < 0.0006 or dCarbon > 0.1:
      return "Movers.Test() MoverNH2Flip linked: Bad hinge motion: "+str(dCarbon)

    dCarbon = (fixed[5] - _lvec3(f.xyz)).length()
    if dCarbon < 0.0004 or dCarbon > 0.1:
      return "Movers.Test() MoverNH2Flip linked: Bad pivot motion: "+str(dCarbon)

    dCarbon = (fixed[6] - _lvec3(ln.xyz)).length()
    if dCarbon < 0.0002 or dCarbon > 0.1:
      return "Movers.Test() MoverNH2Flip linked: Bad linker motion: "+str(dCarbon)

    # Hydrogens come after all linkers
    dHydrogen = (fixed[7] - _lvec3(fh1.xyz)).length()
    if dHydrogen < 0.0004 or dHydrogen > 0.1:
      return "Movers.Test() MoverNH2Flip linked: Bad pivot hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixed[8] - _lvec3(fh2.xyz)).length()
    if dHydrogen < 0.0004 or dHydrogen > 0.1:
      return "Movers.Test() MoverNH2Flip linked: Bad pivot hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixed[9] - _lvec3(lnh1.xyz)).length()
    if dHydrogen < 0.0002 or dHydrogen > 0.1:
      return "Movers.Test() MoverNH2Flip linked: Bad linker hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixed[10] - _lvec3(lnh2.xyz)).length()
    if dHydrogen < 0.0002 or dHydrogen > 0.1:
      return "Movers.Test() MoverNH2Flip linked: Bad linker hydrogen motion: "+str(dHydrogen)

    # Ensure that the Hydrogens moved along with their parent in the flip.
    dHydrogen = (fixed[0] - fixed[2]).length()
    oldDHydrogen = (_lvec3(h1.xyz)-_lvec3(n.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverNH2Flip linked: Bad nitrogen-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (fixed[1] - fixed[2]).length()
    oldDHydrogen = (_lvec3(h2.xyz)-_lvec3(n.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverNH2Flip linked: Bad nitrogen-hydrogen motion: "+str(dHydrogen-oldDHydrogen)


  except Exception as e:
    return "Movers.Test() MoverNH2Flip linked: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  # Test the MoverHistidineFlip class.
  try:
    # Test behavior with offsets for each atom so that we test the generic case.
    # Construct a MoverHistidineFlip that has the ring atoms other than CG located (non-physically)
    # slightly in the +Y direction out of the X-Z plane, with ring atoms and Hydrogens
    # not actually in a ring but with the same topology.
    # They are bonded to a partner (hinge/CG) and friend (pivot/CB) and alpha carbon that are on the Z axis.
    # Non-physically, atoms are 1 unit apart.
    p = pdb.hierarchy.atom()
    p.element = "C"
    p.xyz = [ 0.0, 0.0, 0.0 ]

    f = pdb.hierarchy.atom()
    f.element = "C"
    f.xyz = [ 0.0, 0.0,-1.0 ]

    fh1 = pdb.hierarchy.atom()
    fh1.element = "H"
    fh1.xyz = [ 1.0, 0.0,-1.0 ]

    fh2 = pdb.hierarchy.atom()
    fh2.element = "H"
    fh2.xyz = [-1.0, 0.0,-1.0 ]

    ca = pdb.hierarchy.atom()
    ca.element = "C"
    ca.xyz = [ 0.0, 0.0,-2.0 ]

    # Put the four ring atoms and their Hydrogens in a lattice with H stickout out in X.
    nd1 = pdb.hierarchy.atom()
    nd1.element = "N"
    nd1.xyz = _lvec3([ 1.0, 0.0, 1.0]) + _lvec3([0,0.01,0]) + _lvec3([-0.003, 0.002, 0.003])

    nd1h = pdb.hierarchy.atom()
    nd1h.element = "H"
    nd1h.xyz = [ 2.0, 0.01, 1.0 ]

    ce1 = pdb.hierarchy.atom()
    ce1.element = "C"
    ce1.xyz = _lvec3([ 1.0, 0.0, 2.0]) + _lvec3([0,0.01,0]) + _lvec3([-0.008, 0.001, 0.008])

    ce1h = pdb.hierarchy.atom()
    ce1h.element = "H"
    ce1h.xyz = [ 2.0, 0.01, 2.0 ]

    cd2 = pdb.hierarchy.atom()
    cd2.element = "C"
    cd2.xyz = _lvec3([-1.0, 0.0, 1.0]) + _lvec3([0,0.01,0]) + _lvec3([ 0.007,-0.001, 0.007])

    cd2h = pdb.hierarchy.atom()
    cd2h.element = "H"
    cd2h.xyz = [-2.0, 0.01, 1.0 ]

    ne2 = pdb.hierarchy.atom()
    ne2.element = "N"
    ne2.xyz = _lvec3([-1.0, 0.0, 2.0]) + _lvec3([0,0.01,0]) + _lvec3([-0.002, 0.005, 0.004])

    ne2h = pdb.hierarchy.atom()
    ne2h.element = "H"
    ne2h.xyz = [-2.0, 0.01, 2.0 ]

    # Build the hierarchy so we can reset the i_seq values.
    ag = pdb.hierarchy.atom_group()
    ag.append_atom(nd1)
    ag.append_atom(nd1h)
    ag.append_atom(ce1)
    ag.append_atom(ce1h)
    ag.append_atom(cd2)
    ag.append_atom(cd2h)
    ag.append_atom(ne2)
    ag.append_atom(ne2h)
    ag.append_atom(p)
    ag.append_atom(f)
    ag.append_atom(fh1)
    ag.append_atom(fh2)
    ag.append_atom(ca)
    rg = pdb.hierarchy.residue_group()
    rg.append_atom_group(ag)
    c = pdb.hierarchy.chain()
    c.append_residue_group(rg)
    m = pdb.hierarchy.model()
    m.append_chain(c)
    m.atoms().reset_i_seq()

    bondedNeighborLists = {}
    bondedNeighborLists[nd1] = [ nd1h, p, ce1 ]
    bondedNeighborLists[nd1h] = [ nd1 ]
    bondedNeighborLists[ce1] = [ ce1h, nd1, ne2 ]
    bondedNeighborLists[ce1h] = [ ce1 ]
    bondedNeighborLists[cd2] = [ cd2h, p, ne2 ]
    bondedNeighborLists[cd2h] = [ cd2 ]
    bondedNeighborLists[ne2] = [ ne2h, cd2, ce1 ]
    bondedNeighborLists[ne2h] = [ ne2 ]
    bondedNeighborLists[p] = [ nd1, cd2, f ]
    bondedNeighborLists[f] = [ p, ca, fh1, fh2 ]
    bondedNeighborLists[fh1] = [ f ]
    bondedNeighborLists[fh2] = [ f ]
    bondedNeighborLists[ca] = [ f ]

    mover = MoverHistidineFlip(ne2, bondedNeighborLists)
    fixed = mover.FixUp(1).newPositions

    # Ensure that the coarse-flip results meet the expections:
    # 1) N and C atoms are flipped in pairs
    # 2) H remain at the same distance from the new locations.

    coarse = mover.CoarsePositions()
    if len(coarse.positions) != 2:
      return "Movers.Test() MoverHistidineFlip: Did not find two locations: "+str(len(coarse.positions))
    newPos = coarse.positions[1]
    dist = (_lvec3(newPos[0]) - _lvec3(ce1.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverHistidineFlip: NE2 moved incorrectly: "+str(dist)
    dist = (_lvec3(newPos[2]) - _lvec3(ne2.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverHistidineFlip: CE1 moved incorrectly: "+str(dist)
    dist = (_lvec3(newPos[4]) - _lvec3(cd2.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverHistidineFlip: ND1 moved incorrectly: "+str(dist)
    dist = (_lvec3(newPos[6]) - _lvec3(nd1.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverHistidineFlip: CD2 moved incorrectly: "+str(dist)

    dHydrogen = (_lvec3(newPos[0]) - _lvec3(newPos[1])).length()
    oldDHydrogen = (_lvec3(ne2h.xyz)-_lvec3(ne2.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHistidineFlip: Bad coarse NE2 hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (_lvec3(newPos[2]) - _lvec3(newPos[3])).length()
    oldDHydrogen = (_lvec3(ce1h.xyz)-_lvec3(ce1.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHistidineFlip: Bad coarse CE1 hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (_lvec3(newPos[4]) - _lvec3(newPos[5])).length()
    oldDHydrogen = (_lvec3(nd1h.xyz)-_lvec3(nd1.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHistidineFlip: Bad coarse ND1 hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (_lvec3(newPos[6]) - _lvec3(newPos[7])).length()
    oldDHydrogen = (_lvec3(cd2h.xyz)-_lvec3(cd2.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHistidineFlip: Bad coarse ND1 hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    # Ensure that the FixUp results meet the specifications:
    # 1) New CE1 on the line from the alpha carbon to the old NE2
    # 2) New plane of CE1, NE2, Alpha Carbon matches old plane, but flipped
    # 3) Carbons and pivot Hydrogens move slightly due to rigid-body motion

    newCE1Dir = (fixed[2] - _lvec3(ca.xyz)).normalize()
    oldNE2Dir = (_rvec3(ne2.xyz) - _rvec3(ca.xyz)).normalize()
    if (newCE1Dir * oldNE2Dir)[0] < 0.9999:
      return "Movers.Test() MoverHistidineFlip: Bad CE1 alignment: "+str((newCE1Dir * oldNE2Dir)[0])

    newNE2Dir = (fixed[0] - _lvec3(ca.xyz)).normalize()
    oldCE1Dir = (_rvec3(ce1.xyz) - _rvec3(ca.xyz)).normalize()
    newNormal = (scitbx.matrix.cross_product_matrix(_lvec3(newNE2Dir)) * _rvec3(newCE1Dir)).normalize()
    oldNormal = (scitbx.matrix.cross_product_matrix(_lvec3(oldNE2Dir)) * _rvec3(oldCE1Dir)).normalize()
    dot = (_lvec3(newNormal) * _rvec3(oldNormal))[0]
    if dot > -0.99999:
      return "Movers.Test() MoverHistidineFlip: Bad plane alignment: "+str(dot)

    dCarbon = (fixed[8] - _lvec3(p.xyz)).length()
    if dCarbon < 0.001 or dCarbon > 0.1:
      return "Movers.Test() MoverHistidineFlip: Bad hinge motion: "+str(dCarbon)

    dCarbon = (fixed[9] - _lvec3(f.xyz)).length()
    if dCarbon < 0.0005 or dCarbon > 0.1:
      return "Movers.Test() MoverHistidineFlip: Bad pivot motion: "+str(dCarbon)

    dHydrogen = (fixed[10] - _lvec3(fh1.xyz)).length()
    if dHydrogen < 0.0005 or dHydrogen > 0.1:
      return "Movers.Test() MoverHistidineFlip: Bad pivot hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixed[11] - _lvec3(fh2.xyz)).length()
    if dHydrogen < 0.0005 or dHydrogen > 0.1:
      return "Movers.Test() MoverHistidineFlip: Bad pivot hydrogen motion: "+str(dHydrogen)

    # Ensure that the Hydrogens moved along with their parent in the flip.
    dHydrogen = (fixed[0] - fixed[1]).length()
    oldDHydrogen = (_lvec3(ne2h.xyz)-_lvec3(ne2.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHistidineFlip: Bad NE2-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (fixed[2] - fixed[3]).length()
    oldDHydrogen = (_lvec3(ce1h.xyz)-_lvec3(ce1.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHistidineFlip: Bad CE1-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (fixed[4] - fixed[5]).length()
    oldDHydrogen = (_lvec3(nd1h.xyz)-_lvec3(nd1.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHistidineFlip: Bad ND1-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (fixed[6] - fixed[7]).length()
    oldDHydrogen = (_lvec3(cd2h.xyz)-_lvec3(cd2.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHistidineFlip: Bad CD2-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

  except Exception as e:
    return "Movers.Test() MoverHistidineFlip: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  return ""

##################################################################################
# Internal helper functions to make things that are compatible with vec3_double so
# that we can do math on them.  We need a left-hand and right-hand one so that
# we can make both versions for multiplication.
def _rvec3 (xyz) :
  return scitbx.matrix.rec(xyz, (3,1))
def _lvec3 (xyz) :
  return scitbx.matrix.rec(xyz, (1,3))

##################################################################################
# Internal helper functions for angle manipulation.
def _rotateOppositeFriend(atom, axis, partner, friends):
  '''Rotate the atom to point away from one of the friends.  This means placing it
     within the plane perpendicular to the axis of rotation through the point on that axis
     that it is closest to at the same distance it starts out from that axis and in a
     direction that is the negation of the vector from the partner to the friend projected
     into the plane perpendicular to the axis and then normalized.
     :param atom: iotbx.pdb.hierarchy.atom to be moved.
     :param partner: iotbx.pdb.hierarchy.atom atom that is bonded to the neighbor of
     the atom.
     :param friends: flex iotbx.pdb.hierarchy.atom list of atoms bonded to the partner
     besides the neighbor.
     :param axis: flex array with two scitbx::vec3<double> points, the first
     of which is the origin and the second is a vector pointing in the direction
     of the axis of rotation.  Positive rotations will be right-handed rotations
     around this axis.
     :returns the new location for the atom.
  '''
  normal = _rvec3(axis[1])
  friend = friends[0]
  friendFromPartner = _lvec3(friend.xyz) - _lvec3(partner.xyz)
  alongAxisComponent = (friendFromPartner * normal)
  friendFromPartner = _rvec3(friend.xyz) - _rvec3(partner.xyz)
  inPlane = friendFromPartner - normal*alongAxisComponent
  normalizedOffset = -inPlane.normalize()

  d = -_lvec3(atom.xyz)*normal
  t = - (d + (_lvec3(axis[0])) * normal)
  nearPoint = _rvec3(axis[0]) + normal*t
  distFromNearPoint = (_rvec3(atom.xyz)-nearPoint).length()

  return nearPoint + distFromNearPoint * normalizedOffset

def _rotateAroundAxis(atom, axis, degrees):
  '''Rotate the atom about the specified axis by the specified number of degrees.
     :param atom: iotbx.pdb.hierarchy.atom or scitbx::vec3<double> or
     scitbx.matrix.rec(xyz, (3,1)) to be moved.
     :param axis: flex array with two scitbx::vec3<double> points, the first
     of which is the origin and the second is a vector pointing in the direction
     of the axis of rotation.  Positive rotations will be right-handed rotations
     around this axis.
     :param degrees: How much to rotate the atom around the axis.
     Positive rotation is right-handed around the axis.
     :returns the new location for the atom.
  '''
  # If we have an atom, get its position.  Otherwise, we have a position passed
  # in.
  try:
    pos = atom.xyz
  except:
    pos = _rvec3(atom)

  # Project the atom position onto the axis, finding its closest point on the axis.
  # The position lies on a plane whose normal points along the axis vector.  The
  # point we seek is the intersection of the axis with this plane.
  # The plane equation will be the normalized axis direction vector and the offset
  # from the origin such that N.point + d = 0 defines the plane.
  # Solve for the time at which the line's ray crosses the plane and then solve for
  # the location along the line at that time.  t = - (d + (lineOrigin * planeNormal)) /
  # (lineDirection * planeNormal).  Because the line direction and normal are the
  # same, the divisor is 1.
  normal = _lvec3(axis[1]).normalize()
  d = -normal*pos
  t = - (d + (normal * _rvec3(axis[0])))
  nearPoint = _lvec3(axis[0]) + t * normal

  # Find the vector from the closest point towards the atom, which is its offset
  offset = _lvec3(pos) - nearPoint

  # Rotate the offset vector around the axis by the specified angle.  Add the new
  # offset to the closest point. Store this as the new location for this atom and angle.
  newOffset = offset.rotate_around_origin(_lvec3(axis[1]), degrees*math.pi/180)
  return nearPoint + newOffset

def _rotateHingeDock(movableAtoms, hingeIndex, firstDockIndex, secondDockIndex, alphaCarbon):
  '''Perform the three-step rotate-hinge-dock calculation described in 
     Protein Science Vol 27:293-315 and implemented in
     FlipMemo::RotHingeDock_Flip() in the original Reduce C++ code.
     :param movableAtoms: flex array of iotbx.pdb.hierarchy.atom objects that
     are the original locations of atoms that can be moved, including the Carbon
     atom that is bonded to the alpha carbon but not the alpha carbon itself.
     The order of atoms is as follows: The first ones are in the set that is rotated
     around the vector from the hinge atom to the pivot atom.  Then comes the hinge
     atom, then the pivot atom.  Then all of the atoms between the pivot and the
     Alpha Carbon (if any) including any Hydrogens bound to the pivot atom or the
     linking atoms, excluding the Alpha Carbon itself.
     :param hingeIndex: Index in the movableAtoms array that tells which is
     the atom to hinge around (this will be a Carbon atom).
     :param firstDockIndex: Index in the movableAtoms array that tells which is
     the atom to dock first (this will be the Oxygen atom for Asn and Gln,
     CE1 for His).
     :param secondDockIndex: Index in the movableAtoms array that tells which is
     the atom to dock second (this will be the Nitrogen atom for Asn and Gln,
     NE2 for His).
     :param alphaCarbon: iotbx.pdb.hierarchy.atom that is the alpha carbon.
     The original rotation is around the bond from this atom to the original
     carbon atom.
     :returns A flex array of scitbx.matrix.rec(xyz, (1,3)) entries giving the new locations.
  '''
  #########################
  # A) Rotate the movable atoms that are listed before the hinge atom around the
  #    hinge-pivot vector by 180 degrees.
  # B) Hinge the movable atoms listed before the hings around the hinge to place them back into the
  #    plane that they the two docked atoms originally located in using the axis of least
  #    rotation that passed through the hinge.
  # C) Rotate all of the movable atoms around the alphaCarbon, aligning the docked
  #    atoms with their original locations as much as possible.  This is
  #    done in two steps: 1) Rotating around the shortest axis to make the
  #    first docked atom lie on the original vector from the alphaCarbon to the original second
  #    docked atom.
  #    2) Rotating around the alphaCarbon-to-new-firstDockIndex vector to align the plane
  #    containing the alphaCarbon, new and old secondDockIndex with the plane
  #    containing the alphaCarbon, the original firstDockIndex and the original secondDockIndex
  #    (making the dihedral angle new secondDockIndex, alphaCarbon, old secondDockIndex, old
  #    firstDockIndex be 0).

  # Find the special atoms
  first = movableAtoms[firstDockIndex]
  second = movableAtoms[secondDockIndex]
  hingeAtom = movableAtoms[hingeIndex]
  pivotAtom = movableAtoms[hingeIndex+1]

  # Construct a list of output locations, initializing with the input locations.
  movable = []
  for a in movableAtoms:
    movable.append(a.xyz)

  # A) Rotate the movable atoms that come before the hinge atom by 180 degrees
  # around the pivot-to-hinge vector
  normal = (_rvec3(pivotAtom.xyz) - _rvec3(hingeAtom.xyz)).normalize()
  axis = flex.vec3_double([hingeAtom.xyz, normal])
  for i in range(hingeIndex):
    movable[i] = _rotateAroundAxis(movable[i], axis, 180)

  # B) Hinge the movable atoms around the hinge.
  # Solve for the planes containing the old and new atoms and then the vector
  # that these planes intersect at (cross product).

  cToO = _lvec3(first.xyz) - _lvec3(hingeAtom.xyz)
  nToO = _rvec3(second.xyz) - _rvec3(hingeAtom.xyz)
  oldNormal = (scitbx.matrix.cross_product_matrix(cToO) * nToO).normalize()

  # Flip the order here because we've rotated 180 degrees
  newNToO = _lvec3(movable[secondDockIndex]) - _lvec3(movable[hingeIndex])
  newCToO = _rvec3(movable[firstDockIndex]) - _rvec3(movable[hingeIndex])
  newNormal = (scitbx.matrix.cross_product_matrix(newNToO) * newCToO).normalize()

  # If we don't need to rotate, we'll get a zero-length vector
  hinge = scitbx.matrix.cross_product_matrix(_lvec3(oldNormal))*_rvec3(newNormal)
  if hinge.length() > 0:
    hinge = hinge.normalize()
    axis = flex.vec3_double([hingeAtom.xyz, hinge])
    degrees = 180/math.pi * math.acos((_lvec3(oldNormal)*_rvec3(newNormal))[0])
    for i in range(hingeIndex):
      # Rotate in the opposite direction, taking the new back to the old.
      movable[i] = _rotateAroundAxis(_rvec3(movable[i]), axis, -degrees)

  # C) Rotate the atoms around the alphaCarbon
  #  1) firstDockIndex to alphaCarbon<-->original secondDockIndex line
  aToNewO = _lvec3(movable[firstDockIndex]) - _lvec3(alphaCarbon.xyz)
  aToOldN = _rvec3(second.xyz) - _rvec3(alphaCarbon.xyz)
  # If we don't need to rotate, we'll get a zero-length vector
  normal = scitbx.matrix.cross_product_matrix(aToNewO) * aToOldN
  if normal.length() > 0:
    normal = normal.normalize()
    axis = flex.vec3_double([alphaCarbon.xyz, normal])
    degrees = 180/math.pi * math.acos((_lvec3(aToOldN.normalize())*_rvec3(aToNewO.normalize()))[0])
    for i in range(len(movable)):
      # Rotate in the opposite direction, taking the new back to the old.
      movable[i] = _rotateAroundAxis(_rvec3(movable[i]), axis, -degrees)

  #  2) firstDockIndex to the proper plane.
  sites = flex.vec3_double([ movable[secondDockIndex], alphaCarbon.xyz, second.xyz, first.xyz ])
  degrees = scitbx.math.dihedral_angle(sites=sites, deg=True)
  hinge = _rvec3(alphaCarbon.xyz) - _rvec3(second.xyz)
  if hinge.length() > 0:
    hinge = hinge.normalize()
    axis = flex.vec3_double([alphaCarbon.xyz, hinge])
    for i in range(len(movable)):
      # Rotate in the opposite direction, taking the new back to the old.
      movable[i] = _rotateAroundAxis(_rvec3(movable[i]), axis, -degrees)

  return movable

##################################################################################
# External helper functions.

def getBondedNeighborLists(atoms, bondProxies):
  """Produce a dictionary with one entry for each atom that contains a list of all of
    the atoms that are bonded to it.
    External helper function to produce a dictionary of lists that lists all bonded
    neighbors for each atom in a set of atoms.
    :param atoms: Flex array of atoms (could be obtained using model.get_atoms() if there
    are no chains with multiple conformations, must be a subset of the atoms including
    all in the base conformation and in a particular conformation otherwise).
    :param bondProxies: Flex array of bond proxies for the atoms.  This could be obtained
    using model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart =
    model.get_sites_cart())[0] if the model has only a single conformation.  Otherwise,
    it should be a flex array of atom positions for the atoms that are in the first argument.
    :returns a dictionary with one entry for each atom that contains a list of all of
    the atoms that are bonded to it.
  """
  atomDict = {}
  for a in atoms:
    atomDict[a.i_seq] = a
  bondedNeighbors = {}
  for a in atoms:
    bondedNeighbors[a] = []
  for bp in bondProxies:
    bondedNeighbors[atomDict[bp.i_seqs[0]]].append(atomDict[bp.i_seqs[1]])
    bondedNeighbors[atomDict[bp.i_seqs[1]]].append(atomDict[bp.i_seqs[0]])
  return bondedNeighbors

##################################################################################
# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  ret = Test()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
