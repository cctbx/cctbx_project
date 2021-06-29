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
import scitbx.matrix
from scitbx.array_family import flex
from iotbx import pdb
import traceback

##################################################################################
def getBondedNeighborLists(atoms, bondProxies):
  """Produce a dictionary with one entry for each atom that contains a list of all of
    the atoms that are bonded to it.
    External helper function to produce a dictionary of lists that lists all bonded
    neighbors for each atom in a set of atoms.
    :param atoms: Flex aray of atoms (could be obtained using model.get_atoms()).
    :param bondProxies: Flex array of bond proxies for the atoms (could be obtained
    using model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart =
    model.get_sites_cart()) ).
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
#  - type PositionReturn: ( flex<atom> atoms, flex<flex<vec3>> positions, flex<float> preferenceEnergies,
#                           reduceOptions = None)
#     The reduceOptions is a Phil option subset.  The relevant options for each Mover
#       type below is specified in its __init__ method.
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
  def __init__(self, atom, reduceOptions = None):
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
  def __init__(self, atoms, axis, coarseRange, coarseStepDegrees = None,
                doFineRotations = True, preferenceFunction = None,
                reduceOptions = None):
    """Constructs a Rotator, to be called by a derived class or test program but
       not usually user code.
       A base class for types of Movers that rotate one or more atoms around a single
       axis.  The derived class must determine the initial position for each atom, and
       a preferred-energies function that maps from an angle that is 0 at the initial
       location to an energy.  The range of coarse rotations is also specified, along
       with a flag telling whether fine rotations are allowed.  The coarse step size
       can also be specified, overriding the value from reduceOptions.CoarseStepDegrees.
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
       :param coarseStepDegrees: With the default value of None, the number of
       degrees per coarse step will be determined  by the
       reduceOptions.CoarseStepDegrees or the default value of 30 if neither specifies
       the value.
       To test only two fixed orientations, doFineRotations should be set
       to False, coarseRange to 180, and coarseStepDegrees to 180.
       :param doFineRotations: Specifies whether fine rotations are allowed for
       this instance.  To test only two fixed orientations, this should be set
       to False, coarseRange to 180, and coarseStepDegrees to 180.
       :param preferenceFunction: A function that takes a floating-point
       argument in degrees of rotation about the specified axis and produces
       a floating-point energy value that is multiplied by the value of
       reduceOptions.PreferredOrientationScale and then added to the energy
       returned for this orientation by CoarsePositions() and FinePositions().
       For the default value of None, no addition is performed.
       :param reduceOptions: 
        The reduceOptions is a Phil option subset.  The relevant options for _MoverRotator
          are: CoarseStepDegrees, FineStepDegrees, PreferredOrientationScale.
    """
    self._atoms = atoms
    self._axis = axis
    self._coarseRange = coarseRange
    self._doFineRotations = doFineRotations
    self._preferenceFunction = preferenceFunction
    self._reduceOptions = reduceOptions

    # Determine the coarse step size in degrees.  If it was specified in the
    # constructor, use that.  Otherwise, if it was specified in
    # reduceOptions.CoarseStepDegrees use that.  Otherwise, use a default of
    # 30 degrees.
    self._coarseStepDegrees = 30.0
    try:
      if reduceOptions.CoarseStepDegrees is not None:
        self._coarseStepDegrees = reduceOptions.CoarseStepDegrees
    except:
      pass
    if coarseStepDegrees is not None:
      self._coarseStepDegrees = coarseStepDegrees

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

    # Determine the fine step size from reduceOptions.FineStepDegrees if it
    # exists; otherwise, default to 1 degrees.
    self._fineStepDegrees = 1.0
    try:
      if reduceOptions.FineStepDegrees is not None:
        self._fineStepDegrees = reduceOptions.FineStepDegrees
    except:
      pass

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

    # Determine the preference scale.  The default of 1 can be overridden by
    # a value in reduceOptions.PreferredOrientationScale.
    self._preferredOrientationScale = 1.0
    try:
      if reduceOptions.PreferredOrientationScale is not None:
        self._preferredOrientationScale = reduceOptions.PreferredOrientationScale
    except:
      pass

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
  def __init__(self, atom, bondedNeighborLists, coarseStepDegrees = None, reduceOptions = None):
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
       :param coarseStepDegrees: With the default value of None, the number of
       degrees per coarse step will be determined  by the
       reduceOptions.CoarseStepDegrees or the default value of 30 if neither specifies
       the value.
       :param reduceOptions: 
        The reduceOptions is a Phil option subset.  The relevant options for
          MoverSingleHydrogenRotator are: CoarseStepDegrees, FineStepDegrees, PreferredOrientationScale.
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
    _MoverRotator.__init__(self, atoms, axis, 180, coarseStepDegrees,
      preferenceFunction = preferenceFunction, reduceOptions = reduceOptions)

##################################################################################
class MoverNH3Rotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, coarseStepDegrees = None, reduceOptions = None):
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
       :param coarseStepDegrees: With the default value of None, the number of
       degrees per coarse step will be determined  by the
       reduceOptions.CoarseStepDegrees or the default value of 30 if neither specifies
       the value.
       :param reduceOptions: 
        The reduceOptions is a Phil option subset.  The relevant options for
          MoverNH3Rotator are: CoarseStepDegrees, FineStepDegrees, PreferredOrientationScale.
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
      preferenceFunction = preferenceFunction, reduceOptions = reduceOptions)

##################################################################################
class MoverAromaticMethylRotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, reduceOptions = None):
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
       :param reduceOptions: 
        The reduceOptions is a Phil option subset.  The relevant options for
          MoverAromaticMethylRotator are: None.
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
    _MoverRotator.__init__(self, hydrogens, axis, 180, 180, doFineRotations = False,
      reduceOptions = reduceOptions)

##################################################################################
class MoverTetrahedralMethylRotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, coarseStepDegrees = None, reduceOptions = None):
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
       :param coarseStepDegrees: With the default value of None, the number of
       degrees per coarse step will be determined  by the
       :param reduceOptions: 
        The reduceOptions is a Phil option subset.  The relevant options for
          MoverTetrahedralMethylRotator are: CoarseStepDegrees, FineStepDegrees, PreferredOrientationScale.
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
    _MoverRotator.__init__(self, hydrogens, axis, 180, preferenceFunction = preferenceFunction,
      reduceOptions = reduceOptions)

##################################################################################
# @todo Define each type of Mover

##################################################################################
# Test function to verify that all Movers behave properly.

def Test():
  """Test function for all classes provided above.
  :returns Empty string on success, string describing the problem on failure.
  :returns Empty string on success, string describing the problem on failure.
  """

  # Construct a Class that will behave like the Reduce Phil data structure so that
  # we can specify the probe radius.
  class FakePhil:
    pass
  fakePhil = FakePhil()

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
    rot = _MoverRotator(atoms,axis, 180, 90.0, True, prefFunc)

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

    # Test None preference scale factor and a different nonzero value.
    fakePhil.PreferredOrientationScale = None
    rot = _MoverRotator(atoms,axis, 180, 90.0, True, prefFunc, reduceOptions = fakePhil)
    coarse = rot.CoarsePositions()
    for p in coarse.preferenceEnergies:
      if p != 1:
        return "Movers.Test() _MoverRotator None preference: Expected preference energy = 1, got "+str(p)
    fakePhil.PreferredOrientationScale = 2
    rot = _MoverRotator(atoms,axis, 180, 90.0, True, prefFunc, reduceOptions = fakePhil)
    coarse = rot.CoarsePositions()
    for p in coarse.preferenceEnergies:
      if p != 2:
        return "Movers.Test() _MoverRotator Phil-scaled preference: Expected preference energy = 2, got "+str(p)

    # Test None preference function and a sinusoidal one.
    rot = _MoverRotator(atoms,axis, 180, 90.0, True, None)
    coarse = rot.CoarsePositions()
    for p in coarse.preferenceEnergies:
      if p != 0:
        return "Movers.Test() _MoverRotator None preference function: Expected preference energy = 0, got "+str(p)
    def prefFunc2(x):
      return math.cos(x * math.pi / 180)
    rot = _MoverRotator(atoms,axis, 180, 180, True, prefFunc2)
    coarse = rot.CoarsePositions()
    expected = [1, -1]
    for i,p  in enumerate(coarse.preferenceEnergies):
      val = expected[i]
      if p != val:
        return "Movers.Test() _MoverRotator Sinusoidal preference function: Expected preference energy = "+str(val)+", got "+str(p)

    # Test coarseStepDegrees default behavior and setting via reduceOptions.
    rot = _MoverRotator(atoms,axis, 180)
    coarse = rot.CoarsePositions()
    if len(coarse.positions) != 12:
      return "Movers.Test() _MoverRotator Default coarse step: Expected 12, got "+str(len(coarse.positions))
    fakePhil.CoarseStepDegrees = 15
    rot = _MoverRotator(atoms,axis, 180, reduceOptions = fakePhil)
    coarse = rot.CoarsePositions()
    if len(coarse.positions) != 24:
      return "Movers.Test() _MoverRotator reduceOptions coarse step: Expected 24, got "+str(len(coarse.positions))

    # Test fineStepDegrees setting via reduceOptions.
    fakePhil.CoarseStepDegrees = None
    fakePhil.FineStepDegrees = 1
    rot = _MoverRotator(atoms,axis, 180, reduceOptions = fakePhil)
    fine = rot.FinePositions(0)
    # +/- 15 degrees in 1-degree steps, but we don't do the +15 because it will be handled by the next
    # rotation up.
    if len(fine.positions) != 29:
      return "Movers.Test() _MoverRotator reduceOptions fine step: Expected 29, got "+str(len(fine.positions))

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

    mover = MoverSingleHydrogenRotator(h, bondedNeighborLists)

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
    # Construct a MoverNH3Rotator that has ony hydrogen start out at 45 degrees around Z and the
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
    # Construct a MoverAromaticMethylRotator that has ony hydrogen start out at 45 degrees around Z and the
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
    # Construct a MoverTetrahedralMethylRotator that has ony hydrogen start out at 45 degrees around Z and the
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

  # @todo Test other Mover subclasses

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
     :param atom: iotbx.pdb.hierarchy.atom to be moved.
     :param axis: flex array with two scitbx::vec3<double> points, the first
     of which is the origin and the second is a vector pointing in the direction
     of the axis of rotation.  Positive rotations will be right-handed rotations
     around this axis.
     :param degrees: How much to rotate the atom around the axis.
     Positive rotation is right-handed around the axis.
     :returns the new location for the atom.
  '''
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
  d = -normal*atom.xyz
  t = - (d + (normal * _rvec3(axis[0])))
  nearPoint = _lvec3(axis[0]) + t * normal

  # Find the vector from the closest point towards the atom, which is its offset
  offset = _lvec3(atom.xyz) - nearPoint

  # Rotate the offset vector around the axis by the specified angle.  Add the new
  # offset to the closest point. Store this as the new location for this atom and angle.
  newOffset = offset.rotate_around_origin(_lvec3(axis[1]), degrees*math.pi/180)
  return nearPoint + newOffset

 # If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  ret = Test()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
