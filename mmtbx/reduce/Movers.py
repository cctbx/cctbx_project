##################################################################################
#                Copyright 2021-2023 Richardson Lab at Duke University
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

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import
import math
import scitbx.matrix, scitbx.math
from scitbx.array_family import flex
from iotbx import pdb
import mmtbx_probe_ext as probe
import traceback
from mmtbx.probe.Helpers import rvec3, lvec3, dihedralChoicesForRotatableHydrogens
from mmtbx_reduce_ext import RotateAtomDegreesAroundAxisDir, FindPosesFor

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
#  - type PositionReturn: (
#       flex<atom> atoms,
#       flex<flex<vec3>> positions,
#       flex<flex<probe.ExtraAtomInfo>> extraInfos,
#       flex<flex<bool>> deleteMes,
#       std::vector<double> preferenceEnergies
#    )
#       The atoms element has a list of all of the atoms to be adjusted.
#       The other elements each have a list of entries, where there is one entry
#     for each set of positions.  For the elements that have a double list, the
#     outer list is per entry and the inner list is per atom in the atoms element,
#     each with a corresponding list index.
#       The positions element has the new location of each atom in each set of positions.
#     This array may be shorter in length than the number of atoms because
#     some Movers do not need to change the position for all atoms (for flips, all atoms
#     are involved in fixup but not moved during optimization).  The index in
#     this array will match the index in the atoms array so the earliest atoms will be
#     changed if a subset is present.
#       The extraInfos element has the new ExtraAtomInfo for each atom in each set of positions.
#     This array may be shorter in length than the number of atoms (and may be empty) because
#     some Movers do not need to change the information for any or all atoms.  The index in
#     this array will match the index in the atoms array so the earliest atoms will be
#     changed if a subset is present.
#       The deleteMes element tells whether each atom in each set of positions should be
#     deleted.  This means that it should be ignored in all calculations and also should be
#     deleted from the model if this configuration is chosen.  This array may be shorter in
#     length than the number of atoms (and may be empty) because some Movers do not need to
#     change the information for any or all atoms.  The index in this array will match the
#     index in the atoms array so the earliest atoms will be deleted if a subset is present.
#       The preferenceEnergies entry holds an additional bias term that should be added to
#     the Probe score for each set of positions before comparing them with each other.
#  - PositionReturn CoarsePositions()
#  - PositionReturn FinePositions(coarseIndex)
#     The coarseIndex indicates the index (0 for the first) of the relevant coarse
#       orientation.
#     The list of atoms returned here must match the list returned for this index
#       by CoarsePositions().
#     The return values are the same as for CoarsePositions and they list potential
#       fine motions around the particular coarse position (not including the position
#       itself).  This function can be used by optimizers that wish to do heavy-weight
#       operations at a coarse resolution and then lightweight operations at a finer
#       scale; other optimizers may choose to ask for all of the fine positions and
#       append them to the coarse positions and globally optimize.
#     Note: Some Movers will return empty arrays.
#  - type FixUpReturn: (
#       flex<atom> atoms,
#       flex<vec3> positions,
#       flex<probe.ExtraAtomInfo> extraInfos,
#       flex<bool> deleteMes
#    )
#       These are lists with one entry per atom telling the same information as the
#     PositionReturn does, but only for a single set of positions.
#  - FixUpReturn FixUp(coarseIndex)
#     The coarseIndex indicates the index (0 for the first) of the relevant coarse
#       orientation that was finally chosen.
#     Note: This function may ask to modify atoms other than the ones reported in the
#       coarse and fine position functions.  This is to enable optimizing bonds
#       for Flippers.  None of these depend on fine index.
#     Note: Some Movers will return empty arrays.  This indicates that no fix-up is
#       needed and the atoms are in their final positions.
#
#  - str PoseDescription: (
#       int coarseIndex,   # Coarse position index
#       int fineIndex      # Fine position index
#       bool fixedUp       # Did we do fixup on Movers that have this?
#    )
#     Returns a human-readible description of the state.  This must have the same
#       number of words in all cases for all Movers to make things easier for a
#       program to parse.
#
# The caller is responsible for moving the specified atoms to their positions,
# modifying the ExtraAtomInfo, and deleting/ignoring them before dong any calculations
# with them.  After selecting the coarse or fine adjustment, they must call FixUp()
# and apply any final changes that it returns.
#
# The InteractionGraph.py script provides functions for determining which pairs of
# Movers have overlaps between movable atoms.

##################################################################################
import boost_adaptbx.boost.python as bp
bp.import_ext("mmtbx_reduce_ext")
from mmtbx_reduce_ext import PositionReturn
#class PositionReturn(object):
#  # Return type from CoarsePosition() and FinePosition() calls.
#  def __init__(self, atoms, positions, extraInfos, deleteMes, preferenceEnergies):
#    self.atoms = atoms
#    self.positions = positions
#    self.extraInfos = extraInfos
#    self.deleteMes = deleteMes
#    self.preferenceEnergies = preferenceEnergies

##################################################################################
class FixUpReturn(object):
  # Return type from FixUp() calls.
  def __init__(self, atoms, positions, extraInfos, deleteMes):
    self.atoms = atoms
    self.positions = positions
    self.extraInfos = extraInfos
    self.deleteMes = deleteMes

##################################################################################
class MoverNull(object):
  '''A trivial Mover that returns a single result atom at a single location.
     Useful as a simple and fast test case for programs that use Movers.
     It also serves as a basic example of how to develop a new Mover.
  '''
  def __init__(self, atom, extraAtomInfoMap):
    """Constructs a MoverNull.
       :param atom: Single atom to be moved.
       :param extraAtomInfoMap: probe.ExtraAtomInfoMap that can be used to look
       up the information for atoms whose values need to be changed.  Can be
       obtained by calling mmtbx.probe.Helpers.getExtraAtomInfo().
    """
    self._atom = atom
    # Make a copy of the extra information so that we can be sure it is not modified
    # elsewhere before we return it.
    self._extraAtomInfo = probe.ExtraAtomInfo(extraAtomInfoMap.getMappingFor(self._atom))
    self._coarsePositions = PositionReturn(
        # Single array of atoms
        [ self._atom ],
        # List of lists of positions.  Must be a list rather than a tuple
        [ [ [ self._atom.xyz[0], self._atom.xyz[1], self._atom.xyz[2] ] ] ],
        # Array of arrays of ExtraAtomInfos
        # Return a copy of our data so that someone won't modify it outside of us
        [ [ probe.ExtraAtomInfo(self._extraAtomInfo) ] ],
        # Array of array of DeleteMes
        [ [ False ] ],
        # Single array of preference energies
        [ 0.0 ]
    )

  def CoarsePositions(self):
    # returns: The original atom at its original position.
    return self._coarsePositions
  def FinePositions(self, coarseIndex):
    # returns: No fine positions for any coarse position.
    return PositionReturn([], [], [], [], [])
  def FixUp(self, coarseIndex):
    # No fixups for any coarse index.
    return FixUpReturn([], [], [], [])
  def PoseDescription(self, coarseIndex, fineIndex, fixedUp):
    if coarseIndex >= len(self.CoarsePositions().positions) or fineIndex is not None and (
        fineIndex > 0 and fineIndex >= len(self.FinePositions(0).positions)):
      return "Unrecognized state . ."
    else:
      return "Original location . ."

##################################################################################
class _MoverRotator(object):
  def __init__(self, atoms, axis, dihedral, offset, coarseRange, coarseStepDegrees = 15.0,
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
       :param dihedral: Angle of rotatation in degrees.  Specifies the clockwise
       rotation of the conventional first hydrogen (lowest numbered) in the atoms
       list compared to the atom in the conventional branch three neighbors away.
       The initial rotation for the atoms rotates them back to 0.
       :param offset: Angle of rotation in degrees.  Specifies the desired clockwise
       rotation of the conventional first hydrogen (lowest numbered) in the atoms
       list away from the dihedral.  This is specified separately because we need to
       know it for reporting purposes.  It is added to the negative of the dihedral
       to determine the initial rotation to make for all of the atoms.  It specifies
       the starting orientation (effective 0 degrees) that the coarse and fine steps
       are taken from.
       :param coarseRange: Range in degrees that the atoms can rotate around the axis
       from the starting position.  The range is considered closed on the negative
       end and open on the positive: [-coarseRange..coarseRange).  For example a
       range of 180 will rotate all the way around, from -180 (inclusive) to just
       slightly less than +180.
       :param coarseStepDegrees: The coarse step to take.
       To test only two orientations, doFineRotations should be set
       to False, coarseRange to 180, and coarseStepDegrees to 180.
       :param doFineRotations: Specifies whether fine rotations are allowed for
       this instance.  To test only two orientations, this should be set
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
    self._offset = offset

    # Rotate the atoms so that the conventional atom is at the expected initial
    # orientation before other rotations are made.  This rotates them around the axis by an angle that
    # is the sum of the offset and the dihedral, placing the conventional atom at the
    # specified offset from a 0 dihedral angle.  The other atoms maintain their relative rotations
    # with the conventional atom.
    for atom in atoms:
      atom.xyz = RotateAtomDegreesAroundAxisDir(axis[0], axis[1], atom, offset + dihedral)

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

    # Coarse positions
    self._coarsePositions = self._computeCoarsePositions()

    # Fine positions, one list per coarse index
    self._finePositions = self._computeFinePositions()

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
    return FindPosesFor(angles, self._atoms, self._axis[0], self._axis[1])

  def _computeCoarsePositions(self):
    return PositionReturn(self._atoms,
      self._posesFor(self._coarseAngles),
      [ [] ]  * len(self._coarseAngles),
      [ [] ]  * len(self._coarseAngles),
      self._preferencesFor(self._coarseAngles, self._preferredOrientationScale))

  def CoarsePositions(self):

    # Return the atoms, coarse-angle poses, and coarse-angle preferences
    return self._coarsePositions

  def _computeFinePositions(self):
    ret = []
    for coarseIndex in range(len(self._coarsePositions.positions)):

      if not self._doFineRotations:
        # No fine positions for any coarse position.
        ret.append(PositionReturn([], [], [], [], []))
        continue

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
      ret.append(PositionReturn(self._atoms,
        self._posesFor(angles),
        [ [] ] * len(angles),
        [ [] ] * len(angles),
        self._preferencesFor(angles, self._preferredOrientationScale))
        )
      continue

    return ret

  def FinePositions(self, coarseIndex):
    return self._finePositions[coarseIndex]

  def FixUp(self, coarseIndex):
    # No fixups for any coarse index.
    return FixUpReturn([], [], [], [])

  def PoseDescription(self, coarseIndex, fineIndex, fixedUp):
    if coarseIndex >= len(self.CoarsePositions().positions) or fineIndex is not None and (
        fineIndex > 0 and fineIndex >= len(self.FinePositions(0).positions)):
      return "Unrecognized state . ."
    else:
      fineOffset = 0
      if fineIndex is not None and fineIndex >= 0:
        fineOffset = self._fineAngles[fineIndex]
      angle = self._offset + self._coarseAngles[coarseIndex] + fineOffset
      while angle > 180: angle -= 360
      while angle < -180: angle += 360
      return "Angle {:.1f} deg .".format(angle)

##################################################################################
class MoverSingleHydrogenRotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, extraAtomInfoMap,
                hParameters, potentialAcceptors = [],
                potentialTouches = [],
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
       bonded atoms.  Can be obtained by calling probe.Helpers.getBondedNeighborLists().
       :param extraAtomInfoMap: probe.ExtraAtomInfoMap that can be used to look
       up the information for atoms.  Can be obtained by calling
       mmtbx.probe.Helpers.getExtraAtomInfo().
       :param hParameters: List indexed by sequence ID that stores the riding
       coefficients for hydrogens that have associated dihedral angles.  This can be
       obtained by calling model.setup_riding_h_manager() and then model.get_riding_h_manager().
       :param potentialAcceptors: A flex array of atoms that are nearby potential acceptors.
       Coarse orientations are added that aim the hydrogen in the direction of these potential
       partners.
       :param potentialTouches: A flex array of atoms that are nearby potential touches/clashes.
       :param coarseStepDegrees: The coarse step to take.
       :param fineStepDegrees: The fine step to take.
       :param preferredOrientationScale: How much to scale the preferred-orientation
       energy by before adding it to the total.
    """
    # @todo All callers and tests

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
    # NOTE: We replace the preferenceFunction with None below to match original reduce
    # behavior. We leave the math in here so that if we later decide to have a preference
    # it will be correctly determined based on the number of neighbors.
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
    # Original Reduce did not have a preference, so we remove that here.
    preferenceFunction = None

    # Determine the axis to rotate around, which starts at the partner atom and points at the neighbor.
    normal = (rvec3(neighbor.xyz) - rvec3(partner.xyz)).normalize()
    axis = flex.vec3_double([partner.xyz, normal])

    # Make a list that contains the hydrogen and its bonded neighbor so that its neighbor will
    # be included in energy calculations.  The neighbor will be rotating about an axis that
    # includes it, so will not be moved.
    atoms = [ atom, neighbor ]

    # Move the atom so that it is in one of the preferred locations.  The preferred location depends on
    # whether there are two or three friends, it wants to be in the plane if there are two of them and it
    # wants to be between two edges if there are three.  It turns out that rotating it to point away from
    # the conventional-branch friend works in both of those cases.
    conventionalH, conventionalFriend = dihedralChoicesForRotatableHydrogens(atoms,
      hParameters, friends)
    sites = [ conventionalH.xyz, partner.xyz, neighbor.xyz, conventionalFriend.xyz ]
    dihedral = scitbx.math.dihedral_angle(sites=sites, deg=True)
    offset = 180

    # Construct our parent class, which will do all of the actual work based on our inputs.
    _MoverRotator.__init__(self, atoms, axis, dihedral, offset, 180, coarseStepDegrees = coarseStepDegrees,
      fineStepDegrees = fineStepDegrees, preferenceFunction = preferenceFunction,
      preferredOrientationScale = preferredOrientationScale)

    # Now add orientations that point in the direction of the potential acceptors.
    # Then select from the original angles the one that has the best contact with
    # any touching atoms that is at least the coarse step size away from pointing
    # towards an acceptor.
    # We replace the original coarse angles with this set of 1+ angles to reduce the
    # number of angles to check and to ensure that we always check potential acceptors.

    ###########################
    # Helper utility function to sort atoms consistently from run to run so that we get
    # the same ordering on coarse angles.
    def atomID(a):
      # Return the ID of the atom, which includes its chain, residue name,
      # residue number, atom name, and alternate separated by spaces. This
      # is used to sort the atoms. This must work in the case where we have
      # test atoms that are not completely fleshed out.
      try:
        return ( a.parent().parent().parent().id + a.parent().resname +
          str(a.parent().parent().resseq_as_int()) + a.name + a.parent().altloc )
      except Exception:
        return ""
    #
    ###########################

    # Compute the dihedral angle from the Hydrogen to each potential acceptor through
    # the partner and neighbor.  This is the amount to rotate the hydrogen by.
    # Sort these so that we get the same order each time the program is run.
    acceptorAngles = []
    for a in sorted(potentialAcceptors, key=lambda x:atomID(x)):
      sites = [ atom.xyz, partner.xyz, neighbor.xyz, a.xyz ]
      degrees = scitbx.math.dihedral_angle(sites=sites, deg=True)
      # The degrees can be None in degenerate cases.  We avoid adding an entry in that case.
      # The atom will have been placed in the correct initial orientation by the _MoverRotator
      # constructor, so our dihedral measurement here is with respect to that new location.
      # This means that we merely have to add the degrees offset to its current rotation.
      if degrees is not None:
        acceptorAngles.append(degrees)

    # Find the coarse angle that has the least best contact with potential touches
    # which may also be one of the acceptors (for a weak hydrogen bond, the score
    # can be better for a touch than for an overlap).
    # This is the one whose gap is closest to 0.
    bestTouchAngle = 0
    bestTouchGap = 1e100
    ra = extraAtomInfoMap.getMappingFor(atom).vdwRadius
    for i, ang in enumerate(self._coarseAngles):
      # Find minimum gap with clashing atoms at this angle. This number is
      # negative when there is a clash. It reports the atom that we're most
      # in contact with at this angle.
      minGap = 1e100
      for pt in potentialTouches:
        rt = extraAtomInfoMap.getMappingFor(pt).vdwRadius
        # Measure from the first atom's position (the Hydrogen) at this position to the potential touch
        distance = (rvec3(self._coarsePositions.positions[i][0]) - rvec3(pt.xyz)).length()
        gap = distance - (ra + rt)
        if gap < minGap:
          minGap = gap
      # Find the minimum gap distance that is closest to zero, either
      # above or below zero. This is the one with the best just-touch
      # value.
      if abs(minGap) < bestTouchGap:
        bestTouchGap = abs(minGap)
        bestTouchAngle = ang

    # Check every coarse angle to see whether it is at least 45 degrees away from
    # any of the acceptor angles or the best touch angle.  If it is, then we add
    # it to the list of coarse angles to try. This ensures that we try at least
    # sparsely in all directions.

    sofar = [bestTouchAngle]
    sofar.extend(acceptorAngles)
    for ang in self._coarseAngles:
      minAng = 360
      for a in sofar:
        diff = abs(a - ang)
        if diff < minAng:
          minAng = diff
      if minAng >= 45:
        sofar.append(ang)

    # Replace the coarse angles with the ones that we have found.
    self._coarseAngles = sofar

    # Recompute the coarse and fine positions given the new angles we want to test
    self._coarsePositions = self._computeCoarsePositions()
    self._finePositions = self._computeFinePositions()

##################################################################################
class MoverNH3Rotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, hParameters, coarseStepDegrees = 15.0,
                fineStepDegrees = 1.0, preferredOrientationScale = 1.0):
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
       bonded atoms.  Can be obtained by calling probe.Helpers.getBondedNeighborLists().
       the value.
       :param hParameters: List indexed by sequence ID that stores the riding
       coefficients for hydrogens that have associated dihedral angles.  This can be
       obtained by calling model.setup_riding_h_manager() and then model.get_riding_h_manager().
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
      if a.element_is_hydrogen():
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
    if len(friends) < 3:
      raise ValueError("MoverNH3Rotator(): Partner does not have at least three bonded friends")

    preferenceFunction = None

    # If we have exactly three neighbors, then we set a preference function that prefers
    # the staggered orientation. If there are more, we just leave things alone.
    if len(friends) == 3:

      # Set the preference function to prefer 120-degree rotations away from the starting location.
      # @todo Consider parameterizing the magic constant of 0.1 for the preference magnitude
      # (for example, we might use preferredOrientationScale for this and not set a separate one?)
      def preferenceFunction(degrees): return 0.1 + 0.1 * math.cos(degrees * (math.pi/180) * (360/120))

    # Determine the axis to rotate around, which starts at the partner and points at the neighbor.
    normal = (rvec3(neighbor.xyz) - rvec3(partner.xyz)).normalize()
    axis = flex.vec3_double([partner.xyz, normal])

    # Move the Hydrogens so that they are in one of the preferred locations by rotating the
    # conventional one of them to point away from the conventional one of the friends.
    conventionalH, conventionalFriend = dihedralChoicesForRotatableHydrogens(hydrogens,
      hParameters, friends)
    sites = [ conventionalH.xyz, partner.xyz, neighbor.xyz, conventionalFriend.xyz ]
    dihedral = scitbx.math.dihedral_angle(sites=sites, deg=True)
    offset = 180

    # Make a list that contains the hydrogens and their bonded neighbor so that the neighbor will
    # be included in energy calculations.  The neighbor will be rotating about an axis that
    # includes it, so will not be moved.
    atoms = [ neighbor ]
    atoms.extend(hydrogens)

    # Construct our parent class, which will do all of the actual work based on our inputs.
    _MoverRotator.__init__(self, atoms, axis, dihedral, offset, 60, coarseStepDegrees,
      fineStepDegrees = fineStepDegrees, preferenceFunction = preferenceFunction,
      preferredOrientationScale = preferredOrientationScale)

##################################################################################
class MoverAromaticMethylRotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, hParameters):
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
       bonded atoms.  Can be obtained by calling probe.Helpers.getBondedNeighborLists().
       :param hParameters: List indexed by sequence ID that stores the riding
       coefficients for hydrogens that have associated dihedral angles.  This can be
       obtained by calling model.setup_riding_h_manager() and then model.get_riding_h_manager().
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
      if a.element_is_hydrogen():
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
    normal = (rvec3(neighbor.xyz) - rvec3(partner.xyz)).normalize()
    axis = flex.vec3_double([partner.xyz, normal])

    # Move the Hydrogens so that they are in one of the preferred locations by rotating the
    # conventional one of them to point away from the conventional one of the friends plus 90 degrees.
    conventionalH, conventionalFriend = dihedralChoicesForRotatableHydrogens(hydrogens,
      hParameters, friends)
    sites = [ conventionalH.xyz, partner.xyz, neighbor.xyz, conventionalFriend.xyz ]
    dihedral = scitbx.math.dihedral_angle(sites=sites, deg=True)
    offset = 180 + 90

    # Make a list that contains the hydrogens and their bonded neighbor so that the neighbor will
    # be included in energy calculations.  The neighbor will be rotating about an axis that
    # includes it, so will not be moved.
    atoms = [ neighbor ]
    atoms.extend(hydrogens)

    # Construct our parent class, which will do all of the actual work based on our inputs.
    # We have a coarse step size of 180 degrees and a range of 180 degrees and do not
    # allow fine rotations.
    _MoverRotator.__init__(self, atoms, axis, dihedral, offset, 180, 180, doFineRotations = False)

##################################################################################
class MoverTetrahedralMethylRotator(_MoverRotator):
  def __init__(self, atom, bondedNeighborLists, hParameters, coarseStepDegrees = 15.0,
                  fineStepDegrees = 1.0, preferredOrientationScale = 1.0):
    """ A Mover that rotates three Hydrogens around an axis from their bonded Carbon neighbor
       to the single bonded partner of its partner.  This is designed for use with tetrahedral
       partners whose partner-partner atoms are bonded to three friends but can also be used
       with Methyl's attached to partners with a single other bond like the S in MET.
       The starting orientation has the Hydrogens pointing between the side of the tetrahedron.
       It can rotate to any angle to optimize for hydrogen bonds.
       Note: Reduce does not normally rotate these groups because it is not worth the computational
       cost to do so (and because this can mask mis-placed atoms by forming spurious hydrogen
       bonds), but it will construct them so that they will be aligned staggered to the
       tetrahedron.
       Construct these in Reduce so that they will be staggered but do not optimize them.
       :param atom: Carbon atom bonded to the three Hydrogens that will be rotated.
       It must be bonded to three Hydrogens and a single other
       atom, and the other atom must be bonded to three other atoms.  NOTE: As a side
       effect, the Hydrogens are immediately rotated to lie staggered.
       :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
       structure that the atom from the first parameter interacts with that lists all of the
       bonded atoms.  Can be obtained by calling probe.Helpers.getBondedNeighborLists().
       :param hParameters: List indexed by sequence ID that stores the riding
       coefficients for hydrogens that have associated dihedral angles.  This can be
       obtained by calling model.setup_riding_h_manager() and then model.get_riding_h_manager().
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
      if a.element_is_hydrogen():
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
    if len(friends) != 1 and len(friends) != 3:
      raise ValueError("MoverTetrahedralMethylRotator(): Partner does not have one or three bonded friends")

    # Determine the axis to rotate around, which starts at the partner and points at the neighbor.
    normal = (rvec3(neighbor.xyz) - rvec3(partner.xyz)).normalize()
    axis = flex.vec3_double([partner.xyz, normal])

    # Set the preference function to like 120-degree rotations away from the starting location.
    # @todo Consider parameterizing the magic constant of 0.1 for the preference magnitude
    def preferenceFunction(degrees): return 0.1 + 0.1 * math.cos(degrees * (math.pi/180) * (360/120))

    # Move the Hydrogens so that they are in one of the preferred locations by rotating the
    # conventional one of them to point away from the conventional one of the friends.
    conventionalH, conventionalFriend = dihedralChoicesForRotatableHydrogens(hydrogens,
      hParameters, friends)
    sites = [ conventionalH.xyz, partner.xyz, neighbor.xyz, conventionalFriend.xyz ]
    dihedral = scitbx.math.dihedral_angle(sites=sites, deg=True)
    offset = 180

    # Make a list that contains the hydrogens and their bonded neighbor so that the neighbor will
    # be included in energy calculations.  The neighbor will be rotating about an axis that
    # includes it, so will not be moved.
    atoms = [ neighbor ]
    atoms.extend(hydrogens)

    # Construct our parent class, which will do all of the actual work based on our inputs.
    _MoverRotator.__init__(self, atoms, axis, dihedral, offset, 180, fineStepDegrees = fineStepDegrees,
      preferenceFunction = preferenceFunction,
      preferredOrientationScale = preferredOrientationScale)

##################################################################################
class MoverAmideFlip(object):
  def __init__(self, nh2Atom, caAtomName, bondedNeighborLists, nonFlipPreference):
    """Constructs a Mover that will handle flipping an NH2 with an O, both of which
       are attached to the same Carbon atom (and each of which has no other bonds).
       This Mover will move the hydrogens so that they are located at +/-120 degrees from the
       Carbon-Oxygen bond in the plane of the Nitrogen, Carbon, and Oxygen and at the same
       distance from the Nitrogen as at least one of them started out.
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
       bonded atoms.  Can be obtained by calling probe.Helpers.getBondedNeighborLists().
       :param nonFlipPreference: Score amount by which the original orientation is preferred
       over the flipped orientation.  This will bias things so that the flipped orientation
       is not chosen unless it is this much better than the original.
    """

    # Store away our constructor arguments that we need for later.
    self._nonFlipPreference = nonFlipPreference

    # Verify that we've been run on a valid structure and get a list of all of the
    # atoms up to and including the pivot atom.
    if nh2Atom.element != "N":
      raise ValueError("MoverAmideFlip(): nh2Atom is not a Nitrogen")
    partners = bondedNeighborLists[nh2Atom]
    if len(partners) != 3:
      raise ValueError("MoverAmideFlip(): nh2Atom does not have three bonded neighbors")
    nh2Hydrogens = []
    hinge = None
    for a in partners:
      if a.element_is_hydrogen():
        nh2Hydrogens.append(a)
      else:
        hinge = a
    if len(nh2Hydrogens) != 2:
      raise ValueError("MoverAmideFlip(): nh2Atom does not have two bonded hydrogens")
    if hinge is None:
      raise ValueError("MoverAmideFlip(): nh2Atom does not have bonded (hinge) Carbon friend")

    bonded = bondedNeighborLists[hinge]
    oyxgen = None
    pivot = None
    for b in bonded:
      if b.element == "O":
        oxygen = b
      elif b.element == "C":
        pivot = b
    if pivot is None:
      raise ValueError("MoverAmideFlip(): Hinge does not have bonded (pivot) Carbon friend")
    if oxygen is None:
      raise ValueError("MoverAmideFlip(): Hinge does not have bonded oxygen friend")
    if len(bondedNeighborLists[oxygen]) != 1:
      raise ValueError("MoverAmideFlip(): Oxygen has more than one bonded neighbor")
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
        raise ValueError("MoverAmideFlip(): Linker chain has an element with other than four bonds")
      for b in bonded:
        if b.i_seq == prevID:
          continue
        elif b.element == "C":
          link = b
        elif b.element_is_hydrogen():
          linkerHydrogens.append(b)
      if link is None:
        raise ValueError("MoverAmideFlip(): Did not find Carbon in linker chain step")
      if link.name.strip().upper() == caAtomName.strip().upper():
        caAtom = link
        foundAlpha = True
      else:
        linkers.append(link)
        prevID = cur.i_seq
        cur = link
    if len(linkerHydrogens) != 2*(len(linkers)+1):
      raise ValueError("MoverAmideFlip(): Linker carbons do not have two Hydrogens each: "+
        str(len(linkers))+","+str(len(linkerHydrogens)))
    self._atoms.extend(linkers)
    self._atoms.extend(linkerHydrogens)

    #########################
    # Compute the original positions for the Hydrogens such that they are at the same distance from
    # the Nitrogen as one of them started out and located at +/-120 degrees from the
    # Carbon-Nitrogen bond in the plane of the Nitrogen, Carbon, and Oxygen. The first one should
    # point back towards the mainchain and the second more towards the Oxygen.
    # @todo This assumes that all placement is done like Hydrogenate does, where the first-listed
    # atom is the lower-numbered one.
    cToN = lvec3(nh2Atom.xyz) - lvec3(hinge.xyz)
    cToO = rvec3(oxygen.xyz) - rvec3(hinge.xyz)

    # Normal to the plane containing Nitrogen, Carbon, and Oxygen
    normal = lvec3(scitbx.matrix.cross_product_matrix(cToN) * cToO).normalize()

    hBond0Len = (rvec3(nh2Hydrogens[0].xyz) - rvec3(nh2Atom.xyz)).length()
    hBond1Len = (rvec3(nh2Hydrogens[1].xyz) - rvec3(nh2Atom.xyz)).length()
    nh2Hydrogens[0].xyz = lvec3(nh2Atom.xyz) + ((-cToN.normalize()) * hBond0Len).rotate_around_origin(normal, 120 * math.pi/180)
    nh2Hydrogens[1].xyz = lvec3(nh2Atom.xyz) + ((-cToN.normalize()) * hBond1Len).rotate_around_origin(normal,-120 * math.pi/180)

    #########################
    # Compute the new positions for the Hydrogens such that they are at the same distance from
    # the Oxygen as one of them is from the Nitrogen and located at +/-120 degrees from the
    # Carbon-Oxygen bond in the plane of the Nitrogen, Carbon, and Oxygen.
    cToO = lvec3(oxygen.xyz) - lvec3(hinge.xyz)
    cToN = rvec3(nh2Atom.xyz) - rvec3(hinge.xyz)

    # Normal to the plane containing Nitrogen, Carbon, and Oxygen
    normal = lvec3(scitbx.matrix.cross_product_matrix(cToO) * cToN).normalize()

    hBond0Len = (rvec3(nh2Hydrogens[0].xyz) - rvec3(nh2Atom.xyz)).length()
    hBond1Len = (rvec3(nh2Hydrogens[1].xyz) - rvec3(nh2Atom.xyz)).length()
    newH0 = lvec3(oxygen.xyz) + ((-cToO.normalize()) * hBond0Len).rotate_around_origin(normal, 120 * math.pi/180)
    newH1 = lvec3(oxygen.xyz) + ((-cToO.normalize()) * hBond1Len).rotate_around_origin(normal,-120 * math.pi/180)

    #########################
    # Compute the list of positions for all of the atoms. This consists of the original
    # location and the flipped location where we swap the locations of the two heavy atoms
    # and bring the Hydrogens along for the ride.

    startPos = []
    for a in self._atoms:
      startPos.append(a.xyz)

    newPos = startPos[:]
    newPos[0] = newH0
    newPos[1] = newH1
    newPos[2] = oxygen.xyz
    newPos[3] = nh2Atom.xyz

    # Only consider the first 5 atoms when optimizing, the four that move and the one they may shield
    self._coarsePositions = [ startPos[:5], newPos[:5] ]

    #########################
    # Compute the list of Fixup returns.
    hingeIndex = 4
    firstDockIndex = 3
    secondDockIndex = 2
    movable = _rotateHingeDock(self._atoms, hingeIndex, firstDockIndex, secondDockIndex, caAtom)

    # No fix-up for coarse position 0, do the above adjustment for position 1
    self._fixUpPositions = [ [], movable ]

  def CoarsePositions(self):
    # returns: The two possible coarse positions with an energy penalty of -nonFlipPreference for the flipped.
    # The -nonFlipPreference penalty is to prevent uncertain flips from happening -- unless the
    # score is this much better we leave it alone.
    return PositionReturn(self._atoms, self._coarsePositions,
      [ [], [] ],
      [ [], [] ],
      [0.0, -self._nonFlipPreference])

  def FinePositions(self, coarseIndex):
    # returns: No fine positions for any coarse position.
    return PositionReturn([], [], [], [], [])

  def FixUp(self, coarseIndex):
    # Return the appropriate fixup
    return FixUpReturn(self._atoms, self._fixUpPositions[coarseIndex], [], [])

  def PoseDescription(self, coarseIndex, fineIndex, fixedUp):
    if coarseIndex == 1:
      if fixedUp:
        fString = 'AnglesAdjusted'
      else:
        fString = 'AnglesNotAdjusted'
    else:
      fString = '.'
    if coarseIndex >= len(self.CoarsePositions().positions) or fineIndex is not None and (
        fineIndex > 0 and fineIndex >= len(self.FinePositions(0).positions)):
      return "Unrecognized state ."
    elif coarseIndex == 0:
      return "Unflipped . . {}".format(fString)
    else:
      return "Flipped . . {}".format(fString)

##################################################################################
class MoverHisFlip(object):
  def __init__(self, ne2Atom, bondedNeighborLists, extraAtomInfoMap, nonFlipPreference,
               enabledFlipStates = 3, enableFixup = True):
    """Constructs a Mover that will handle flipping a Histidine ring.
       This Mover uses a simple swap of the center positions of the heavy atoms (with
       repositioning of the Hydrogens to lie in the same directions)
       for its testing, but during FixUp it adjusts the bond lengths and angles for
       additional atoms per Protein Science Vol 27:293-315.
       :param ne2Atom: NE2 atom within the Histidine ring.
       :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
       structure that the atom from the first parameter interacts with that lists all of the
       bonded atoms.  Can be obtained by calling probe.Helpers.getBondedNeighborLists().
       :param extraAtomInfoMap: probe.ExtraAtomInfoMap that can be used to look
       up the information for atoms whose values need to be changed.  Can be
       obtained by calling mmtbx.probe.Helpers.getExtraAtomInfo().
       :param nonFlipPreference: Score amount by which the original orientation is preferred
       over the flipped orientation.  This will bias things so that the flipped orientation
       :param enabledFlipStates: Which flip states are enabled? A value of 3 enables both,
       a value of 1 enables non-flipped and a value of 2 enables flipped. If only one state
       is enabled.
       :param enableFixup: This allows the constructor to disable fixup by making it the same
       as not fixed up. This is used by the FlipKin generation code.
    """

    # Store away our constructor arguments that we need for later.
    self._nonFlipPreference = nonFlipPreference
    self._enabledFlipStates = enabledFlipStates
    self._enableFixup = enableFixup

    # Verify that we've been run on a valid structure and get a list of all of the
    # atoms up to and including the pivot atom.
    if ne2Atom.element != "N":
      raise ValueError("MoverHisFlip(): ne2Atom is not a Nitrogen")
    partners = bondedNeighborLists[ne2Atom]
    if len(partners) < 3:
      raise ValueError("MoverHisFlip(): NE2 does not have three bonded neighbors")
    hydrogens = []
    carbons = []
    for a in partners:
      if a.element_is_hydrogen():
        hydrogens.append(a)
      elif a.element == "C":
        carbons.append(a)
    if len(hydrogens) != 1:
      raise ValueError("MoverHisFlip(): NE2 does not have one bonded hydrogen (probably ionically bound)")
    if len(carbons) != 2:
      raise ValueError("MoverHisFlip(): NE2 does not have two bonded carbons")
    ne2HAtom = hydrogens[0]

    # Determine if the first Carbon is CE1, which is bonded to two Nitrogens (zero Carbons).  If not,
    # then the second one must be.  Fill in CD2 and CE1 based on this information, then we
    # can check them and continue parsing.
    cTest = carbons[0]
    bonded = bondedNeighborLists[cTest]
    if len(bonded) != 3:
      raise ValueError("MoverHisFlip(): NE2 neighbor does not have three bonded neighbors")
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
      raise ValueError("MoverHisFlip(): CE1 does not have three bonded neighbors")
    ce1HAtom = None
    for b in bonded:
      if b.element_is_hydrogen():
        ce1HAtom = b
    if ce1HAtom is None:
      raise ValueError("MoverHisFlip(): Could not find Hydrogen attached to CE1")

    # Find the Hydrogen on CD2.
    bonded = bondedNeighborLists[cd2Atom]
    if len(bonded) != 3:
      raise ValueError("MoverHisFlip(): CD2 does not have three bonded neighbors")
    cd2HAtom = None
    for b in bonded:
      if b.element_is_hydrogen():
        cd2HAtom = b
    if cd2HAtom is None:
      raise ValueError("MoverHisFlip(): CD2 does not have a bonded hydrogen")

    # Find CG on the other side of CD2.
    cgAtom = None
    bonded = bondedNeighborLists[cd2Atom]
    if len(bonded) < 2:
      raise ValueError("MoverHisFlip(): CD2 does not have at least two bonded neighbors")
    for b in bonded:
      if b.i_seq != cd2Atom.i_seq and b.element == "C":
        cgAtom = b
    if cgAtom is None:
      raise ValueError("MoverHisFlip(): Could not find CG")

    # Find CB and ND1 on the other side of CG
    cbAtom = None
    nd1Atom = None
    bonded = bondedNeighborLists[cgAtom]
    if len(bonded) != 3:
      raise ValueError("MoverHisFlip(): CG does not have three bonded neighbors")
    for b in bonded:
      if b.i_seq == cd2Atom.i_seq: # We already know about CD2
        continue
      elif b.element == "N":
        nd1Atom = b
      elif b.element == "C":
        cbAtom = b
    if nd1Atom is None:
      raise ValueError("MoverHisFlip(): Could not find ND1")
    if cbAtom is None:
      raise ValueError("MoverHisFlip(): Could not find CB")

    # Find the Hydrogen attached to ND1
    partners = bondedNeighborLists[nd1Atom]
    if len(partners) < 3:
      raise ValueError("MoverHisFlip(): ND1 does not have three bonded neighbors")
    hydrogens = []
    carbons = []
    for a in partners:
      if a.element_is_hydrogen():
        hydrogens.append(a)
      elif a.element == "C":
        carbons.append(a)
    if len(hydrogens) != 1:
      raise ValueError("MoverHisFlip(): ND1 does not have one bonded hydrogen (probably ionically bound)")
    if len(carbons) != 2:
      raise ValueError("MoverHisFlip(): ND1 does not have two bonded carbons")
    nd1HAtom = hydrogens[0]

    # Find CA on the other side of CB and find CBs Hydrogens
    caAtom = None
    cbHydrogens = []
    bonded = bondedNeighborLists[cbAtom]
    if len(bonded) != 4:
      raise ValueError("MoverHisFlip(): CB does not have four bonded neighbors, has "+str(len(bonded)))
    for b in bonded:
      if b.i_seq == cgAtom.i_seq:
        continue
      elif b.element == "C":
        caAtom = b
      elif b.element_is_hydrogen():
        cbHydrogens.append(b)
    if caAtom is None:
      raise ValueError("MoverHisFlip(): Could not find CA")
    if len(cbHydrogens) != 2:
      raise ValueError("MoverHisFlip(): Could not find Hydrogens on CB")

    self._atoms = [ ne2Atom, ne2HAtom, ce1Atom, ce1HAtom, nd1Atom, nd1HAtom, cd2Atom, cd2HAtom,
      cgAtom, cbAtom, cbHydrogens[0], cbHydrogens[1] ]

    #########################
    # Compute the new positions for the Hydrogens such that they are at the same distance from
    # their swapped parent atoms and in the direction of the Hydrogens from the original atoms at
    # each location.  We swap ND1 with CD2 and NE2 with CE1.
    nd1HVec = lvec3(nd1HAtom.xyz) - lvec3(nd1Atom.xyz)
    ne2HVec = lvec3(ne2HAtom.xyz) - lvec3(ne2Atom.xyz)
    ce1HVec = lvec3(ce1HAtom.xyz) - lvec3(ce1Atom.xyz)
    cd2HVec = lvec3(cd2HAtom.xyz) - lvec3(cd2Atom.xyz)

    nd1HNew = lvec3(cd2Atom.xyz) + nd1HVec.length() * cd2HVec.normalize()
    cd2HNew = lvec3(nd1Atom.xyz) + cd2HVec.length() * nd1HVec.normalize()
    ce1HNew = lvec3(ne2Atom.xyz) + ce1HVec.length() * ne2HVec.normalize()
    ne2HNew = lvec3(ce1Atom.xyz) + ne2HVec.length() * ce1HVec.normalize()

    #########################
    # There are eight possible states for the flipped Histidine.  The first four
    # use the original orientation and the second four use the swapped or fixed-up
    # orientation.  The swapped orientation is used for the coarse tests and if
    # one is accepted, then the fixed-up orienations are used in place of the swapped.
    #   For each location, we have four cases of Hydrogen
    # placement; both (as found), first N Hydrogen removed, second N Hydrogen
    # removed, and both Hydrogens removed.  When a Hydrogen is removed, the associated
    # N is turned into an Acceptor; when it is present, the N is not an acceptor.

    #########################
    # Compute the list of positions for all of the atoms. This consists of the original
    # location and the flipped location where we swap the locations of the two pairs of heavy atoms
    # and bring the Hydrogens along for the ride.
    # We only put the options in place for the enabled flip states.

    startPos = []
    for a in self._atoms:
      startPos.append(a.xyz)

    newPos = startPos[:]
    newPos[0] = ce1Atom.xyz   # ne2 moved to this location
    newPos[1] = ne2HNew
    newPos[2] = ne2Atom.xyz   # ce1 moved to this location
    newPos[3] = ce1HNew
    newPos[4] = cd2Atom.xyz   # nd1 moved to this location
    newPos[5] = nd1HNew
    newPos[6] = nd1Atom.xyz   # cd2 moved to this location
    newPos[7] = cd2HNew

    self._coarsePositions = []
    if self._enabledFlipStates & 1:
      for i in range(4):
        # Only move the first 9 atoms when optimizing, the ones that move and the one they may shield.
        self._coarsePositions.append(startPos[:9])
    if self._enabledFlipStates & 2:
      for i in range(4):
        # Only move the first 9 atoms when optimizing, the ones that move and the one they may shield.
        self._coarsePositions.append(newPos[:9])

    #########################
    # Compute the list of Fixup returns.
    hingeIndex = 8
    firstDockIndex = 0
    secondDockIndex = 2
    fixedUp = _rotateHingeDock(self._atoms, hingeIndex, firstDockIndex, secondDockIndex, caAtom)

    # No fix-up for coarse positions 0-3, do the above adjustment for position3 4-7
    self._fixUpPositions = [ ]
    if self._enabledFlipStates & 1:
      for i in range(4):
        self._fixUpPositions.append([])
    if self._enabledFlipStates & 2:
      for i in range(4):
        if enableFixup:
          self._fixUpPositions.append(fixedUp)
        else:
          self._fixUpPositions.append([])

    #########################
    # Compute the ExtraAtomInfo and deleteMe values.  They are all as provided and False
    # for the initial configuration (0th and 4th).  The first Hydrogen and
    # Nitrogen are removed and adjusted for cases 1&3, 5&7, the second for 2&3 and 6&7.
    # We make copies of each by constructing new ones so we can independently change them.
    self._extras = []
    self._deleteMes = []
    for i in range(len(self._coarsePositions)):
      # Copy the initial values
      extras = []
      deleteMes = []
      for a in self._atoms:
        extras.append(probe.ExtraAtomInfo(extraAtomInfoMap.getMappingFor(a)))
        deleteMes.append(False)

      # Replace any that need it for this configuration.
      if i % 4 == 1 or i % 4 == 3: # Remove the Hydrogen from NE2
        extras[0].isAcceptor = True
        deleteMes[1] = True

      if i % 4 == 2 or i % 4 == 3: # Remove the Hydrogen from ND1
        extras[4].isAcceptor = True
        deleteMes[5] = True

      # Append to the lists.
      self._extras.append(extras)
      self._deleteMes.append(deleteMes)

    #########################
    # Compute the preference energies.
    # There is an energy penalty of -nonFlipPreference
    # for the flipped orientations, and a penalty of -0.05 for keeping both Hydrogens;
    # The doubly-deprotenated case (both Hydrogens removed) has a penalty of -1.0.
    # The -nonFlipPreference penalty is to prevent uncertain flips from happening -- unless the
    # score is this much better we leave it alone.
    self._preferenceEnergies = []
    if self._enabledFlipStates & 1:
      self._preferenceEnergies.extend([ 0.0 - 0.05,  0.0,  0.0,  0.0 - 1.0])
    if self._enabledFlipStates & 2:
      self._preferenceEnergies.extend([
        -self._nonFlipPreference - 0.05,
        -self._nonFlipPreference,
        -self._nonFlipPreference,
        -self._nonFlipPreference - 1.0])

  def CoarsePositions(self):
    # returns: The potential coarse positions.
    return PositionReturn(self._atoms, self._coarsePositions,
      self._extras, self._deleteMes, self._preferenceEnergies)

  def FinePositions(self, coarseIndex):
    # returns: No fine positions for any coarse position.
    return PositionReturn([], [], [], [], [])

  def FixUp(self, coarseIndex):
    # Return the appropriate fixup
    return FixUpReturn(self._atoms, self._fixUpPositions[coarseIndex],
      self._extras[coarseIndex], self._deleteMes[coarseIndex])

  def PoseDescription(self, coarseIndex, fineIndex, fixedUp):
    if self._enabledFlipStates == 3:
      if coarseIndex >= len(self.CoarsePositions().positions) or fineIndex is not None and (
          fineIndex > 0 and fineIndex >= len(self.FinePositions(0).positions)):
        return "Unrecognized state . ."
      elif coarseIndex < 4:
        ret = "Unflipped"
      else:
        ret = "Flipped"
    elif self._enabledFlipStates == 1:
      ret = "Unflipped"
    elif self._enabledFlipStates == 2:
      ret = "Flipped"
    else:
      return "Unrecognized flip states ."

    if coarseIndex % 4 == 0 or coarseIndex % 4 == 1:
      ret += " HD1Placed"
    else:
      ret += " HD1NotPlaced"
    if coarseIndex % 4 == 0 or coarseIndex % 4 == 2:
      ret += " HE2Placed"
    else:
      ret += " HE2NotPlaced"
    # First, check to see if we are in a state where fixup might occur; this is when we
    # have been locked into the flipped state or when we are allowed to choose either state
    # and we're in the uppper half of the coarse positions. If not, we just print '.'
    # because fixup was never an option. If we are, check whether fixup is enabled and if
    # so and we've been told above to fix up, then report angles adjusted.
    if (self._enabledFlipStates == 2) or (
       (self._enabledFlipStates == 3) and (coarseIndex >= len(self._coarsePositions) / 2)):
      if self._enableFixup and fixedUp:
        ret += ' AnglesAdjusted'
      else:
        ret += ' AnglesNotAdjusted'
    else:
      ret += ' .'
    return ret

##################################################################################
# Internal helper functions for angle manipulation.
def _rotateOppositeFriend(atom, axis, partner, friend):
  '''Rotate the atom to point away from the friend.  This means placing it
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
  normal = rvec3(axis[1])
  friendFromPartner = lvec3(friend.xyz) - lvec3(partner.xyz)
  alongAxisComponent = (friendFromPartner * normal)
  friendFromPartner = rvec3(friend.xyz) - rvec3(partner.xyz)
  inPlane = friendFromPartner - normal*alongAxisComponent
  normalizedOffset = -inPlane.normalize()

  d = -lvec3(atom.xyz)*normal
  t = - (d + (lvec3(axis[0])) * normal)
  nearPoint = rvec3(axis[0]) + normal*t
  distFromNearPoint = (rvec3(atom.xyz)-nearPoint).length()

  return nearPoint + distFromNearPoint * normalizedOffset

from mmtbx_reduce_ext import RotatePointDegreesAroundAxisDir
def _rotateAroundAxis(atom, axis, degrees):
  '''Rotate the atom about the specified axis by the specified number of degrees.
     :param atom: iotbx.pdb.hierarchy.atom or scitbx::vec3<double> or
     scitbx.matrix.rec(xyz, (3,1)) to be moved.
     :param axis: flex array with two scitbx::vec3<double> points, the first
     of which is the origin in space around which to rotate and the second is
     a vector pointing from the origin in the direction of the axis of rotation.
     Positive rotations will be right-handed rotations around this axis.
     :param degrees: How much to rotate the atom around the axis.
     Positive rotation is right-handed around the axis.
     :returns the new location for the atom.
  '''
  # If we have an atom, get its position.  Otherwise, we have a position passed
  # in.
  try:
    pos = atom.xyz
  except Exception:
    pos = rvec3(atom)

  # The axis of rotation for this function is specified as the two ends of the axis.
  # The axis passed in has the point around which to rotate and the direction vector
  # from the origin, so we need to add those together to get the other end of the axis.
  # (This is done in the C++ code by the RotatePointDegreesAroundAxisDir function.)
  return lvec3(RotatePointDegreesAroundAxisDir(axis[0], axis[1], pos, degrees))
  #return lvec3(scitbx.matrix.rotate_point_around_axis(
  #    axis_point_1 = axis[0], axis_point_2 = rvec3(axis[0]) + rvec3(axis[1]),
  #    point = pos, angle = degrees, deg = True))

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
  normal = (rvec3(pivotAtom.xyz) - rvec3(hingeAtom.xyz)).normalize()
  axis = flex.vec3_double([hingeAtom.xyz, normal])
  for i in range(hingeIndex):
    movable[i] = _rotateAroundAxis(movable[i], axis, 180)

  # B) Hinge the movable atoms around the hinge.
  # Solve for the planes containing the old and new atoms and then the vector
  # that these planes intersect at (cross product).

  cToO = lvec3(first.xyz) - lvec3(hingeAtom.xyz)
  nToO = rvec3(second.xyz) - rvec3(hingeAtom.xyz)
  oldNormal = (scitbx.matrix.cross_product_matrix(cToO) * nToO).normalize()

  # Flip the order here because we've rotated 180 degrees
  newNToO = lvec3(movable[secondDockIndex]) - lvec3(movable[hingeIndex])
  newCToO = rvec3(movable[firstDockIndex]) - rvec3(movable[hingeIndex])
  newNormal = (scitbx.matrix.cross_product_matrix(newNToO) * newCToO).normalize()

  # If we don't need to rotate, we'll get a zero-length vector
  hinge = scitbx.matrix.cross_product_matrix(lvec3(oldNormal))*rvec3(newNormal)
  if hinge.length() > 0:
    hinge = hinge.normalize()
    axis = flex.vec3_double([hingeAtom.xyz, hinge])
    degrees = 180/math.pi * math.acos((lvec3(oldNormal)*rvec3(newNormal))[0])
    for i in range(hingeIndex):
      # Rotate in the opposite direction, taking the new back to the old.
      movable[i] = _rotateAroundAxis(rvec3(movable[i]), axis, -degrees)

  # C) Rotate the atoms around the alphaCarbon
  #  1) firstDockIndex to alphaCarbon<-->original secondDockIndex line
  aToNewO = lvec3(movable[firstDockIndex]) - lvec3(alphaCarbon.xyz)
  aToOldN = rvec3(second.xyz) - rvec3(alphaCarbon.xyz)
  # If we don't need to rotate, we'll get a zero-length vector
  normal = scitbx.matrix.cross_product_matrix(aToNewO) * aToOldN
  if normal.length() > 0:
    normal = normal.normalize()
    axis = flex.vec3_double([alphaCarbon.xyz, normal])
    degrees = 180/math.pi * math.acos((lvec3(aToOldN.normalize())*rvec3(aToNewO.normalize()))[0])
    for i in range(len(movable)):
      movable[i] = _rotateAroundAxis(rvec3(movable[i]), axis, degrees)

  #  2) firstDockIndex to the proper plane.
  sites = flex.vec3_double([ movable[secondDockIndex], alphaCarbon.xyz, second.xyz, first.xyz ])
  degrees = scitbx.math.dihedral_angle(sites=sites, deg=True)
  hinge = rvec3(alphaCarbon.xyz) - rvec3(second.xyz)
  if hinge.length() > 0:
    hinge = hinge.normalize()
    axis = flex.vec3_double([alphaCarbon.xyz, hinge])
    for i in range(len(movable)):
      # Rotate in the opposite direction, taking the new back to the old.
      movable[i] = _rotateAroundAxis(rvec3(movable[i]), axis, -degrees)

  return movable

##################################################################################
# Test function and its helpers to verify that all Movers behave properly.

def _StableUnderAtomMotion(mover, atom):
  cI = mover.CoarsePositions().positions
  fI = mover.FinePositions(0).positions
  atom.xyz = ( atom.xyz[0]+10, atom.xyz[1], atom.xyz[2] )
  cC = mover.CoarsePositions().positions
  fC = mover.FinePositions(0).positions
  for i,p in enumerate(cI):
    for j, atom in enumerate(p):
      if cI[i][j] != cC[i][j]:
        return False
  for i,p in enumerate(fI):
    for j, atom in enumerate(p):
      if fI[i][j] != fC[i][j]:
        return False
  return True

def Test():
  """Test function for all classes provided above.
  :returns Empty string on success, string describing the problem on failure.
  :returns Empty string on success, string describing the problem on failure.
  """

  # Test the behavior of probe.ExtraAtomInfo to ensure that when we make a copy
  # we can independently modify it.
  info = probe.ExtraAtomInfo()
  info2 = probe.ExtraAtomInfo(info)
  info2.isAcceptor = True
  if info == info2:
    return "Movers.Test() Got unexpected behavior when modifying copy-constructed probe.ExtraAtomInfo() "

  # Test the MoverNull class.
  try:
    atom = pdb.hierarchy.atom()
    atom.name = "C"
    atoms = [atom]
    extras = [probe.ExtraAtomInfo()]
    extrasMap = probe.ExtraAtomInfoMap(atoms, extras)
    m = MoverNull(atom, extrasMap)
    coarse = m.CoarsePositions()
    if len(coarse.atoms) != 1:
      return "Movers.Test() MoverNull: Expected 1 atom for CoarsePositions, got "+str(len(coarse.atoms))
    fine = m.FinePositions(0)
    if len(fine.atoms) != 0:
      return "Movers.Test() MoverNull: Expected 0 atoms for FinePositions, got "+str(len(fine.atoms))
    fixUp = m.FixUp(0)
    if len(fixUp.atoms) != 0:
      return "Movers.Test() MoverNull: Expected 0 atoms for FixUp, got "+str(len(fixUp.atoms))
    if m.PoseDescription(0,0, False) != "Original location . .":
      return "Movers.Test() MoverNull: Unexpected results for PoseDescription, got "+m.PoseDescription(0,0, False)
    if m.PoseDescription(1,0, False) != "Unrecognized state . .":
      return "Movers.Test() MoverNull: Unexpected results for PoseDescription, got "+m.PoseDescription(1,0, False)
    if m.PoseDescription(0,1, False) != "Unrecognized state . .":
      return "Movers.Test() MoverNull: Unexpected results for PoseDescription, got "+m.PoseDescription(0,1, False)

    # Verify that the coarse and fine results don't change when the atom position is moved after
    # the Mover has been constructed.
    if not _StableUnderAtomMotion(m, atom):
      return "Movers.Test() MoverNull: Positions not stable when atom moved"

  except Exception as e:
    return "Movers.Test() MoverNull: Exception during test of MoverNull: "+str(e)+"\n"+traceback.format_exc()

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
      atom.name = "H"
      a.xyz = coords[i]
      atoms.append(a)
    axis = flex.vec3_double([ [ 0.0, 0.0,-2.0], [ 0.0, 0.0, 1.0] ])
    def prefFunc(x):
      return 1.0
    rot = _MoverRotator(atoms,axis, 0, 0, 180, 90.0, True, preferenceFunction = prefFunc)

    # See if the results of each of the functions are what we expect in terms of sizes and locations
    # of atoms and preferences.  We'll use the default fine step size.
    # The first coarse rotation should be by -90 degrees, moving the first atom to (0, -1, 1)
    coarse = rot.CoarsePositions()
    if len(coarse.atoms) != 3:
      return "Movers.Test() _MoverRotator basic: Expected 3 atoms for CoarsePositions, got "+str(len(coarse.atoms))
    atom0pos1 = coarse.positions[1][0]
    if (lvec3(atom0pos1) - lvec3([0,-1,1])).length() > 1e-5:
      return "Movers.Test() _MoverRotator basic: Expected location = (0,-1,1), got "+str(atom0pos1)

    # The first fine rotation (index 0) around the second coarse index (index 1) should be to -91 degrees,
    # moving the first atom to the appropriately rotated location around the Z axis
    rad = -91 / 180 * math.pi
    x = math.cos(rad)
    y = math.sin(rad)
    z = 1
    fine = rot.FinePositions(1)
    atom0pos1 = fine.positions[0][0]
    if (lvec3(atom0pos1) - lvec3([x,y,z])).length() > 1e-5:
      return "Movers.Test() _MoverRotator basic: Expected fine location = "+str([x,y,z])+", got "+str(atom0pos1)

    # The preference function should always return 1.
    for p in fine.preferenceEnergies:
      if p != 1:
        return "Movers.Test() _MoverRotator basic: Expected preference energy = 1, got "+str(p)

    # Test different preference scale value.
    rot = _MoverRotator(atoms,axis, 0, 0, 180, 90.0, True, preferenceFunction = prefFunc, preferredOrientationScale = 2)
    coarse = rot.CoarsePositions()
    for p in coarse.preferenceEnergies:
      if p != 2:
        return "Movers.Test() _MoverRotator Scaled preference: Expected preference energy = 2, got "+str(p)

    # Test None preference function and a sinusoidal one.
    rot = _MoverRotator(atoms,axis, 0, 0, 180, 90.0, True, preferenceFunction = None)
    coarse = rot.CoarsePositions()
    for p in coarse.preferenceEnergies:
      if p != 0:
        return "Movers.Test() _MoverRotator None preference function: Expected preference energy = 0, got "+str(p)
    def prefFunc2(x):
      return math.cos(x * math.pi / 180)
    rot = _MoverRotator(atoms,axis, 0, 0, 180, 180, True, preferenceFunction = prefFunc2)
    coarse = rot.CoarsePositions()
    expected = [1, -1]
    if len(coarse.preferenceEnergies) != len(expected):
      return "Movers.Test() _MoverRotator Sinusoidal preference function: Unexpected preference length:  "+str(len(coarse.preferenceEnergies))
    for i,p  in enumerate(coarse.preferenceEnergies):
      val = expected[i]
      if p != val:
        return "Movers.Test() _MoverRotator Sinusoidal preference function: Expected preference energy = "+str(val)+", got "+str(p)

    # Test coarseStepDegrees default behavior.
    rot = _MoverRotator(atoms,axis, 0, 0, 180)
    coarse = rot.CoarsePositions()
    if len(coarse.positions) != 24:
      return "Movers.Test() _MoverRotator Default coarse step: Expected 24, got "+str(len(coarse.positions))

    # Test doFineRotations = False and 180 degree coarseStepDegrees.
    rot = _MoverRotator(atoms,axis, 0, 0, 180, 180, False)
    coarse = rot.CoarsePositions()
    if len(coarse.positions) != 2:
      return "Movers.Test() _MoverRotator 180 coarse steps: Expected 2, got "+str(len(coarse.positions))
    fine = rot.FinePositions(0)
    if len(fine.positions) != 0:
      return "Movers.Test() _MoverRotator 180 coarse steps: Expected 0, got "+str(len(fine.positions))

    # Test fineStepDegrees setting.
    rot = _MoverRotator(atoms,axis, 0, 0, 180, fineStepDegrees = 2)
    fine = rot.FinePositions(0)
    # +/- 7.5 degrees in 2-degree steps, but we wouldn't do the +7.5 because it will be handled by the next
    # rotation up. So we get +/- 2, 4, and 6
    if len(fine.positions) != 6:
      return "Movers.Test() _MoverRotator setting fine step: Expected 6, got "+str(len(fine.positions))

    # Test the PoseDescription
    if rot.PoseDescription(1,1, False) != "Angle -13.0 deg .":
      return "Movers.Test() _MoverRotator: Unexpected results for PoseDescription, got "+rot.PoseDescription(1,1, False)

    # Test setting an offset and counterbalancing dihedral.
    rot = _MoverRotator(atoms,axis, -10, 10, 180, fineStepDegrees = 2)
    if rot.PoseDescription(1,1, False) != "Angle -3.0 deg .":
      return "Movers.Test() _MoverRotator: Unexpected results for offset, got "+rot.PoseDescription(1,1, False)

    # Verify that the coarse and fine results don't change when the atom position is moved after
    # the Mover has been constructed.
    if not _StableUnderAtomMotion(rot, atoms[0]):
      return "Movers.Test() _MoverRotator: Positions not stable when atom moved"

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
    h.name = "H"
    h.xyz = [ 1.0, 1.0, 1.0 ]

    n = pdb.hierarchy.atom()
    n.name = "C"
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.name = "C"
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.name = "C1"
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.name = "C2"
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

    # Prepare our extraAtomInfoMap
    atoms = c.atoms()
    extras = []
    for a in atoms:
      extras.append(probe.ExtraAtomInfo())
    extrasMap = probe.ExtraAtomInfoMap(atoms, extras)

    # Add a non-bonded potential acceptor atom at 13 degrees rotation towards the Y axis from
    # the X axis.
    acc = pdb.hierarchy.atom()
    acc.name = "C"
    acc.xyz = [ math.cos(13*math.pi/180), math.sin(13*math.pi/180), 1.0 ]

    # Construct a stand-in riding-coefficients-producing hParams object that will provide an
    # appropriate dihedral atom for the hydrogen to use.
    class Item:
      def __init__(self, paramN, paramA2):
        self.n = paramN
        self.a2 = paramA2
    item = Item(0, f1.i_seq)
    hParams = {}
    hParams[h.i_seq] = item

    mover = MoverSingleHydrogenRotator(h, bondedNeighborLists, extrasMap, hParams, [acc])

    # Check for hydrogen rotated into -X plane at a distance of sqrt(2) from Z axis.
    # It should have been rotated 180 degrees from f1 because f1 is the conventional branch based on its name.
    if h.xyz[2] != 1 or abs(-h.xyz[0]-math.sqrt(2)) > 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator pair: bad H placement"

    ''' Preference function has been removed from this class
    # Check fitness function preferring 0 and 180 rotations
    zero = mover._preferenceFunction(0)
    ninety = mover._preferenceFunction(90)
    oneEighty = mover._preferenceFunction(180)
    if abs(zero - oneEighty) > 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator pair: bad preference function"
    if zero - ninety < 1e-5:
      return "Movers.Test() MoverSingleHydrogenRotator pair: bad preference function"
    '''

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

    # Verify that the coarse and fine results don't change when the atom position is moved after
    # the Mover has been constructed.
    if not _StableUnderAtomMotion(mover, h):
      return "Movers.Test() MoverSingleHydrogenRotator: Positions not stable when atom moved"

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
    h.name = "H"
    h.xyz = [ 1.0, 1.0, 1.0 ]

    n = pdb.hierarchy.atom()
    n.name = "C"
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.name = "C"
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.name = "C1"
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.name = "C2"
    f2.xyz = _rotateAroundAxis(f1, axis, -120)

    f3 = pdb.hierarchy.atom()
    f3.name = "C3"
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

    # Prepare our extraAtomInfoMap
    atoms = c.atoms()
    extras = []
    for a in atoms:
      extras.append(probe.ExtraAtomInfo())
    extrasMap = probe.ExtraAtomInfoMap(atoms, extras)

    # Construct a stand-in riding-coefficients-producing hParams object that will provide an
    # appropriate dihedral atom for the hydrogen to use.
    class Item:
      def __init__(self, paramN, paramA2):
        self.n = paramN
        self.a2 = paramA2
    item = Item(0, f1.i_seq)
    hParams = {}
    hParams[h.i_seq] = item

    mover = MoverSingleHydrogenRotator(h, bondedNeighborLists, extrasMap, hParams)

    # Check for a hydrogen on the -X axis at a distance of sqrt(2) from Z axis,
    # it should have picked f1 as the conventional friend to be opposite to based on its name.
    angle = 0
    if not (h.xyz[2] == 1 and h.xyz[0]+math.sqrt(2) < 1e-5):
      return "Movers.Test() MoverSingleHydrogenRotator triple: bad H placement"

    ''' Preference function has been removed from this class
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
    '''

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
    h1.name = "H1"
    h1.xyz = [ 1.0, 1.0, 1.0 ]

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.name = "H2"
    h2.xyz = _rotateAroundAxis(h1, axis, -120)

    h3 = pdb.hierarchy.atom()
    h3.element = "H"
    h3.name = "H3"
    h3.xyz = _rotateAroundAxis(h1, axis, 120)

    n = pdb.hierarchy.atom()
    n.element = "N"
    n.name = "N"
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.name = "C"
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.name = "C1"
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.name = "C2"
    f2.xyz = _rotateAroundAxis(f1, axis, -120)

    f3 = pdb.hierarchy.atom()
    f3.name = "C3"
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

    # Construct a stand-in riding-coefficients-producing hParams object that will provide an
    # appropriate dihedral atom for the appropriate hydrogen to use.
    class Item:
      def __init__(self, paramN, paramA2):
        self.n = paramN
        self.a2 = paramA2
    item = Item(0, f1.i_seq)
    hParams = {}
    hParams[h3.i_seq] = item

    mover = MoverNH3Rotator(n, bondedNeighborLists, hParams)

    # Check for the third hydrogen on the -X axis at a distance of sqrt(2) from Z the axis.
    if not (h3.xyz[2] == 1 and h3.xyz[0]+math.sqrt(2) < 1e-5):
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

    # Verify that the coarse and fine results don't change when the atom position is moved after
    # the Mover has been constructed.
    if not _StableUnderAtomMotion(mover, h1):
      return "Movers.Test() MoverNH3Rotator: Positions not stable when atom moved"

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
    h1.name = "H1"
    h1.xyz = [ 1.0, 1.0, 1.0 ]

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.name = "H2"
    h2.xyz = _rotateAroundAxis(h1, axis, -120)

    h3 = pdb.hierarchy.atom()
    h3.element = "H"
    h3.name = "H3"
    h3.xyz = _rotateAroundAxis(h1, axis, 120)

    n = pdb.hierarchy.atom()
    n.element = "C"
    n.name = "C"
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.name = "C"
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.name = "C1"
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.name = "C2"
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

    # Construct a stand-in riding-coefficients-producing hParams object that will provide an
    # appropriate dihedral atom for the appropriate hydrogen to use.
    class Item:
      def __init__(self, paramN, paramA2):
        self.n = paramN
        self.a2 = paramA2
    item = Item(0, f1.i_seq)
    hParams = {}
    hParams[h3.i_seq] = item

    mover = MoverAromaticMethylRotator(n, bondedNeighborLists, hParams)

    # Check for the third hydrogen on the +Y axis at a distance of sqrt(2) from Z the axis.
    if not (h3.xyz[2] == 1 and abs(h3.xyz[1]-math.sqrt(2)) < 1e-5):
      return "Movers.Test() MoverAromaticMethylRotator basic: bad H placement"

    # Check that we get two coarse and no fine orientations
    coarse = mover.CoarsePositions().positions
    if len(coarse) != 2:
      return "Movers.Test() MoverAromaticMethylRotator basic: bad coarse count: "+str(len(coarse))
    fine = mover.FinePositions(0).positions
    if len(fine) != 0:
      return "Movers.Test() MoverAromaticMethylRotator basic: bad fine count: "+str(len(fine))

    # Verify that the coarse and fine results don't change when the atom position is moved after
    # the Mover has been constructed.
    if not _StableUnderAtomMotion(mover, h1):
      return "Movers.Test() MoverNH3Rotator: Positions not stable when atom moved"

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
    h1.name = "H1"
    h1.xyz = [ 1.0, 1.0, 1.0 ]

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.name = "H2"
    h2.xyz = _rotateAroundAxis(h1, axis, -120)

    h3 = pdb.hierarchy.atom()
    h3.element = "H"
    h3.name = "H3"
    h3.xyz = _rotateAroundAxis(h1, axis, 120)

    n = pdb.hierarchy.atom()
    n.element = "C"
    n.name = "C"
    n.xyz = [ 0.0, 0.0, 0.0 ]

    p = pdb.hierarchy.atom()
    p.name = "C"
    p.xyz = [ 0.0, 0.0,-1.0 ]

    f1 = pdb.hierarchy.atom()
    f1.name = "C1"
    f1.xyz = [ 1.0, 0.0,-2.0 ]

    f2 = pdb.hierarchy.atom()
    f2.name = "C2"
    f2.xyz = _rotateAroundAxis(f1, axis, -120)

    f3 = pdb.hierarchy.atom()
    f3.name = "C3"
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

    # Construct a stand-in riding-coefficients-producing hParams object that will provide an
    # appropriate dihedral atom for the appropriate hydrogen to use.
    class Item:
      def __init__(self, paramN, paramA2):
        self.n = paramN
        self.a2 = paramA2
    item = Item(0, f1.i_seq)
    hParams = {}
    hParams[h3.i_seq] = item

    mover = MoverTetrahedralMethylRotator(n, bondedNeighborLists, hParams)

    # Check for the third hydrogen on the -X axis at a distance of sqrt(2) from Z the axis.
    if not (h3.xyz[2] == 1 and abs(-h3.xyz[0]-math.sqrt(2)) < 1e-5):
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

    # Verify that the coarse and fine results don't change when the atom position is moved after
    # the Mover has been constructed.
    if not _StableUnderAtomMotion(mover, h1):
      return "Movers.Test() MoverTetrahedralMethylRotator: Positions not stable when atom moved"

  except Exception as e:
    return "Movers.Test() MoverTetrahedralMethylRotator basic: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  # Test the MoverAmideFlip class with no linker (similar to Asn).
  try:
    # Test behavior with offsets for each atom so that we test the generic case.
    # Construct a MoverAmideFlip that has the N, H's and Oxygen located (non-physically)
    # slightly in the +Y direction out of the X-Z plane, with the Hydrogens 120 around the
    # same offset axis.
    # They are bonded to a Carbon and friend (pivot) and alpha carbon that are on the Z axis.
    # Non-physically, atoms are 1 unit apart.
    p = pdb.hierarchy.atom()
    p.name = "C"
    p.xyz = [ 0.0, 0.0, 0.0 ]

    f = pdb.hierarchy.atom()
    f.element = "C"
    f.name = "C"
    f.xyz = [ 0.0, 0.0,-1.0 ]

    fh1 = pdb.hierarchy.atom()
    fh1.element = "H"
    fh1.name = "H"
    fh1.xyz = [ 1.0, 0.0,-1.0 ]

    fh2 = pdb.hierarchy.atom()
    fh2.element = "H"
    fh2.name = "H"
    fh2.xyz = [-1.0, 0.0,-1.0 ]

    ca = pdb.hierarchy.atom()
    ca.element = "C"
    ca.name = "CA"
    ca.xyz = [ 0.0, 0.0,-2.0 ]

    # Nitrogen and Oxygen are +/-120 degrees from carbon-carbon bond
    axis = flex.vec3_double([ [0,0,0], [0,1,0] ])
    n = pdb.hierarchy.atom()
    n.element = "N"
    n.name = "N"
    n.xyz = _rotateAroundAxis(f, axis,-120) + lvec3([0,0.01,0]) + lvec3([ 0.002, 0.003,-0.004])

    o = pdb.hierarchy.atom()
    o.element = "O"
    o.name = "O"
    o.xyz = _rotateAroundAxis(f, axis, 120) + lvec3([0,0.01,0]) + lvec3([-0.003, 0.002, 0.003])

    # Hydrogens are +/-120 degrees from nitrogen-carbon bond
    axis = flex.vec3_double([ n.xyz, [0,1,0] ])
    h1 = pdb.hierarchy.atom()
    h1.element = "H"
    h1.name = "H1"
    h1.xyz = _rotateAroundAxis(p, axis,-120) + lvec3([0,0.01,0]) + lvec3([-0.008, 0.001, 0.008])

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.name = "H2"
    h2.xyz = _rotateAroundAxis(p, axis, 120) + lvec3([0,0.01,0]) + lvec3([ 0.007,-0.001, 0.007])

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

    mover = MoverAmideFlip(n, ca.name, bondedNeighborLists, 0.5)

    # Ensure that the hydrogens have been rotated to have a 0 and 180 dihedral with
    # the Oxygen atom.
    sites = [ h1.xyz, n.xyz, p.xyz, o.xyz ]
    degrees = scitbx.math.dihedral_angle(sites=sites, deg=True)
    while degrees > 180:
        degrees -= 360
    while degrees <= -180:
        degrees += 360
    if abs(degrees) > 1e-5 and abs(180 - degrees) > 1e-5:
      return "Movers.Test() MoverAmideFlip: h1 dihedral not 0 or 180: {:.2f}".format(degrees)
    sites = [ h2.xyz, n.xyz, p.xyz, o.xyz ]
    degrees = scitbx.math.dihedral_angle(sites=sites, deg=True)
    while degrees > 180:
        degrees -= 360
    while degrees <= -180:
        degrees += 360
    if abs(degrees) > 1e-5 and abs(180 - degrees) > 1e-5:
      return "Movers.Test() MoverAmideFlip: h2 dihedral not 0 or 180: {:.2f}".format(degrees)

    # Ensure that the coarse-flip results meet the expections:
    # 1) N and O are flipped in position
    # 2) H remain at the same distance from the new N.

    coarse = mover.CoarsePositions()
    if coarse.preferenceEnergies[0] <= coarse.preferenceEnergies[1]:
      return "Movers.Test() MoverAmideFlip: Original orientation not preferred"
    if len(coarse.positions) != 2:
      return "Movers.Test() MoverAmideFlip basic: Did not find two locations: "+str(len(coarse.positions))
    newPos = coarse.positions[1]
    dist = (lvec3(newPos[2]) - lvec3(o.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverAmideFlip basic: Nitrogen moved incorrectly: "+str(dist)
    dist = (lvec3(newPos[3]) - lvec3(n.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverAmideFlip basic: Oxygen moved incorrectly: "+str(dist)

    dHydrogen = (lvec3(newPos[0]) - lvec3(newPos[2])).length()
    oldDHydrogen = (lvec3(h1.xyz)-lvec3(n.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverAmideFlip basic: Bad coarse hydrogen1 motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (lvec3(newPos[1]) - lvec3(newPos[2])).length()
    oldDHydrogen = (lvec3(h2.xyz)-lvec3(n.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverAmideFlip basic: Bad coarse hydrogen2 motion: "+str(dHydrogen-oldDHydrogen)

    # Ensure that the Fixup results meet the specifications:
    # 1) New Oxygen on the line from the alpha carbon to the old Nitrogen
    # 2) New plane of Oxygen, Nitrogen, Alpha Carbon matches old plane, but flipped
    # 3) Carbons and pivot Hydrogens move slightly due to rigid-body motion

    fixedUp = mover.FixUp(1).positions
    newODir = (fixedUp[3] - lvec3(f.xyz)).normalize()
    oldNDir = (rvec3(n.xyz) - rvec3(f.xyz)).normalize()
    if (newODir * oldNDir)[0] < 0.9999:
      return "Movers.Test() MoverAmideFlip basic: Bad oxygen alignment: "+str((newODir * oldNDir)[0])

    newNDir = (fixedUp[2] - lvec3(f.xyz)).normalize()
    oldODir = (rvec3(o.xyz) - rvec3(f.xyz)).normalize()
    newNormal = (scitbx.matrix.cross_product_matrix(lvec3(newNDir)) * rvec3(newODir)).normalize()
    oldNormal = (scitbx.matrix.cross_product_matrix(lvec3(oldNDir)) * rvec3(oldODir)).normalize()
    dot = (lvec3(newNormal) * rvec3(oldNormal))[0]
    if dot > -0.99999:
      return "Movers.Test() MoverAmideFlip basic: Bad plane alignment: "+str(dot)

    dCarbon = (fixedUp[4] - lvec3(p.xyz)).length()
    if dCarbon < 0.001 or dCarbon > 0.1:
      return "Movers.Test() MoverAmideFlip basic: Bad hinge motion: "+str(dCarbon)

    dCarbon = (fixedUp[5] - lvec3(f.xyz)).length()
    if dCarbon < 0.0005 or dCarbon > 0.1:
      return "Movers.Test() MoverAmideFlip basic: Bad pivot motion: "+str(dCarbon)

    dHydrogen = (fixedUp[6] - lvec3(fh1.xyz)).length()
    if dHydrogen < 0.0005 or dHydrogen > 0.1:
      return "Movers.Test() MoverAmideFlip basic: Bad pivot hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixedUp[7] - lvec3(fh2.xyz)).length()
    if dHydrogen < 0.0005 or dHydrogen > 0.1:
      return "Movers.Test() MoverAmideFlip basic: Bad pivot hydrogen motion: "+str(dHydrogen)

    # Test the PoseDescription
    if mover.PoseDescription(1, 0, True) != "Flipped . . AnglesAdjusted":
      return "Movers.Test() MoverAmideFlip basic: Unexpected results for PoseDescription, got "+mover.PoseDescription(1,1, True)

    # Verify that the coarse and fine results don't change when the atom position is moved after
    # the Mover has been constructed.
    if not _StableUnderAtomMotion(mover, n):
      return "Movers.Test() MoverAmideFlip: Positions not stable when atom moved"

  except Exception as e:
    return "Movers.Test() MoverAmideFlip basic: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  # Test the MoverAmideFlip class with a linker (similar to Gln).
  try:
    # Test behavior with offsets for each atom so that we test the generic case.
    # Construct a MoverAmideFlip that has the N, H's and Oxygen located (non-physically)
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
    n.xyz = _rotateAroundAxis(f, axis,-120) + lvec3([0,0.01,0]) + lvec3([ 0.002, 0.003,-0.004])

    o = pdb.hierarchy.atom()
    o.element = "O"
    o.xyz = _rotateAroundAxis(f, axis, 120) + lvec3([0,0.01,0]) + lvec3([-0.003, 0.002, 0.003])

    # Hydrogens are +/-120 degrees from nitrogen-carbon bond
    axis = flex.vec3_double([ n.xyz, [0,1,0] ])
    h1 = pdb.hierarchy.atom()
    h1.element = "H"
    h1.xyz = _rotateAroundAxis(p, axis,-120) + lvec3([0,0.01,0]) + lvec3([-0.008, 0.001, 0.008])

    h2 = pdb.hierarchy.atom()
    h2.element = "H"
    h2.xyz = _rotateAroundAxis(p, axis, 120) + lvec3([0,0.01,0]) + lvec3([ 0.007,-0.001, 0.007])

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

    mover = MoverAmideFlip(n, ca.name, bondedNeighborLists, 0.5)
    fixedUp = mover.FixUp(1).positions

    # Ensure that the results meet the specifications:
    # 1) New Oxygen on the line from the alpha carbon to the old Nitrogen
    # 2) New plane of Oxygen, Nitrogen, Alpha Carbon matches old plane, but flipped
    # 3) Pivot and linker Carbons and Hydrogens move slightly due to rigid-body motion

    newODir = (fixedUp[3] - lvec3(ca.xyz)).normalize()
    oldNDir = (rvec3(n.xyz) - rvec3(ca.xyz)).normalize()
    if (newODir * oldNDir)[0] < 0.9999:
      return "Movers.Test() MoverAmideFlip linked: Bad oxygen alignment: "+str((newODir * oldNDir)[0])

    newNDir = (fixedUp[2] - lvec3(ca.xyz)).normalize()
    oldODir = (rvec3(o.xyz) - rvec3(ca.xyz)).normalize()
    newNormal = (scitbx.matrix.cross_product_matrix(lvec3(newNDir)) * rvec3(newODir)).normalize()
    oldNormal = (scitbx.matrix.cross_product_matrix(lvec3(oldNDir)) * rvec3(oldODir)).normalize()
    dot = (lvec3(newNormal) * rvec3(oldNormal))[0]
    if dot > -0.99999:
      return "Movers.Test() MoverAmideFlip linked: Bad plane alignment: "+str(dot)

    dCarbon = (fixedUp[4] - lvec3(p.xyz)).length()
    if dCarbon < 0.0006 or dCarbon > 0.1:
      return "Movers.Test() MoverAmideFlip linked: Bad hinge motion: "+str(dCarbon)

    dCarbon = (fixedUp[5] - lvec3(f.xyz)).length()
    if dCarbon < 0.0004 or dCarbon > 0.1:
      return "Movers.Test() MoverAmideFlip linked: Bad pivot motion: "+str(dCarbon)

    dCarbon = (fixedUp[6] - lvec3(ln.xyz)).length()
    if dCarbon < 0.0002 or dCarbon > 0.1:
      return "Movers.Test() MoverAmideFlip linked: Bad linker motion: "+str(dCarbon)

    # Hydrogens come after all linkers
    dHydrogen = (fixedUp[7] - lvec3(fh1.xyz)).length()
    if dHydrogen < 0.0004 or dHydrogen > 0.1:
      return "Movers.Test() MoverAmideFlip linked: Bad pivot hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixedUp[8] - lvec3(fh2.xyz)).length()
    if dHydrogen < 0.0004 or dHydrogen > 0.1:
      return "Movers.Test() MoverAmideFlip linked: Bad pivot hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixedUp[9] - lvec3(lnh1.xyz)).length()
    if dHydrogen < 0.0002 or dHydrogen > 0.1:
      return "Movers.Test() MoverAmideFlip linked: Bad linker hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixedUp[10] - lvec3(lnh2.xyz)).length()
    if dHydrogen < 0.0002 or dHydrogen > 0.1:
      return "Movers.Test() MoverAmideFlip linked: Bad linker hydrogen motion: "+str(dHydrogen)

    # Ensure that the Hydrogens moved along with their parent in the flip.
    dHydrogen = (fixedUp[0] - fixedUp[2]).length()
    oldDHydrogen = (lvec3(h1.xyz)-lvec3(n.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverAmideFlip linked: Bad nitrogen-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (fixedUp[1] - fixedUp[2]).length()
    oldDHydrogen = (lvec3(h2.xyz)-lvec3(n.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverAmideFlip linked: Bad nitrogen-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    # Test the PoseDescription
    if mover.PoseDescription(0, 0, True) != "Unflipped . . .":
      return "Movers.Test() MoverAmideFlip: Unexpected results for PoseDescription 0T, got "+mover.PoseDescription(0,0, True)
    if mover.PoseDescription(1, 0, False) != "Flipped . . AnglesNotAdjusted":
      return "Movers.Test() MoverAmideFlip: Unexpected results for PoseDescription 1F, got "+mover.PoseDescription(1,0, False)
    if mover.PoseDescription(1, 0, True) != "Flipped . . AnglesAdjusted":
      return "Movers.Test() MoverAmideFlip: Unexpected results for PoseDescription 1T, got "+mover.PoseDescription(1,0, True)

  except Exception as e:
    return "Movers.Test() MoverAmideFlip linked: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  # Test the MoverHisFlip class.
  try:
    # Test behavior with offsets for each atom so that we test the generic case.
    # Construct a MoverHisFlip that has the ring atoms other than CG located (non-physically)
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
    nd1.xyz = lvec3([ 1.0, 0.0, 1.0]) + lvec3([0,0.01,0]) + lvec3([-0.003, 0.002, 0.003])

    nd1h = pdb.hierarchy.atom()
    nd1h.element = "H"
    nd1h.xyz = [ 2.0, 0.01, 1.0 ]

    ce1 = pdb.hierarchy.atom()
    ce1.element = "C"
    ce1.xyz = lvec3([ 1.0, 0.0, 2.1]) + lvec3([0,0.01,0]) + lvec3([-0.008, 0.001, 0.008])

    ce1h = pdb.hierarchy.atom()
    ce1h.element = "H"
    ce1h.xyz = [ 2.0, 0.01, 2.0 ]

    cd2 = pdb.hierarchy.atom()
    cd2.element = "C"
    cd2.xyz = lvec3([-1.0, 0.0, 1.0]) + lvec3([0,0.01,0]) + lvec3([ 0.007,-0.001, 0.007])

    cd2h = pdb.hierarchy.atom()
    cd2h.element = "H"
    cd2h.xyz = [-2.0, 0.01, 1.0 ]

    ne2 = pdb.hierarchy.atom()
    ne2.element = "N"
    ne2.xyz = lvec3([-1.0, 0.0, 2.0]) + lvec3([0,0.01,0]) + lvec3([-0.002, 0.005, 0.004])

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

    # Prepare our bonded-neighbor information
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

    # Prepare our extraAtomInfoMap
    atoms = c.atoms()
    extras = []
    for a in atoms:
      extras.append(probe.ExtraAtomInfo())
    extrasMap = probe.ExtraAtomInfoMap(atoms, extras)

    mover = MoverHisFlip(ne2, bondedNeighborLists, extrasMap, 0.5)

    # Ensure that the coarse-flip results meet the expections (spot check 4th position):
    # 1) N and C atoms are flipped in pairs
    # 2) H remain at the same distance from the new locations.

    coarse = mover.CoarsePositions()
    if len(coarse.positions) != 8:
      return "Movers.Test() MoverHisFlip: Did not find 8 locations: found "+str(len(coarse.positions))
    if coarse.preferenceEnergies[0] <= coarse.preferenceEnergies[4]:
      return "Movers.Test() MoverHisFlip: Original orientation not preferred"
    if coarse.preferenceEnergies[0] >= coarse.preferenceEnergies[1]:
      return "Movers.Test() MoverHisFlip: Hydrogen removal not preferred"
    if coarse.preferenceEnergies[3] >= coarse.preferenceEnergies[0]:
      return "Movers.Test() MoverHisFlip: Both Hydrogen removal preferred"
    newPos = coarse.positions[4]
    dist = (lvec3(newPos[0]) - lvec3(ce1.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverHisFlip: NE2 moved incorrectly: "+str(dist)
    dist = (lvec3(newPos[2]) - lvec3(ne2.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverHisFlip: CE1 moved incorrectly: "+str(dist)
    dist = (lvec3(newPos[4]) - lvec3(cd2.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverHisFlip: ND1 moved incorrectly: "+str(dist)
    dist = (lvec3(newPos[6]) - lvec3(nd1.xyz)).length()
    if dist > 0.01:
      return "Movers.Test() MoverHisFlip: CD2 moved incorrectly: "+str(dist)

    dHydrogen = (lvec3(newPos[0]) - lvec3(newPos[1])).length()
    oldDHydrogen = (lvec3(ne2h.xyz)-lvec3(ne2.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHisFlip: Bad coarse NE2 hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (lvec3(newPos[2]) - lvec3(newPos[3])).length()
    oldDHydrogen = (lvec3(ce1h.xyz)-lvec3(ce1.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHisFlip: Bad coarse CE1 hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (lvec3(newPos[4]) - lvec3(newPos[5])).length()
    oldDHydrogen = (lvec3(nd1h.xyz)-lvec3(nd1.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHisFlip: Bad coarse ND1 hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (lvec3(newPos[6]) - lvec3(newPos[7])).length()
    oldDHydrogen = (lvec3(cd2h.xyz)-lvec3(cd2.xyz)).length()
    if abs(dHydrogen - oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHisFlip: Bad coarse ND1 hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    # Ensure that the FixUp results meet the specifications (spot check 4th position):
    # 1) New CE1 on the line from the alpha carbon to the old NE2
    # 2) The line from the new NE2 to the alpha carbon lies on the
    #    line from the original CE1 to the alpha carbon
    # 3) New plane of CE1, NE2, Alpha Carbon matches old plane, but flipped
    # 4) Carbons and pivot Hydrogens move slightly due to rigid-body motion

    fixedUp = mover.FixUp(4).positions
    newCE1Dir = (fixedUp[2] - lvec3(ca.xyz)).normalize()
    oldNE2Dir = (rvec3(ne2.xyz) - rvec3(ca.xyz)).normalize()
    if (newCE1Dir * oldNE2Dir)[0] < 0.9999:
      return "Movers.Test() MoverHisFlip: Bad CE1 alignment: "+str((newCE1Dir * oldNE2Dir)[0])

    newNE2Dir = (fixedUp[0] - lvec3(ca.xyz)).normalize()
    oldCE1Dir = (rvec3(ce1.xyz) - rvec3(ca.xyz)).normalize()
    if (newNE2Dir * oldCE1Dir)[0] < 0.9999:
      return "Movers.Test() MoverHisFlip: Bad NE2 alignment: "+str((newNE2Dir * oldCE1Dir)[0])

    newNormal = (scitbx.matrix.cross_product_matrix(lvec3(newNE2Dir)) * rvec3(newCE1Dir)).normalize()
    oldNormal = (scitbx.matrix.cross_product_matrix(lvec3(oldNE2Dir)) * rvec3(oldCE1Dir)).normalize()
    dot = (lvec3(newNormal) * rvec3(oldNormal))[0]
    if dot > -0.99999:
      return "Movers.Test() MoverHisFlip: Bad plane alignment: "+str(dot)

    dCarbon = (fixedUp[8] - lvec3(p.xyz)).length()
    if dCarbon < 0.001 or dCarbon > 0.1:
      return "Movers.Test() MoverHisFlip: Bad hinge motion: "+str(dCarbon)

    dCarbon = (fixedUp[9] - lvec3(f.xyz)).length()
    if dCarbon < 0.0005 or dCarbon > 0.1:
      return "Movers.Test() MoverHisFlip: Bad pivot motion: "+str(dCarbon)

    dHydrogen = (fixedUp[10] - lvec3(fh1.xyz)).length()
    if dHydrogen < 0.0005 or dHydrogen > 0.1:
      return "Movers.Test() MoverHisFlip: Bad pivot hydrogen motion: "+str(dHydrogen)

    dHydrogen = (fixedUp[11] - lvec3(fh2.xyz)).length()
    if dHydrogen < 0.0005 or dHydrogen > 0.1:
      return "Movers.Test() MoverHisFlip: Bad pivot hydrogen motion: "+str(dHydrogen)

    # Ensure that the Hydrogens moved along with their parent in the flip.
    dHydrogen = (fixedUp[0] - fixedUp[1]).length()
    oldDHydrogen = (lvec3(ne2h.xyz)-lvec3(ne2.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHisFlip: Bad NE2-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (fixedUp[2] - fixedUp[3]).length()
    oldDHydrogen = (lvec3(ce1h.xyz)-lvec3(ce1.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHisFlip: Bad CE1-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (fixedUp[4] - fixedUp[5]).length()
    oldDHydrogen = (lvec3(nd1h.xyz)-lvec3(nd1.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHisFlip: Bad ND1-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    dHydrogen = (fixedUp[6] - fixedUp[7]).length()
    oldDHydrogen = (lvec3(cd2h.xyz)-lvec3(cd2.xyz)).length()
    if abs(dHydrogen-oldDHydrogen) > 0.0001:
      return "Movers.Test() MoverHisFlip: Bad CD2-hydrogen motion: "+str(dHydrogen-oldDHydrogen)

    # Ensure that the extraInfo and deleteMe information matches expectations:
    # 1) First hydrogen removed and first N is an acceptor on % = 1
    # 2) Second hydrogen removed and second N is an acceptor on % = 2
    for i in range(6):
      coarseNE2Accept = coarse.extraInfos[i][0].isAcceptor
      coarseNE2HDelete = coarse.deleteMes[i][1]
      fixedUpNE2Accept = mover.FixUp(i).extraInfos[0].isAcceptor
      fixedUpNE2HDelete = mover.FixUp(i).deleteMes[1]
      coarseND1Accept = coarse.extraInfos[i][4].isAcceptor
      coarseND1HDelete = coarse.deleteMes[i][5]
      fixedUpND1Accept = mover.FixUp(i).extraInfos[4].isAcceptor
      fixedUpND1HDelete = mover.FixUp(i).deleteMes[5]
      if i % 4 == 1 or i % 4 == 3:
        if not coarseNE2Accept:
          return "Movers.Test() MoverHisFlip: No NE2 acceptor, pos "+str(i)
        if not coarseNE2HDelete:
          return "Movers.Test() MoverHisFlip: Missing NE2 hydrygen deletion, pos "+str(i)
        if not fixedUpNE2Accept:
          return "Movers.Test() MoverHisFlip: No fixup NE2 acceptor, pos "+str(i)
        if not fixedUpNE2HDelete:
          return "Movers.Test() MoverHisFlip: Missing fixup NE2 hydrygen deletion, pos "+str(i)
      else:
        if coarseNE2Accept:
          return "Movers.Test() MoverHisFlip: Unexpected NE2 acceptor, pos "+str(i)
        if coarseNE2HDelete:
          return "Movers.Test() MoverHisFlip: Unexpected NE2 hydrygen deletion, pos "+str(i)
        if fixedUpNE2Accept:
          return "Movers.Test() MoverHisFlip: Unexpected fixup NE2 acceptor, pos "+str(i)
        if fixedUpNE2HDelete:
          return "Movers.Test() MoverHisFlip: Unexpected fixup NE2 hydrygen deletion, pos "+str(i)
      if i % 4 == 2 or i % 4 == 3:
        if not coarseND1Accept:
          return "Movers.Test() MoverHisFlip: No ND1 acceptor, pos "+str(i)
        if not coarseND1HDelete:
          return "Movers.Test() MoverHisFlip: Missing ND1 hydrygen deletion, pos "+str(i)
        if not fixedUpND1Accept:
          return "Movers.Test() MoverHisFlip: No fixup ND1 acceptor, pos "+str(i)
        if not fixedUpND1HDelete:
          return "Movers.Test() MoverHisFlip: Missing fixup ND1 hydrygen deletion, pos "+str(i)
      else:
        if coarseND1Accept:
          return "Movers.Test() MoverHisFlip: Unexpected ND1 acceptor, pos "+str(i)
        if coarseND1HDelete:
          return "Movers.Test() MoverHisFlip: Unexpected ND1 hydrygen deletion, pos "+str(i)
        if fixedUpND1Accept:
          return "Movers.Test() MoverHisFlip: fixup Unexpected ND1 acceptor, pos "+str(i)
        if fixedUpND1HDelete:
          return "Movers.Test() MoverHisFlip: fixup Unexpected ND1 hydrygen deletion, pos "+str(i)

    # Test the PoseDescription
    if mover.PoseDescription(1,0, False) != "Unflipped HD1Placed HE2NotPlaced .":
      return "Movers.Test() MoverHisFlip: Unexpected results for PoseDescription 1, got "+mover.PoseDescription(1,0, False)
    if mover.PoseDescription(2,0, False) != "Unflipped HD1NotPlaced HE2Placed .":
      return "Movers.Test() MoverHisFlip: Unexpected results for PoseDescription 2, got "+mover.PoseDescription(2,0, False)
    if mover.PoseDescription(3,0, False) != "Unflipped HD1NotPlaced HE2NotPlaced .":
      return "Movers.Test() MoverHisFlip: Unexpected results for PoseDescription 3, got "+mover.PoseDescription(3,0, False)
    if mover.PoseDescription(4,0, True) != "Flipped HD1Placed HE2Placed AnglesAdjusted":
      return "Movers.Test() MoverHisFlip: Unexpected results for PoseDescription 4, got "+mover.PoseDescription(4,0, True)

    # Verify that the coarse and fine results don't change when the atom position is moved after
    # the Mover has been constructed.
    if not _StableUnderAtomMotion(mover, ne2):
      return "Movers.Test() MoverHisFlip: Positions not stable when atom moved"

    # Try locking down the state to only non-flipped and then only flipped and make sure we get
    # the expected number of states and behavior.
    mover = MoverHisFlip(ne2, bondedNeighborLists, extrasMap, 0.5, 1, True)
    if len(mover.CoarsePositions().positions) != 4:
      return "Movers.Test() MoverHisFlip: Unexpected position count for Unflipped, got "+len(mover.CoarsePositions().positions)
    if 'Unflipped' not in mover.PoseDescription(0, 0, False):
      return "Movers.Test() MoverHisFlip: Unexpected PoseDescription for Unflipped, got "+mover.PoseDescription(0, 0, False)
    if '.' != mover.PoseDescription(0, 0, True).split()[-1]:
      return "Movers.Test() MoverHisFlip: Unexpected angle description Unflipped, got "+mover.PoseDescription(0, 0, False)
    mover = MoverHisFlip(ne2, bondedNeighborLists, extrasMap, 0.5, 2, True)
    if len(mover.CoarsePositions().positions) != 4:
      return "Movers.Test() MoverHisFlip: Unexpected position count for Flipped, got "+len(mover.CoarsePositions().positions)
    if 'Flipped' not in mover.PoseDescription(0, 0, True):
      return "Movers.Test() MoverHisFlip: Unexpected PoseDescription for Flipped, got "+mover.PoseDescription(0, 0, False)
    if 'AnglesAdjusted' != mover.PoseDescription(0, 0, True).split()[-1]:
      return "Movers.Test() MoverHisFlip: Unexpected angle description Flipped, got "+mover.PoseDescription(0, 0, True)
    if 'AnglesNotAdjusted' != mover.PoseDescription(0, 0, False).split()[-1]:
      return "Movers.Test() MoverHisFlip: Unexpected angle description Flipped, got "+mover.PoseDescription(0, 0, False)

  except Exception as e:
    return "Movers.Test() MoverHisFlip: Exception during test: "+str(e)+"\n"+traceback.format_exc()

  return ""

##################################################################################
# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  ret = Test()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
