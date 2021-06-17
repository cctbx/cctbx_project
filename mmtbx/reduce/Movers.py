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

##################################################################################
# This is a set of classes that implement Reduce's "Movers".  These are sets of
# atoms that have more than one potential set of locations.
#  - They may have a-priori preferences for some locations over others.
#  - They may have both coarse locations that they test and a set of fine
#    locations around each coarse location.
#  - They may have final tune-up behaviors once their final locations are selected.
#
# There are two basic types of Movers:
#  - Rotater: One or more Hydrogen atoms that have a set of orientations spinning
#    around a common axis.
#  - Flipper: A structure that has pairs of atoms that have similar densities, such
#    that flipping them across a center axis produces similar density fits.
#    They sometimes have optional hydrogen placements on some of the atoms.
#
# All Movers have the following methods, parameters, and return types:
#  - type PositionReturn: ( flex<atom> atoms, flex<flex<vec3>> positions, flex<float> preferenceEnergies )
#  - PositionReturn CoarsePositions(reduceOptions)
#     The reduceOptions is a Phil option subset.  The relevant options for Movers
#       are: CoarseStepDegrees, FineStepDegrees, PreferredOrientationScale
#     The first return lists all of the hierarchy atoms that move.
#     The second return has a vector that has some number of positions, each of
#       which is a vector the same length as the first return
#       with a vec3 position for each atom.
#     The third return has a vector of the same length as the number of positions in
#       each element of the third return; it indicates the relative favorability of
#       each possible location and should be added to the total energy by an
#       optimizer to break ties in otherwise-equivalent orientations.
#  - PositionReturn FinePositions(coarseIndex, reduceOptions)
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
#  - FixUpReturn FixUp(coarseIndex, reduceOptions)
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
# adjustment, they must call FixUp() and apply any final changes.
#
# The InteractionGraph.py script provides functions for determining which pairs of
# Movers have overlaps between movable atoms.
#

##################################################################################
# Return type from CoarsePosition() and FinePosition() calls.
class PositionReturn:
  def __init__(self, atoms, positions, preferenceEnergies):
    self.atoms = atoms
    self.positions = positions
    self.preferenceEnergies = preferenceEnergies

##################################################################################
# Return type from FixUp() calls.
class FixUpReturn:
  def __init__(self, atoms, newPositions):
    self.atoms = atoms
    self.newPositions = newPositions

##################################################################################
# A trivial Mover that returns a single result atom at a single location.
# Useful as a simple and fast test case for programs that use Movers.
# It also serves as a basic example of how to develop a new Mover.
class MoverNull:
  def __init__(self, atom):
    self._atom = atom
  def CoarsePositions(self, reduceOptions):
    # returns: The original atom at its original position.
    return PositionReturn([ self._atom ],
        [ [ [ self._atom.xyz[0], self._atom.xyz[1], self._atom.xyz[2] ] ] ],
        [ 0.0 ])
  def FinePositions(self, coarseIndex, reduceOptions):
    # returns: No fine positions for any coarse position.
    return PositionReturn([], [], [])
  def FixUp(self, coarseIndex, reduceOptions):
    # No fixups for any coarse index.
    return FixUpReturn([], [])

##################################################################################
# A base class for types of Movers that rotate one or more atoms around a single
# axis.  The derived class must determine the initial position for each atom, and
# a preferred-energies function that maps from an angle that is 0 at the initial
# location to an energy.  The range of coarse rotations is also specified, along
# with a flag telling whether fine rotations are allowed.  The coarse step size
# can also be specified, overriding the value from reduceOptions.CoarseStepDegrees.
class MoverRotater:
  def __init__(self, atoms, axis, coarseRange, coarseStepDegrees = None,
                doFineRotations = True, preferenceFunction = None):
    """Constructs a Rotator, to be called by a derived class or test program but
       not usually user code.
       :param atoms: flex array of atoms that will be rotated.  These should be
       placed in one of the preferred orientations if there are any.
       : param axis: flex array with two scitbx::vec3<double> points, the first
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
       reduceOptions.CoarseStepDegrees entry passed to CoarsePositions().
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
    """
    self._atoms = atoms
    self._axis = axis
    self._coarseRange = coarseRange
    self._doFineRotations = doFineRotations
    self._preferenceFunction = preferenceFunction

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
    curStep = CoarseStepDegrees
    while curStep <= self._coarseRange:
      # The interval is closed on the left, so we place
      self._coarseAngles.append(-curStep)
      # We only place on the right if we're strictly less than because
      # the interval is open on the right.
      if curStep < self._coarseRange:
        self._coarseAngles.append(curStep)

  def CoarsePositions(self, reduceOptions):

    # Produce a list of lists, where each entry has the positions of each atom
    # rotated about the specified axis by the specified amount.  There is an entry
    # for each orientation in the coarse angles list.
    poses = []
    for a in self._coarseAngles:
      @todo

    # Produce a list of preference function values, one per angle.  If we don't
    # have a preference function, the value is 0 for all of them.
    preferences = []
    val = 0.0
    for a in self._coarseAngles:
      if self._preferenceFunction is not None:
        val = self._preferenceFunction(a)
      preferences.append(val)
    @todo Multiply by the reduce option value

    # returns: The original atom at its original position.
    return PositionReturn([ self._atoms ], poses, preferences)

  def FinePositions(self, coarseIndex, reduceOptions):
    if not self._doFineRotations:
      # returns: No fine positions for any coarse position.
      return PositionReturn([], [], [])
    else:
      @todo

  def FixUp(self, coarseIndex, reduceOptions):
    # No fixups for any coarse index.
    return FixUpReturn([], [])


# @todo Define each type of Mover

def Test():
  """Test function for all classes provided above.
  returns: Empty string on success, string describing the problem on failure.
  """
  # @todo


