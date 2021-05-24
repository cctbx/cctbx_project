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
#  - { flex<atom> atoms, flex<flex<vec3>> positions, flex<float> preferenceEnergies } CoarsePositions(reduceOptions)
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
#  - { flex<atom> atoms, flex<flex<vec3>> positions, flex<float> preferenceEnergies } FinePositions(coarseIndex, reduceOptions)
#     The coarseIndex indicates the index (0 for the first) of the relevant coarse
#       orientation.
#     The return values are the same as for CoarsePositions and they list potential
#       fine motions around the particular coarse position (not including the position
#       itself).  This function can be used by optimizers that wish to do heavy-weight
#       operations at a coarse resolution and then lightweight operations at a finer
#       scale; other optimizers may choose to ask for all of the fine positions and
#       append them to the coarse positions and globally optimize.
#     Note: Some Movers will return empty arrays.
#  - { flex<atom> atoms, flex<vec3> newPositions } FixUp(coarseIndex, reduceOptions)
#     The coarseIndex indicates the index (0 for the first) of the relevant coarse
#       orientation that was finally chosen.
#     The first return lists the atoms that should be repositioned.
#     The second return lists the new locations for each atom.
#     Note: This function may modify atoms other than the ones it reported in its
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
# Movers have overlaps between movable atoms in both of them.
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
class NullMover:
  def __init__(self, atom):
    self._atom = atom
  def CoarsePositions(self, reduceOptions):
    return PositionReturn([ self._atom ],
        [ [ [ self._atom.xyz[0], self._atom.xyz[1], self._atom.xyz[2] ] ] ],
        [ 0.0 ])
  def FinePositions(self, coarseIndex, reduceOptions):
    return PositionReturn([], [], [])
  def FixUp(self, coarseIndex, reduceOptions):
    return FixUpReturn([], [])


# @todo Define each type of Mover

def Test():
  """Test function for all functions provided above.
  returns: Empty string on success, string describing the problem on failure.
  """
  # @todo


