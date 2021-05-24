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
# This is a set of functions that implement Reduce's Interaction Graphs.  These are
# functions that produce Boost graphs of Movers, enabling easy determination
# of Cliques (connected components of this graph).
from boost_adaptbx import graph
from boost_adaptbx.graph import connected_component_algorithm as cca
import Movers

def AABBOverlap(box1, box2):
  """Helper function that tells whether two axis-aligned bounding boxes overlap.
  :param box1: list of three ranges, for X, Y, and Z, that indicate one axis-aligned
  bounding box.
  :param box2: list of three ranges, for X, Y, and Z, that indicate another axis-aligned
  bounding box.
  :returns True if the boxes overlap, False if not.
  """

  return ( (box1[0][0] <= box2[0][1] and box1[0][1] >= box2[0][0]) and
           (box1[1][0] <= box2[1][1] and box1[1][1] >= box2[1][0]) and
           (box1[2][0] <= box2[2][1] and box1[2][1] >= box2[2][0]) )

def InteractionGraphAABB(movers, extraAtomInfo, reduceOptions):
  """Uses the overlap of the axis-aligned bounding boxes (AABBs) of all possible
  positions of all movable atoms in the set of movers passed in to construct the
  graph of which might overlap across all possible orientations of each.  The use
  of bounding boxes makes this an overestimate but its time complexity is linear
  in the product of the number of movers times the number of atoms times the number
  of possible positions for each and quadratic in the number of Movers.

  :param movers: flex array of movers to add to the graph.
  :param extaAtomInfo: flex array of Probe.ExtraAtomInfo classes that have a vdwRadius
  property.  Warning: The i_seq values from the atoms in the Movers are used to look
  up directly in this vector so they must not have changed (due to structure modification)
  since the extaInfo vector or the atom structure used by the Movers were generated.
  The positions of individual atoms can have been moved but atoms cannot have been
  removed and re-added to the structure.
  :param reduceOptions: a Phil option subset.  The relevant option is probeRadius.
  If it is not set, the default value of 0.25 will be used.
  :returns An undirected Boost graph whose nodes are Movers and whose edges
  indicate which Movers might overlap in any of their states.
  """

  # Determine the probe radius to use.  If it is specified in the reduceOptions,
  # then use its value.  Otherwise, use the default value of 0.25
  try:
    pr = reduceOptions.probeRadius
  except:
    pr = 0.25

  # Add all of the Movers as nodes in the graph
  # Compute the axis-aligned bounding box for each Mover
  ret = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "Mover",
        edge_type = "set",
        )
  AABBs = []
  for m in movers:
    ret.add_vertex(m)

    # Find all possible positions, coarse and fine.
    coarses = m.CoarsePositions(reduceOptions)
    atoms = coarses.atoms
    coarsePositions = coarses.positions
    total = coarsePositions.copy()
    for c in len(coarsePositions):
      total.extend(m.FinePositions(c, reduceOptions).positions)

    # Find the range of positions of all atoms in X, Y, and Z
    xRange = [ 1e10, -1e10 ]
    yrange = [ 1e10, -1e10 ]
    zRange = [ 1e10, -1e10 ]
    for pos in total:
      for i, atomLoc in enumerate(pos):
        # Find the radius of the atom, which is used to extend it in all directions
        # so that we catch all potential overlaps.
        r = ExtraAtomInfo[atoms[i].i_seq].vdwRadius

        x = atomLoc[0]
        xRange[0] = min(xRange[0], x - r)
        xRange[1] = max(xRange[1], x + r)
            
        y = atomLoc[1]
        yRange[0] = min(yRange[0], y - r)
        yRange[1] = max(yRange[1], y + r)

        z = atomLoc[2]
        zRange[0] = min(zRange[0], z - r)
        zRange[1] = max(zRange[1], z + r)

    # Dilate the bounding box by the radius of the probe.
    xRange = [ xRange[0] - pr, xRange[1] + pr ]
    yRange = [ yRange[0] - pr, yRange[1] + pr ]
    zRange = [ zRange[0] - pr, zRange[1] + pr ]

    # Store the bounding boxes for this Mover
    AABBs.append( [xRange, yRange, zRange] )

  # For each pair of Movers whose bounding boxes overlap, add an
  # edge to the graph
  for i in range(len(movers)-1):
    for j in range(i+1, len(movers)):
      if AABBOverlap(AABBs[i], AABBs[j]):
        ret.add_edge( vertex1 = movers[i], vertex2 = movers[j])

  return ret

#######################################################################################################
# Test code and objects below here

from iotbx import pdb
import math
import mmtbx_probe_ext as probe

def Test():
  """Test function for all functions provided above.
  returns: Empty string on success, string describing the problem on failure.
  """
  # @todo

  # Construct a Class that will behave like the Reduce Phil data structure so that
  # we can specify the probe radius.
  class FakePhil:
    pass
  fakePhil = FakePhil()

  # Construct a class that will behave like a Mover but which returns a single result
  # atom at a single location.  We'll use this rather than actual Movers as a simple and
  # fast test case.
  class FakeMover():
    def __init__(self, atom):
      self._atom = atom
    def CoarsePositions(reduceOptions):
      return ( [ self._atom ], [ [ self.atom.data.xyz[0], self.atom.data.xyz[1], self.atom.data.xyz[2] ] ], [ 0.0 ] )
    def FinePositions(coarseIndex, reduceOptions):
      return ( [], [], [] )
    def FixUp(coarseIndex, reduceOptions):
      return ( [], [] )

  # Construct a set of Mover/atoms that will be used to test the routines.  They will all be part of the
  # same residue and they will all have the same unit radius in the extraAtomInfo associated with them.
  # There will be a set of five along the X axis, with one pair overlapping slightly and the others
  # spaced 0.45 units apart so that they will overlap when using a probe radius of 0.25 (diameter 0.5).
  # There will be another one that is obliquely located away from the first such that it will overlap
  # in a bounding-box test but not in a true atom-comparison test for a probe with radius 0.25.  There
  # will be a final one 10 units above the origin.
  rad = 1.0
  probeRad = 0.25
  locs = [ [0.0, 0.0, 0.0], [1.9, 0.0, 0.0] ]
  for i in range(1,4):
    loc = [1.9 + 2.1*i, 0.0, 0.0]
    locs.append(loc)
  delta = 2*rad + 2*probeRad + 0.1
  dist = - delta * math.cos(math.pi/4)
  dist = - delta * math.sin(math.pi/4)
  locs.append([dist, dist, 0.0])
  locs.append([0.0, 0.0, 10.0])

  name = " H  ";
  ag = pdb.hierarchy.atom_group()
  ag.resname = "LYS"
  atoms = pdb.hierarchy.af_shared_atom()
  extras = []
  movers = []
  baseAtom = pdb.hierarchy.atom()
  for i in range(len(locs)):
    a = pdb.hierarchy.atom(parent = ag, other=baseAtom)
    a.name = name
    a.xyz = locs[i]
    atoms.append(a)
    e = probe.ExtraAtomInfo(rad)
    extras.append(e)
    movers.append(FakeMover(a))
  # Fix the sequence numbers, which are otherwise all 0
  atoms.reset_i_seq()

  # Generate a table of parameters and expected results.  The first entry in each row is the
  # function to call (the type of interactiong graph).  The second is the probe radius.  The third
  # is the expected number of connected components.  The fourth is the size of the largest connected
  # component.
  _expectedCases = [
    [ InteractionGraphAABB, 0.0, 5, 2 ],
    [ InteractionGraphAABB, probeRad, 2, 6 ],
    [ InteractionGraphAABB, 100, 1, 7 ]
  ]

  # Specify the probe radius and run the test.  Compare the results to what we expect.
  for e in _expectedCases:
    fakePhil.probeRadius = e[1]
    graph = e[0](movers, extras, fakePhil)

    # Find the connected components of the graph and compare their counts and maximum size to
    # what is expected.
    # @todo


  # Specify a probe radius that is so large it will make all of the Movers overlap.
  fakePhil.probeRadius = 10000

# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  ret = Test()
  if len(ret) == 0:
    print('Success!')

  assert (len(ret) == 0)
