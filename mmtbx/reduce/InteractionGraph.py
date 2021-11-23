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

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import

from boost_adaptbx import graph
import mmtbx.reduce.Movers as Movers
import mmtbx_probe_ext as probeExt

def InteractionGraphAllPairs(movers, extraAtomInfoMap, probeRadius = 0.25):
  """Tests for overlap of all possible positions of all movable atoms between each
  pair of Movers in the set of Movers passed in to construct the
  graph of which overlap across all possible orientations of each.

  :param movers: flex array of movers to add to the graph.  Note that this list must
  not be modified after the graph has been constructed because that will change the
  index of its elements, making the graph point to the wrong elements (or to elements
  that no longer exist).
  :param extraAtomInfoMap: probe.ExtraAtomInfoMap that can be used to look
  up the information for atoms whose values need to be changed.  Can be
  obtained by calling mmtbx.probe.Helpers.getExtraAtomInfo().
  :param probeRadius: Radius of the probe to use to determine neighbor contact.
  If it is not set, the default value of 0.25 will be used.
  :returns (1) An undirected Boost graph whose nodes are indices into the movers list
  and whose edges indicate which Movers might overlap in any of their states.  Note that
  the mover list must not be modified after the graph has been constructed because
  that will change the index of its elements, making the graph point to the wrong
  elements (or to elements that no longer exist). (2) A dictionary with atoms as the
  key that returns the set of Movers that the atom interacts with; each has at least
  the Mover that it is a part of and may contain additional ones when they overlap.
  """

  # Run the AABB test to get a superset of the list of pairs that we need to check for
  # overlap.  If we try to brute-force all of the Movers against all of the others, it
  # takes too long.
  myGraph = _InteractionGraphAABB(movers, extraAtomInfoMap, probeRadius)

  # Dictionary looked up by atom that returns the set of Movers that atom interacts
  # with.
  atomMoverSets = {}

  # Dictionaries of list of atoms per mover and dictionary of list of positions per atom per mover.
  # Each of these is indexed the same way that movers is, so finding the index of the
  # mover gets the same index for them.
  atoms = {}
  positions = {}
  for m in movers:

    # Find all possible positions, coarse and fine, for each atom
    # in this mover.
    coarses = m.CoarsePositions()
    coarsePositions = coarses.positions
    total = coarsePositions[:]
    for c in range(len(coarsePositions)):
      total.extend(m.FinePositions(c).positions)

    # Add the atoms and positions into our dictionaries
    atoms[m] = coarses.atoms
    positions[m] = total
    for a in coarses.atoms:
      atomMoverSets[a] = {m}

  # For each pair of movers that are connected by an edge in the graph produced
  # by the AABB algorithm to see if they actually overlap.  If not, remove that edge.
  for e in myGraph.edges():
    sourceMover = myGraph.vertex_label( myGraph.source(e) )
    targetMover = myGraph.vertex_label( myGraph.target(e) )
    if not _PairsOverlap(sourceMover, atoms[sourceMover], positions[sourceMover],
        targetMover, atoms[targetMover], positions[targetMover],
        extraAtomInfoMap, probeRadius,
        atomMoverSets):
      myGraph.remove_edge( e )

  return myGraph, atomMoverSets

#######################################################################################################
# Internal helper functions defined here

def _InteractionGraphAABB(movers, extraAtomInfoMap, probeRadius = 0.25):
  """Uses the overlap of the axis-aligned bounding boxes (AABBs) of all possible
  positions of all movable atoms in the set of movers passed in to construct the
  graph of which might overlap across all possible orientations of each.  The use
  of bounding boxes makes this an overestimate but its time complexity is linear
  in the product of the number of movers times the number of atoms times the number
  of possible positions for each and quadratic in the number of Movers.

  :param movers: flex array of movers to add to the graph.  Note that this list must
  not be modified after the graph has been constructed because that will change the
  index of its elements, making the graph point to the wrong elements (or to elements
  that no longer exist).
  :param extraAtomInfoMap: probe.ExtraAtomInfoMap that can be used to look
  up the information for atoms whose values need to be changed.  Can be
  obtained by calling mmtbx.probe.Helpers.getExtraAtomInfo().
  :param probeRadius: Radius of the probe to use to determine neighbor contact.
  If it is not set, the default value of 0.25 will be used.
  :returns An undirected Boost graph whose nodes are indices into the movers list
  and whose edges indicate which Movers might overlap in any of their states.  Note that
  the mover list must not be modified after the graph has been constructed because
  that will change the index of its elements, making the graph point to the wrong
  elements (or to elements that no longer exist).
  """

  pr = probeRadius

  # Add all of the Movers as nodes in the graph
  # Compute the axis-aligned bounding box for each Mover
  ret = graph.adjacency_list(
        vertex_type = "list",   # List so that deletions do not invalidate iterators and descriptors
        )
  AABBs = []
  verts = []
  for m in movers:
    verts.append(ret.add_vertex(m))

    # Find all possible positions, coarse and fine.
    coarses = m.CoarsePositions()
    atoms = coarses.atoms
    coarsePositions = coarses.positions
    total = coarsePositions[:]
    for c in range(len(coarsePositions)):
      total.extend(m.FinePositions(c).positions)

    # Find the range of positions of all atoms in X, Y, and Z
    xRange = [ 1e10, -1e10 ]
    yRange = [ 1e10, -1e10 ]
    zRange = [ 1e10, -1e10 ]
    for pos in total:
      for i, atomLoc in enumerate(pos):
        # Find the radius of the atom, which is used to extend it in all directions
        # so that we catch all potential overlaps.
        r = extraAtomInfoMap.getMappingFor(atoms[i]).vdwRadius

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
    # Because we're dilating each box by this radius, we're properly
    # checking to twice the probe radius between two Movers.
    xRange = [ xRange[0] - pr, xRange[1] + pr ]
    yRange = [ yRange[0] - pr, yRange[1] + pr ]
    zRange = [ zRange[0] - pr, zRange[1] + pr ]

    # Store the bounding boxes for this Mover
    AABBs.append( [xRange, yRange, zRange] )

  # For each pair of Movers whose bounding boxes overlap, add an
  # edge to the graph.  We add them based on their indices.
  for i in range(len(movers)-1):
    for j in range(i+1, len(movers)):
      if _AABBOverlap(AABBs[i], AABBs[j]):
        ret.add_edge( vertex1 = verts[i], vertex2 = verts[j])

  return ret

def _AABBOverlap(box1, box2):
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

def _PairsOverlap(mover1, atoms1, positions1,
                  mover2, atoms2, positions2, extraAtomInfoMap, probeRad, atomMoverSets):
  """Helper function that tells whether any pair of atoms from two Movers overlap.
  :param mover1: The first Mover
  :param atoms1: Atom list for the first Mover
  :param positions1: probe.PositionReturn.positions holding possible positions for each.
  :param mover2: The first Mover
  :param atoms2: Atom list for the second Mover
  :param positions2: probe.PositionReturn.positions holding possible positions for each.
  :param extraAtomInfoMap: probe.ExtraAtomInfoMap that can be used to look
  up the information for atoms whose values need to be changed.  Can be
  obtained by calling mmtbx.probe.Helpers.getExtraAtomInfo().
  :param ProbeRad: Probe radius
  :param atomMoverSets: Parameter that is modified in place to record all Movers that
  a particular atom interacts with.  An entry is created whenever there is overlap
  with an atom in another Mover.
  :returns True if a pair of atoms with one from each overlap, False if not.
  """

  for i1, p1 in enumerate(positions1):
    for ai1 in range(len(p1)):
      r1 = extraAtomInfoMap.getMappingFor(atoms1[ai1]).vdwRadius
      for i2, p2 in enumerate(positions2):
        for ai2 in range(len(p2)):
          r2 = extraAtomInfoMap.getMappingFor(atoms2[ai2]).vdwRadius
          dx = p1[ai1][0] - p2[ai2][0]
          dy = p1[ai1][1] - p2[ai2][1]
          dz = p1[ai1][2] - p2[ai2][2]
          dSquared = dx*dx + dy*dy + dz*dz
          limit = r1 + r2 + 2*probeRad
          limitSquared = limit*limit
          if dSquared <= limitSquared:
            # Add the opposite Mover to each atom; they interact
            atomMoverSets[atoms1[ai1]].add(mover2)
            atomMoverSets[atoms2[ai2]].add(mover1)
            return True
  return False

#######################################################################################################
# Test code and objects below here

from iotbx import pdb
import math
import mmtbx_probe_ext as probe
from boost_adaptbx.graph import connected_component_algorithm as cca

def Test():
  """Test function for all functions provided above.
  returns: Empty string on success, string describing the problem on failure.
  """

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
  delta = 2*rad + 2*probeRad - 0.1
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
    extrasMap = probeExt.ExtraAtomInfoMap(atoms, extras)
    movers.append(Movers.MoverNull(a, extrasMap))
  # Fix the sequence numbers, which are otherwise all 0
  atoms.reset_i_seq()

  # Generate a table of parameters and expected results.  The first entry in each row is
  # the probe radius.  The second is the expected number of connected components.
  # The third is the size of the largest connected component.
  _expectedCases = [
    [ 0.0, 5, 3 ],
    [ probeRad, 2, 6 ],
    [ 100, 1, 7 ]
  ]

  # Specify the probe radius and run the test.  Compare the results to what we expect.
  for i, e in enumerate(_expectedCases):
    probeRadius = e[0]
    g = _InteractionGraphAABB(movers, extrasMap, probeRadius)

    # Find the connected components of the graph and compare their counts and maximum size to
    # what is expected.
    components = cca.connected_components( graph = g )
    if len(components) != e[1]:
      return "AABB Expected "+str(e[1])+" components, found "+str(len(components))+" for case "+str(i)
    maxLen = -1
    for c in components:
      if len(c) > maxLen:
        maxLen = len(c)
    if maxLen != e[2]:
      return "AABB Expected max sized component of "+str(e[2])+", found "+str(maxLen)+" for case "+str(i)

  # Generate a table of parameters and expected results.  The first entry in each row is
  # the probe radius.  The second is the expected number of connected components.
  # The third is the size of the largest connected component.
  # The fourth (not present in the AABB table above) is the set of expected sizes of
  # atomMoverSets across all atoms; not one per atom but across all atoms what answers are
  # expected.  The easiest to explain is the 100-radius entry, which should have all atoms interacting
  # with all Movers so the only answer across all atoms is 7.  The 0-radius case has only one pair
  # of overlaps, so only up to 2 Movers per atom.  The middle case has some Movers overlapping with
  # two neighbors, so up to 3 Movers associated with a given atom.
  _expectedCases = [
    # One of the pairs actually does not overlap for the all-pairs test.  Other conditions are the same
    # as the AABB tests.
    [ 0.0, 6, 2, {1,2} ],
    [ probeRad, 2, 6, {1,2,3} ],
    [ 100, 1, 7, {7} ]
  ]

  # Specify the probe radius and run the test.  Compare the results to what we expect.
  for i, e in enumerate(_expectedCases):
    probeRadius = e[0]
    g, am = InteractionGraphAllPairs(movers, extrasMap, probeRadius)

    # Find the connected components of the graph and compare their counts and maximum size to
    # what is expected.
    components = cca.connected_components( graph = g )
    if len(components) != e[1]:
      return "Expected "+str(e[1])+" components, found "+str(len(components))+" for case "+str(i)
    maxLen = -1
    for c in components:
      if len(c) > maxLen:
        maxLen = len(c)
    if maxLen != e[2]:
      return "Expected max sized component of "+str(e[2])+", found "+str(maxLen)+" for case "+str(i)

    # Check atom/Mover overlaps by finding the set of lengths that are present accross all atoms.
    lengths = set()
    for a in atoms:
      lengths.add(len(am[a]))
    if lengths != e[3]:
      return "Expected set of overlap counts "+str(e[3])+", found "+str(lengths)+" for case "+str(i)

  return ""

# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  ret = Test()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
