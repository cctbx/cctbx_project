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
# functions that produce NetworkX graphs sets of Movers, enabling easy determination
# of Cliques (connected components of this graph).
import networkx as nx
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

def InteractionGraphAABB(movers, reduceOptions):
  """Uses the overlap of the axis-aligned bounding boxes (AABBs) of all possible
  positions of all movable atoms in the set of movers passed in to construct the
  graph of which might overlap across all possible orientations of each.  The use
  of bounding boxes makes this an overestimate but its time complexity is linear
  in the product of the number of movers times the number of atoms times the number
  of possible positions for each and quadratic in the number of Movers.

  :param movers: flex array of movers to add to the graph.
  :param reduceOptions: a Phil option subset.  The relevant option is probeRadius.
  :returns An undirected NetworkX graph whose nodes are Movers and whose edges
  indicate which Movers might overlap in any of their states.
  """

  # Add all of the Movers as nodes in the graph
  # Compute the axis-aligned bounding box for each Mover
  ret = nx.Graph()
  AABBs = []
  for m in movers:
    ret.add_node(m)

    # Find all possible positions, coarse and fine.
    atoms, coarse, favorable = m.CoarsePositions(reduceOptions)
    total = coarse.copy()
    for c in len(coarse):
      total.extend(m.FinePositions(c, reduceOptions)[1])

    # Find the range of positions of all atoms in X, Y, and Z
    xRange = [ 1e10, -1e10 ]
    yrange = [ 1e10, -1e10 ]
    zRange = [ 1e10, -1e10 ]
    for pos in total:
      for atomLoc in pos:
         x = atomLoc[0]
         xRange[0] = min(xRange[0], x)
         xRange[1] = max(xRange[1], x)
            
         y = atomLoc[1]
         yRange[0] = min(yRange[0], y)
         yRange[1] = max(yRange[1], y)

         z = atomLoc[2]
         zRange[0] = min(zRange[0], z)
         zRange[1] = max(zRange[1], z)

    # Dilate the bounding box by the radius of the probe if it is
    # specified in the parameter set.
    try:
      r = reduceOptions.probeRadius
      if r is not None:
        xRange = [ xRange[0] - r, xRange[1] + r ]
        yRange = [ yRange[0] - r, yRange[1] + r ]
        zRange = [ zRange[0] - r, zRange[1] + r ]

    # Store the bounding boxes for this Mover
    AABBs.append( [xRange, yRange, zRange] )


  # For each pair of Movers whose bounding boxes overlap, add an
  # edge to the graph
  for i in range(len(movers)-1):
    for j in range(i+1, len(movers)):
      if AABBOverlap(AABBs[i], AABBs[j]):
        ret.add_edge(movers[i], movers[j])

  return ret

def Test():
  """Test function for all functions provided above.
  returns: Empty string on success, string describing the problem on failure.
  """
  # @todo