// Copyright(c) 2023, Richardson Lab at Duke
// Licensed under the Apache 2 license
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissionsand
// limitations under the License.

#include "Optimizers.h"
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/subgraph.hpp>

typedef scitbx::vec3<double> vec3;

static std::vector< std::vector<unsigned> > generateAllStates(std::vector<unsigned> const &numStates)
{
  std::vector< std::vector<unsigned> > ret;

  // Initialize the states with all zeros -- the first state.
  std::vector<unsigned> curStateValues(numStates.size(), 0);
  unsigned curState = 0;

  // Cycle the curStateValues[] list through all possible states for each element,
  // incremementing each until it rolls over and then jumping up to the next.
  // This is similar to doing + 1 arithmetic with carry on a multi - digit number.
  // Increment the state.We do this by increasing the current element until it reaches its
  // number of values, then we bump it and all of its neighbors to the right back to 0 and, if
  // we're not at the left end, bump the next one up.  If we are at the left end, we're done.
  // When done, return.

  while (true) {
    ret.push_back(curStateValues);

    // Go to the next state value, if there is one.
    curStateValues[curState]++;
    bool rippled = false;
    while (curStateValues[curState] == numStates[curState]) {
      // Clear all values to the right, and ours, because we're rolling over
      for (unsigned i = 0; i <= curState; i++) {
        curStateValues[i] = 0;
      }
      // If we're the left-most state, we're done
      if (curState+1 >= numStates.size()) {
        return ret;
      } else {
        curState++;
        curStateValues[curState]++;
        rippled = true;
      }
    }
    // If we rippled, bump back to the right-most column and start counting there again.
    if (rippled) {
      curState = 0;
    }
  }
}

static std::string setMoverState(molprobity
  ::reduce::PositionReturn & positionReturn,
  unsigned index,
  molprobity::probe::SpatialQuery& spatialQuery,
  molprobity::probe::ExtraAtomInfoMap& extraAtomInfoMap,
  boost::python::object& deleteMes,
  int verbosity)
{
  std::string ret;

  // Move the atoms to their new positions, updating the spatial query structure
  // by removing the old and adding the new location.
  for (size_t i = 0; i < positionReturn.atoms.size(); i++) {
    iotbx::pdb::hierarchy::atom &a = positionReturn.atoms[i];

    spatialQuery.remove(a);
    // Overwrite the location of the atom
    for (size_t j = 0; j < 3; j++) {
      a.data->xyz[j] = positionReturn.positions[index][i][j];
    }
    spatialQuery.add(a);
  }

  // Update the extraAtomInfo associated with each atom.
  // Note that there may be fewer entries than atoms, but they all correspond
  // so we can look up the atom by index in the atoms array just like above.
  for (size_t i = 0; i < positionReturn.extraInfos[index].size(); i++) {
    iotbx::pdb::hierarchy::atom& a = positionReturn.atoms[i];
    extraAtomInfoMap.setMappingFor(a, positionReturn.extraInfos[index][i]);
  }

  // Manage the deletion status of each atom, including ensuring
  // consistency with the spatial - query structure.
  // Note that there may be fewer entries than atoms, but they all correspond
  // so we can look up the atom by index in the atoms array just like above.
  for (size_t i = 0; i < positionReturn.deleteMes[index].size(); i++) {
    iotbx::pdb::hierarchy::atom& a = positionReturn.atoms[i];
    bool doDelete = boost::python::extract<bool>(deleteMes.contains(a));
    if (doDelete) {
      spatialQuery.remove(a);
      deleteMes.attr("add")(a);
      if (verbosity >= 10) {
        ret += "          Deleting atom\n";
      }
    } else {
      spatialQuery.add(a);
      deleteMes.attr("discard")(a);
      if (verbosity >= 10) {
        ret += "           Ensuring deletable atom is present\n";
      }
    }
  }

  return ret;
}

// Score all atoms that have not been marked for deletion, calling the Python object's
// scoring function (which may or may not use caching to do the scoring).
static double scorePosition(boost::python::object& self,
  molprobity::reduce::PositionReturn& states, size_t index)
{
  double ret = 0;
  for (size_t a = 0; a < states.atoms.size(); a++) {
    // There may not be as many deleteMes as there are atoms, so we need to check for that.
    if ((a >= states.deleteMes[index].size()) || !states.deleteMes[index][a]) {
      boost::python::object result = self.attr("_scoreAtom")(states.atoms[a]);
      double resultValue = boost::python::extract<double>(result);
      ret += resultValue;
    }
  }
  return ret;
}

namespace molprobity {
  namespace reduce {

boost::python::tuple OptimizeCliqueCoarseBruteForceC(
  boost::python::object &self,
  int verbosity,
  double preferenceMagnitude,
  scitbx::af::shared<boost::python::object> movers,
  molprobity::probe::SpatialQuery &spatialQuery,
  molprobity::probe::ExtraAtomInfoMap &extraAtomInfoMap,
  boost::python::object &deleteMes,
  boost::python::dict &coarseLocations,
  boost::python::dict &highScores
)
{
  // Information to pass back about what we did, if verbosity is high enough.
  std::string infoString;

  if (verbosity >= 3) {
    std::ostringstream oss;
    oss << "   Optimizing clique of size " << movers.size() << " using brute force\n";
    infoString = oss.str();
  }

  // Fill in vectors with the states of each mover.
  scitbx::af::shared<molprobity::reduce::PositionReturn> states;
  for (scitbx::af::shared<boost::python::object>::const_iterator it = movers.begin(); it != movers.end(); ++it) {
    boost::python::object const&m = *it;
    states.push_back(boost::python::extract<molprobity::reduce::PositionReturn>(m.attr("CoarsePositions")()));
  }

  // Keep track of the best score and the state where we found it.
  // It starts out empty, letting us know to fill it.
  double bestScore = -1e100;
  std::vector<unsigned> bestState;

  // Find the length of each state and record it into a vector
  std::vector<unsigned> numStates;
  for (scitbx::af::shared<molprobity::reduce::PositionReturn>::const_iterator it = states.begin();
       it != states.end(); ++it) {
    molprobity::reduce::PositionReturn const& s = *it;
    numStates.push_back(s.positions.size());
  }

  // Generate a vector of state combinations, storing the index of the selected
  // element from each state. We will iterate over these to find the best state.
  /// @todo Consider pulling this generation into here as a loop rather than using up
  // memory to store this.
  std::vector< std::vector<unsigned> > allStates = generateAllStates(numStates);

  // Cycle through all states and find the one with the largest score.
  for (std::vector< std::vector<unsigned> >::const_iterator it = allStates.begin(); it != allStates.end(); ++it) {
    std::vector<unsigned> const& curStateValues = *it;

    // Set all movers to match the state list
    for (unsigned m = 0; m < movers.size(); m++) {
      // Only change this mover if it is different from the last time
      boost::python::extract<unsigned> value(coarseLocations.get(movers[m]));
      if (value != curStateValues[m]) {
        // Set the mover to this state
        infoString += setMoverState(states[m], curStateValues[m], spatialQuery, extraAtomInfoMap, deleteMes,
          verbosity);
        coarseLocations[movers[m]] = curStateValues[m];
      }
    }

    // Compute the score over all atoms in all Movers.
    double score = 0;
    for (size_t m = 0; m < movers.size(); m++) {

      score += preferenceMagnitude * states[m].preferenceEnergies[curStateValues[m]];

      // Score all atoms that have not been marked for deletion, calling the Python object's
      // scoring function (which may or may not use caching to do the scoring).
      score += scorePosition(self, states[m], curStateValues[m]);
    }
    if (verbosity >= 5) {
      std::ostringstream oss;
      oss << "    Score is " << score << " at [";
      infoString += oss.str();
      for (unsigned i = 0; i < curStateValues.size(); i++) {
        std::ostringstream oss2;
        oss2 << curStateValues[i];
        infoString += oss2.str();
        if (i < curStateValues.size() - 1) {
          infoString += ", ";
        }
      }
      infoString += "]\n";
    }

    // See if it is the best score.  If so, update the best.
    if ((score > bestScore) || (bestState.size() == 0)) {
      if (verbosity >= 4) {
        std::ostringstream oss;
        oss << "    New best score is " << score << " at [";
        infoString += oss.str();
        for (unsigned i = 0; i < curStateValues.size(); i++) {
          std::ostringstream oss2;
          oss2 << curStateValues[i];
          infoString += oss2.str();
          if (i < curStateValues.size() - 1) {
            infoString += ", ";
          }
        }
        infoString += "]\n";
      }
      bestScore = score;
      bestState = curStateValues;
    }
  }

  // Put each Mover into its state in the best configuration and compute its high-score value.
  // Store the individual scores for these Movers in the best config for use in later fine - motion
  // processing.
  for (size_t m = 0; m < movers.size(); m++) {
    infoString += setMoverState(states[m], bestState[m], spatialQuery, extraAtomInfoMap, deleteMes,
      verbosity);
    coarseLocations[movers[m]] = bestState[m];
    double score = preferenceMagnitude * states[m].preferenceEnergies[bestState[m]];
    score += scorePosition(self, states[m], bestState[m]);
    highScores[movers[m]] = score;
    if (verbosity >= 3) {
      std::ostringstream oss;
      oss << "    Setting Mover in clique to coarse orientation " << bestState[m]
        << ", max score = " << score << "\n";
      infoString += oss.str();
    }
  }

  // Return the result
  return boost::python::make_tuple(bestScore, infoString);
}

// Must use vector style for second (vertex) entry for connected_components to work
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::python::object*>
  CliqueGraph;

/// @brief Return a vector of vectors of integers representing all combinations of m integers
///        from 0 to n-1, inclusive.
std::vector< std::vector<int> > nChooseM(int n, int m) {
  std::vector<std::vector<int> > result;

  std::vector<int> indices(m);
  for (int i = 0; i < m; ++i) {
    indices[i] = i + 1;
  }

  while (true) {
    std::vector<int> currentCombination;
    for (size_t it = 0; it < indices.size(); ++it) {
      // Use zero-based indexing for the results
      currentCombination.push_back(indices[it] - 1);
    }
    result.push_back(currentCombination);

    int i = m - 1;
    while (i >= 0 && indices[i] == n - m + i + 1) {
      --i;
    }

    if (i < 0) {
      break;
    }

    ++indices[i];
    for (int j = i + 1; j < m; ++j) {
      indices[j] = indices[j - 1] + 1;
    }
  }

  return result;
}

/// @brief Find the subset of a clique graph that contains a given set of movers and edges between them.
static CliqueGraph subsetGraph(CliqueGraph const& graph, std::vector<boost::python::object*>& keepMovers)
{
  CliqueGraph ret;

  // We keep a from labels to the vertex in the new graph that points to
  // that label so that we can construct edges in the new graph.
  std::map<boost::python::object*, CliqueGraph::vertex_descriptor> vertexMap;

  // Construct a subgraph that consists only of vertices to be kept and
  // edges both of whose ends are on these vertices.
  for (std::vector<boost::python::object*>::iterator it = keepMovers.begin(); it != keepMovers.end(); ++it) {
    // Add a vertex for this mover, and keep track of its descriptor.
    boost::python::object* mover = *it;
    vertexMap[mover] = boost::add_vertex(mover, ret);
  }
  boost::iterator_range<CliqueGraph::edge_iterator> edges = boost::edges(graph);
  for (CliqueGraph::edge_iterator e = edges.begin(); e != edges.end(); e++) {
    CliqueGraph::vertex_descriptor v1 = boost::source(*e, graph);
    CliqueGraph::vertex_descriptor v2 = boost::target(*e, graph);
    boost::python::object* mover1 = graph[v1];
    boost::python::object* mover2 = graph[v2];
    if ( (std::find(keepMovers.begin(), keepMovers.end(), mover1) != keepMovers.end()) &&
         (std::find(keepMovers.begin(), keepMovers.end(), mover2) != keepMovers.end())) {
      // Both ends of this edge are in the subset, so add it to the subset graph.
      boost::add_edge(vertexMap[mover1], vertexMap[mover2], ret);
    }
  }

  return ret;
}

/// @brief Find one of the smallest vertex cuts in a clique graph.
/// @param [in] graph: The clique graph to find the vertex cut in.
/// @param [out] cutMovers: The movers that correspond to the vertices
/// that are removed.
/// @param [out] cutGraph: The graph that stores the subset of clique
/// whose vertices have been removed.
static void findVertexCut(CliqueGraph const& graph,
  std::vector<boost::python::object*>& cutMovers, CliqueGraph& cutGraph)
{
  // Check all vertex cut sizes from 1 to 2 less than the number of vertices(we must
  // have at least 2 vertices left to have a disconnected graph).
  size_t numVerts = boost::num_vertices(graph);
  for (size_t n = 0; n < numVerts; n++) {
    // Iterate over all sets of vertices of size n that might be removed
    std::vector< std::vector<int> > vertexSets = nChooseM(numVerts, n);
    for (std::vector< std::vector<int> >::const_iterator it = vertexSets.begin(); it != vertexSets.end(); ++it) {
      std::vector<int> const& vertexSet = *it;

      // Get a list of the vertices to keep by removing the ones in the removal set.
      std::vector<boost::python::object*> keepMovers;
      for (size_t i = 0; i < numVerts; i++) {
        if (std::find(vertexSet.begin(), vertexSet.end(), i) == vertexSet.end()) {
          keepMovers.push_back(graph[i]);
        }
      }

      // Get the subset of the graph that contains only the vertices to keep.
      CliqueGraph currentGraph = subsetGraph(graph, keepMovers);

      // Check if the graph is disconnected
      std::vector<int> componentMap(boost::num_vertices(currentGraph));
      int numComponents = boost::connected_components(currentGraph, &componentMap[0]);
      if (numComponents > 1) {
        // The graph is disconnected, so we have found a vertex cut.
        // Make a list of the movers that correspond to the vertices that were removed.
        for (size_t i = 0; i < vertexSet.size(); i++) {
          cutMovers.push_back(graph[vertexSet[i]]);
        }
        cutGraph = currentGraph;
        return;
      }
    }
  }

  // We didn't find an answer. Return a complete copy of the graph and an empty list of movers.
  cutGraph = graph;
  cutMovers.clear();
}


boost::python::tuple OptimizeCliqueCoarseVertexCutC(
  boost::python::object& self,
  int verbosity,
  double preferenceMagnitude,
  scitbx::af::shared<boost::python::object> movers,
  scitbx::af::versa<int, scitbx::af::flex_grid<> > &interactions,
  molprobity::probe::SpatialQuery& spatialQuery,
  molprobity::probe::ExtraAtomInfoMap& extraAtomInfoMap,
  boost::python::object& deleteMes,
  boost::python::dict& coarseLocations,
  boost::python::dict& highScores
)
{
  // Information to pass back about what we did, if verbosity is high enough.
  std::string infoString;

  if (verbosity >= 3) {
    std::ostringstream oss;
    oss << "   Optimizing clique of size " << movers.size() << " using recursion\n";
    infoString += oss.str();
  }

  // If we've gotten down to a clique of size 2, we terminate recursion and call our parent's method
  // because we can never split this into two connected components.
  if (movers.size() <= 2) {
    if (verbosity >= 3) {
      std::ostringstream oss;
      oss << "   Recursion terminated at clique of size " << movers.size() << "\n";
      infoString += oss.str();
    }
    boost::python::tuple ret = OptimizeCliqueCoarseBruteForceC(self, verbosity, preferenceMagnitude, movers,
      spatialQuery, extraAtomInfoMap, deleteMes, coarseLocations, highScores);
    double doubleValue = boost::python::extract<double>(ret[0]);
    std::string stringValue = boost::python::extract<std::string>(ret[1]);
    return boost::python::make_tuple(doubleValue, infoString + stringValue);
  }

  // Map from pointers to Movers (Python objects) to PositionReturn objects.
  // This must be a map because we're going to deal with subsets of Movers so
  // the indices will change.
  std::map<boost::python::object*, molprobity::reduce::PositionReturn> states;
  for (scitbx::af::shared<boost::python::object>::iterator it = movers.begin(); it != movers.end(); ++it) {
    boost::python::object& m = *it;
    states[&m] = boost::python::extract<molprobity::reduce::PositionReturn>(m.attr("CoarsePositions")());
  }

  // Construct a graph of the movers and their interactions.
  size_t nInteractions = interactions.accessor().all()[0];
  size_t nIndices = interactions.accessor().all()[1];
  if ((nInteractions > 0) && (nIndices != 2)) {
    infoString += "ERROR: OptimizeCliqueCoarseVertexCutC(): Internal error: invalid array size\n";
    return boost::python::make_tuple(-1e100, infoString);
  }
  CliqueGraph clique;
  for (scitbx::af::shared<boost::python::object>::iterator it = movers.begin(); it != movers.end(); ++it) {
    boost::python::object& m = *it;
    boost::add_vertex(&m, clique);
  }
  for (size_t i = 0; i < nInteractions; i++) {
    boost::add_edge(boost::vertex(interactions(i, 0),clique), boost::vertex(interactions(i, 1), clique), clique);
  }

  // Find a vertex cut for the graph we were given.
  std::vector<boost::python::object*> cutMovers;
  CliqueGraph cutGraph;
  findVertexCut(clique, cutMovers, cutGraph);
  if (cutMovers.size() == 0) {
    // No vertex cut found, so we call the brute-force method to solve this subgraph.
    if (verbosity >= 3) {
      std::ostringstream oss;
      oss << "   No vertex cut for clique of size " << movers.size() << ", calling parent\n";
      infoString += oss.str();
    }
    boost::python::tuple ret = OptimizeCliqueCoarseBruteForceC(self, verbosity, preferenceMagnitude, movers,
           spatialQuery, extraAtomInfoMap, deleteMes, coarseLocations, highScores);
    double doubleValue = boost::python::extract<double>(ret[0]);
    std::string stringValue = boost::python::extract<std::string>(ret[1]);
    return boost::python::make_tuple(doubleValue, infoString + stringValue);
  }

  // Run through all states of the vertex cut and for each recursively find the best score for all of
  // the connected components, followed by the score for the vertex cut. Keep track of the best state
  // and score across all of them and set back to that at the end.
  if (verbosity >= 3) {
    std::ostringstream oss;
    oss << "   Found vertex cut of size " << cutMovers.size() << "\n";
    infoString += oss.str();
  }

  // Keep track of the best score and the state where we found it. We can use a
  // vector here because we'll always be dealing with the total state and looking it
  // up each time. It starts out empty, letting us know to fill it.
  std::vector<unsigned> bestState;
  double bestScore = -1e100;

  // Find the number of options for each Mover and record it into a vector
  std::vector<unsigned> numStates;
  for (std::vector<boost::python::object*>::iterator it = cutMovers.begin(); it != cutMovers.end(); ++it) {
    boost::python::object* m = *it;
    numStates.push_back(states[m].positions.size());
  }

  // Generate a vector of state combinations, storing the index of the selected
  // element from each state. We will iterate over these to find the best state.
  /// @todo Consider pulling this generation into here as a loop rather than using up
  // memory to store this.
  std::vector< std::vector<unsigned> > allStates = generateAllStates(numStates);
  for (std::vector< std::vector<unsigned> >::const_iterator it = allStates.begin(); it != allStates.end(); ++it) {
    std::vector<unsigned> const& curStateValues = *it;

    // Set all cutMovers to match the state list
    for (unsigned m = 0; m < cutMovers.size(); m++) {
      // Only change this mover if it is different from the last time
      boost::python::extract<unsigned> value(coarseLocations.get(*cutMovers[m]));
      if (value != curStateValues[m]) {
        // Set the mover to this state
        infoString += setMoverState(states[cutMovers[m]], curStateValues[m],
          spatialQuery, extraAtomInfoMap, deleteMes, verbosity);
        coarseLocations[*cutMovers[m]] = curStateValues[m];
      }
    }

    // Recursively compute the best score across all connected components in the cutGraph.
    // This will leave each subgraph in its best state for this set of cutMovers states.
    double score = 0;
    std::vector<int> componentMap(boost::num_vertices(cutGraph));
    int numComponents = boost::connected_components(cutGraph, &componentMap[0]);
    // The components are tagged in the componentMap with integers from 0 to numComponents-1.
    for (int i = 0; i < numComponents; i++) {

      std::vector<boost::python::object*> subMovers;
      scitbx::af::shared<boost::python::object> afSubMovers;
      for (size_t m = 0; m < boost::num_vertices(cutGraph); m++) {
        if (componentMap[m] == i) {
          subMovers.push_back(cutGraph[m]);
          afSubMovers.push_back(*cutGraph[m]);
        }
      }
      CliqueGraph subGraph = subsetGraph(cutGraph, subMovers);

      // Fill in the interactions (edges) for this subgraph, which we'll need for the recursive call.
      scitbx::af::versa<int, scitbx::af::flex_grid<> > subInteractions(scitbx::af::flex_grid<>(boost::num_edges(subGraph), 2));
      boost::iterator_range<CliqueGraph::edge_iterator> edgeIterator = boost::edges(subGraph);
      size_t e = 0;
      for (CliqueGraph::edge_iterator edgeIt = edgeIterator.begin(); edgeIt != edgeIterator.end(); ++edgeIt) {
        subInteractions(e, 0) = boost::source(*edgeIt, subGraph);
        subInteractions(e, 1) = boost::target(*edgeIt, subGraph);
        e++;
      }

      // Recursively call this function to find the best score for this subgraph.
      boost::python::tuple ret = OptimizeCliqueCoarseVertexCutC(self, verbosity, preferenceMagnitude, afSubMovers,
        subInteractions, spatialQuery, extraAtomInfoMap, deleteMes, coarseLocations, highScores);
      double doubleValue = boost::python::extract<double>(ret[0]);
      std::string stringValue = boost::python::extract<std::string>(ret[1]);
      score += doubleValue;
      infoString += stringValue;
    }

    // Add the score over all atoms in the vertex-cut Movers and see if it is the best.  If so,
    // update the best.
    for (size_t i = 0; i < cutMovers.size(); i++) {
      score += preferenceMagnitude * states[cutMovers[i]].preferenceEnergies[curStateValues[i]];
      score += scorePosition(self, states[cutMovers[i]], curStateValues[i]);
    }
    if (verbosity >= 5) {
      std::ostringstream oss;
      oss << "    Cut score is " << score << " at[";
      infoString += oss.str();
      for (size_t i = 0; i < curStateValues.size(); i++) {
        std::ostringstream oss2;
        oss2 << curStateValues[i];
        infoString += oss2.str();
        if (i < curStateValues.size() - 1) {
          infoString += ", ";
        }
      }
      infoString += "]\n";
    }
    if ((score > bestScore) || (bestState.size() == 0)) {
      if (verbosity >= 4) {
        std::ostringstream oss;
        oss << "    New best score is " << score << " at [";
        infoString += oss.str();
        for (unsigned i = 0; i < curStateValues.size(); i++) {
          std::ostringstream oss2;
          oss2 << curStateValues[i];
          infoString += oss2.str();
          if (i < curStateValues.size() - 1) {
            infoString += ", ";
          }
        }
        infoString += "]\n";
      }
      bestScore = score;
      // Get the current state for all Movers in the Clique, not just the vertex-cut Movers
      bestState.clear();
      for (scitbx::af::shared<boost::python::object>::iterator it = movers.begin(); it != movers.end(); ++it) {
        boost::python::object& m = *it;
        boost::python::extract<unsigned> value(coarseLocations.get(m));
        bestState.push_back(value);
      }
    }
  }

  // Put each Mover in the entire Clique into its best state and compute its high-score value.
  // Compute the best individual scores for these Movers for use in later fine - motion
  // processing.
  for (size_t m = 0; m < movers.size(); m++) {
    infoString += setMoverState(states[&movers[m]], bestState[m], spatialQuery, extraAtomInfoMap, deleteMes,
      verbosity);
    coarseLocations[movers[m]] = bestState[m];
    double score = preferenceMagnitude * states[&movers[m]].preferenceEnergies[bestState[m]];
    score += scorePosition(self, states[&movers[m]], bestState[m]);
    highScores[movers[m]] = score;
    if (verbosity >= 3) {
      std::ostringstream oss;
      oss << "    Setting Mover in clique to coarse orientation " << bestState[m]
        << ", max score = " << score << "\n";
      infoString += oss.str();
    }
  }

  // Return the result
  return boost::python::make_tuple(bestScore, infoString);
}

// Helper functions to let us initialize vectors conveniently in C++98
static std::vector<unsigned> init_uVec(unsigned const* values, size_t size)
{
  std::vector<unsigned> ret;
  for (size_t i = 0; i < size; i++) {
    ret.push_back(values[i]);
  }
  return ret;
}

static std::vector<int> init_intVec(int const* values, size_t size)
{
  std::vector<int> ret;
  for (size_t i = 0; i < size; i++) {
    ret.push_back(values[i]);
  }
  return ret;
}

std::string Optimizers_test()
{
  // Test generateAllStates()
  unsigned ns[] = { 2, 3, 5 };
  std::vector<unsigned> numStates = init_uVec(ns, sizeof(ns)/sizeof(unsigned));
  size_t product = 1;
  for (size_t it = 0; it < numStates.size(); it++) {
    product *= numStates[it];
  }
  std::vector<unsigned> zeroes(numStates.size(), 0);
  std::vector<unsigned> lastState;
  for (size_t it = 0; it < numStates.size(); it++) {
    lastState.push_back(numStates[it] - 1);
  }
  std::vector< std::vector<unsigned> > allStates = generateAllStates(numStates);
  if (allStates.size() != product) {
    return "mmtbx_reduce_ext.Optimizers_test(): Number of states not as expected";
  }
  if (allStates.front() != zeroes) {
    return "mmtbx_reduce_ext.Optimizers_test(): First states not as expected";
  }
  if (allStates.back() != lastState) {
    return "mmtbx_reduce_ext.Optimizers_test(): Last states not as expected";
  }

  // setMoverState() is tested as part of testing OptimizeCliqueCoarseBruteForce()
  // OptimizeCliqueCoarseBruteForce() is tested from Python

  // Test nChooseM()
  std::vector < std::vector<int> > choices = nChooseM(5, 3);
  if (choices.size() != 10) {
    return "mmtbx_reduce_ext.Optimizers_test(): nChooseM(5, 3) failed";
  }
  choices = nChooseM(10, 1);
  if ((choices.size() != 10) || (0 != choices[0][0])) {
    return "mmtbx_reduce_ext.Optimizers_test(): nChooseM(10, 1) failed";
  }
  choices = nChooseM(3, 2);
  if (choices.size() != 3) {
    return "mmtbx_reduce_ext.Optimizers_test(): nChooseM(3, 2) failed";
  }

  // Test subsetGraph()
  {
    std::vector<boost::python::object> objs(10);
    std::vector< std::vector<boost::python::object*> > verts;
    /* = {
          { &objs[0],& objs[1]},
          { &objs[0],& objs[1],& objs[2] },
          { &objs[0],& objs[1],& objs[2],& objs[3] },
          { &objs[0],& objs[1],& objs[2],& objs[3] }
         }; */

    std::vector<boost::python::object*> s0;
    s0.push_back(&objs[0]); s0.push_back(&objs[1]);
    verts.push_back(s0);
    std::vector<boost::python::object*> s1;
    s1.push_back(&objs[0]); s1.push_back(&objs[1]); s1.push_back(&objs[2]);
    verts.push_back(s1);
    std::vector<boost::python::object*> s2;
    s2.push_back(&objs[0]); s2.push_back(&objs[1]); s2.push_back(&objs[2]); s2.push_back(&objs[3]);
    verts.push_back(s2);
    std::vector<boost::python::object*> s3;
    s3.push_back(&objs[0]); s3.push_back(&objs[1]); s3.push_back(&objs[2]); s3.push_back(&objs[3]);
    verts.push_back(s3);
    std::vector< std::vector< std::vector<unsigned> > > edges;
    /* = {
       { { 0, 1 } },
       { { 0, 1 }, { 1, 2 } },
       { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } },
       { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } }
     }; */
    unsigned e01[] = { 0, 1 };
    unsigned e12[] = { 1, 2 };
    unsigned e23[] = { 2, 3 };
    unsigned e30[] = { 3, 0 };
    std::vector< std::vector<unsigned> > e0;
    e0.push_back(init_uVec(e01, sizeof(e01) / sizeof(unsigned)));
    edges.push_back(e0);
    std::vector< std::vector<unsigned> > e1;
    e1.push_back(init_uVec(e01, sizeof(e01) / sizeof(unsigned)));
    e1.push_back(init_uVec(e12, sizeof(e12) / sizeof(unsigned)));
    edges.push_back(e1);
    std::vector< std::vector<unsigned> > e2;
    e2.push_back(init_uVec(e01, sizeof(e01) / sizeof(unsigned)));
    e2.push_back(init_uVec(e12, sizeof(e12) / sizeof(unsigned)));
    e2.push_back(init_uVec(e23, sizeof(e23) / sizeof(unsigned)));
    e2.push_back(init_uVec(e30, sizeof(e30) / sizeof(unsigned)));
    edges.push_back(e2);
    std::vector< std::vector<unsigned> > e3;
    e3.push_back(init_uVec(e01, sizeof(e01) / sizeof(unsigned)));
    e3.push_back(init_uVec(e12, sizeof(e12) / sizeof(unsigned)));
    e3.push_back(init_uVec(e23, sizeof(e23) / sizeof(unsigned)));
    e3.push_back(init_uVec(e30, sizeof(e30) / sizeof(unsigned)));
    edges.push_back(e3);
    unsigned ki0[] = { 0 };
    unsigned ki1[] = { 0, 1 };
    unsigned ki2[] = { 0, 1, 2, 3 };
    unsigned ki3[] = { 0, 2 };
    std::vector< std::vector<unsigned> > keepIndices;
    keepIndices.push_back(init_uVec(ki0, sizeof(ki0) / sizeof(unsigned)));
    keepIndices.push_back(init_uVec(ki1, sizeof(ki1) / sizeof(unsigned)));
    keepIndices.push_back(init_uVec(ki2, sizeof(ki2) / sizeof(unsigned)));
    keepIndices.push_back(init_uVec(ki3, sizeof(ki3) / sizeof(unsigned)));
    int evc[] = { 1, 2, 4, 2 };
    std::vector<int> expectedVertexCounts = init_intVec(evc, sizeof(evc)/sizeof(int));
    int eec[] = { 0, 1, 4, 0 };
    std::vector<int> expectedEdgeCounts = init_intVec(eec, sizeof(eec) / sizeof(int));

    for (size_t i = 0; i < verts.size(); i++) {
      CliqueGraph g;
      for (std::vector<boost::python::object*>::iterator it = verts[i].begin(); it != verts[i].end(); ++it) {
        boost::add_vertex(*it, g);
      }
      for (std::vector< std::vector<unsigned> >::const_iterator it = edges[i].begin(); it != edges[i].end(); ++it) {
        boost::add_edge(boost::vertex((*it)[0], g), boost::vertex((*it)[1], g), g);
      }
      std::vector<boost::python::object*> keepers;
      for (size_t v = 0; v < keepIndices[i].size(); v++) {
        keepers.push_back(g[keepIndices[i][v]]);
      }
      CliqueGraph subset = subsetGraph(g, keepers);
      if (boost::num_vertices(subset) != expectedVertexCounts[i]) {
        return "mmtbx_reduce_ext.Optimizers_test(): subsetGraph() failed with unexpected vertex counts";
      }
      if (boost::num_edges(subset) != expectedEdgeCounts[i]) {
        return "mmtbx_reduce_ext.Optimizers_test(): subsetGraph() failed with unexpected edge counts";
      }
    }
  }

  // Test findVertexCut()
  {
    std::vector<boost::python::object> objs(10);
    std::vector< std::vector<boost::python::object*> > verts;
    /* = {
       { &objs[0], &objs[1]},
       { &objs[0], &objs[1], &objs[2] },
       { &objs[0], &objs[1], &objs[2], &objs[3] }
     }; */
    std::vector<boost::python::object*> s0;
    s0.push_back(&objs[0]); s0.push_back(&objs[1]);
    verts.push_back(s0);
    std::vector<boost::python::object*> s1;
    s1.push_back(&objs[0]); s1.push_back(&objs[1]); s1.push_back(&objs[2]);
    verts.push_back(s1);
    std::vector<boost::python::object*> s2;
    s2.push_back(&objs[0]); s2.push_back(&objs[1]); s2.push_back(&objs[2]); s2.push_back(&objs[3]);
    verts.push_back(s2);
    std::vector< std::vector< std::vector<unsigned> > > edges;
    /*= {
      { { 0, 1 } },
      { { 0, 1 }, { 1, 2 } },
      { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } }
    };*/
    unsigned e01[] = { 0, 1 };
    unsigned e12[] = { 1, 2 };
    unsigned e23[] = { 2, 3 };
    unsigned e30[] = { 3, 0 };
    std::vector< std::vector<unsigned> > e0;
    e0.push_back(init_uVec(e01, sizeof(e01) / sizeof(unsigned)));
    edges.push_back(e0);
    std::vector< std::vector<unsigned> > e1;
    e1.push_back(init_uVec(e01, sizeof(e01) / sizeof(unsigned)));
    e1.push_back(init_uVec(e12, sizeof(e12) / sizeof(unsigned)));
    edges.push_back(e1);
    std::vector< std::vector<unsigned> > e2;
    e2.push_back(init_uVec(e01, sizeof(e01) / sizeof(unsigned)));
    e2.push_back(init_uVec(e12, sizeof(e12) / sizeof(unsigned)));
    e2.push_back(init_uVec(e23, sizeof(e23) / sizeof(unsigned)));
    e2.push_back(init_uVec(e30, sizeof(e30) / sizeof(unsigned)));
    edges.push_back(e2);
    int evc[] = { 0, 1, 2 };
    std::vector<int> expectedCutCounts = init_intVec(evc, sizeof(evc) / sizeof(int));
    int eec[] = { 2, 2, 2 };
    std::vector<int> expectedVertexCounts = init_intVec(eec, sizeof(eec) / sizeof(int));
    for (size_t i = 0; i < verts.size(); i++) {
      CliqueGraph g;
      for (std::vector<boost::python::object*>::iterator it = verts[i].begin(); it != verts[i].end(); ++it) {
        boost::add_vertex(*it, g);
      }
      for (std::vector< std::vector<unsigned> >::const_iterator it = edges[i].begin(); it != edges[i].end(); ++it) {
        boost::add_edge(boost::vertex((*it)[0], g), boost::vertex((*it)[1], g), g);
      }
      std::vector<boost::python::object*> cut;
      CliqueGraph cutGraph;
      findVertexCut(g, cut, cutGraph);
      if (boost::num_vertices(cutGraph) != expectedVertexCounts[i]) {
        return "mmtbx_reduce_ext.Optimizers_test(): findVertexCut() failed with unexpected vertex counts";
      }
      if (cut.size() != expectedCutCounts[i]) {
        return "mmtbx_reduce_ext.Optimizers_test(): findVertexCut() failed with unexpected cut sizes";
      }
    }
  }

  // OptimizeCliqueCoarseVertexCutC() is tested by the Python code.

  /// @todo

  // All tests passed.
  return "";
}

  } // end namespace reduce
} // end namespace molprobity
