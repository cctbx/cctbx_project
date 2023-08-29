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
#include "../probe/DotSpheres.h"
#include <list>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/subgraph.hpp>

static std::string roundToTwoDigits(double d)
{
  std::ostringstream oss;
  oss.precision(2);
  oss << std::fixed << d;
  return oss.str();
}

static std::string stripWhitespace(const std::string& input) {
  std::string result = input;

  // Remove leading whitespace
  result.erase(result.begin(), std::find_if(result.begin(), result.end(), [](unsigned char ch) {
    return !std::isspace(ch);
    }));

  // Remove trailing whitespace
  result.erase(std::find_if(result.rbegin(), result.rend(), [](unsigned char ch) {
    return !std::isspace(ch);
    }).base(), result.end());

  return result;
}

static std::string toUpperCase(const std::string& input) {
  std::string result = input;
  std::transform(result.begin(), result.end(), result.begin(), ::toupper);
  return result;
}

static std::string fromSmall(const char* buf, unsigned maxLen)
{
  std::string ret;

  // Scan through the buffer. As long as we do not find null termination, continue
  // copying into the string.
  for (unsigned i = 0; i < maxLen; i++) {
    if (buf[i] == '\0') {
      break;
    }
    ret += buf[i];
  }
  return ret;
}

static std::string resNameAndID(iotbx::pdb::hierarchy::atom const& a)
{
  std::string chainID = a.parent().get().parent().get().parent().get().data->id;
  std::string resName = toUpperCase(stripWhitespace(a.parent().get().data->resname));
  std::string resID = stripWhitespace(fromSmall(a.parent().get().parent().get().data->resseq.elems, 4));
  std::string altLoc = stripWhitespace(fromSmall(a.parent().get().data->altloc.elems,1));
  // Don't print the code if it is a space (blank).
  std::string insertionCode = stripWhitespace(fromSmall(a.parent().get().parent().get().data->icode.elems,1));
  return "chain " + chainID + " " + altLoc + resName + " " + resID + insertionCode;
}

static std::string describeMover(boost::python::object const& mover,
  iotbx::pdb::hierarchy::atom const& atom)
{
  // Extract the type name between the quotes and then the last class name from
  // the module path.
  boost::python::object objectType = mover.attr("__class__");
  std::string typeName =
    boost::python::extract<std::string>(objectType.attr("__name__"));

  return typeName + ' ' + resNameAndID(atom);
}

namespace molprobity {
  namespace reduce {

std::vector< std::vector<unsigned> > OptimizerC::generateAllStates(
  std::vector<unsigned> const &numStates)
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

std::vector< std::vector<int> > OptimizerC::nChooseM(int n, int m)
{
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

OptimizerC::CliqueGraph OptimizerC::subsetGraph(CliqueGraph const& graph,
  std::vector<boost::python::object*>& keepMovers)
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
    if ((std::find(keepMovers.begin(), keepMovers.end(), mover1) != keepMovers.end()) &&
      (std::find(keepMovers.begin(), keepMovers.end(), mover2) != keepMovers.end())) {
      // Both ends of this edge are in the subset, so add it to the subset graph.
      boost::add_edge(vertexMap[mover1], vertexMap[mover2], ret);
    }
  }

  return ret;
}

void OptimizerC::findVertexCut(CliqueGraph const& graph,
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

OptimizerC::OptimizerC(boost::python::object& self, int verbosity, double preferenceMagnitude,
  double maxVDWRadius,
  double minOccupancy,
  double probeRadius,
  double probeDensity,
  boost::python::dict& exclude,
  boost::python::dict& dotSpheres,
  boost::python::dict& atomMoverLists,
  molprobity::probe::SpatialQuery& spatialQuery,
  molprobity::probe::ExtraAtomInfoMap& extraAtomInfoMap,
  boost::python::object& deleteMes,
  boost::python::dict& coarseLocations,
  boost::python::dict& fineLocations,
  boost::python::dict& highScores)
  : m_self(self)
  , m_verbosity(verbosity)
  , m_preferenceMagnitude(preferenceMagnitude)
  , m_maxVDWRadius(maxVDWRadius)
  , m_minOccupancy(minOccupancy)
  , m_probeRadius(probeRadius)
  , m_probeDensity(probeDensity)
  , m_exclude(exclude)
  , m_dotSpheres(dotSpheres)
  , m_atomMoverLists(atomMoverLists)
  , m_spatialQuery(spatialQuery)
  , m_extraAtomInfoMap(extraAtomInfoMap)
  , m_deleteMes(deleteMes)
  , m_coarseLocations(coarseLocations)
  , m_fineLocations(fineLocations)
  , m_highScores(highScores)
  , m_cachedScores(0)
  , m_calculatedScores(0)
{
  // Look up the self._dotScorer object and store a pointer to it.
  boost::python::object dotScorerObj = m_self.attr("_dotScorer");
  molprobity::probe::DotScorer& temp = boost::python::extract<molprobity::probe::DotScorer&>(dotScorerObj);
  m_dotScorer = &temp;

  /// @todo We can look up a lot of the parameters above from the self object.
}

double OptimizerC::scoreAtom(iotbx::pdb::hierarchy::atom const& a)
{
  // Keep track of calculated scores whether or not we're doing caching
  m_calculatedScores++;

  // Find the maximum radius of the pair of atoms, to limit our neightbor search.
  double maxRadiusWithoutProbe = m_maxVDWRadius + m_extraAtomInfoMap.getMappingFor(a).getVdwRadius();

  // Find the excluded atoms for this atom. This is a dictionary looked up by i_seq that has a list of atoms.
  /// @todo We'd like to do this without a copy but we can't get a reference.
  /// @todo Consider building a C++ map for all of the atoms in the current clique and using it.
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude =
    boost::python::extract<scitbx::af::shared<iotbx::pdb::hierarchy::atom> >(m_exclude[a.data->i_seq]);

  // Find the dots for this atom. This is a dictionary looked up by i_seq that returns a DotSphere.
  /// @todo Consider building a C++ map for all of the atoms in the current clique and using it.
  molprobity::probe::DotSphere& ds =
    boost::python::extract<molprobity::probe::DotSphere&>(m_dotSpheres[a.data->i_seq]);

  // Score the dots for this atom.
  return m_dotScorer->score_dots(a, m_minOccupancy, m_spatialQuery, maxRadiusWithoutProbe,
    m_probeRadius, exclude, ds.dots(), m_probeDensity, false).totalScore();
}

double OptimizerC::scoreAtomCached(iotbx::pdb::hierarchy::atom const& a)
{
  // Find the vector that is the state of all Movers related to this atom.
  /// @todo We may want to cache this vector as well.
  std::vector<unsigned> state;
  boost::python::list listObj = boost::python::extract<boost::python::list>(m_atomMoverLists[a.data->i_seq])();
  for (int i = 0; i < boost::python::len(listObj); ++i) {
    // Find each mover.
    boost::python::object moverObj = listObj[i];
    // Go through and look up the state of each mover.
    boost::python::extract<unsigned> value(m_coarseLocations.get(moverObj));
    state.push_back(value);
  }

  // See if we have a cached score for this state.  If so, use it. If not, calculate it and store it.
  if ((*m_scoreCacheMap)[a.data->i_seq].find(state) != (*m_scoreCacheMap)[a.data->i_seq].end()) {
    m_cachedScores++;
    return (*m_scoreCacheMap)[a.data->i_seq][state];
  } else {
    double score = scoreAtom(a);
    (*m_scoreCacheMap)[a.data->i_seq][state] = score;
    return score;
  }
}

double OptimizerC::scorePosition(molprobity::reduce::PositionReturn& states, size_t index)
{
  double ret = 0;
  for (size_t a = 0; a < states.atoms.size(); a++) {
    // There may not be as many deleteMes as there are atoms, so we need to check for that.
    if ((a >= states.deleteMes[index].size()) || !states.deleteMes[index][a]) {
      if (m_scoreCacheMap) {
        ret += scoreAtomCached(states.atoms[a]);
      } else {
        ret += scoreAtom(states.atoms[a]);
      }
    }
  }
  return ret;
}

std::string OptimizerC::setMoverState(molprobity::reduce::PositionReturn& positionReturn,
  unsigned index)
{
  std::string ret;

  // Move the atoms to their new positions, updating the spatial query structure
  // by removing the old and adding the new location.
  for (size_t i = 0; i < positionReturn.atoms.size(); i++) {
    iotbx::pdb::hierarchy::atom& a = positionReturn.atoms[i];

    m_spatialQuery.remove(a);
    // Overwrite the location of the atom
    for (size_t j = 0; j < 3; j++) {
      a.data->xyz[j] = positionReturn.positions[index][i][j];
    }
    m_spatialQuery.add(a);
  }

  // Update the extraAtomInfo associated with each atom.
  // Note that there may be fewer entries than atoms, but they all correspond
  // so we can look up the atom by index in the atoms array just like above.
  for (size_t i = 0; i < positionReturn.extraInfos[index].size(); i++) {
    m_extraAtomInfoMap.setMappingFor(positionReturn.atoms[i],
      positionReturn.extraInfos[index][i]);
  }

  // Manage the deletion status of each atom, including ensuring
  // consistency with the spatial - query structure.
  // Note that there may be fewer entries than atoms, but they all correspond
  // so we can look up the atom by index in the atoms array just like above.
  for (size_t i = 0; i < positionReturn.deleteMes[index].size(); i++) {
    iotbx::pdb::hierarchy::atom& a = positionReturn.atoms[i];
    bool doDelete = positionReturn.deleteMes[index][i];
    if (doDelete) {
      m_spatialQuery.remove(a);
      m_deleteMes.attr("add")(a);
      if (m_verbosity >= 10) {
        ret += "          Deleting atom\n";
      }
    } else {
      m_spatialQuery.add(a);
      m_deleteMes.attr("discard")(a);
      if (m_verbosity >= 10) {
        ret += "          Ensuring deletable atom is present\n";
      }
    }
  }

  return ret;
}

std::string OptimizerC::Initialize(scitbx::af::shared<boost::python::object> movers)
{
  std::string infoString;

  for (size_t i = 0; i < movers.size(); i++) {
    boost::python::object const& mover = movers[i];

    molprobity::reduce::PositionReturn coarse =
      boost::python::extract<molprobity::reduce::PositionReturn>(mover.attr("CoarsePositions")());
    double score = m_preferenceMagnitude * coarse.preferenceEnergies[0];
    setMoverState(coarse, 0);
    score += scorePosition(coarse, 0);
    m_coarseLocations[mover] = 0;
    m_highScores[mover] = score;
  }

  return infoString;
}

std::pair<double, std::string> OptimizerC::OptimizeCliqueCoarseBruteForce(
  std::map<boost::python::object*, molprobity::reduce::PositionReturn>& states,
  CliqueGraph& clique)
{
  // Information to pass back about what we did, if verbosity is high enough.
  std::string infoString;

  // Get the list of movers from the graph
  /// @todo Speed this up by just getting a vertex iterator and using it?
  std::vector<boost::python::object*> movers;
  for (size_t i = 0; i < boost::num_vertices(clique); i++) {
    movers.push_back(clique[i]);
  }

  if (m_verbosity >= 3) {
    std::ostringstream oss;
    oss << "   Optimizing clique of size " << movers.size() << " using brute force\n";
    infoString = oss.str();
  }

  // Keep track of the best score and the state where we found it.
  // It starts out empty, letting us know to fill it.
  double bestScore = -1e100;
  std::vector<unsigned> bestState;

  // Find the number of options for each Mover and record it into a vector
  std::vector<unsigned> numStates;
  for (std::vector<boost::python::object*>::iterator it = movers.begin(); it != movers.end(); ++it) {
    boost::python::object* m = *it;
    numStates.push_back(states[m].positions.size());
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
      boost::python::extract<unsigned> value(m_coarseLocations.get(*movers[m]));
      if (value != curStateValues[m]) {
        // Set the mover to this state
        infoString += setMoverState(states[movers[m]], curStateValues[m]);
        m_coarseLocations[*movers[m]] = curStateValues[m];
      }
    }

    // Compute the score over all atoms in all Movers.
    double score = 0;
    for (size_t m = 0; m < movers.size(); m++) {

      score += m_preferenceMagnitude * states[movers[m]].preferenceEnergies[curStateValues[m]];

      // Score all atoms that have not been marked for deletion, calling the Python object's
      // scoring function (which may or may not use caching to do the scoring).
      score += scorePosition(states[movers[m]], curStateValues[m]);
    }
    if (m_verbosity >= 5) {
      std::ostringstream oss;
      oss << "     Score is " << roundToTwoDigits(score) << " at [";
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
      if (m_verbosity >= 4) {
        std::ostringstream oss;
        oss << "    New best score is " << roundToTwoDigits(score) << " at [";
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
  double ret = 0.0;
  for (size_t m = 0; m < movers.size(); m++) {
    infoString += setMoverState(states[movers[m]], bestState[m]);
    m_coarseLocations[*movers[m]] = bestState[m];
    double myScore = m_preferenceMagnitude * states[movers[m]].preferenceEnergies[bestState[m]];
    myScore += scorePosition(states[movers[m]], bestState[m]);
    m_highScores[*movers[m]] = myScore;
    ret += myScore;
    if (m_verbosity >= 3) {
      std::ostringstream oss;
      oss << "   Setting Mover in clique to coarse orientation " << bestState[m]
        << ", max score = " << roundToTwoDigits(myScore) << "\n";
      infoString += oss.str();
    }
  }

  // Return the result
  return std::pair<double, std::string>(ret, infoString);
}

std::pair<double, std::string> OptimizerC::OptimizeCliqueCoarseVertexCut(
  std::map<boost::python::object*, molprobity::reduce::PositionReturn>& states,
  CliqueGraph& clique)
{
  // Information to pass back about what we did, if verbosity is high enough.
  std::string infoString;

  // Get the list of movers from the graph
  /// @todo Speed this up by just getting a vertex iterator and using it?
  std::vector<boost::python::object*> movers;
  for (size_t i = 0; i < boost::num_vertices(clique); i++) {
    movers.push_back(clique[i]);
  }

  if (m_verbosity >= 3) {
    std::ostringstream oss;
    oss << "   Optimizing clique of size " << movers.size() << " using recursion\n";
    infoString += oss.str();
  }

  // If we've gotten down to a clique of size 2, we terminate recursion and call the brute-force
  // because we can never split this into two connected components.
  if (movers.size() <= 2) {
    if (m_verbosity >= 3) {
      std::ostringstream oss;
      oss << "   Recursion terminated at clique of size " << movers.size() << "\n";
      infoString += oss.str();
    }
    std::pair<double, std::string> ret = OptimizeCliqueCoarseBruteForce(states, clique);
    return std::pair<double, std::string>(ret.first, infoString + ret.second);
  }

  // Find a vertex cut for the clique we were given.
  std::vector<boost::python::object*> cutMovers;
  CliqueGraph cutGraph;
  findVertexCut(clique, cutMovers, cutGraph);
  if (cutMovers.size() == 0) {
    // No vertex cut found, so we call the brute-force method to solve this subgraph.
    if (m_verbosity >= 3) {
      std::ostringstream oss;
      oss << "   No vertex cut for clique of size " << movers.size() << ", calling parent\n";
      infoString += oss.str();
    }
    std::pair<double, std::string> ret = OptimizeCliqueCoarseBruteForce(states, clique);
    return std::pair<double, std::string>(ret.first, infoString + ret.second);
  }

  // Run through all states of the vertex cut and for each recursively find the best score for all of
  // the connected components, followed by the score for the vertex cut. Keep track of the best state
  // and score across all of them and set back to that at the end.
  if (m_verbosity >= 3) {
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
      boost::python::extract<unsigned> value(m_coarseLocations.get(*cutMovers[m]));
      if (value != curStateValues[m]) {
        // Set the mover to this state
        infoString += setMoverState(states[cutMovers[m]], curStateValues[m]);
        m_coarseLocations[*cutMovers[m]] = curStateValues[m];
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
      for (size_t m = 0; m < boost::num_vertices(cutGraph); m++) {
        if (componentMap[m] == i) {
          subMovers.push_back(cutGraph[m]);
        }
      }
      CliqueGraph subGraph = subsetGraph(cutGraph, subMovers);

      // Recursively call this function to find the best score for this subgraph.
      std::pair<double, std::string> ret = OptimizeCliqueCoarseVertexCut(states, subGraph);
      score += ret.first;
      infoString += ret.second;
    }

    // Add the score over all atoms in the vertex-cut Movers and see if it is the best.  If so,
    // update the best.
    for (size_t i = 0; i < cutMovers.size(); i++) {
      score += m_preferenceMagnitude * states[cutMovers[i]].preferenceEnergies[curStateValues[i]];
      score += scorePosition(states[cutMovers[i]], curStateValues[i]);
    }
    if (m_verbosity >= 5) {
      std::ostringstream oss;
      oss << "     Cut score is " << roundToTwoDigits(score) << " at[";
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
      if (m_verbosity >= 4) {
        std::ostringstream oss;
        oss << "    New best score is " << roundToTwoDigits(score) << " at [";
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
      for (std::vector<boost::python::object*>::iterator it = movers.begin(); it != movers.end(); ++it) {
        boost::python::object* m = *it;
        boost::python::extract<unsigned> value(m_coarseLocations.get(*m));
        bestState.push_back(value);
      }
    }
  }

  // Put each Mover in the entire Clique into its best state and compute its high-score value.
  // Compute the best individual scores for these Movers for use in later fine - motion
  // processing.
  double ret = 0.0;
  for (size_t m = 0; m < movers.size(); m++) {
    infoString += setMoverState(states[movers[m]], bestState[m]);
    m_coarseLocations[*movers[m]] = bestState[m];
    double score = m_preferenceMagnitude * states[movers[m]].preferenceEnergies[bestState[m]];
    score += scorePosition(states[movers[m]], bestState[m]);
    m_highScores[*movers[m]] = score;
    ret += score;
    if (m_verbosity >= 3) {
      std::ostringstream oss;
      oss << "   Setting Mover in clique to coarse orientation " << bestState[m]
        << ", max score = " << roundToTwoDigits(score) << "\n";
      infoString += oss.str();
    }
  }

  // Return the result
  return std::pair<double, std::string>(ret, infoString);
}

boost::python::tuple OptimizerC::OptimizeSingleMoverCoarse(boost::python::object const& mover)
{
  // Information to pass back about what we did, if verbosity is high enough.
  std::string infoString;

  molprobity::reduce::PositionReturn coarse =
    boost::python::extract<molprobity::reduce::PositionReturn>(mover.attr("CoarsePositions")());
  std::vector<double> scores = coarse.preferenceEnergies;
  for (size_t i = 0; i < scores.size(); i++) {
    scores[i] *= m_preferenceMagnitude;
  }
  for (size_t i = 0; i < coarse.positions.size(); i++) {
    infoString += setMoverState(coarse, i);
    scores[i] += scorePosition(coarse, i);
    if (m_verbosity >= 5) {
      std::ostringstream oss;
      oss << "     Single Mover " << describeMover(mover, coarse.atoms[0]) << " score at orientation " << i
        << " = " << roundToTwoDigits(scores[i]) << "\n";
      infoString += oss.str();
    }
  }

  // Find the maximum score, keeping track of the index of the maximum score.
  double maxScore = scores[0];
  unsigned maxIndex = 0;
  for (size_t i = 1; i < scores.size(); i++) {
    if (scores[i] > maxScore) {
      maxScore = scores[i];
      maxIndex = i;
    }
  }

  // Put the Mover into its final position (which may be the same as its initial position).
  if (m_verbosity >= 3) {
    std::ostringstream oss;
    oss << "   Setting single Mover to coarse orientation " << maxIndex
      << ", max score = " << roundToTwoDigits(maxScore)
      << " (initial score " << roundToTwoDigits(scores[0])
      << ")\n";
    infoString += oss.str();
  }
  infoString += setMoverState(coarse, maxIndex);
  m_coarseLocations[mover] = maxIndex;

  // Record and return the best score for this Mover.
  m_highScores[mover] = maxScore;
  return boost::python::make_tuple(maxScore, infoString);
}

boost::python::tuple OptimizerC::OptimizeSingleMoverFine(boost::python::object const& mover)
{
  std::string infoString;

  // Record this in case we need to put it back.
  molprobity::reduce::PositionReturn coarse =
    boost::python::extract<molprobity::reduce::PositionReturn>(mover.attr("CoarsePositions")());

  boost::python::extract<double> initialScore(m_highScores.get(mover));
  double maxScore = initialScore;
  boost::python::extract<unsigned> cl(m_coarseLocations.get(mover));
  unsigned coarseLoc = cl;
  molprobity::reduce::PositionReturn fine =
    boost::python::extract<molprobity::reduce::PositionReturn>(mover.attr("FinePositions")(coarseLoc));
  if (fine.positions.size() > 0) {

    std::vector<double> scores = fine.preferenceEnergies;
    for (size_t i = 0; i < scores.size(); i++) {
      scores[i] *= m_preferenceMagnitude;
    }

    for (size_t i = 0; i < fine.positions.size(); i++) {
      infoString += setMoverState(fine, i);
      scores[i] += scorePosition(fine, i);
      if (m_verbosity >= 5) {
        std::ostringstream oss;
        oss << "     Mover score at fine orientation " << i << " = " << roundToTwoDigits(scores[i]) << "\n";
        infoString += oss.str();
      }
    }

    // Find the maximum score, keeping track of the index of the maximum score.
    maxScore = scores[0];
    unsigned maxIndex = 0;
    for (size_t i = 1; i < scores.size(); i++) {
      if (scores[i] > maxScore) {
        maxScore = scores[i];
        maxIndex = i;
      }
    }

    // Put the Mover into its final position (which may be back to its initial position)
    // and update the high score.
    if (maxScore > m_highScores[mover]) {
      m_fineLocations[mover] = maxIndex;
      if (m_verbosity >= 3) {
        std::ostringstream oss;
        oss << "   Setting Mover to fine orientation " << maxIndex
          << ", max score = " << roundToTwoDigits(maxScore)
          << " (coarse score " << roundToTwoDigits(initialScore)
          << ")\n";
        infoString += oss.str();
      }
      setMoverState(fine, maxIndex);

      // Record the high score for this Mover.
      m_highScores[mover] = maxScore;
    } else {
      // Put us back into the initial coarse location and don't change the high score
      setMoverState(coarse, coarseLoc);
      if (m_verbosity >= 3) {
        infoString += "   Leaving Mover at coarse orientation\n";
      }
      // Leave the m_fineLocations at its original None value.
    }
  }

  // Return the result
  return boost::python::make_tuple(maxScore, infoString);
}

boost::python::tuple OptimizerC::OptimizeCliqueCoarse(
  scitbx::af::shared<boost::python::object> movers,
  scitbx::af::versa<int, scitbx::af::flex_grid<> >& interactions)
{
  // Information to pass back about what we did, if verbosity is high enough.
  std::string infoString;

  if (m_verbosity >= 3) {
    std::ostringstream oss;
    oss << "   Optimizing clique of size " << movers.size()
      << " using atom-score cache\n";
    infoString += oss.str();
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
    boost::add_edge(boost::vertex(interactions(i, 0), clique), boost::vertex(interactions(i, 1), clique), clique);
  }

  // Construct the ScoreCacheMap for the clique before calling and then remove it when done.
  // This will mean that we only use score caching, and only on our new map, for this clique.
  m_scoreCacheMap = new ScoreCacheMap();
  std::pair<double, std::string> ret = OptimizeCliqueCoarseVertexCut(states, clique);
  delete m_scoreCacheMap;
  m_scoreCacheMap = nullptr;
  infoString += ret.second;

  // Format and return the result
  return boost::python::make_tuple(ret.first, infoString);
}

// Test-generation helper functions to let us initialize vectors conveniently in C++98
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

std::string OptimizerC::Test()
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
    return "mmtbx_reduce_ext.OptimizerC::Test(): Number of states not as expected";
  }
  if (allStates.front() != zeroes) {
    return "mmtbx_reduce_ext.OptimizerC::Test(): First states not as expected";
  }
  if (allStates.back() != lastState) {
    return "mmtbx_reduce_ext.OptimizerC::Test(): Last states not as expected";
  }

  // setMoverState() is tested as part of testing OptimizeCliqueCoarseBruteForce()
  // OptimizeCliqueCoarseBruteForce() is tested from Python

  // Test nChooseM()
  std::vector < std::vector<int> > choices = nChooseM(5, 3);
  if (choices.size() != 10) {
    return "mmtbx_reduce_ext.OptimizerC::Test(): nChooseM(5, 3) failed";
  }
  choices = nChooseM(10, 1);
  if ((choices.size() != 10) || (0 != choices[0][0])) {
    return "mmtbx_reduce_ext.OptimizerC::Test(): nChooseM(10, 1) failed";
  }
  choices = nChooseM(3, 2);
  if (choices.size() != 3) {
    return "mmtbx_reduce_ext.OptimizerC::Test(): nChooseM(3, 2) failed";
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
        return "mmtbx_reduce_ext.OptimizerC::Test(): subsetGraph() failed with unexpected vertex counts";
      }
      if (boost::num_edges(subset) != expectedEdgeCounts[i]) {
        return "mmtbx_reduce_ext.OptimizerC::Test(): subsetGraph() failed with unexpected edge counts";
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
        return "mmtbx_reduce_ext.OptimizerC::Test(): findVertexCut() failed with unexpected vertex counts";
      }
      if (cut.size() != expectedCutCounts[i]) {
        return "mmtbx_reduce_ext.OptimizerC::Test(): findVertexCut() failed with unexpected cut sizes";
      }
    }
  }

  // OptimizeCliqueCoarseVertexCutC() is tested by the Python code.

  /// @todo

  // All tests passed.
  return "";
}

std::string Optimizers_test()
{
  return OptimizerC::Test();
}

  } // end namespace reduce
} // end namespace molprobity
