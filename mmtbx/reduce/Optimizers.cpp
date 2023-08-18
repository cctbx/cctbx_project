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

#include "InteractionGraph.h"
#include "PositionReturn.h"
#include "../probe/SpatialQuery.h"
#include "../probe/Scoring.h"
#include <vector>
#include <string>
#include <chrono>

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
  // Fill in a vector with the states of each mover.
  scitbx::af::shared<molprobity::reduce::PositionReturn> states;
  for (auto m : movers) {
    states.push_back(boost::python::extract<molprobity::reduce::PositionReturn>(m.attr("CoarsePositions")()));
  }

  // Keep track of the best score and the state where we found it.
  double bestScore = -1e100;
  std::vector<unsigned> bestState;

  std::string infoString;

  // Find the length of each state and record it into a vector
  std::vector<unsigned> numStates;
  for (auto const& s : states) {
    numStates.push_back(s.positions.size());
  }

  // Generate a vector of state combinations, storing the index of the selected
  // element from each state. We will iterate over these to find the best state.
  /// @todo Consider pulling this generation into here as a loop rather than using up
  // memory to store this.
  std::vector< std::vector<unsigned> > allStates = generateAllStates(numStates);

  // Cycle through all states and find the one with the largest score.
  for (auto const& curStateValues : allStates) {

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
      if (verbosity >= 5) {
        infoString += "      score is " + std::to_string(score) + " at [";
        for (unsigned i = 0; i < curStateValues.size(); i++) {
          infoString += std::to_string(curStateValues[i]);
          if (i < curStateValues.size() - 1) {
            infoString += ", ";
          }
        }
        infoString += "]\n";
      }
    }

    // See if it is the best score.  If so, update the best.
    if ((score > bestScore) || (bestState.size() == 0)) {
      if (verbosity >= 4) {
        infoString += "    New best score is " + std::to_string(score) + " at [";
        for (unsigned i = 0; i < curStateValues.size(); i++) {
          infoString += std::to_string(curStateValues[i]);
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
      infoString += "    Setting Mover in clique to coarse orientation " + std::to_string(bestState[m])
        + ", max score = " + std::to_string(score) + "\n";
    }
  }

  // Return the result
  return boost::python::make_tuple(bestScore, infoString);
}

std::string Optimizers_test()
{
  // Test generateAllStates()
  std::vector<unsigned> numStates = { 2, 3, 5 };
  size_t product = 1;
  for (auto n : numStates) {
    product *= n;
  }
  std::vector<unsigned> zeroes(numStates.size(), 0);
  std::vector<unsigned> lastState;
  for (auto n : numStates) {
    lastState.push_back(n - 1);
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

  /// @todo

  // All tests passed.
  return "";
}

  } // end namespace reduce
} // end namespace molprobity
