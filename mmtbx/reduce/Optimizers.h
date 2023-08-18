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

// Contains functions and classes to support the Optimizers.py classes.
// They are linked into the mmtbx_reduce_ext library by boost_python/reduce_bpl.cpp.

#pragma once

#include <string>
#include <utility>
#include <boost/python.hpp>
#include "../probe/Scoring.h"
#include "PositionReturn.h"

namespace molprobity {
  namespace reduce {

    /** @brief Function to perform brute-force optimization on a clique of Movers.
        @param [in] verbosity: Controls how much information is added to the string.
        @param [in] preferenceMagnitude: Multiples the preference energies, so that we
                can scale down their importance if we want.
        @param [in] movers: A list of Movers to jointly optimize.
        @param [in] states: List of responses to calling CoarsePositions() on
                each of the movers, with each entry corresponding to its mover. These are
                easier to fill in via Python than to determine here.
        @param [inOut] spatialQuery: Spatial-query structure telling which atoms are where
        @param [inOut] extraAtomInfoMap: Map containing extra information about each atom.
        @param [inOut] deleteMes: Set of atoms to be deleted, passed as a Python object.
        @param [out] coarseLocations: Dictionary looked up by mover that records the
                coarse locations of the movers once they are optimized. These are modified
                to point to the result.
        @param [out] highScores: Dictionary looked up by mover with the scores at best locations.
        @return A tuple, where the first is the score at the best position for all movers
                and the second is a string describing what was done, which may be empty
                if verbosity is too small.
    */
    boost::python::tuple OptimizeCliqueCoarseBruteForceC(
      boost::python::object &self,
      int verbosity,
      double preferenceMagnitude,
      scitbx::af::shared<boost::python::object> movers,
      scitbx::af::shared<molprobity::reduce::PositionReturn> states,
      molprobity::probe::SpatialQuery &spatialQuery,
      molprobity::probe::ExtraAtomInfoMap &extraAtomInfoMap,
      boost::python::object &deleteMes,
      boost::python::dict &coarseLocations,
      boost::python::dict &highScores
    );

    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string Optimizers_test();
  }
}
