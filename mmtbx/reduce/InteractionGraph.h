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

#pragma once

#include <string>
#include <boost/python.hpp>
#include "../../boost_adaptbx/graph/graph_type.hpp"
#include "../probe/Scoring.h"

namespace molprobity {
  namespace reduce {

    /** @brief Function to report whether any atoms overlap between two Movers.

        Helper function that tells whether any pair of atoms from two Movers overlap.
        @param mover1 : The first Mover
        @param atoms1 : Atom list for the first Mover
        @param positions1 : probe.PositionReturn.positions holding possible positions for each.
        @param mover2 : The second Mover
        @param atoms2 : Atom list for the second Mover
        @param positions2 : probe.PositionReturn.positions holding possible positions for each.
        @param extraAtomInfoMap : probe.ExtraAtomInfoMap that can be used to look
          up the information for atoms whose values need to be changed.Can be
          obtained by calling mmtbx.probe.Helpers.getExtraAtomInfo().
        @param ProbeRad : Probe radius
        @param atomMoverSets : Parameter that is modified in place to record all Movers that
          a particular atom interacts with.An entry is created whenever there is overlap
          with an atom in another Mover.
        @return True if a pair of atoms with one from each overlap, False if not.
    */
    bool PairsOverlap(
      boost::python::object mover1,
      scitbx::af::shared<iotbx::pdb::hierarchy::atom>  atoms1,
      scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > positions1,
      boost::python::object mover2,
      scitbx::af::shared<iotbx::pdb::hierarchy::atom>  atoms2,
      scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > positions2,
      molprobity::probe::ExtraAtomInfoMap extraAtomInfoMap,
      double probeRad,
      boost::python::dict atomMoverSets
    );

    typedef boost::adjacency_list <
      boost::setS, boost::listS, boost::undirectedS,
      boost::property<boost::vertex_name_t, boost::python::object>,
      boost::property<boost::edge_weight_t, double> > PythonObjectList;
    /** @brief For each pair of movers that are connected by an edge in the graph produced
        by the AABB algorithm to see if they actually overlap. If not, remove that edge.
        @todo Need to pass in the other objects needed to call PairsOverlap.
    */
    void RemoveFalseEdges(PythonObjectList& graph);


    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string InteractionGraph_test();
  }
}
