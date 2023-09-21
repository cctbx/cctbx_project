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

// Contains functions and classes to support the InteractionGraph.py classes.
// They are linked into the mmtbx_reduce_ext library by boost_python/reduce_bpl.cpp.

#pragma once

#include <string>
#include <vector>
#include <boost/python.hpp>
#include "../../boost_adaptbx/graph/graph_type.hpp"
#include "../probe/Scoring.h"

namespace molprobity {
  namespace reduce {

    /** @brief Class to keep track of the Movers that each atom interacts with. */
    class AtomMoverLists {
    public:

      /// @brief Default constructor.
      AtomMoverLists();

      /// @brief Add a Mover to the list for an atom.
      /// @param i_seq : Index of the atom
      /// @param mover : Mover to add
      void AddAtomMoverEntry(unsigned i_seq, boost::python::object mover);

      /// @brief Clear the list of Movers for all atoms.
      void Clear();

      /// @brief Get the list of Movers for an atom.
      /// @param i_seq : Index of the atom
      /// @return List of Movers for the atom. Throws a Python exception if it beyond those added.
      std::vector< boost::python::object > const& GetAtomMoverList(unsigned i_seq) const;

    protected:
      /// Vector indexed by i_seq of the atom. Each entry is a list of Movers that the atom interacts with.
      std::vector< std::vector< boost::python::object > > m_atomMoverLists;
    };

    /** @brief Function to report whether any atoms overlap between two Movers.

        Helper function that tells whether any pair of atoms from two Movers overlap.
        @param mover1 : The first Mover
        @param atoms1 : Atom list for the first Mover
        @param positions1 : probe.PositionReturn.positions holding possible positions for each.
          This is indexed first by position, then by atom in that position.
        @param mover2 : The second Mover
        @param atoms2 : Atom list for the second Mover
        @param positions2 : probe.PositionReturn.positions holding possible positions for each.
          This is indexed first by position, then by atom in that position.
        @param extraAtomInfoMap : probe.ExtraAtomInfoMap that can be used to look
          up the information for atoms whose values need to be changed. Can be
          obtained by calling mmtbx.probe.Helpers.getExtraAtomInfo().
        @param ProbeRad : Probe radius
        @param atomMoverSets : Parameter that is modified in place to record all Movers that
          a particular atom interacts with. An entry is created whenever there is overlap
          with an atom in another Mover. Indexed by i_seq of the atom. Each entry includes
          a list of Movers that the atom interacts with.
        @return True if a pair of atoms with one from each overlap, False if not.
    */
    bool PairsOverlap(
      boost::python::object const &mover1,
      boost::python::object const &mover2,
      molprobity::probe::ExtraAtomInfoMap const &extraAtomInfoMap,
      double probeRad,
      AtomMoverLists &atomMoverLists
    );

    /** @brief Function to determine which Movers have overlapping bounding boxes.
    * @param movers : List of Movers to check
    * @param extraAtomInfoMap : probe.ExtraAtomInfoMap that can be used to look
    * up the information for atoms whose values need to be changed. Can be
    * obtained by calling mmtbx.probe.Helpers.getExtraAtomInfo().
    * @param probeRad : Probe radius
    * @return Vector of vectors of indices of movers that overlap. Each inner vector has
    * two values.
    */
    scitbx::af::shared< scitbx::af::shared<int> > FindOverlappingMoversAABB(
      scitbx::af::shared<boost::python::object> const &movers,
      molprobity::probe::ExtraAtomInfoMap const &extraAtomInfoMap,
      double probeRad);

    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string InteractionGraph_test();
  }
}
