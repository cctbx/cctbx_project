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
#include <vector>
#include "../probe/Scoring.h"
#include <iotbx/pdb/hierarchy.h>
#include <scitbx/boost_python/container_conversions.h>

namespace molprobity {
  namespace reduce {

    /// @brief Structure to hold atom-behavior information returned from Mover methods.
    class PositionReturn {
    public:

      /// @brief Default constructor
      PositionReturn() {};

      /// @brief Constructor with all arguments
      ///
      /// The arguments must all be const references so that the automagic wrapping
      /// works correctly with them.
      /// @param p_atoms A list of all of the atoms to be adjusted.
      /// @param p_positions The positions element has the new location of each atom in each set
      ///     of positions.
      /// @param p_extraInfos The extraInfos element has the new ExtraAtomInfo for each atom in
      ///     each set of positions that tells new information for a subset of the atoms.
      ///     If there are no changes, this list is empty. It can be smaller than the
      ///     p_atoms list, only effecting the first n atoms.
      /// @param p_deleteMes The deleteMes element has a list of booleans for each atom in each
      ///     set of positions that tells whether the atom should be deleted in this set.
      ///     If there are no deletions, this list is empty. It can be smaller than the
      ///     p_atoms list, only effecting the first n atoms.
      /// @param p_preferenceEnergies The preferenceEnergies element has an energy for
      ///     each set of positions.
      /// @return A PositionReturn object with the given values.
      PositionReturn(
        scitbx::af::shared<iotbx::pdb::hierarchy::atom>  const& p_atoms
        , scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> >  const& p_positions
        , scitbx::af::shared< scitbx::af::shared<molprobity::probe::ExtraAtomInfo> > const& p_extraInfos
        , scitbx::af::shared< scitbx::af::shared<bool> > const& p_deleteMes
        , scitbx::af::shared<double>  const& p_preferenceEnergies
      ) : atoms(p_atoms)
        , positions(p_positions)
        , extraInfos(p_extraInfos)
        , deleteMes(p_deleteMes)
        /// @todo preferenceEnergies(p_preferenceEnergies)
      {
        /// @todo Remove once we are using af::shared rather than std::vector.
        preferenceEnergies.resize(p_preferenceEnergies.size());
        for (size_t i = 0; i < p_preferenceEnergies.size(); ++i) {
          preferenceEnergies[i] = p_preferenceEnergies[i];
        }
      };

      /// A list of all of the atoms to be adjusted.
      scitbx::af::shared<iotbx::pdb::hierarchy::atom>  atoms;
      /// The positions element has the new location of each atom in each set of positions.
      scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> >  positions;
      /// The extraInfos element has the new ExtraAtomInfo for each atom in each set of
      /// positions. This array may be shorter in length than the number of atoms
      /// (and may be empty) because some Movers do not need to change the information
      /// for any or all atoms.The index in this array will match the index in the atoms
      /// array so the earliest atoms will be changed if a subset is present.
      scitbx::af::shared< scitbx::af::shared<molprobity::probe::ExtraAtomInfo> > extraInfos;
      /// The deleteMes element tells whether each atom in each set of positions
      /// should be deleted. This means that it should be ignored in all calculations
      /// and also should be deleted from the model if this configuration is chosen.
      /// This array may be shorter in length than the number of atoms (and may be empty)
      /// because some Movers do not need to change the information for any or all
      /// atoms.The index in this array will match the index in the atoms array so the
      /// earliest atoms will be deleted if a subset is present.
      scitbx::af::shared< scitbx::af::shared<bool> >  deleteMes;
      /// The preferenceEnergies entry holds an additional bias term that should be
      /// added to the Probe score for each set of positions before comparing them
      /// with each other.
      /// @todo We use a std::vector here to avoid double-wrapping the af::shared:double.
      std::vector<double>  preferenceEnergies;
    };

    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string PositionReturn_test();
  }
}
