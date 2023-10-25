// Copyright(c) 2021, Richardson Lab at Duke
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

#include <vector>
#include <boost/array.hpp>
#include <set>
#include <algorithm>
#include "Common.h"
#include <iotbx/pdb/hierarchy.h>

namespace molprobity {
  namespace probe {

    /// @brief Spatial query acceleration object that tells which atoms are close to a Point.
    class SpatialQuery {
    public:
      /// @brief Construct with a spatial extent and bin size but to atoms.
      /// @param [in] lowerBounds Point indicating the lower bounds of each axis.
      /// @param [in] upperBounds Point indicating the upper bounds of each axis.
      ///         Each axis must be larger in value than the corresponding lowerBounds
      ///         axis or it will be set to the lowerBounds axis.
      /// @param [in] binSize The bin size in Angstroms for all three axes.  If this is
      ///         less than or equal to zero it will be set to 1.
      SpatialQuery(Point lowerBounds, Point upperBounds, Point binSize);

      /// @brief Construct with a hierarchy, filling in all atoms.
      ///
      /// This is a helper constructor that generates a filled-in structure
      /// whose extent matches that of the atom vector passed in, and that is
      /// pre-filled with all of its atoms.  The bin sizes are such that there
      /// at most 50 bins per axis.  This is equivalent to calling the bounds-based
      /// constructor and then calling add() on all of the atoms.
      /// @param [in] atoms Vector of atoms used to determine spatial extent.
      ///             the bin sizes are 3 Angstroms on each axis but with a maximum of
      ///             50 steps along each axis (125,000 total bins).  All atoms in the
      ///             model are added.
      ///             Note: For most cases, these atoms should all be from the same
      ///             conformation of a hierarchy, although it is of course possible to
      ///             use the parent() methods to chase up the hierarchy and verify that
      ///             the atom is in a specific group.
      SpatialQuery(scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &atoms);

      /// @brief Add an atom to the query object
      /// @param [in] a Atom to be added to the query object.
      /// @return True if the atom was added, false if it was not (because it was already there)
      bool add(iotbx::pdb::hierarchy::atom a);

      /// @brief Remove an atom from the query object
      ///
      /// This routine will normally be used when an atom has moved from one location to another;
      /// the old atom location will be removed and then the new atom will be added.
      /// @param [in] a Atom to be added to the query object.
      /// @return True if the atom was removed, false if it was not (because it was not there)
      bool remove(iotbx::pdb::hierarchy::atom a);

      /// @brief Locate atoms within a range of distances from a Point
      /// @param [in] p Point to measure distances from
      /// @param [in] min_distance Minimum distance from the Point to the atom.  When looking
      ///     for neighbors to an atom, setting this larger than 0 will cause the atom itself
      ///     not to be returned.  This can also be used to look for potentially-bonded atoms
      ///     by setting the minimum to the closest the two atoms should be, slightly less than
      ///     the sum of their radii.
      /// @param [in] max_distance Maximum distance from the point to the atom.  Specifies the
      ///     furthest an atom can be and be considered a neighbor.
      /// @return Vector of atoms that are within the specified distance from the Point. NOTE:
      ///     These are not guaranteed to be in any particular order from run to run, so if
      ///     you need them to be in the same order, you should sort them by i_seq.
      scitbx::af::shared<iotbx::pdb::hierarchy::atom> neighbors(
        Point const& p, double min_distance, double max_distance);

      //===========================================================================
      // Seldom-used methods below here.

      /// @brief Test method to verify that the class is behaving as intended.
      /// @return Empty string on success, string telling what went wrong on failure.
      static std::string test();

    protected:
      /// @brief the meat of the constructor, so it can be called by all constructors.
      void initialize(Point lowerBounds, Point upperBounds, Point binSize);

      /// @brief Default size of a grid bin in X, Y, and Z.
      ///
      /// The default size is set in the .cpp file to be of a size that would be expected
      /// to contain only a few atoms.  If the molecule is so large that this would produce
      /// too many bins along any axis, the size is made larger to avoid using too much
      /// memory.  This is the bin size chosen when it is not specified in the constructor.
      /// There is a constructor that lets the caller specify the bin size, in which case
      /// this is ignored.
      static const Coord DEFAULT_BIN_SIZE;

      Point   m_lowerBounds;            ///< Location of the lower corner of the grid in all dimensions
      boost::array<size_t, 3> m_gridSize; ///< Number of grid points in each axis
      Point   m_binSize;                ///< Width of a bin in each of the 3 directions

      /// We need the less-than operator to be defined on our atom type so that we can
      /// insert it into a set.
      struct atom_less : public std::binary_function<iotbx::pdb::hierarchy::atom, iotbx::pdb::hierarchy::atom, bool> {
        bool operator()(const iotbx::pdb::hierarchy::atom& lhs, const iotbx::pdb::hierarchy::atom& rhs) const
        {
          return lhs.data.get() < rhs.data.get();
        }
      };

      /// Grid that stores sets of atoms within each spatial location.
      /// X coordinate varies fastest, then Y, then Z.  Use the grid_point() method to
      /// get a reference to the vector where a specified Point is located.  Use the
      /// grid_index() method to get an index to that vector in the grid.
      typedef std::set<iotbx::pdb::hierarchy::atom, atom_less> GridPoint;
      std::vector<GridPoint> m_grid;

      /// @brief Return the index of the grid element that this point falls in.
      /// @param [in] p Point to find grid index for.  For points that lie outside
      ///         the grid, the edge element closest to the point is returned.
      /// @return Index of the grid element containing or closest to p.
      size_t  grid_index(Point const& p) const {
        boost::array<size_t, 3> xyz;
        for (size_t i = 0; i < 3; i++) {
          if (p[i] < m_lowerBounds[i]) { xyz[i] = 0; }
          else { xyz[i] = static_cast<size_t>(floor((p[i] - m_lowerBounds[i]) / m_binSize[i])); }
          if (xyz[i] >= m_gridSize[i]) { xyz[i] = m_gridSize[i] - 1; }
        }
        return xyz[0] + m_gridSize[0] * (xyz[1] + m_gridSize[1] * (xyz[2]));
      }

      /// @brief Return the grid element containing or closest to a Point.
      GridPoint& grid_point(Point const& p) { return m_grid[grid_index(p)]; }
    };

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string SpatialQuery_test();
  }
}
