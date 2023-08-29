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

#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <boost/python.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "../probe/Scoring.h"
#include "../probe/SpatialQuery.h"
#include "PositionReturn.h"

namespace molprobity {
  namespace reduce {

    class OptimizerC {
    public:
      /** @brief Constructor.
          @param [in] self: The Python object that constructed us. /// @todo Pass all values and remove this.
          @param [in] verbosity: Controls how much information is added to the string.
          @param [in] preferenceMagnitude: Multiples the preference energies, so that we
                  can scale down their importance if we want.
          @param [in] minOccupancy: The minimum occupancy for an atom to be considered.
          @param [in] probeRadius: The radius of the probe sphere, in A.
          @param [in] probeDensity: The density of the probe sphere, in A^-3.
          @param [in] exclude: Dictionary of atoms to exclude from collisions, looked up by i_seq.
          @param [in] dotSpheres: Dictionary of dot spheres, looked up by i_seq.
          @param [in] atomMoverLists: Dictionary of list of movers, looked up by i_seq.
          @param [inOut] spatialQuery: Spatial-query structure telling which atoms are where
          @param [inOut] extraAtomInfoMap: Map containing extra information about each atom.
          @param [inOut] deleteMes: Set of atoms to be deleted, passed as a Python object.
          @param [out] coarseLocations: Dictionary looked up by mover that records the
                  coarse locations of the movers once they are optimized. These are modified
                  to point to the result.
          @param [out] highScores: Dictionary looked up by mover with the scores at best locations.
      */
      OptimizerC(
        boost::python::object& self,
        int verbosity,
        double preferenceMagnitude,
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
        boost::python::dict& highScores);

      /** @brief Function to perform fast optimization on a singleton Mover.
          @param [in] mover: The Mover to optimize.
          @return A tuple, where the first is the score at the best position for all movers
                  and the second is a string describing what was done, which may be empty
                  if verbosity is too small.
      */
      boost::python::tuple OptimizeSingleMoverCoarse(boost::python::object const &mover);

      /** @brief Function to perform fast optimization on a clique of Movers.
      * 
      *   This is the main function that is called from Python. It unwraps the
      *   Python objects to C++ and rewraps the return to Python, enabling us to
      *   run faster in the recursive calls.
          @param [in] movers: A list of Movers to jointly optimize.
          @param [in] interactions: A list of edges between movers, as a 2D array where the first
                  index is the number of edges and the second is 2. It stores the index
                  into the movers array of the mover at each end of the edge. This lists
                  the pairs of movers that interact with each other.
          @return A tuple, where the first is the score at the best position for all movers
                  and the second is a string describing what was done, which may be empty
                  if verbosity is too small.
      */
      boost::python::tuple OptimizeCliqueCoarse(
        scitbx::af::shared<boost::python::object> movers,
        scitbx::af::versa<int, scitbx::af::flex_grid<> >& interactions);

      /// @brief Returns the number of calculated atom scores within cliques
      size_t GetNumCalculatedAtoms() const { return m_calculatedScores; }

      /// @brief Returns the number of cached atom scores within cliques
      size_t GetNumCachedAtoms() const { return m_cachedScores; }

      /// @brief Test function that checks our privcate static methods.
      ///
      /// Does not do a complete behavior test of the class -- that must be done
      /// from Python code that can generate models and run the methods.
      /// @return Empty string if all tests pass, otherwise a string describing the failure.
      static std::string Test();

    protected:
      boost::python::object m_self;           //< Make a copy so it will persist
      int m_verbosity;
      double m_maxVDWRadius;
      double m_preferenceMagnitude;
      double m_minOccupancy;
      double m_probeRadius;
      double m_probeDensity;
      boost::python::dict m_exclude;          //< Make a copy so it will persist
      boost::python::dict m_dotSpheres;       //< Make a copy so it will persist
      boost::python::dict m_atomMoverLists;    //< Make a copy so it will persist
      molprobity::probe::SpatialQuery& m_spatialQuery;
      molprobity::probe::ExtraAtomInfoMap& m_extraAtomInfoMap;
      boost::python::object m_deleteMes;      //< Make a copy so it will persist
      boost::python::dict m_coarseLocations;  //< Make a copy so it will persist
      boost::python::dict m_highScores;       //< Make a copy so it will persist
      molprobity::probe::DotScorer *m_dotScorer = nullptr;  //< Pointer to the DotScorer object

      /// Caches scores for atoms that have already been calculated based on the
      /// values of the Movers that they depend on.
      /// This is the map stored for a particular atom, based on its relevant Mover positions.
      typedef std::map< std::vector<unsigned>, double > ScoreCache;
      /// This is the map stored for all atoms, looked up by i_seq.
      typedef std::map< unsigned, ScoreCache > ScoreCacheMap;

      /// tracks how many cached vs. calculated scores we have.
      size_t m_cachedScores;
      size_t m_calculatedScores;

      /// This is a pointer to the ScoreCacheMap for the current clique, if it exists.
      /// When it is a nullptr, we do not use caching. This is set to non-nullptr
      /// when we are optimizing a clique, and set back to nullptr when we are done.
      ScoreCacheMap *m_scoreCacheMap = nullptr;

      /// @brief Score an atom, not using the cache.
      double scoreAtom(iotbx::pdb::hierarchy::atom const& a);

      /// @brief Score an atom, using the cache.
      ///
      /// The m_scoreCacheMap must be non-nullptr when this is called.
      double scoreAtomCached(iotbx::pdb::hierarchy::atom const& a);

      /// @brief Score all atoms that have not been marked for deletion, calling the Python object's
      /// scoring function (which may or may not use caching to do the scoring).
      double scorePosition(molprobity::reduce::PositionReturn& states, size_t index);

      std::string setMoverState(molprobity::reduce::PositionReturn& positionReturn, unsigned index);

      // Must use vector style for second (vertex) entry for connected_components to work
      typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::python::object*>
        CliqueGraph;

      static std::vector< std::vector<unsigned> > generateAllStates(std::vector<unsigned> const& numStates);

      /// @brief Return a vector of vectors of integers representing all combinations of m integers
      ///        from 0 to n-1, inclusive.
      static std::vector< std::vector<int> > nChooseM(int n, int m);

      /// @brief Find the subset of a clique graph that contains a given set of movers and edges between them.
      static CliqueGraph subsetGraph(CliqueGraph const& graph, std::vector<boost::python::object*>& keepMovers);

      /// @brief Find one of the smallest vertex cuts in a clique graph.
      /// @param [in] graph: The clique graph to find the vertex cut in.
      /// @param [out] cutMovers: The movers that correspond to the vertices
      /// that are removed.
      /// @param [out] cutGraph: The graph that stores the subset of clique
      /// whose vertices have been removed.
      static void findVertexCut(CliqueGraph const& graph,
        std::vector<boost::python::object*>& cutMovers, CliqueGraph& cutGraph);

      /** @brief Function to perform brute-force optimization on a clique of Movers.
          @param [in] movers: A list of Movers to jointly optimize.
          @return A tuple, where the first is the score at the best position for all movers
                  and the second is a string describing what was done, which may be empty
                  if verbosity is too small.
      */
      std::pair<double, std::string> OptimizeCliqueCoarseBruteForce(
        std::map<boost::python::object*, molprobity::reduce::PositionReturn>& states,
        CliqueGraph& clique);

      /** @brief Function to perform vertex-cut optimization on a clique of Movers.
          @param [in] movers: A list of Movers to jointly optimize.
          @param [in] interactions: A list of edges between movers, as a 2D array where the first
                  index is the number of edges and the second is 2. It stores the index
                  into the movers array of the mover at each end of the edge. This lists
                  the pairs of movers that interact with each other.
          @return A tuple, where the first is the score at the best position for all movers
                  and the second is a string describing what was done, which may be empty
                  if verbosity is too small.
      */
      std::pair<double, std::string> OptimizeCliqueCoarseVertexCut(
        std::map<boost::python::object*, molprobity::reduce::PositionReturn> &states,
        CliqueGraph &clique);
    };

    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    ///
    /// Does not do a complete behavior test of the class -- that must be done
    /// from Python code that can generate models and run the methods.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string Optimizers_test();
  }
}
