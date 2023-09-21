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
#include "../probe/Common.h"
#include "../probe/Scoring.h"
#include "../probe/DotSpheres.h"
#include "../probe/SpatialQuery.h"
#include "PositionReturn.h"
#include "InteractionGraph.h"

namespace molprobity {
  namespace reduce {

    class OptimizerC {
    public:
      /** @brief Constructor.
          @param [in] verbosity: Controls how much information is added to the string.
          @param [in] preferenceMagnitude: Multiples the preference energies, so that we
                  can scale down their importance if we want.
          @param [in] minOccupancy: The minimum occupancy for an atom to be considered.
          @param [in] probeRadius: The radius of the probe sphere, in A.
          @param [in] probeDensity: The density of the probe sphere, in A^-3.
          @param [in] atoms: List of atoms.
          @param [in] exclude: Dictionary of atoms to exclude from collisions, looked up by i_seq.
          @param [in] dotScorer: Dot scorer to use.
          @param [in] dotSphereCache: Dot sphere cache to use to generate spheres for atoms.
          @param [in] atomMoverLists: Class with lists of movers that interact with an atom, looked up by i_seq.
          @param [inOut] spatialQuery: Spatial-query structure telling which atoms are where
          @param [inOut] extraAtomInfoMap: Map containing extra information about each atom.
          @param [inOut] deleteMes: Set of atoms to be deleted, passed as a Python object.
      */
      OptimizerC(
        int verbosity,
        double preferenceMagnitude,
        double maxVDWRadius,
        double minOccupancy,
        double probeRadius,
        double probeDensity,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& atoms,
        boost::python::dict& exclude,
        boost::python::object& dotScorer,
        boost::python::object& dotSphereCache,
        AtomMoverLists& atomMoverLists,       //< Defined in InteractionGraph.h
        molprobity::probe::SpatialQuery& spatialQuery,
        molprobity::probe::ExtraAtomInfoMap& extraAtomInfoMap,
        boost::python::object& deleteMes);

      /** @brief Initialize the Movers to their 0 coarse states and record initial scores.
          @param [in] movers: A list of Movers that will be processed.
          @return: A string describing what was done, which may be empty if verbosity is too small.
      */
      std::string Initialize(scitbx::af::shared<boost::python::object> const &movers);

      /** @brief Function to perform coarse optimization on a singleton Mover.
          @param [in] mover: The Mover to optimize.
          @return A tuple, where the first is the score at the best position for all movers
                  and the second is a string describing what was done, which may be empty
                  if verbosity is too small.
          side effect: Changes the value of GetHighScore(mover) to the score at the coarse position
          selected if one is selected.
      */
      boost::python::tuple OptimizeSingleMoverCoarse(boost::python::object const &mover);

      /** @brief Function to perform fine optimization on a singleton Mover.

          Find the score for the Mover in all fine orientations by moving each atom into the
          specified position and summing the scores over all of them.  Determine the best
          orientation by selecting the highest scorer.
          Add the preference energy to the sum for each orientation scaled by our preference
          @param [in] mover: The Mover to optimize.
          @return A tuple, where the first is the score at the best position for all movers
                  and the second is a string describing what was done, which may be empty
                  if verbosity is too small.
          side effect: Changes the value of GetHighScore(mover) to the score at the fine position
          selected if one is selected.
      */
      boost::python::tuple OptimizeSingleMoverFine(boost::python::object const& mover);

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

      /// @brief Returns the coarse location of the specified mover.
      unsigned GetCoarseLocation(boost::python::object const& mover) {
        return m_coarseLocations[mover.ptr()];
      }

      /// @brief Returns the fine location of the specified mover, -1 if none.
      int GetFineLocation(boost::python::object const& mover) {
        return m_fineLocations[mover.ptr()];
      }

      /// @brief Returns the high score of the specified mover, 0 if none.
      double GetHighScore(boost::python::object const& mover) {
        return m_highScores[mover.ptr()];
      }

      /// @brief Returns the dots for a specified atom i_seq.
      scitbx::af::shared<molprobity::probe::Point> GetDots(unsigned atom_i_seq) {
        return m_dotSpheres[atom_i_seq];
      }

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
      int m_verbosity;
      double m_maxVDWRadius;
      double m_preferenceMagnitude;
      double m_minOccupancy;
      double m_probeRadius;
      double m_probeDensity;
      scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& m_atoms;
      boost::python::dict m_exclude;          //< Make a copy so it will persist
      molprobity::probe::DotScorer& m_dotScorer;
      molprobity::probe::DotSphereCache& m_dotSphereCache;
      AtomMoverLists &m_atomMoverLists;
      molprobity::probe::SpatialQuery& m_spatialQuery;
      molprobity::probe::ExtraAtomInfoMap& m_extraAtomInfoMap;
      boost::python::object m_deleteMes;      //< Make a copy so it will persist
      std::map<PyObject*, unsigned> m_coarseLocations;
      std::map<PyObject*, int> m_fineLocations;
      std::map<PyObject*, double> m_highScores;

      //=========================================================================
      // Section dealing with reducing the number of dots that must be checked
      // for each atom.

      /// Stores the full list of dots for each atom, indexed by i_seq. This may also
      /// include dots that are inside excluded atoms.
      std::map<unsigned, scitbx::af::shared<molprobity::probe::Point> > m_dotSpheres;

      /// Stores the list of dots for each atom that are not inside excluded atoms,
      /// indexed by both i_seq and the coarse position of the Mover that the atom
      /// is part of. This is cleared during Initialize() and filled in lazily as they
      /// are encountered within scoreAtom().
      /// The set of dots that are excluded varies as this atom's Mover changes coarse
      /// location but it does not depend on other Movers.
      std::map< std::pair<unsigned, unsigned>, scitbx::af::shared<molprobity::probe::Point> > m_excludedDots;

      //=========================================================================
      // Section dealing with caching of computed atom scores

      /// Cached scores for atoms that have already been calculated based on the
      /// coarse positions of the Movers that they depend on.
      /// This is the map stored for a particular atom, based on its relevant Mover positions.
      typedef std::map< std::vector<unsigned>, double > ScoreCache;
      /// This is the type of a map that stores cached scores for all atoms, looked up by i_seq.
      typedef std::map< unsigned, ScoreCache > ScoreCacheMap;

      /// This is a pointer to the ScoreCacheMap for the current clique, if it exists.
      /// When it is a nullptr, we do not use caching. This is set to non-nullptr
      /// when we are optimizing a clique, and set back to nullptr when we are done.
      ScoreCacheMap *m_scoreCacheMap;

      /// Tracks how many cached vs. calculated scores we have over the course of our calculations.
      /// These are used for reporting purposes and are accumulated over the entire life of
      /// the OptimizeC object.
      size_t m_cachedScores;
      size_t m_calculatedScores;

      /// @brief Score an atom, not using the atom-score cache.
      /// @param [in] a: The atom to score.
      /// @param [in] locationIndex: The index for the location of the Mover that the
      ///         atom is part of. This must be different for every location, with a
      ///         separate range for coarse and fine locations. It is used to cache
      ///         dots.
      /// @return The score for the atom.
      double scoreAtom(iotbx::pdb::hierarchy::atom const& a, unsigned locationIndex);

      /// @brief Score an atom, using the cache.
      ///
      /// The m_scoreCacheMap must be non-nullptr when this is called.
      /// @param [in] a: The atom to score.
      /// @param [in] locationIndex: The index for the location of the Mover that the
      ///         atom is part of. This must be different for every location, with a
      ///         separate range for coarse and fine locations. It is used to cache
      ///         dots.
      /// @return The score for the atom.
      double scoreAtomCached(iotbx::pdb::hierarchy::atom const& a, unsigned locationIndex);

      /// @brief Score all atoms that have not been marked for deletion, calling appropriate
      /// scoring function (use caching to do the scoring if we have an atom-score cache).
      /// @param [in] states: The set of positions that we are scoring (may be coarse, may be
      ///         fine).
      /// @param [in] index: The index of the state that we are scoring.
      /// @param [in] dotCacheOffset: This is added to the index of the state to get
      ///         the parameter to pass to the atom-scoring code. This lets us make
      ///         separate ranges for the coarse and fine scoring.
      double scorePosition(molprobity::reduce::PositionReturn& states, unsigned index,
        unsigned dotCacheOffset);

      //=========================================================================
      // Section dealing with optimization of a clique

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
      ///
      /// This finds a set of Movers that, if removed, would disconnect the graph.
      /// It starts looking for a cut of size 1, then 2, etc. until it finds one
      /// or finds that the graph cannot be cut.
      /// @param [in] graph: The clique graph to find the vertex cut in.
      /// @param [out] cutMovers: The Movers that correspond can be removed to cause a cut.
      /// @param [out] cutGraph: The graph that stores the subset of clique
      /// whose vertices have not been removed. If there is more than one Mover
      /// in cutMovers, then this graph will have two or more connected components.
      /// If there are no Movers, then this will be the original graph.
      static void findVertexCut(CliqueGraph const& graph,
        std::vector<boost::python::object*>& cutMovers, CliqueGraph& cutGraph);

      /// @brief Function to perform brute-force optimization on a clique of Movers.
      ///
      /// This function will try all combinations of states for the Movers in the clique
      /// and return the best one. It will be called when a subgraph cannot be cut into
      /// more than one connected component.
      /// @param [in] states: A map from Movers to the potential states for each. This
      ///         holds the result of calling CoarsePositions() on each Mover.
      /// @param [in] clique: The clique graph to optimize.
      /// @return A tuple, where the first is the score at the best position for all movers
      ///         and the second is a string describing what was done, which may be empty
      ///         if verbosity is too small.
      std::pair<double, std::string> OptimizeCliqueCoarseBruteForce(
        std::map<boost::python::object*, molprobity::reduce::PositionReturn>& states,
        CliqueGraph& clique);

      /// @brief Function to perform vertex-cut optimization on a clique of Movers.
      ///
      /// This attempts to find a vertex cut in the clique graph, and if it does, it
      /// optimizes each connected component separately. If it cannot find a vertex
      /// cut, it calls OptimizeCliqueCoarseBruteForce.
      ///
      /// When it finds a vertex cut, it will recursively call itself on each connected
      /// component and jointly optimize the cut vertices, then combine the results.
      /// @param [in] states: A map from Movers to the potential states for each. This
      ///         holds the result of calling CoarsePositions() on each Mover.
      /// @param [in] clique: The clique graph to optimize.
      /// @return A tuple, where the first is the score at the best position for all movers
      ///         and the second is a string describing what was done, which may be empty
      ///         if verbosity is too small.
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
