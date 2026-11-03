// Licensed under the Apache 2 license
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once
#include "Scoring.h"
#include "SpatialQuery.h"
#include "DotSpheres.h"
#include <iotbx/pdb/hierarchy.h>
#include <scitbx/array_family/shared.h>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <cmath>
#include <vector>

namespace molprobity {
  namespace probe {

    /// The atoms bonded to a within N hops that are in a compatible conformation and close
    /// enough to touch, not including a (see Helpers.getAtomsWithinNBonds for the arguments).
    /// neighbors(x, out) appends the atoms bonded to x; compatible(a, n) compares conformations.
    template <class Neighbors, class Compatible>
    std::vector<iotbx::pdb::hierarchy::atom> atomsWithinNBonds(
      iotbx::pdb::hierarchy::atom const& a, Neighbors const& neighbors, Compatible const& compatible,
      ExtraAtomInfoMap const& info, double probeRad, int N, double nonHydrogenN)
    {
      typedef iotbx::pdb::hierarchy::atom atom;
      Point const& aLoc = a.data->xyz;
      double aRad = info.getMappingFor(a).getVdwRadius();
      bool atomIsHydrogen = a.element_is_hydrogen();
      // Only ever allow a Phantom Hydrogen to be bonded to its parent Oxygen.
      if (info.getMappingFor(a).getIsDummyHydrogen()) { N = 1; }
      std::vector<atom> found(1, a);   // a handful of atoms, so linear membership tests
      std::vector<atom> bonded;
      for (int i = 0; i < N; i++) {
        size_t levelSize = found.size();
        for (size_t k = 0; k < levelSize; k++) {
          bonded.clear();
          neighbors(found[k], bonded);
          for (size_t b = 0; b < bonded.size(); b++) {
            atom const& n = bonded[b];
            bool seen = false;
            for (size_t f = 0; f < found.size() && !seen; f++) { seen = (found[f].data == n.data); }
            if (seen) { continue; }
            if ((i < nonHydrogenN || atomIsHydrogen || n.element_is_hydrogen()) && compatible(a, n)) {
              Point const& nLoc = n.data->xyz;
              double nRad = info.getMappingFor(n).getVdwRadius();
              double dx = nLoc[0] - aLoc[0], dy = nLoc[1] - aLoc[1], dz = nLoc[2] - aLoc[2];
              if (std::sqrt(dx*dx + dy*dy + dz*dz) <= aRad + nRad + 2 * probeRad) {
                found.push_back(n);
              }
            }
          }
        }
      }
      found.erase(found.begin());
      return found;
    }

    /// atomsWithinNBonds() reading bondedNeighborLists, a Python mapping from each atom to the
    /// list of atoms bonded to it; same arguments as Helpers.getAtomsWithinNBonds.
    boost::python::list getAtomsWithinNBonds(iotbx::pdb::hierarchy::atom const& a,
      boost::python::object const& bondedNeighborLists, ExtraAtomInfoMap const& info,
      double probeRad, int N, double nonHydrogenN);

    /// @brief probe2's per-atom dot loop, with contacts summarized per (source, cause) atom pair.
    ///
    /// Filters the atoms near each source atom and scores its dots with
    /// DotScorer::check_dots(). For each pair it keeps the dot count per InteractionType,
    /// the minimum gap, and the minimum-gap Bump dot (the line clashscore2 builds a clash
    /// from); with set_keep_dots() it also keeps every reported dot, for probe2's output.
    /// Candidate atoms are passed to check_dots in i_seq order, so ties between atoms at
    /// equal gap resolve the same way every run.
    class ContactDB {
    public:
      enum { N_TYPES = 7 };
      enum AtomFlag { IGNORE = 1, MAINCHAIN = 2, HET = 4, WATER = 8, TARGET = 16, SIDECHAIN = 32 };

      struct PairRecord {
        unsigned cause;
        unsigned counts[N_TYPES];
        double minGap;
        double typeMinGap[N_TYPES];   ///< Smallest gap per InteractionType (1e9 if none)
        bool hasBump;
        double bumpGap;     ///< Smallest gap among this pair's Bump dots (first one on ties)
        Point bumpLoc;      ///< Location of that dot
      };

      /// One reported dot, with the fields of probe2's DotInfo.
      struct DotRecord {
        unsigned src;       ///< i_seq
        unsigned cause;     ///< i_seq
        int itype;          ///< DotScorer::InteractionType
        int overlapType;    ///< DotScorer::OverlapType
        char ptmaster;      ///< ' ', 'M', 'S', 'P' or 'O'
        double gap;
        Point loc;
        Point spike;
      };

      /// The grids are referenced, not copied; the Python binding keeps them alive.
      ContactDB(ExtraAtomInfoMap const& info, DotScorer const& scorer, SpatialQuery& grid,
                SpatialQuery& phantoms, DotSphereCache const& dotCache);

      /// Per-atom properties indexed by i_seq. altloc: 0 means "" or " " (compatible with all).
      void set_atom_data(scitbx::af::shared<int> flags, scitbx::af::shared<int> model,
                         scitbx::af::shared<int> chain, scitbx::af::shared<int> atomGroup,
                         scitbx::af::shared<int> altloc);

      /// Bonded neighbors in CSR form: neighbors of i_seq i are targets[offsets[i]..offsets[i+1]).
      void set_bonds(scitbx::af::shared<int> offsets,
                     scitbx::af::shared<iotbx::pdb::hierarchy::atom> targets);

      void set_params(double probeRadius, double maxRadius, double minOccupancy,
                      int excludedBondChainLength, bool includeMainchainMainchain,
                      bool includeWaterWater, double overlapScale, bool separateWorseClashes,
                      bool reportHBonds, bool reportClashes, bool reportVdws, bool onlyBadClashes);

      /// Also keep every reported dot (off by default; ~70 bytes per dot).
      void set_keep_dots(bool keep) { m_keepDots = keep; }

      /// Score every source atom, in the given order.
      void build(scitbx::af::shared<iotbx::pdb::hierarchy::atom> sources);

      /// Flat arrays, one entry per pair in source order: (source i_seq, cause i_seq, min gap,
      /// N_TYPES counts per pair, has-Bump flag, Bump gap, Bump location, N_TYPES minimum
      /// gaps per pair).
      boost::python::tuple pairs() const;

      /// Flat arrays, one entry per kept dot in the order they were scored: (source i_seq,
      /// cause i_seq, InteractionType, OverlapType, gap, location, spike end, ptmaster
      /// characters as one string).
      boost::python::tuple dots() const;

    private:
      unsigned idx(iotbx::pdb::hierarchy::atom const& a) const;
      bool compatible(iotbx::pdb::hierarchy::atom const& a, iotbx::pdb::hierarchy::atom const& b) const;
      std::vector<iotbx::pdb::hierarchy::atom> excluded_for(iotbx::pdb::hierarchy::atom const& a,
                                                            unsigned N, unsigned nonHydrogenN) const;
      void score_atom(iotbx::pdb::hierarchy::atom const& src, std::vector<PairRecord>& out);
      void keep_dot(unsigned si, int sf, unsigned ci, Point const& sXYZ, double srcVdw,
                    int itype, DotScorer::CheckDotResult const& res);

      ExtraAtomInfoMap m_info;
      DotScorer m_scorer;
      SpatialQuery* m_grid;
      SpatialQuery* m_phantoms;
      DotSphereCache m_dotCache;
      std::vector<int> m_flags, m_model, m_chain, m_group, m_altloc, m_bondOffsets;
      std::vector<iotbx::pdb::hierarchy::atom> m_bondTargets, m_sources;
      std::vector<std::vector<PairRecord> > m_records;
      std::vector<DotRecord> m_dots;
      double m_probeRadius, m_maxRadius, m_minOccupancy, m_overlapScale;
      int m_excludedBondChainLength;
      bool m_includeMainchainMainchain, m_includeWaterWater, m_separateWorseClashes,
           m_reportHBonds, m_reportClashes, m_reportVdws, m_onlyBadClashes, m_keepDots;
    };

  }
}
