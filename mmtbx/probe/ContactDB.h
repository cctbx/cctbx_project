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
#include <boost/python/tuple.hpp>
#include <vector>

namespace molprobity {
  namespace probe {

    /// @brief Contacts summarized per (source, cause) atom pair instead of per dot.
    ///
    /// Implements the per-atom loop of probe2's _generate_interaction_dots() with the
    /// same filters and the same DotScorer::check_dots(), but keeps, for each pair, the
    /// dot count per InteractionType, the minimum gap, and the minimum-gap Bump dot (the
    /// line clashscore2 builds a clash from). Candidate atoms are passed to check_dots in
    /// i_seq order, so ties between atoms at equal gap resolve the same way every run.
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

      /// One reported dot, as probe2.py's _save_dot() records it.
      struct DotRecord {
        unsigned src;       ///< Index into the build() list
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

      /// Flat arrays, one entry per pair in source order: (source index into the build()
      /// list, cause i_seq, min gap, N_TYPES counts per pair, has-Bump flag, Bump gap,
      /// Bump location, N_TYPES minimum gaps per pair).
      boost::python::tuple pairs() const;

      /// Flat arrays, one entry per kept dot in the order probe2.py saves them: (source
      /// index into the build() list, cause i_seq, InteractionType, OverlapType, gap,
      /// location, spike end, ptmaster characters as one string).
      boost::python::tuple dots() const;

    private:
      unsigned idx(iotbx::pdb::hierarchy::atom const& a) const;
      bool compatible(iotbx::pdb::hierarchy::atom const& a, iotbx::pdb::hierarchy::atom const& b) const;
      std::vector<iotbx::pdb::hierarchy::atom> excluded_for(iotbx::pdb::hierarchy::atom const& a,
                                                            unsigned N, unsigned nonHydrogenN) const;
      void score_atom(iotbx::pdb::hierarchy::atom const& src, unsigned srcIndex,
                      std::vector<PairRecord>& out);
      void keep_dot(unsigned srcIndex, int sf, unsigned ci, Point const& sXYZ, double srcVdw,
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
