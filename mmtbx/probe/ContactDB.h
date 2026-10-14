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
      enum AtomFlag { IGNORE = 1, MAINCHAIN = 2, HET = 4, WATER = 8, TARGET = 16 };

      struct PairRecord {
        unsigned cause;
        unsigned counts[N_TYPES];
        double minGap;
        double typeMinGap[N_TYPES];   ///< Smallest gap per InteractionType (1e9 if none)
        bool hasBump;
        double bumpGap;     ///< Smallest gap among this pair's Bump dots (first one on ties)
        Point bumpLoc;      ///< Location of that dot
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

      /// Score every source atom, in the given order.
      void build(scitbx::af::shared<iotbx::pdb::hierarchy::atom> sources);

      /// Flat arrays, one entry per pair in source order: (source index into the build()
      /// list, cause i_seq, min gap, N_TYPES counts per pair, has-Bump flag, Bump gap,
      /// Bump location, N_TYPES minimum gaps per pair).
      boost::python::tuple pairs() const;

    private:
      unsigned idx(iotbx::pdb::hierarchy::atom const& a) const;
      bool compatible(iotbx::pdb::hierarchy::atom const& a, iotbx::pdb::hierarchy::atom const& b) const;
      std::vector<iotbx::pdb::hierarchy::atom> excluded_for(iotbx::pdb::hierarchy::atom const& a,
                                                            unsigned N, unsigned nonHydrogenN) const;
      void score_atom(iotbx::pdb::hierarchy::atom const& src, std::vector<PairRecord>& out);

      ExtraAtomInfoMap m_info;
      DotScorer m_scorer;
      SpatialQuery* m_grid;
      SpatialQuery* m_phantoms;
      DotSphereCache m_dotCache;
      std::vector<int> m_flags, m_model, m_chain, m_group, m_altloc, m_bondOffsets;
      std::vector<iotbx::pdb::hierarchy::atom> m_bondTargets, m_sources;
      std::vector<std::vector<PairRecord> > m_records;
      double m_probeRadius, m_maxRadius, m_minOccupancy, m_overlapScale;
      int m_excludedBondChainLength;
      bool m_includeMainchainMainchain, m_includeWaterWater, m_separateWorseClashes,
           m_reportHBonds, m_reportClashes, m_reportVdws, m_onlyBadClashes;
    };

  }
}
