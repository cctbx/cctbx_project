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

#include "ContactDB.h"
#include <boost/python.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// The distance tests below must round exactly as the probe2.py code they replace.
#ifdef __clang__
#pragma clang fp contract(off)
#endif

namespace molprobity {
  namespace probe {

    using iotbx::pdb::hierarchy::atom;
    namespace af = scitbx::af;

    namespace {
      bool by_iseq(atom const& a, atom const& b) { return a.data->i_seq < b.data->i_seq; }

      bool contains(std::vector<atom> const& v, unsigned i_seq) {
        for (size_t k = 0; k < v.size(); k++) { if (v[k].data->i_seq == i_seq) { return true; } }
        return false;
      }
    }

    ContactDB::ContactDB(ExtraAtomInfoMap const& info, DotScorer const& scorer, SpatialQuery& grid,
                         SpatialQuery& phantoms, DotSphereCache const& dotCache)
      : m_info(info), m_scorer(scorer), m_grid(&grid), m_phantoms(&phantoms), m_dotCache(dotCache)
      , m_probeRadius(0.25), m_maxRadius(0), m_minOccupancy(0), m_overlapScale(0.5)
      , m_excludedBondChainLength(4), m_includeMainchainMainchain(true), m_includeWaterWater(false)
      , m_separateWorseClashes(false), m_reportHBonds(true), m_reportClashes(true)
      , m_reportVdws(true), m_onlyBadClashes(false)
    {}

    void ContactDB::set_atom_data(af::shared<int> flags, af::shared<int> model, af::shared<int> chain,
                                  af::shared<int> atomGroup, af::shared<int> altloc)
    {
      m_flags.assign(flags.begin(), flags.end());
      m_model.assign(model.begin(), model.end());
      m_chain.assign(chain.begin(), chain.end());
      m_group.assign(atomGroup.begin(), atomGroup.end());
      m_altloc.assign(altloc.begin(), altloc.end());
    }

    void ContactDB::set_bonds(af::shared<int> offsets, af::shared<atom> targets)
    {
      m_bondOffsets.assign(offsets.begin(), offsets.end());
      m_bondTargets.assign(targets.begin(), targets.end());
    }

    void ContactDB::set_params(double probeRadius, double maxRadius, double minOccupancy,
                               int excludedBondChainLength, bool includeMainchainMainchain,
                               bool includeWaterWater, double overlapScale, bool separateWorseClashes,
                               bool reportHBonds, bool reportClashes, bool reportVdws, bool onlyBadClashes)
    {
      m_probeRadius = probeRadius; m_maxRadius = maxRadius; m_minOccupancy = minOccupancy;
      m_excludedBondChainLength = excludedBondChainLength;
      m_includeMainchainMainchain = includeMainchainMainchain; m_includeWaterWater = includeWaterWater;
      m_overlapScale = overlapScale; m_separateWorseClashes = separateWorseClashes;
      m_reportHBonds = reportHBonds; m_reportClashes = reportClashes; m_reportVdws = reportVdws;
      m_onlyBadClashes = onlyBadClashes;
    }

    void ContactDB::build(af::shared<atom> sources)
    {
      m_sources.assign(sources.begin(), sources.end());
      m_records.assign(m_sources.size(), std::vector<PairRecord>());
      for (size_t i = 0; i < m_sources.size(); i++) {
        score_atom(m_sources[i], m_records[i]);
      }
    }

    boost::python::tuple ContactDB::pairs() const
    {
      af::shared<int> src, cause, counts;
      af::shared<double> gap, bumpGap, typeGap;
      af::shared<bool> hasBump;
      af::shared<Point> bumpLoc;
      for (size_t i = 0; i < m_sources.size(); i++) {
        for (size_t j = 0; j < m_records[i].size(); j++) {
          PairRecord const& r = m_records[i][j];
          src.push_back(static_cast<int>(i));
          cause.push_back(static_cast<int>(r.cause));
          gap.push_back(r.minGap);
          for (int t = 0; t < N_TYPES; t++) {
            counts.push_back(static_cast<int>(r.counts[t]));
            typeGap.push_back(r.typeMinGap[t]);
          }
          hasBump.push_back(r.hasBump);
          bumpGap.push_back(r.bumpGap);
          bumpLoc.push_back(r.bumpLoc);
        }
      }
      return boost::python::make_tuple(src, cause, gap, counts, hasBump, bumpGap, bumpLoc, typeGap);
    }

    unsigned ContactDB::idx(atom const& a) const
    {
      unsigned i = a.data->i_seq;
      if (i >= m_flags.size()) { throw std::out_of_range("ContactDB: atom i_seq beyond per-atom arrays"); }
      return i;
    }

    bool ContactDB::compatible(atom const& a, atom const& b) const
    {
      int x = m_altloc[idx(a)], y = m_altloc[idx(b)];
      return x == 0 || y == 0 || x == y;
    }

    // Port of mmtbx.probe.Helpers.getAtomsWithinNBonds.
    std::vector<atom> ContactDB::excluded_for(atom const& a, unsigned N, unsigned nonHydrogenN) const
    {
      Point const& aLoc = a.data->xyz;
      double aRad = m_info.getMappingFor(a).getVdwRadius();
      bool atomIsHydrogen = a.element_is_hydrogen();
      if (m_info.getMappingFor(a).getIsDummyHydrogen()) { N = 1; }
      std::vector<atom> found(1, a);   // a handful of atoms, so linear membership tests
      for (unsigned i = 0; i < N; i++) {
        size_t levelSize = found.size();
        for (size_t k = 0; k < levelSize; k++) {
          unsigned s = idx(found[k]);
          if (s + 1 >= m_bondOffsets.size()) { throw std::out_of_range("ContactDB: bond offsets too short"); }
          for (int b = m_bondOffsets[s]; b < m_bondOffsets[s + 1]; b++) {
            atom const& n = m_bondTargets[b];
            if (contains(found, n.data->i_seq)) { continue; }
            if (i < nonHydrogenN || atomIsHydrogen || n.element_is_hydrogen()) {
              if (compatible(a, n)) {
                Point const& nLoc = n.data->xyz;
                double nRad = m_info.getMappingFor(n).getVdwRadius();
                double dx = nLoc[0] - aLoc[0], dy = nLoc[1] - aLoc[1], dz = nLoc[2] - aLoc[2];
                double d = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (d <= aRad + nRad + 2 * m_probeRadius) {
                  found.push_back(n);
                }
              }
            }
          }
        }
      }
      found.erase(found.begin());
      return found;
    }

    // Port of probe2.Program._generate_interaction_dots for one source atom.
    void ContactDB::score_atom(atom const& src, std::vector<PairRecord>& out)
    {
      unsigned si = idx(src);
      int sf = m_flags[si];
      if ((sf & IGNORE) || src.data->occ < m_minOccupancy) { return; }
      af::shared<atom> nearby = m_grid->neighbors(src.data->xyz, 0.00001, m_maxRadius);
      double srcVdw = m_info.getMappingFor(src).getVdwRadius();
      Point const& sXYZ = src.data->xyz;

      std::vector<atom> atomSet;
      for (size_t k = 0; k < nearby.size(); k++) {
        atom const& n = nearby[k];
        unsigned ni = idx(n);
        int nf = m_flags[ni];
        Point const& nXYZ = n.data->xyz;
        double dx = nXYZ[0] - sXYZ[0], dy = nXYZ[1] - sXYZ[1], dz = nXYZ[2] - sXYZ[2];
        double d_sq = dx*dx + dy*dy + dz*dz;
        double cutoff = m_info.getMappingFor(n).getVdwRadius() + srcVdw + 2 * m_probeRadius;
        if (!(d_sq <= cutoff * cutoff)) { continue; }
        if (!m_includeMainchainMainchain && (sf & MAINCHAIN) && (nf & MAINCHAIN)
            && !((sf & HET) || (nf & HET)) && m_chain[si] == m_chain[ni]) { continue; }
        if (nf & IGNORE) { continue; }
        else if (n.data->occ < m_minOccupancy) { continue; }
        else if ((sf & WATER) && (nf & WATER)) {
          if (!m_includeWaterWater) { continue; }
          else if (m_group[si] == m_group[ni]) { continue; }
        }
        else if (!compatible(src, n)) { continue; }
        else if (m_model[si] != m_model[ni]) { continue; }
        atomSet.push_back(n);
      }
      if (atomSet.empty()) { return; }

      std::vector<atom> excluded = excluded_for(src, m_excludedBondChainLength, 3);
      if (m_info.getMappingFor(src).getIsDummyHydrogen()) {
        af::shared<atom> nearbyPhantoms = m_phantoms->neighbors(src.data->xyz, 0.00001, m_maxRadius);
        for (size_t k = 0; k < atomSet.size(); k++) {
          if (!m_info.getMappingFor(atomSet[k]).getIsAcceptor()) { excluded.push_back(atomSet[k]); }
        }
        excluded.insert(excluded.end(), nearbyPhantoms.begin(), nearbyPhantoms.end());
      }

      std::sort(atomSet.begin(), atomSet.end(), by_iseq);
      af::shared<atom> interacting, exclude;
      for (size_t k = 0; k < excluded.size(); k++) { exclude.push_back(excluded[k]); }
      for (size_t k = 0; k < atomSet.size(); k++) {
        if (contains(excluded, atomSet[k].data->i_seq)) { continue; }
        if (!interacting.empty() && interacting.back().data->i_seq == atomSet[k].data->i_seq) { continue; }
        interacting.push_back(atomSet[k]);
      }

      af::shared<Point> const& dots = m_dotCache.get_sphere(srcVdw).dots();
      af::shared<DotScorer::CheckDotResult> results =
        m_scorer.check_dots(src, dots, m_probeRadius, interacting, exclude, m_overlapScale);
      bool srcWild = m_altloc[si] == 0;
      for (size_t k = 0; k < results.size(); k++) {
        DotScorer::CheckDotResult const& res = results[k];
        unsigned ci = idx(res.cause);
        if (!(m_flags[ci] & TARGET)) { continue; }
        if (res.overlapType != DotScorer::HydrogenBond) {
          if (srcWild || m_altloc[ci] == 0) {
            if (src.data->occ + res.cause.data->occ <= 1) { continue; }
          }
        }
        DotScorer::InteractionType itype = m_scorer.interaction_type(res.overlapType, res.gap, m_separateWorseClashes);
        if (itype == DotScorer::Invalid) { continue; }
        bool show;
        if (!m_onlyBadClashes) {
          show = (m_reportHBonds && res.overlapType == DotScorer::HydrogenBond)
              || (m_reportClashes && res.overlapType == DotScorer::Clash)
              || (m_reportVdws && res.overlapType == DotScorer::NoOverlap);
        } else {
          show = m_reportClashes && (itype == DotScorer::Bump || itype == DotScorer::BadBump);
        }
        if (!show) { continue; }
        PairRecord* rec = 0;
        for (size_t r = 0; r < out.size(); r++) {
          if (out[r].cause == ci) { rec = &out[r]; break; }
        }
        if (!rec) {
          PairRecord fresh;
          fresh.cause = ci; fresh.minGap = 1e9; fresh.hasBump = false; fresh.bumpGap = 1e9;
          std::fill(fresh.counts, fresh.counts + N_TYPES, 0u);
          std::fill(fresh.typeMinGap, fresh.typeMinGap + N_TYPES, 1e9);
          out.push_back(fresh);
          rec = &out.back();
        }
        rec->counts[itype]++;
        if (res.gap < rec->minGap) { rec->minGap = res.gap; }
        if (res.gap < rec->typeMinGap[itype]) { rec->typeMinGap[itype] = res.gap; }
        if (itype == DotScorer::Bump && (!rec->hasBump || res.gap < rec->bumpGap)) {
          rec->hasBump = true;
          rec->bumpGap = res.gap;
          // Same arithmetic as probe2.py: loc = srcXYZ + dotOffset, per component.
          rec->bumpLoc = Point(sXYZ[0] + res.dotOffset[0], sXYZ[1] + res.dotOffset[1],
                               sXYZ[2] + res.dotOffset[2]);
        }
      }
    }

  }
}
