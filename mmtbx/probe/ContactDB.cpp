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
#include <boost/python/stl_iterator.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

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

      // Bonds as ContactDB::set_bonds() stores them, indexed by i_seq.
      struct TableNeighbors {
        std::vector<int> const& offsets;
        std::vector<atom> const& targets;
        TableNeighbors(std::vector<int> const& o, std::vector<atom> const& t) : offsets(o), targets(t) {}
        void operator()(atom const& x, std::vector<atom>& out) const {
          unsigned s = x.data->i_seq;
          if (s + 1 >= offsets.size()) { throw std::out_of_range("ContactDB: bond offsets too short"); }
          out.insert(out.end(), targets.begin() + offsets[s], targets.begin() + offsets[s + 1]);
        }
      };

      struct TableCompatible {
        std::vector<int> const& altloc;
        explicit TableCompatible(std::vector<int> const& a) : altloc(a) {}
        bool operator()(atom const& a, atom const& b) const {
          if (a.data->i_seq >= altloc.size() || b.data->i_seq >= altloc.size()) {
            throw std::out_of_range("ContactDB: atom i_seq beyond per-atom arrays");
          }
          int x = altloc[a.data->i_seq], y = altloc[b.data->i_seq];
          return x == 0 || y == 0 || x == y;
        }
      };

      struct MappingNeighbors {
        boost::python::object const& lists;
        explicit MappingNeighbors(boost::python::object const& l) : lists(l) {}
        void operator()(atom const& x, std::vector<atom>& out) const {
          boost::python::object nbrs = lists[boost::python::object(x)];
          out.insert(out.end(), boost::python::stl_input_iterator<atom>(nbrs),
                     boost::python::stl_input_iterator<atom>());
        }
      };

      // Helpers.compatibleConformations
      struct AltlocCompatible {
        static bool blank(atom const& a) {
          boost::shared_ptr<iotbx::pdb::hierarchy::atom_group_data> p = a.parent_ptr();
          return !p || p->altloc.empty() || p->altloc == " ";
        }
        bool operator()(atom const& a, atom const& b) const {
          return blank(a) || blank(b) || a.parent_ptr()->altloc == b.parent_ptr()->altloc;
        }
      };
    }

    boost::python::list getAtomsWithinNBonds(atom const& a, boost::python::object const& bondedNeighborLists,
                                             ExtraAtomInfoMap const& info, double probeRad, int N,
                                             double nonHydrogenN)
    {
      std::vector<atom> found = atomsWithinNBonds(a, MappingNeighbors(bondedNeighborLists),
                                                  AltlocCompatible(), info, probeRad, N, nonHydrogenN);
      boost::python::list ret;
      for (size_t i = 0; i < found.size(); i++) { ret.append(found[i]); }
      return ret;
    }

    ContactDB::ContactDB(ExtraAtomInfoMap const& info, DotScorer const& scorer, SpatialQuery& grid,
                         SpatialQuery& phantoms, DotSphereCache const& dotCache)
      : m_info(info), m_scorer(scorer), m_grid(&grid), m_phantoms(&phantoms), m_dotCache(dotCache)
      , m_probeRadius(0.25), m_maxRadius(0), m_minOccupancy(0), m_overlapScale(0.5)
      , m_excludedBondChainLength(4), m_includeMainchainMainchain(true), m_includeWaterWater(false)
      , m_separateWorseClashes(false), m_reportHBonds(true), m_reportClashes(true)
      , m_reportVdws(true), m_onlyBadClashes(false), m_keepDots(false)
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
      m_dots.clear();
      for (size_t i = 0; i < m_sources.size(); i++) {
        score_atom(m_sources[i], static_cast<unsigned>(i), m_records[i]);
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

    boost::python::tuple ContactDB::dots() const
    {
      af::shared<int> src, cause, itype, overlapType;
      af::shared<double> gap;
      af::shared<Point> loc, spike;
      std::string ptmaster;
      ptmaster.reserve(m_dots.size());
      for (size_t i = 0; i < m_dots.size(); i++) {
        DotRecord const& d = m_dots[i];
        src.push_back(static_cast<int>(d.src));
        cause.push_back(static_cast<int>(d.cause));
        itype.push_back(d.itype);
        overlapType.push_back(d.overlapType);
        gap.push_back(d.gap);
        loc.push_back(d.loc);
        spike.push_back(d.spike);
        ptmaster.push_back(d.ptmaster);
      }
      return boost::python::make_tuple(src, cause, itype, overlapType, gap, loc, spike, ptmaster);
    }

    void ContactDB::keep_dot(unsigned srcIndex, int sf, unsigned ci, Point const& sXYZ, double srcVdw,
                             int itype, DotScorer::CheckDotResult const& res)
    {
      // Main/side chain class of the pair, then the dot and spike-end locations.
      int cf = m_flags[ci];
      bool srcMC = sf & MAINCHAIN, srcSC = sf & SIDECHAIN, causeMC = cf & MAINCHAIN,
           causeSC = cf & SIDECHAIN, noHet = !(sf & HET) && !(cf & HET);
      char ptmaster = ' ';
      if (srcMC && causeMC) { if (noHet) { ptmaster = 'M'; } }
      else if (srcSC && causeSC) { if (noHet) { ptmaster = 'S'; } }
      else if ((srcMC && causeSC) || (srcSC && causeMC)) { if (noHet) { ptmaster = 'P'; } }
      else { ptmaster = 'O'; }

      Point const& d = res.dotOffset;
      DotRecord rec;
      rec.src = srcIndex; rec.cause = ci; rec.itype = itype;
      rec.overlapType = res.overlapType; rec.ptmaster = ptmaster; rec.gap = res.gap;
      rec.loc = Point(sXYZ[0] + d[0], sXYZ[1] + d[1], sXYZ[2] + d[2]);
      double dvLen = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      if (dvLen > 0) {
        double spikeScale = (srcVdw - res.overlap) / dvLen;
        rec.spike = Point(sXYZ[0] + d[0]*spikeScale, sXYZ[1] + d[1]*spikeScale, sXYZ[2] + d[2]*spikeScale);
      } else {
        rec.spike = rec.loc;
      }
      m_dots.push_back(rec);
    }

    unsigned ContactDB::idx(atom const& a) const
    {
      unsigned i = a.data->i_seq;
      if (i >= m_flags.size()) { throw std::out_of_range("ContactDB: atom i_seq beyond per-atom arrays"); }
      return i;
    }

    bool ContactDB::compatible(atom const& a, atom const& b) const
    {
      return TableCompatible(m_altloc)(a, b);
    }

    std::vector<atom> ContactDB::excluded_for(atom const& a, unsigned N, unsigned nonHydrogenN) const
    {
      return atomsWithinNBonds(a, TableNeighbors(m_bondOffsets, m_bondTargets), TableCompatible(m_altloc),
                               m_info, m_probeRadius, static_cast<int>(N), nonHydrogenN);
    }

    void ContactDB::score_atom(atom const& src, unsigned srcIndex, std::vector<PairRecord>& out)
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
        if (m_keepDots) { keep_dot(srcIndex, sf, ci, sXYZ, srcVdw, itype, res); }
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
          // Placed as keep_dot() places dots.
          rec->bumpLoc = Point(sXYZ[0] + res.dotOffset[0], sXYZ[1] + res.dotOffset[1],
                               sXYZ[2] + res.dotOffset[2]);
        }
      }
    }

  }
}
