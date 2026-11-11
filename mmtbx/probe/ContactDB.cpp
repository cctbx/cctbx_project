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
      , m_extremaDots(false), m_extremaMinOverlap(0.4)
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
          src.push_back(static_cast<int>(idx(m_sources[i])));
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

    void ContactDB::keep_dot(unsigned si, int sf, unsigned ci, Point const& sXYZ, double srcVdw,
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
      rec.src = si; rec.cause = ci; rec.itype = itype;
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

    // ---- extrema_dots: dots at the exact worst-overlap points of close pairs --------------
    //
    // With set_extrema_dots(true), score_atom() also scores a few dots that are placed
    // analytically for each close pair (source s, target t), after the dot sphere.
    //
    // Why: the worst gap reported for a pair is the minimum over the dots kept for it. On a
    // finite dot sphere that minimum is above the true one (less overlap), by an amount that
    // depends on where the dots happen to fall, so it also changes when the model is rotated,
    // and pairs close to the clash cutoff are counted or not depending on orientation.
    //
    // Which pairs: those whose centre overlap r_s + r_t - d_t can reach the clash cutoff.
    //
    // Geometry. A dot is the surface point p = c_s + r_s n of the source, n a unit vector. With
    // u_t the unit vector from c_s to t's centre c_t, and d_t their distance, the distance from
    // the dot to c_t is
    //
    //     rho^2 = |p - c_t|^2 = r_s^2 + d_t^2 - 2 r_s d_t (n.u_t),
    //
    // and the gap to t is rho - r_t. So the gap depends on n only through y = n.u_t, and it gets
    // more negative as y grows: the worst gap for (s, t) is reached at the largest n.u_t among
    // the dots kept for t. Conversely, y and rho determine each other:
    //
    //     y = Y(rho) = (r_s^2 + d_t^2 - rho^2) / (2 r_s d_t).                          (1)
    //
    // The per-dot rules of DotScorer::check_dot(), written as regions of the unit sphere:
    //   - the dot touches t: rho < r_t, i.e. n.u_t > y0 (a cap around u_t);
    //   - it is outside every excluded atom m: |p - c_m| >= r_m, i.e.
    //     n.u_m <= lam_m = (r_s^2 + d_m^2 - r_m^2) / (2 r_s d_m). A cap around u_m is removed;
    //     its edge is a circle, the "exclusion circle";
    //   - t is its nearest interacting atom: for every other atom a that could take the dot,
    //     gap(a) >= gap(t), i.e. |p - c_a| >= rho + (r_a - r_t). The edge of that region (the
    //     "border" with a) is the curve |p - c_a| - |p - c_t| = r_a - r_t on the sphere, in
    //     general not a circle. At fixed rho, i.e. on a fixed latitude circle about u_t, it is a
    //     cap condition again: n.u_a <= A(rho), with
    //
    //     A(rho) = (r_s^2 + d_a^2 - (rho + r_a - r_t)^2) / (2 r_s d_a).                 (2)
    //
    // Where the maximum is. We maximise the linear function n.u_t over a region of the sphere
    // bounded by exclusion circles and borders. At the maximum either no constraint is active
    // (then n = u_t), or one boundary is active and the latitude circle about u_t is tangent to
    // it (the point of that boundary nearest u_t), or two boundaries meet there (a corner). So
    // the worst kept gap is reached at one of these points:
    //   1. the point facing t, n = u_t;
    //   2. the point of an exclusion circle nearest u_t. It lies on the great circle through u_t
    //      and u_m: n = lam_m u_m + sqrt(1 - lam_m^2) e, where e is the unit vector along
    //      u_t - (u_t.u_m) u_m. Closed form;
    //   3. a crossing of two exclusion circles: n.u_1 = lam_1, n.u_2 = lam_2, |n| = 1, so
    //      n = x u_1 + y u_2 +/- z (u_1 x u_2). Closed form;
    //   4. the point of a border nearest u_t. By symmetry it lies on the great circle through
    //      u_t and u_a: n = Y u_t + sqrt(1 - Y^2) w, where w is the unit vector along u_a - c u_t
    //      and c = u_a.u_t. Then n.u_a = c Y + s sqrt(1 - Y^2), with s = sqrt(1 - c^2), and
    //      setting it equal to A gives s^2 (1 - Y^2) = (A - c Y)^2: a quartic in rho, since Y
    //      and A are quadratics in rho (1, 2);
    //   5. a crossing of an exclusion circle and a border: n.u_t = Y(rho), n.u_m = lam_m,
    //      n.u_a = A(rho). Solving this 3x3 linear system for n gives n = v0 + v1 rho + v2 rho^2,
    //      and |n|^2 = 1 is a quartic in rho;
    //   6. a crossing of two borders (e.g. both phantom hydrogens of a water): n.u_t = Y(rho),
    //      n.u_a = A(rho), n.u_b = B(rho), solved the same way.
    // The real roots of each quartic are isolated between the roots of its derivatives.
    // Degenerate geometry has no single nearest point, and is sampled instead (kExtScan samples,
    // then bisection):
    //   - an atom centred on the axis (u_a = +-u_t): its border is itself a latitude circle
    //     about u_t;
    //   - an exclusion circle centred on the axis: all its points are equally near;
    //   - three nearly coplanar directions: the 3x3 system is singular.
    //
    // Each point is moved off its defining boundaries by the smallest step for which
    // DotScorer::check_dot() itself keeps it for t, then scored like any other dot; only
    // bumps (Clash) are kept. These are real surface points: they can only reveal overlap
    // that is there.
    namespace {
      const int kExtScan = 512;
      const double kExtPi = 3.14159265358979323846;
      const double kExtStepMin = 1.0e-9, kExtStepMax = 1.0e-3;

      struct ExtNbr {          // one neighbor of the source atom
        double u[3];           // unit vector from the source
        double d, r;           // centre distance and radius
        double lam;            // its sphere meets the source surface at n.u = lam
      };

      double eDot(const double* a, const double* b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
      bool eUnit(double* v)
      {
        double s = std::sqrt(eDot(v, v));
        if (s < 1.0e-12) { return false; }
        v[0] /= s; v[1] /= s; v[2] /= s;
        return true;
      }
      void eCross(const double* a, const double* b, double* c)
      {
        c[0] = a[1]*b[2] - a[2]*b[1]; c[1] = a[2]*b[0] - a[0]*b[2]; c[2] = a[0]*b[1] - a[1]*b[0];
      }
      void eTangent(const double* v, const double* n, double* g)
      {
        double c = eDot(v, n);
        for (int k = 0; k < 3; k++) { g[k] = v[k] - c*n[k]; }
        if (!eUnit(g)) { g[0] = g[1] = g[2] = 0; }
      }

      /// One directed pair: the source, its target t, and the constraints that reach t's cap.
      struct ExtPair {
        DotScorer* scorer;
        iotbx::pdb::hierarchy::atom const* src;
        iotbx::pdb::hierarchy::atom const* tAtom;
        af::shared<iotbx::pdb::hierarchy::atom> const* interacting;
        af::shared<iotbx::pdb::hierarchy::atom> const* exclude;
        double probeRadius, overlapScale, rs, y0;
        ExtNbr t;
        std::vector<ExtNbr> msk, cmp;
        std::vector<Point>* out;

        double rho(const double* n, ExtNbr const& a) const
        {
          double q = rs*rs + a.d*a.d - 2*rs*a.d*eDot(n, a.u);
          return std::sqrt(q > 0 ? q : 0);
        }
        /// gap(a) - gap(t) at n: > 0 when t, not a, is the nearer atom
        double border(const double* n, ExtNbr const& a) const
        {
          return (rho(n, a) - a.r) - (rho(n, t) - t.r);
        }
        /// Does probe2 keep the dot r_s n for t (as a dot with overlap)?
        bool allowed(const double* n) const
        {
          if (eDot(n, t.u) < y0) { return false; }
          Point off(rs*n[0], rs*n[1], rs*n[2]);
          DotScorer::CheckDotResult res =
            scorer->check_dot(*src, off, probeRadius, *interacting, *exclude, overlapScale);
          // A Clash for t (for an H-bond pair: one that is too close, its gap already shifted by
          // the H-bond cutoff). A normal H-bond dot is never a bad overlap, so it does not count.
          return res.overlapType == DotScorer::Clash && res.cause.data == tAtom->data && res.gap < 0;
        }
        /// n0 lies on the boundaries of masks m1, m2 and/or borders a1, a2 (at most two).
        /// Step into the allowed side, along x with g.x = 1 for each inward normal g, by the
        /// smallest step that allowed() accepts, and keep that point.
        void emit(const double* n0, ExtNbr const* m1, ExtNbr const* m2,
                  ExtNbr const* a1, ExtNbr const* a2) const
        {
          double n[3] = {n0[0], n0[1], n0[2]}, p[3], x[3] = {0, 0, 0}, g[4][3], gr[3];
          int ng = 0;
          if (!eUnit(n)) { return; }
          ExtNbr const* m[2] = {m1, m2};
          ExtNbr const* a[2] = {a1, a2};
          for (int j = 0; j < 2; j++) {
            if (m[j]) { eTangent(m[j]->u, n, g[ng]); for (int k = 0; k < 3; k++) { g[ng][k] = -g[ng][k]; } ng++; }
          }
          for (int j = 0; j < 2; j++) {
            if (a[j]) {
              double ka = -rs*a[j]->d / std::max(rho(n, *a[j]), 1e-12), kt = rs*t.d / std::max(rho(n, t), 1e-12);
              for (int k = 0; k < 3; k++) { gr[k] = ka*a[j]->u[k] + kt*t.u[k]; }
              eTangent(gr, n, g[ng]); ng++;
            }
          }
          bool ok = false;
          if (ng == 0) {
            for (int k = 0; k < 3; k++) { p[k] = n[k]; }
            ok = allowed(p);
          } else {
            if (ng == 1) { for (int k = 0; k < 3; k++) { x[k] = g[0][k]; } }
            else {
              double c = eDot(g[0], g[1]);
              if (c < -1.0 + 1e-12) { return; }   // boundaries touch: no area
              for (int k = 0; k < 3; k++) { x[k] = (g[0][k] + g[1][k]) / (1.0 + c); }
            }
            for (double step = kExtStepMin; step <= kExtStepMax && !ok; step *= 2) {
              for (int k = 0; k < 3; k++) { p[k] = n[k] + step*x[k]; }
              ok = eUnit(p) && allowed(p);
            }
          }
          if (ok) { out->push_back(Point(rs*p[0], rs*p[1], rs*p[2])); }
        }

        // ---- polynomials in rho, ascending coefficients
        static double peval(const double* c, int deg, double x)
        {
          double v = 0;
          for (int i = deg; i >= 0; i--) { v = v*x + c[i]; }
          return v;
        }
        /// real roots in [lo, hi], isolated between the roots of the derivative
        static int proots(const double* c, int deg, double lo, double hi, double* out)
        {
          double d[4], cp[6], scale = 0;
          int n = 0;
          while (deg > 0 && c[deg] == 0.0) { deg--; }
          if (deg <= 0 || hi < lo) { return 0; }
          if (deg == 1) {
            double m = -c[0] / c[1];
            if (m >= lo && m <= hi) { out[n++] = m; }
            return n;
          }
          double big = std::max(1.0, std::max(std::fabs(lo), std::fabs(hi)));
          for (int i = 0; i <= deg; i++) { scale += std::fabs(c[i]) * std::pow(big, i); }
          for (int i = 1; i <= deg; i++) { d[i-1] = i*c[i]; }
          int ncp = proots(d, deg - 1, lo, hi, cp + 1);
          cp[0] = lo; cp[ncp + 1] = hi;
          for (int i = 0; i <= ncp; i++) {
            double a = cp[i], b = cp[i+1], fa = peval(c, deg, a), fb = peval(c, deg, b);
            if (std::fabs(fa) <= 1e-13*scale) {
              if (n == 0 || a - out[n-1] > 1e-12) { out[n++] = a; }
            } else if (std::fabs(fb) > 1e-13*scale && (fa < 0) != (fb < 0)) {
              for (int it = 0; it < 200 && b - a > 1e-15; it++) {
                double m = 0.5*(a + b), fm = peval(c, deg, m);
                if ((fm < 0) == (fa < 0)) { a = m; fa = fm; } else { b = m; }
              }
              out[n++] = 0.5*(a + b);
            }
          }
          double fh = peval(c, deg, hi);
          if (std::fabs(fh) <= 1e-13*scale && (n == 0 || hi - out[n-1] > 1e-12)) { out[n++] = hi; }
          return n;
        }
        static void pmul2(const double* a, const double* b, double* o)
        {
          for (int i = 0; i < 5; i++) { o[i] = 0; }
          for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) { o[i+j] += a[i]*b[j]; } }
        }
        /// n.u_t = Y(rho), eq. (1) above, and on the border with a, n.u_a = A(rho), eq. (2);
        /// coefficients of 1, rho, rho^2
        void ypoly(double* y) const
        {
          double Q = 2*rs*t.d;
          y[0] = (rs*rs + t.d*t.d) / Q; y[1] = 0; y[2] = -1.0 / Q;
        }
        void apoly(ExtNbr const& a, double* A) const
        {
          double D = a.r - t.r, Q = 2*rs*a.d;
          A[0] = (rs*rs + a.d*a.d - D*D) / Q; A[1] = -2*D / Q; A[2] = -1.0 / Q;
        }
        void rhoRange(ExtNbr const* a, ExtNbr const* b, double& lo, double& hi) const
        {
          lo = std::fabs(rs - t.d);
          if (a) { lo = std::max(lo, t.r - a->r); }
          if (b) { lo = std::max(lo, t.r - b->r); }
          hi = t.r;
        }
        /// n with n.R_i = p_i(rho) (quadratics): n = v0 + v1 rho + v2 rho^2; |n|^2 = 1 is a quartic.
        static bool solve3(const double* R0, const double* R1, const double* R2,
                           const double* p0, const double* p1, const double* p2, double v[3][3], double* q)
        {
          double C[3][3];
          eCross(R1, R2, C[0]); eCross(R2, R0, C[1]); eCross(R0, R1, C[2]);
          double det = eDot(R0, C[0]);
          if (std::fabs(det) < 1e-6) { return false; }
          const double* p[3] = {p0, p1, p2};
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              v[j][k] = 0;
              for (int i = 0; i < 3; i++) { v[j][k] += p[i][j]*C[i][k]; }
              v[j][k] /= det;
            }
          }
          q[0] = eDot(v[0], v[0]) - 1; q[1] = 2*eDot(v[0], v[1]);
          q[2] = eDot(v[1], v[1]) + 2*eDot(v[0], v[2]); q[3] = 2*eDot(v[1], v[2]); q[4] = eDot(v[2], v[2]);
          return true;
        }
        static int circleCircle(const double* u1, double l1, const double* u2, double l2, double n[2][3])
        {
          double G = eDot(u1, u2), den = 1 - G*G, cc[3], base[3];
          if (den < 1e-14) { return 0; }
          double x = (l1 - l2*G) / den, y = (l2 - l1*G) / den;
          eCross(u1, u2, cc);
          for (int q = 0; q < 3; q++) { base[q] = x*u1[q] + y*u2[q]; }
          double z2 = (1 - eDot(base, base)) / eDot(cc, cc);
          if (z2 < 0) { return 0; }
          for (int q = 0; q < 3; q++) { n[0][q] = base[q] + std::sqrt(z2)*cc[q]; n[1][q] = base[q] - std::sqrt(z2)*cc[q]; }
          return 2;
        }
        bool borderAt(ExtNbr const& a, double r, double n[2][3]) const
        {
          double ra = r + a.r - t.r;
          if (ra < 0) { return false; }
          return circleCircle(t.u, (rs*rs + t.d*t.d - r*r) / (2*rs*t.d),
                              a.u, (rs*rs + a.d*a.d - ra*ra) / (2*rs*a.d), n) == 2;
        }

        // ---- fallbacks for near-degenerate geometry
        /// curve kind 0: great circle cos s u_t + sin s w; kind 1: circle o + rho (cos s f1 + sin s f2)
        void curve(int kind, double s, const double* p, double* n) const
        {
          for (int k = 0; k < 3; k++) {
            n[k] = kind == 0 ? std::cos(s)*p[k] + std::sin(s)*p[3+k]
                             : p[k] + p[9]*(std::cos(s)*p[3+k] + std::sin(s)*p[6+k]);
          }
        }
        void scanCurve(ExtNbr const& a, ExtNbr const* m, int kind, const double* p, double s0, double s1) const
        {
          double n[3];
          curve(kind, s0, p, n);
          double sa = s0, ha = border(n, a);
          for (int i = 1; i <= kExtScan; i++) {
            double sb = s0 + (s1 - s0)*i/kExtScan;
            curve(kind, sb, p, n);
            double hb = border(n, a);
            if ((ha < 0) != (hb < 0)) {
              double lo = sa, hi = sb, hlo = ha;
              for (int it = 0; it < 60; it++) {
                double sm = 0.5*(lo + hi);
                curve(kind, sm, p, n);
                double hm = border(n, a);
                if ((hm < 0) == (hlo < 0)) { lo = sm; hlo = hm; } else { hi = sm; }
              }
              curve(kind, 0.5*(lo + hi), p, n);
              emit(n, m, 0, &a, 0);
            }
            sa = sb; ha = hb;
          }
        }
        void scanBorders(ExtNbr const& a, ExtNbr const& b) const
        {
          double lo, hi, n[2][3];
          rhoRange(&a, 0, lo, hi);
          if (hi <= lo) { return; }
          for (int br = 0; br < 2; br++) {
            bool ok0 = borderAt(a, lo, n);
            double r0 = lo, h0 = ok0 ? border(n[br], b) : 0;
            for (int i = 1; i <= kExtScan; i++) {
              double r1 = lo + (hi - lo)*i/kExtScan;
              bool ok1 = borderAt(a, r1, n);
              double h1 = ok1 ? border(n[br], b) : 0;
              if (ok0 && ok1 && ((h0 < 0) != (h1 < 0))) {
                double ra = r0, rb = r1, hlo = h0;
                for (int it = 0; it < 60; it++) {
                  double rm = 0.5*(ra + rb);
                  if (!borderAt(a, rm, n)) { break; }
                  double hm = border(n[br], b);
                  if ((hm < 0) == (hlo < 0)) { ra = rm; hlo = hm; } else { rb = rm; }
                }
                if (borderAt(a, 0.5*(ra + rb), n)) { emit(n[br], 0, 0, &a, &b); }
              }
              r0 = r1; h0 = h1; ok0 = ok1;
            }
          }
        }

        // ---- the six kinds
        void all() const
        {
          emit(t.u, 0, 0, 0, 0);                                              // 1.
          for (size_t j = 0; j < msk.size(); j++) {                           // 2.
            double e1[3], n[3], rho = std::sqrt(std::max(0.0, 1 - msk[j].lam*msk[j].lam));
            double c = eDot(t.u, msk[j].u);
            for (int k = 0; k < 3; k++) { e1[k] = t.u[k] - c*msk[j].u[k]; }
            if (eUnit(e1)) {
              for (int k = 0; k < 3; k++) { n[k] = msk[j].lam*msk[j].u[k] + rho*e1[k]; }
              emit(n, &msk[j], 0, 0, 0);
            } else {   // circle centred on the axis: every point is equally near
              double f1[3], f2[3], ref[3] = {1, 0, 0};
              if (std::fabs(msk[j].u[0]) > 0.9) { ref[0] = 0; ref[1] = 1; }
              eCross(msk[j].u, ref, f1); eUnit(f1); eCross(msk[j].u, f1, f2);
              for (double s = 0; s < 2*kExtPi - 1e-9; s += kExtPi/8) {
                for (int k = 0; k < 3; k++) { n[k] = msk[j].lam*msk[j].u[k] + rho*(std::cos(s)*f1[k] + std::sin(s)*f2[k]); }
                emit(n, &msk[j], 0, 0, 0);
              }
            }
          }
          for (size_t j = 0; j < msk.size(); j++) {                           // 3.
            for (size_t k = j + 1; k < msk.size(); k++) {
              double nn[2][3];
              if (circleCircle(msk[j].u, msk[j].lam, msk[k].u, msk[k].lam, nn)) {
                emit(nn[0], &msk[j], &msk[k], 0, 0);
                emit(nn[1], &msk[j], &msk[k], 0, 0);
              }
            }
          }
          double y[3], A[3], B[3], v[3][3], q[5], roots[8], n[3], lo, hi;
          ypoly(y);
          for (size_t j = 0; j < cmp.size(); j++) {                           // 4.
            ExtNbr const& a = cmp[j];
            double c = eDot(a.u, t.u), w[3], d1[3], q1[5], q2[5];
            for (int k = 0; k < 3; k++) { w[k] = a.u[k] - c*t.u[k]; }
            if (eUnit(w)) {
              double s = std::sqrt(std::max(0.0, 1 - c*c));
              apoly(a, A);
              for (int k = 0; k < 3; k++) { d1[k] = A[k] - c*y[k]; }
              pmul2(y, y, q1); pmul2(d1, d1, q2);
              for (int k = 0; k < 5; k++) { q[k] = -s*s*q1[k] - q2[k]; }
              q[0] += s*s;
              rhoRange(&a, 0, lo, hi);
              int nr = proots(q, 4, lo, hi, roots);
              for (int i = 0; i < nr; i++) {
                double yy = peval(y, 2, roots[i]), qq = (peval(A, 2, roots[i]) - c*yy) / s;
                for (int k = 0; k < 3; k++) { n[k] = yy*t.u[k] + qq*w[k]; }
                emit(n, 0, 0, &a, 0);
              }
            } else {   // a on the axis: the border is a circle around it, sample its azimuth
              double f1[3], f2[3], ref[3] = {1, 0, 0}, p[6];
              double angt = std::acos(std::max(-1.0, std::min(1.0, y0)));
              if (std::fabs(t.u[0]) > 0.9) { ref[0] = 0; ref[1] = 1; }
              eCross(t.u, ref, f1); eUnit(f1); eCross(t.u, f1, f2);
              for (int k = 0; k < 3; k++) { p[k] = t.u[k]; }
              for (double az = 0; az < kExtPi - 1e-9; az += kExtPi/8) {
                for (int k = 0; k < 3; k++) { p[3+k] = std::cos(az)*f1[k] + std::sin(az)*f2[k]; }
                scanCurve(a, 0, 0, p, -angt, angt);
              }
            }
          }
          for (size_t j = 0; j < msk.size(); j++) {                           // 5.
            ExtNbr const& m = msk[j];
            double lam[3] = {m.lam, 0, 0};
            for (size_t k = 0; k < cmp.size(); k++) {
              ExtNbr const& a = cmp[k];
              apoly(a, A);
              if (solve3(t.u, a.u, m.u, y, A, lam, v, q)) {
                rhoRange(&a, 0, lo, hi);
                int nr = proots(q, 4, lo, hi, roots);
                for (int i = 0; i < nr; i++) {
                  for (int kk = 0; kk < 3; kk++) { n[kk] = v[0][kk] + v[1][kk]*roots[i] + v[2][kk]*roots[i]*roots[i]; }
                  emit(n, &m, 0, &a, 0);
                }
              } else {   // nearly coplanar: bracket along the circle
                double p[10], ref[3] = {1, 0, 0}, c = eDot(t.u, m.u);
                for (int kk = 0; kk < 3; kk++) { p[kk] = m.lam*m.u[kk]; p[3+kk] = t.u[kk] - c*m.u[kk]; }
                if (!eUnit(p + 3)) {
                  if (std::fabs(m.u[0]) > 0.9) { ref[0] = 0; ref[1] = 1; }
                  eCross(m.u, ref, p + 3); eUnit(p + 3);
                }
                eCross(m.u, p + 3, p + 6);
                p[9] = std::sqrt(std::max(0.0, 1 - m.lam*m.lam));
                scanCurve(a, &m, 1, p, -kExtPi, kExtPi);
              }
            }
          }
          for (size_t j = 0; j < cmp.size(); j++) {                           // 6.
            for (size_t k = j + 1; k < cmp.size(); k++) {
              ExtNbr const& a = cmp[j];
              ExtNbr const& b = cmp[k];
              apoly(a, A); apoly(b, B);
              if (solve3(t.u, a.u, b.u, y, A, B, v, q)) {
                rhoRange(&a, &b, lo, hi);
                int nr = proots(q, 4, lo, hi, roots);
                for (int i = 0; i < nr; i++) {
                  for (int kk = 0; kk < 3; kk++) { n[kk] = v[0][kk] + v[1][kk]*roots[i] + v[2][kk]*roots[i]*roots[i]; }
                  emit(n, 0, 0, &a, &b);
                }
              } else {
                scanBorders(a, b);
                scanBorders(b, a);
              }
            }
          }
        }
      };

      ExtNbr makeNbr(Point const& sXYZ, double rs, Point const& xyz, double r)
      {
        ExtNbr e;
        e.u[0] = xyz[0] - sXYZ[0]; e.u[1] = xyz[1] - sXYZ[1]; e.u[2] = xyz[2] - sXYZ[2];
        e.d = std::sqrt(eDot(e.u, e.u));
        e.r = r;
        if (!eUnit(e.u)) { e.d = 0; }
        e.lam = (rs*rs + e.d*e.d - r*r) / (2*rs*(e.d > 0 ? e.d : 1.0));
        return e;
      }
    }

    af::shared<Point> ContactDB::extrema_dots(atom const& src, double srcVdw,
      af::shared<atom> const& interacting, af::shared<atom> const& exclude)
    {
      std::vector<Point> out;
      Point const& sXYZ = src.data->xyz;
      ExtPair P;
      P.scorer = &m_scorer; P.src = &src; P.interacting = &interacting; P.exclude = &exclude;
      P.probeRadius = m_probeRadius; P.overlapScale = m_overlapScale; P.rs = srcVdw; P.out = &out;
      std::vector<ExtNbr> ex, in;
      for (size_t k = 0; k < exclude.size(); k++) {
        ex.push_back(makeNbr(sXYZ, srcVdw, exclude[k].data->xyz, m_info.getMappingFor(exclude[k]).getVdwRadius()));
      }
      for (size_t k = 0; k < interacting.size(); k++) {
        in.push_back(makeNbr(sXYZ, srcVdw, interacting[k].data->xyz, m_info.getMappingFor(interacting[k]).getVdwRadius()));
      }
      for (size_t i = 0; i < interacting.size(); i++) {
        ExtNbr const& t = in[i];
        unsigned ti = idx(interacting[i]);
        if (!(m_flags[ti] & TARGET) || m_info.getMappingFor(interacting[i]).getIsDummyHydrogen() || t.d <= 0) { continue; }
        if (srcVdw + t.r - t.d < m_extremaMinOverlap) { continue; }   // cannot be a clash
        if (t.lam >= 1) { continue; }                                   // t inside src: no overlap
        P.t = t; P.tAtom = &interacting[i]; P.y0 = t.lam;
        double angt = std::acos(std::max(-1.0, t.lam));
        P.msk.clear(); P.cmp.clear();
        bool allMasked = false;
        for (size_t j = 0; j < ex.size() && !allMasked; j++) {
          ExtNbr const& c = ex[j];
          if (c.d <= 0 || c.lam >= 1) { continue; }
          if (c.lam <= -1) { allMasked = true; break; }
          double ang = std::acos(std::max(-1.0, std::min(1.0, eDot(c.u, t.u))));
          if (ang < angt + std::acos(c.lam) + 1e-9) { P.msk.push_back(c); }
        }
        if (allMasked) { continue; }
        for (size_t j = 0; j < in.size(); j++) {
          ExtNbr const& c = in[j];
          if (j == i || c.d <= 0 || c.lam >= 1) { continue; }
          double ang = std::acos(std::max(-1.0, std::min(1.0, eDot(c.u, t.u))));
          if (ang < angt + std::acos(std::max(c.lam, -1.0)) + 1e-9) { P.cmp.push_back(c); }
        }
        P.all();
      }
      af::shared<Point> ret;
      for (size_t k = 0; k < out.size(); k++) { ret.push_back(out[k]); }
      return ret;
    }

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
      auto add = [&](DotScorer::CheckDotResult const& res) {
        unsigned ci = idx(res.cause);
        if (!(m_flags[ci] & TARGET)) { return; }
        if (res.overlapType != DotScorer::HydrogenBond) {
          if (srcWild || m_altloc[ci] == 0) {
            if (src.data->occ + res.cause.data->occ <= 1) { return; }
          }
        }
        DotScorer::InteractionType itype = m_scorer.interaction_type(res.overlapType, res.gap, m_separateWorseClashes);
        if (itype == DotScorer::Invalid) { return; }
        bool show;
        if (!m_onlyBadClashes) {
          show = (m_reportHBonds && res.overlapType == DotScorer::HydrogenBond)
              || (m_reportClashes && res.overlapType == DotScorer::Clash)
              || (m_reportVdws && res.overlapType == DotScorer::NoOverlap);
        } else {
          show = m_reportClashes && (itype == DotScorer::Bump || itype == DotScorer::BadBump);
        }
        if (!show) { return; }
        if (m_keepDots) { keep_dot(si, sf, ci, sXYZ, srcVdw, itype, res); }
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
      };
      for (size_t k = 0; k < results.size(); k++) { add(results[k]); }
      if (m_extremaDots && !m_info.getMappingFor(src).getIsDummyHydrogen()) {
        af::shared<DotScorer::CheckDotResult> extra = m_scorer.check_dots(src,
          extrema_dots(src, srcVdw, interacting, exclude), m_probeRadius, interacting, exclude, m_overlapScale);
        for (size_t k = 0; k < extra.size(); k++) {
          // bumps only (Bump/BadBump, as -EXTREMA in C Probe), not small overlaps
          DotScorer::InteractionType xt = m_scorer.interaction_type(extra[k].overlapType, extra[k].gap,
                                                                    m_separateWorseClashes);
          if (xt == DotScorer::Bump || xt == DotScorer::BadBump) { add(extra[k]); }
        }
      }
    }

  }
}
