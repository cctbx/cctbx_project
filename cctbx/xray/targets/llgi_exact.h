#ifndef CCTBX_XRAY_TARGETS_LLGI_EXACT_H
#define CCTBX_XRAY_TARGETS_LLGI_EXACT_H

#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <limits>
#include <scitbx/constants.h>
#include <scitbx/math/quadrature.h>
#include <cctbx/french_wilson.h>
#include <boost/math/special_functions/bessel.hpp>
#include <algorithm>
#include <cmath>
#include <vector>

namespace cctbx { namespace xray { namespace targets { namespace llgi_exact {

  /*! Exact log-likelihood-gain on intensities (TEPS == 1), for the
      reflections where the Rice approximation with moment-matched
      (Dobs, Eeff) is not accurate enough (no Rice solution, or measurement
      error not small compared with the model error; see class hybrid).

      With J = |E|^2, the observed intensity on the E^2 scale is
      eo_sq +/- sig (Gaussian), the model prior for E is the Rice (or,
      centric, Gaussian) distribution with mean a*ec and variance
      S = 1 - a^2, and the null hypothesis is Wilson. Writing b = 1
      (acentric) or 1/2 (centric) and nu = a*ec,

        LLGI = b*(ln(1/S) - nu^2/S) + M(kappa, lambda) - M(0, b)
        M(kappa, lambda) = ln int_0^inf J^p exp(-(J-eo_sq)^2/(2 sig^2)
                                                 - lambda J) B(kappa sqrt J) dJ
        lambda = b/S,  kappa = 2 b nu/S,
        p = 0, B = I0 (acentric);  p = -1/2, B = cosh (centric).

      Derivatives of M are moments under the "tilted" posterior weights of
      the same integrand:
        dM/dkappa = <sqrt(J) B'/B>,  dM/dlambda = -<J>,
      so the gradient costs nothing beyond the integral itself. <sqrt(J) B'/B> is also the
      posterior expected E along the model phase (the exact map
      coefficient), and dLLGI/dec = (2 b a/S)(<E> - a ec).

      The integrand is unimodal in s = sqrt(J/sig) in every regime (well
      measured, low information, strongly tilted by the Bessel factor at
      high sigmaA and large ec), and smooth (the centric J^(-1/2)
      singularity is absorbed by the substitution), so a single
      Gauss-Legendre rule over the interval where the log-integrand is
      within `drop` of its maximum is used throughout. The interval is
      placed around the mode of the actual integrand, so the rule follows
      the tilt.
  */

  //! 32-point Gauss-Legendre nodes and weights on [-1,1], computed once.
  struct gauss_legendre_32
  {
    af::shared<double> x, w;
    gauss_legendre_32()
    {
      scitbx::math::quadrature::gauss_legendre_engine<double> engine(32);
      x = engine.x();
      w = engine.w();
    }
    static gauss_legendre_32 const& get()
    {
      static const gauss_legendre_32 gl;
      return gl;
    }
  };

  //! ln B(x), x >= 0: ln I0(x) (acentric) or ln cosh(x) (centric).
  inline double ln_b(double x, bool centric)
  {
    if (centric) return x + std::log1p(std::exp(-2. * x)) - std::log(2.);
    if (x < 500.) return std::log(boost::math::cyl_bessel_i(0, x));
    // Asymptotic expansion of I0(x) exp(-x) sqrt(2 pi x); 6 terms are
    // exact to double precision for x >= 500.
    double t = 1. / (8. * x), s = 1., term = 1.;
    for (int k = 1; k <= 6; k++) {
      double c = (2. * k - 1.);
      term *= c * c * t / k;
      s += term;
    }
    return x - 0.5 * std::log(2. * scitbx::constants::pi * x) + std::log(s);
  }

  //! B'(x)/B(x), x >= 0: I1(x)/I0(x) (acentric) or tanh(x) (centric).
  inline double ratio_b(double x, bool centric)
  {
    if (centric) return std::tanh(x);
    if (x < 500.) {
      return boost::math::cyl_bessel_i(1, x) / boost::math::cyl_bessel_i(0, x);
    }
    double t = 1. / (8. * x), s0 = 1., s1 = 1., term0 = 1., term1 = 1.;
    for (int k = 1; k <= 6; k++) {
      double c = (2. * k - 1.);
      term0 *= c * c * t / k;
      term1 *= -(4. - c * c) * t / k;
      s0 += term0;
      s1 += term1;
    }
    return s1 / s0;
  }

  //! Log-integral and tilted-posterior moments of
  //! int_0^inf J^p exp(-(J-eo_sq)^2/(2 sig^2) - lambda J) B(kappa sqrt J) dJ
  struct integral
  {
    double log_z;   // ln of the integral
    double e_j;     // <J>
    double e_sr;    // <sqrt(J) R(kappa sqrt J)>, R = B'/B
    double e_s;     // <sqrt(J)>, the posterior mean amplitude
    double s_mode, s_lo, s_hi; // quadrature interval in s = sqrt(J/sig)

    integral(
      double eo_sq, double sig, double lambda, double kappa, bool centric,
      double drop = 36.)
    {
      CCTBX_ASSERT(sig > 0);
      const double p = centric ? -0.5 : 0.;
      const double q = 2. * p + 1.;   // power of s in the s-integrand
      const double mu = (eo_sq - lambda * sig * sig) / sig;
      const double k = kappa * std::sqrt(sig);
      const double a0 = centric ? 0.5 : 0.25;
      // psi(s) = -(s^2-mu)^2/2 + ln B(k s) + q ln s, s = sqrt(J/sig),
      // evaluated without its constant -mu^2/2 (which is added back
      // analytically to log_z below): for large sig, |mu| is large and
      // forming (s^2-mu)^2 would leave roundoff of order mu^2 * 1e-16 in
      // every node's weight.
      struct psi_t {
        double mu, k, q; bool centric;
        double operator()(double s) const {
          double s2 = s * s;
          double v = s2 * (mu - 0.5 * s2) + ln_b(k * s, centric);
          if (q != 0) v += q * std::log(s);
          return v;
        }
        double d(double s) const {
          double v = -2. * s * (s*s - mu) + k * ratio_b(k * s, centric);
          if (q != 0) v += q / s;
          return v;
        }
      } psi = {mu, k, q, centric};
      // Mode of psi. psi' crosses zero at most once on (0, inf); for the
      // centric case psi'(0) = 0 and the mode is at 0 unless
      // psi''(0) = 2 mu + k^2 > 0.
      s_mode = 0;
      bool interior = (q != 0) || (mu + a0 * k * k > 0);
      if (interior) {
        double hi = std::sqrt(std::max(mu, 0.) + a0 * k * k + 2.) + 1.;
        while (psi.d(hi) > 0) hi *= 2.;
        double lo = (q != 0) ? hi * 1e-12 : 0.;
        if (q != 0) { while (psi.d(lo) < 0) lo *= 1e-3; }
        else {
          // psi'(s) ~ s (2 mu + k^2) near 0: start just above 0
          lo = std::min(1e-8, hi * 1e-8);
        }
        // Illinois (regula falsi) on psi'
        double flo = psi.d(lo), fhi = psi.d(hi);
        int side = 0;
        for (int it = 0; it < 100; it++) {
          double s = (lo * fhi - hi * flo) / (fhi - flo);
          if (!(s > lo && s < hi)) s = 0.5 * (lo + hi);
          double fs = psi.d(s);
          if (fs > 0) {
            lo = s; flo = fs;
            if (side == 1) fhi *= 0.5;
            side = 1;
          }
          else {
            hi = s; fhi = fs;
            if (side == -1) flo *= 0.5;
            side = -1;
          }
          if (hi - lo < 1e-9 * hi) break;
        }
        s_mode = 0.5 * (lo + hi);
      }
      double psi_max = interior ? psi(s_mode) : psi(0.);
      double target = psi_max - drop;
      // Curvature at the mode sets the first guess of the interval.
      double h = 1e-6 * std::max(s_mode, 1e-3);
      double curv = -(psi.d(s_mode + h) - psi.d(std::max(s_mode - h, 0.)))
                    / (s_mode + h - std::max(s_mode - h, 0.));
      double delta = (curv > 1e-12) ? std::sqrt(2. * drop / curv)
                                    : std::pow(2. * drop, 0.25);
      // Upper end: bracket, then a few bisections.
      s_hi = s_mode + delta;
      while (psi(s_hi) > target) s_hi = s_mode + 2. * (s_hi - s_mode);
      {
        double lo = s_mode, hi = s_hi;
        for (int it = 0; it < 12; it++) {
          double m = 0.5 * (lo + hi);
          if (psi(m) > target) lo = m; else hi = m;
        }
        s_hi = hi;
      }
      // Lower end.
      s_lo = 0.;
      if (interior && s_mode > 0) {
        double cand = s_mode - delta;
        while (cand > 0 && psi(cand) > target) {
          cand = s_mode - 2. * (s_mode - cand);
        }
        // psi(0) is -inf for the acentric case (q = 1)
        cand = std::max(cand, 0.);
        if (cand > 0 || q != 0 || psi(0.) <= target) {
          double lo = cand, hi = s_mode;
          for (int it = 0; it < 12; it++) {
            double m = 0.5 * (lo + hi);
            if (psi(m) > target) hi = m; else lo = m;
          }
          s_lo = lo;
        }
      }
      // Gauss-Legendre on [s_lo, s_hi]; J = sig s^2,
      // dJ J^p = sig^(1+p) 2 s^q ds.
      gauss_legendre_32 const& gl = gauss_legendre_32::get();
      std::size_t n = gl.x.size();
      double half = 0.5 * (s_hi - s_lo);
      std::vector<double> lw(n), s(n);
      double lw_max = -1e300;
      for (std::size_t i = 0; i < n; i++) {
        s[i] = s_lo + half * (gl.x[i] + 1.);
        lw[i] = std::log(gl.w[i] * half * 2.) + psi(s[i]);
        lw_max = std::max(lw_max, lw[i]);
      }
      double z = 0, sj = 0, ssr = 0, ss = 0;
      for (std::size_t i = 0; i < n; i++) {
        double wi = std::exp(lw[i] - lw_max);
        double jj = sig * s[i] * s[i];
        double rj = std::sqrt(jj);
        double x = kappa * rj;
        double r = ratio_b(x, centric);
        z += wi;
        sj += wi * jj;
        ssr += wi * rj * r;
        ss += wi * rj;
      }
      // -mu^2/2 - lambda eo_sq + lambda^2 sig^2/2 = -eo_sq^2/(2 sig^2)
      log_z = lw_max + std::log(z)
            - 0.5 * (eo_sq / sig) * (eo_sq / sig)
            + (1. + p) * std::log(sig);
      e_j = sj / z;
      e_sr = ssr / z;
      e_s = ss / z;
    }
  };

  struct result
  {
    double ll;          // log-likelihood gain
    double d_ll_d_ec;   // d LLGI / d ec
    double d_ll_d_a;    // d LLGI / d a
    double e_expected;  // posterior <E> along the model phase
    double e_abs_expected; // posterior <|E|>
  };

  //! Exact LLGI and derivatives for one reflection. eo_sq, sig: observed
  //! E^2 and its standard deviation (sig > 0); ec >= 0; 0 <= a < 1.
  //! null_log_z, if supplied (not NaN), is M(0, b) for this reflection,
  //! which does not depend on ec or a and can be computed once.
  inline
  result
  evaluate(
    double eo_sq, double sig, double ec, double a, bool centric,
    double null_log_z = std::numeric_limits<double>::quiet_NaN())
  {
    CCTBX_ASSERT(sig > 0);
    CCTBX_ASSERT(a >= 0 && a < 1);
    CCTBX_ASSERT(ec >= 0);
    const double b = centric ? 0.5 : 1.;
    const double s = 1. - a * a;
    const double lambda = b / s;
    const double kappa = 2. * b * a * ec / s;
    integral m(eo_sq, sig, lambda, kappa, centric);
    if (null_log_z != null_log_z) {
      null_log_z = integral(eo_sq, sig, b, 0., centric).log_z;
    }
    const double m_k = m.e_sr;
    const double m_l = -m.e_j;
    const double s2 = s * s, a2 = a * a;
    const double p_val = b * (-std::log(s) - a2 * ec * ec / s);
    const double p_e = -2. * b * a2 * ec / s;
    const double p_a = 2. * b * a / s - 2. * b * a * ec * ec / s2;
    const double k_e = 2. * b * a / s;
    const double k_a = 2. * b * ec * (1. + a2) / s2;
    const double l_a = 2. * b * a / s2;
    result r;
    r.ll = p_val + m.log_z - null_log_z;
    r.d_ll_d_ec = p_e + m_k * k_e;
    r.d_ll_d_a = p_a + m_k * k_a + m_l * l_a;
    r.e_expected = m_k;
    r.e_abs_expected = m.e_s;
    return r;
  }

  //! M(0, b): the ec- and a-independent null-hypothesis term of evaluate().
  inline
  double
  null_log_z(double eo_sq, double sig, bool centric)
  {
    return integral(eo_sq, sig, centric ? 0.5 : 1., 0., centric).log_z;
  }

  //! evaluate() over arrays (for Python-side fits and tests). null_log_z
  //! may be empty (computed here) or one value per reflection.
  class evaluate_many
  {
    public:
      af::shared<double> ll, d_ll_d_ec, d_ll_d_a, e_expected;
      af::shared<double> e_abs_expected;

      evaluate_many(
        af::const_ref<double> const& eo_sq,
        af::const_ref<double> const& sig,
        af::const_ref<double> const& ec,
        af::const_ref<double> const& a,
        af::const_ref<bool> const& centric,
        af::const_ref<double> const& null_log_z_ = af::const_ref<double>(0,0))
      {
        std::size_t n = eo_sq.size();
        CCTBX_ASSERT(sig.size() == n && ec.size() == n);
        CCTBX_ASSERT(a.size() == n && centric.size() == n);
        CCTBX_ASSERT(null_log_z_.size() == 0 || null_log_z_.size() == n);
        ll.reserve(n); d_ll_d_ec.reserve(n); d_ll_d_a.reserve(n);
        e_expected.reserve(n);
        for (std::size_t i = 0; i < n; i++) {
          result r = evaluate(eo_sq[i], sig[i], ec[i], a[i], centric[i],
            null_log_z_.size() ? null_log_z_[i]
                               : std::numeric_limits<double>::quiet_NaN());
          ll.push_back(r.ll);
          d_ll_d_ec.push_back(r.d_ll_d_ec);
          d_ll_d_a.push_back(r.d_ll_d_a);
          e_expected.push_back(r.e_expected);
          e_abs_expected.push_back(r.e_abs_expected);
        }
      }
  };

  //! Per-reflection data for the hybrid LLGI. The Rice approximation with
  //! moment-matched (Dobs, Eeff) has variance 1 - D^2 a^2 about D a ec, of
  //! which 1 - D^2 comes from measurement error and D^2 (1 - a^2) from model
  //! error. It is exact as the measurement error vanishes, and its worst
  //! case (an ec far from what the observation implies) grows as the model
  //! error shrinks, so the exact likelihood is used where
  //!   1 - D^2 > rice_kappa (1 - a^2)^2
  //! i.e. where the measurement variance is not small compared with the
  //! square of the model variance, or where force_exact (no Rice solution);
  //! Rice everywhere else. Reflections with
  //! sig_e_obs_sq <= 0 (no intensity error estimate) always use Rice.
  //! rice_kappa <= 0 means exact wherever possible, for fits of sigmaA
  //! itself (a sigmaA-dependent switch would make the fitted objective
  //! discontinuous).
  class hybrid
  {
    public:
      af::shared<double> e_obs_sq, sig_e_obs_sq, null_log_z, dsqr;
      af::shared<bool> force_exact;
      double rice_kappa;

      hybrid() : rice_kappa(0) {}

      hybrid(
        af::const_ref<double> const& e_obs_sq_,
        af::const_ref<double> const& sig_e_obs_sq_,
        af::const_ref<bool> const& force_exact_,
        af::const_ref<bool> const& centric,
        double rice_kappa_)
      :
        e_obs_sq(e_obs_sq_.begin(), e_obs_sq_.end()),
        sig_e_obs_sq(sig_e_obs_sq_.begin(), sig_e_obs_sq_.end()),
        null_log_z(e_obs_sq_.size(), 0.),
        dsqr(e_obs_sq_.size(), 0.),
        force_exact(force_exact_.begin(), force_exact_.end()),
        rice_kappa(rice_kappa_)
      {
        std::size_t n = e_obs_sq.size();
        CCTBX_ASSERT(sig_e_obs_sq.size() == n);
        CCTBX_ASSERT(force_exact.size() == n && centric.size() == n);
        for (std::size_t i = 0; i < n; i++) {
          if (sig_e_obs_sq[i] > 0) {
            null_log_z[i] = llgi_exact::null_log_z(
              e_obs_sq[i], sig_e_obs_sq[i], centric[i]);
            cctbx::rice_moments rm = cctbx::rice_from_intensity(
              e_obs_sq[i], sig_e_obs_sq[i], centric[i]);
            if (rm.valid) dsqr[i] = rm.dsqr;
          }
        }
      }

      std::size_t size() const { return e_obs_sq.size(); }

      //! Subset (e.g. after fmodel.select()), without recomputing.
      hybrid
      select(af::const_ref<bool> const& selection) const
      {
        CCTBX_ASSERT(selection.size() == size());
        hybrid result;
        result.rice_kappa = rice_kappa;
        for (std::size_t i = 0; i < size(); i++) {
          if (!selection[i]) continue;
          result.e_obs_sq.push_back(e_obs_sq[i]);
          result.sig_e_obs_sq.push_back(sig_e_obs_sq[i]);
          result.null_log_z.push_back(null_log_z[i]);
          result.dsqr.push_back(dsqr[i]);
          result.force_exact.push_back(force_exact[i]);
        }
        return result;
      }

      //! Copy with a different rice_kappa (0: exact wherever there is an
      //! intensity error estimate, as needed by the sigmaA and ScatFrac fits).
      hybrid
      with_rice_kappa(double t) const
      {
        hybrid result(*this);
        result.rice_kappa = t;
        return result;
      }

      bool
      use_exact(std::size_t i, double a) const
      {
        if (!(sig_e_obs_sq[i] > 0 && a > 0 && a < 0.999)) return false;
        if (force_exact[i] || rice_kappa <= 0) return true;
        double c = 1. - a * a;
        return 1. - dsqr[i] > rice_kappa * c * c;
      }

      result
      evaluate_at(std::size_t i, double ec, double a, bool centric) const
      {
        return evaluate(e_obs_sq[i], sig_e_obs_sq[i], std::max(ec, 0.), a,
          centric, null_log_z[i]);
      }

      //! Which reflections would use the exact likelihood at these sigmaA.
      af::shared<bool>
      exact_selection(af::const_ref<double> const& sigmaa) const
      {
        CCTBX_ASSERT(sigmaa.size() == size());
        af::shared<bool> result(size(), false);
        for (std::size_t i = 0; i < size(); i++) {
          result[i] = use_exact(i, sigmaa[i]);
        }
        return result;
      }
  };

  //! cctbx::rice_from_intensity over arrays.
  class rice_moments_many
  {
    public:
      af::shared<bool> valid;
      af::shared<double> mu2, mu4, dsqr, eeff;

      rice_moments_many(
        af::const_ref<double> const& eo_sq,
        af::const_ref<double> const& sig,
        af::const_ref<bool> const& centric)
      {
        std::size_t n = eo_sq.size();
        CCTBX_ASSERT(sig.size() == n && centric.size() == n);
        for (std::size_t i = 0; i < n; i++) {
          cctbx::rice_moments r = cctbx::rice_from_intensity(
            eo_sq[i], sig[i], centric[i]);
          valid.push_back(r.valid);
          mu2.push_back(r.mu2);
          mu4.push_back(r.mu4);
          dsqr.push_back(r.dsqr);
          eeff.push_back(r.eeff);
        }
      }
  };

}}}} // namespace cctbx::xray::targets::llgi_exact

#endif // CCTBX_XRAY_TARGETS_LLGI_EXACT_H
