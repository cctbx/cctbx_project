#ifndef CCTBX_XRAY_TARGETS_LLGI_EXACT_H
#define CCTBX_XRAY_TARGETS_LLGI_EXACT_H

#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <limits>
#include <scitbx/constants.h>
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
      and similarly for the second derivatives, so the gradient costs
      nothing beyond the integral itself. <sqrt(J) B'/B> is also the
      posterior expected E along the model phase (the exact map
      coefficient), and dLLGI/dec = (2 b a/S)(<E> - a ec).

      The integrand is unimodal in s = sqrt(J/sig) in every regime (well
      measured, low information, strongly tilted by the Bessel factor at
      high sigmaA and large ec), and smooth (the centric J^(-1/2)
      singularity is absorbed by the substitution), so a single
      Gauss-Legendre rule over the interval where the log-integrand is
      within `drop` of its maximum is used throughout. Node placement
      depends on all the arguments, so the rule adapts to the tilt; fixed
      nodes per (reflection, sigmaA), as in the handoff prototype, fail
      badly there.
  */

  //! Gauss-Legendre nodes/weights on [-1,1], computed once.
  class gauss_legendre
  {
    public:
      std::vector<double> x, w;
      explicit gauss_legendre(std::size_t n)
      : x(n), w(n)
      {
        const double pi = scitbx::constants::pi;
        for (std::size_t i = 0; i < (n+1)/2; i++) {
          double z = std::cos(pi * (i + 0.75) / (n + 0.5));
          double pp = 0;
          for (int it = 0; it < 100; it++) {
            double p1 = 1, p2 = 0;
            for (std::size_t j = 1; j <= n; j++) {
              double p3 = p2; p2 = p1;
              p1 = ((2.*j - 1.) * z * p2 - (j - 1.) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.);
            double z1 = z;
            z = z1 - p1 / pp;
            if (std::abs(z - z1) < 1e-16) break;
          }
          x[i] = -z; x[n-1-i] = z;
          w[i] = w[n-1-i] = 2. / ((1. - z * z) * pp * pp);
        }
      }
      static gauss_legendre const& get32()
      {
        static const gauss_legendre gl(32);
        return gl;
      }
      static gauss_legendre const& get96()
      {
        static const gauss_legendre gl(96);
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

  //! (B'(x)/B(x))/x, smooth at x = 0.
  inline double ratio_b_over_x(double x, bool centric)
  {
    if (x < 1e-4) return centric ? 1. - x * x / 3. : 0.5 - x * x / 16.;
    return ratio_b(x, centric) / x;
  }

  //! Log-integral and tilted-posterior moments of
  //! int_0^inf J^p exp(-(J-eo_sq)^2/(2 sig^2) - lambda J) B(kappa sqrt J) dJ
  struct integral
  {
    double log_z;   // ln of the integral
    double e_j;     // <J>
    double e_jj;    // <J^2>
    double e_sr;    // <sqrt(J) R(kappa sqrt J)>, R = B'/B
    double e_j_sr;  // <J sqrt(J) R>
    double e_j_b2;  // <J B''/B>
    double e_s;     // <sqrt(J)>, the posterior mean amplitude
    double s_mode, s_lo, s_hi; // quadrature interval in s = sqrt(J/sig)

    integral(
      double eo_sq, double sig, double lambda, double kappa, bool centric,
      double drop = 36., bool high_precision = false)
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
      gauss_legendre const& gl = high_precision ? gauss_legendre::get96()
                                                : gauss_legendre::get32();
      std::size_t n = gl.x.size();
      double half = 0.5 * (s_hi - s_lo);
      std::vector<double> lw(n), s(n);
      double lw_max = -1e300;
      for (std::size_t i = 0; i < n; i++) {
        s[i] = s_lo + half * (gl.x[i] + 1.);
        lw[i] = std::log(gl.w[i] * half * 2.) + psi(s[i]);
        lw_max = std::max(lw_max, lw[i]);
      }
      double z = 0, sj = 0, sjj = 0, ssr = 0, sjsr = 0, sjb2 = 0, ss = 0;
      for (std::size_t i = 0; i < n; i++) {
        double wi = std::exp(lw[i] - lw_max);
        double jj = sig * s[i] * s[i];
        double rj = std::sqrt(jj);
        double x = kappa * rj;
        double r = ratio_b(x, centric);
        double b2 = centric ? 1. : 1. - ratio_b_over_x(x, centric);
        z += wi;
        sj += wi * jj;
        sjj += wi * jj * jj;
        ssr += wi * rj * r;
        sjsr += wi * jj * rj * r;
        sjb2 += wi * jj * b2;
        ss += wi * rj;
      }
      // -mu^2/2 - lambda eo_sq + lambda^2 sig^2/2 = -eo_sq^2/(2 sig^2)
      log_z = lw_max + std::log(z)
            - 0.5 * (eo_sq / sig) * (eo_sq / sig)
            + (1. + p) * std::log(sig);
      e_j = sj / z;
      e_jj = sjj / z;
      e_sr = ssr / z;
      e_j_sr = sjsr / z;
      e_j_b2 = sjb2 / z;
      e_s = ss / z;
    }
  };

  struct result
  {
    double ll;          // log-likelihood gain
    double d_ll_d_ec;   // d LLGI / d ec
    double d_ll_d_a;    // d LLGI / d a
    double d2_ll_d_a2;  // d^2 LLGI / d a^2
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
    const double m_kk = m.e_j_b2 - m.e_sr * m.e_sr;
    const double m_ll = m.e_jj - m.e_j * m.e_j;
    const double m_kl = -(m.e_j_sr - m.e_j * m.e_sr);
    const double s2 = s * s, s3 = s2 * s, a2 = a * a;
    const double p_val = b * (-std::log(s) - a2 * ec * ec / s);
    const double p_e = -2. * b * a2 * ec / s;
    const double p_a = 2. * b * a / s - 2. * b * a * ec * ec / s2;
    const double p_aa = b * (2. * (1. + a2) / s2
                             - 2. * ec * ec * (1. + 3. * a2) / s3);
    const double k_e = 2. * b * a / s;
    const double k_a = 2. * b * ec * (1. + a2) / s2;
    const double k_aa = 4. * b * ec * a * (3. + a2) / s3;
    const double l_a = 2. * b * a / s2;
    const double l_aa = 2. * b * (1. + 3. * a2) / s3;
    result r;
    r.ll = p_val + m.log_z - null_log_z;
    r.d_ll_d_ec = p_e + m_k * k_e;
    r.d_ll_d_a = p_a + m_k * k_a + m_l * l_a;
    r.d2_ll_d_a2 = p_aa + m_kk * k_a * k_a + 2. * m_kl * k_a * l_a
                 + m_ll * l_a * l_a + m_k * k_aa + m_l * l_aa;
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

  //! 2nd/4th-moment Rice parameters from the French-Wilson posterior, as in
  //! phasertng's math::rice_from_intensity: <E^2> by the same quadrature
  //! (exact for both acentric and centric), <E^4> = m <E^2> + k sig^2 with
  //! k = 1 (acentric) or 1/2 (centric) and m = eo_sq - k sig^2.
  //! For large sig that identity is a difference of two terms of order
  //! sig^2, and near the edge of the Rice family (D -> 0) the solution
  //! depends on a further near-cancellation, so <E^2> is computed with
  //! 96 nodes (to about machine precision) rather than the 32 used for
  //! the likelihood; phasertng gets the same precision from closed forms.
  struct rice_moments
  {
    bool valid;
    double mu2, mu4, dsqr, eeff;

    rice_moments(double eo_sq, double sig, bool centric)
    : valid(false), mu2(0), mu4(0), dsqr(0), eeff(0)
    {
      if (!(sig > 0)) {
        valid = true;
        mu2 = std::max(eo_sq, 0.);
        mu4 = mu2 * mu2;
        dsqr = 1.;
        eeff = std::sqrt(mu2);
        return;
      }
      const double k = centric ? 0.5 : 1.;
      const double m = eo_sq - k * sig * sig;
      mu2 = integral(eo_sq, sig, k, 0., centric, 36., true).e_j;
      mu4 = m * mu2 + k * sig * sig;
      const double eta = mu2 - 1.;
      double gap, disc;
      if (centric) {
        const double zeta = mu4 - 3.;
        gap = eta * eta + 6. * eta - zeta;   // 2(S^2 - eta^2)
        disc = (3. * eta * eta + 6. * eta - zeta) / 2.;
      }
      else {
        const double zeta = mu4 - 2.;
        gap = eta * eta + 4. * eta - zeta;   // S^2 - eta^2
        disc = 2. * eta * eta + 4. * eta - zeta;
      }
      if (!(disc >= 0.)) return;
      const double sq = std::sqrt(disc);
      double d2 = (eta >= 0.) ? (centric ? gap / 2. : gap) / (eta + sq)
                              : sq - eta;
      if (!(d2 > 0.)) return;
      d2 = std::min(d2, 1.);
      valid = true;
      dsqr = d2;
      eeff = std::sqrt(std::max(1. + eta / d2, 0.));
    }
  };

  //! evaluate() over arrays (for Python-side fits and tests). null_log_z
  //! may be empty (computed here) or one value per reflection.
  class evaluate_many
  {
    public:
      af::shared<double> ll, d_ll_d_ec, d_ll_d_a, d2_ll_d_a2, e_expected;
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
        d2_ll_d_a2.reserve(n); e_expected.reserve(n);
        for (std::size_t i = 0; i < n; i++) {
          result r = evaluate(eo_sq[i], sig[i], ec[i], a[i], centric[i],
            null_log_z_.size() ? null_log_z_[i]
                               : std::numeric_limits<double>::quiet_NaN());
          ll.push_back(r.ll);
          d_ll_d_ec.push_back(r.d_ll_d_ec);
          d_ll_d_a.push_back(r.d_ll_d_a);
          d2_ll_d_a2.push_back(r.d2_ll_d_a2);
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
  //! square of the model variance (hybrid LLGI handoff, revision 2, sec. 5.5,
  //! applied at every sigmaA a), or where force_exact (no Rice solution);
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
            rice_moments rm(e_obs_sq[i], sig_e_obs_sq[i], centric[i]);
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

      //! Fraction of the Rice variance due to measurement error at sigmaA a.
      double
      measurement_fraction(std::size_t i, double a) const
      {
        return (1. - dsqr[i]) / (1. - dsqr[i] * a * a);
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

      //! measurement_fraction() for every reflection (0 where there is no
      //! intensity error estimate), for diagnostics.
      af::shared<double>
      measurement_fractions(af::const_ref<double> const& sigmaa) const
      {
        CCTBX_ASSERT(sigmaa.size() == size());
        af::shared<double> result(size(), 0.);
        for (std::size_t i = 0; i < size(); i++) {
          if (sig_e_obs_sq[i] > 0) {
            result[i] = measurement_fraction(i, std::min(sigmaa[i], 0.999));
          }
        }
        return result;
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

  //! French-Wilson posterior of |E| in units of sqrt(sigma_I): for
  //! J/sigma_I = u >= 0 with posterior u^p exp(-(u - h)^2/2) (p = 0
  //! acentric, -1/2 centric), where h = I/sigma_I - sigma_I/<I>
  //! (acentric) or I/sigma_I - sigma_I/(2<I>) (centric), the French-Wilson
  //! F = <sqrt(u)> sqrt(sigma_I) and SIGF = sd(sqrt(u)) sqrt(sigma_I).
  struct french_wilson_moments
  {
    double mean, sd;
    french_wilson_moments(double h, bool centric)
    {
      integral m(h, 1.0, 0.0, 0.0, centric);
      mean = m.e_s;
      sd = std::sqrt(std::max(m.e_j - m.e_s * m.e_s, 0.0));
    }
  };

  //! Recover (I, sigma_I) from French-Wilson amplitudes. SIGF/F depends on
  //! h alone and decreases monotonically from sqrt(4/pi - 1) (acentric) or
  //! sqrt(pi/2 - 1) (centric) as h -> -infinity, to 1/(2h) for large h, so
  //! h follows from SIGF/F by bisection; then sigma_I = (F/<sqrt(u)>)^2
  //! and I = sigma_I*(h + c*sigma_I/mean_intensity), c = 1 (acentric) or
  //! 1/2 (centric). mean_intensity is the prior <I> the French-Wilson
  //! calculation used (it only enters the last, weak-data term). Where
  //! SIGF/F is at or above the h = h_min limit, h is set to h_min and
  //! prior_dominated is true: such amplitudes carry essentially no
  //! information about I.
  class french_wilson_inverse
  {
    public:
      af::shared<double> i_obs, sig_i_obs, h;
      af::shared<bool> valid, prior_dominated;

      french_wilson_inverse(
        af::const_ref<double> const& f,
        af::const_ref<double> const& sigf,
        af::const_ref<double> const& mean_intensity,
        af::const_ref<bool> const& centric,
        double h_min = -6.0)
      {
        std::size_t n = f.size();
        CCTBX_ASSERT(sigf.size() == n && mean_intensity.size() == n);
        CCTBX_ASSERT(centric.size() == n);
        for (std::size_t i = 0; i < n; i++) {
          bool ok = f[i] > 0 && sigf[i] > 0 && mean_intensity[i] > 0;
          double hh = 0, si = 0, io = 0;
          bool prior = false;
          if (ok) {
            double r = sigf[i] / f[i];
            french_wilson_moments m_lo(h_min, centric[i]);
            if (r >= m_lo.sd / m_lo.mean) {
              hh = h_min;
              prior = true;
            }
            else {
              double lo = h_min, hi = std::max(10.0, 2.0 / r);
              while (true) {
                french_wilson_moments m(hi, centric[i]);
                if (m.sd / m.mean < r) break;
                hi *= 2;
              }
              for (int it = 0; it < 80; it++) {
                double mid = 0.5 * (lo + hi);
                french_wilson_moments m(mid, centric[i]);
                if (m.sd / m.mean > r) lo = mid; else hi = mid;
                if (hi - lo < 1e-10 * std::max(1.0, std::abs(hi))) break;
              }
              hh = 0.5 * (lo + hi);
            }
            french_wilson_moments m(hh, centric[i]);
            double sqrt_si = f[i] / m.mean;
            si = sqrt_si * sqrt_si;
            double c = centric[i] ? 0.5 : 1.0;
            io = si * (hh + c * si / mean_intensity[i]);
          }
          i_obs.push_back(io);
          sig_i_obs.push_back(si);
          h.push_back(hh);
          valid.push_back(ok);
          prior_dominated.push_back(prior);
        }
      }
  };

  //! rice_moments over arrays.
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
          rice_moments r(eo_sq[i], sig[i], centric[i]);
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
