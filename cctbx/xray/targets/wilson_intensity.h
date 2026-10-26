#ifndef CCTBX_XRAY_TARGETS_WILSON_INTENSITY_H
#define CCTBX_XRAY_TARGETS_WILSON_INTENSITY_H

#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/xray/targets/llgi_exact.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/constants.h>
#include <scitbx/math/erf.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace cctbx { namespace xray { namespace targets { namespace wilson_intensity {

  /*! Wilson intensity likelihood with Gaussian measurement error.

      For an observed intensity I with standard deviation sig_i, and
      Sigma = epsilon * SigmaN * TEPS > 0, the true intensity J has the
      Wilson density
        acentric: p(J) = exp(-J/Sigma)/Sigma
        centric:  p(J) = (2 pi Sigma J)^(-1/2) exp(-J/(2 Sigma))
      and
        p(I) = int_0^inf p(J) N(I; J, sig_i^2) dJ.
      The acentric case follows phasertng pIntensity::pIobsAcen_grad_hess_by_dSN
      (erfc for moderate, erfcx with an asymptotic series for large erfc
      argument), recast in logs and in ln Sigma. The centric case uses the
      Gauss-Legendre integral of llgi_exact::integral and its posterior
      moments of J.
      The tail probabilities port phasertng pIntensity (E^2 scale).
  */

  namespace detail {
    const double sqrt2 = 1.4142135623730951;
    const double sqrt_pi = 1.7724538509055159;
    const double sqrt_2pi = 2.5066282746310002;
    inline double sqr(double x) { return x * x; }
    inline double log_half_max()
    {
      return std::log(std::numeric_limits<double>::max() / 2.);
    }
  }

  struct result
  {
    double log_likelihood;
    double d_ll_d_ln_sigma;
    double d2_ll_d_ln_sigma2;
  };

  //! Acentric ln p(I), d/d ln Sigma, d2/d ln Sigma2 (sigma_i > 0, Sigma > 0).
  inline
  result
  evaluate_acentric(double iobs, double sigma_i, double sigma_n)
  {
    using detail::sqr;
    using detail::sqrt2;
    using detail::sqrt_pi;
    using detail::sqrt_2pi;
    CCTBX_ASSERT(sigma_i > 0);
    CCTBX_ASSERT(sigma_n > 0);
    const double crossover1 = 6.;
    const double sn = sigma_n, sn2 = sqr(sn);
    const double si = sigma_i, si2 = sqr(si);
    const double exparg = sqr(iobs / si) / 2.;
    // X: argument of erfc
    const double x = (-iobs * sn + si2) / (sqrt2 * si * sn);
    double ln_p, g_over_p, h_over_p; // g, h: derivatives of p wrt Sigma
    if (x < crossover1) {
      // prob = exp(a) erfc(x) / (2 sn); expterm = exp(-exparg)
      const double a = (si2 - 2. * iobs * sn) / (2. * sn2);
      ln_p = a + std::log(scitbx::math::erfc(x)) - std::log(2. * sn);
      // expterm/prob without underflow
      const double r = std::exp(-exparg - ln_p);
      g_over_p = (r * si / sqrt_2pi - (si2 + sn * (sn - iobs))) / (sn2 * sn);
      h_over_p = ((sqr(si2) + si2 * sn * (5. * sn - 2. * iobs)
                   + sn2 * (sqr(iobs) - 4. * iobs * sn + 2. * sn2))
                  - r * si * (si2 + sn * (4. * sn - iobs)) / sqrt_2pi)
                 / (sn2 * sn2 * sn2);
    }
    else {
      // prob = expterm erfcx(x) / (2 sn)
      const double erfcx = scitbx::math::erfcx(x);
      ln_p = -exparg + std::log(erfcx) - std::log(2. * sn);
      // Asymptotic series in 1/x (see phasertng), terms to order 4
      const double sndx = -si / (sqrt2 * sn);
      const double ix1 = 1. / x;
      const double ix2 = ix1 / x;
      const double ix3 = ix2 / x;
      const double ix4 = ix3 / x;
      const double ix5 = ix4 / x;
      const double ix6 = ix5 / x;
      const double ix7 = ix6 / x;
      const double ix8 = ix7 / x;
      const double ix9 = ix8 / x;
      const double ix10 = ix9 / x;
      const double ix11 = ix10 / x;
      const double gs =
            -(ix1 +      sndx * ix2) / 2.
          +  (ix3 +   3. * sndx * ix4) / 4.
          -  (3. * ix5 +  15. * sndx * ix6) / 8.
          + (15. * ix7 + 105. * sndx * ix8) / 16.
          - (105. * ix9 + 945. * sndx * ix10) / 32.;
      const double sndxsqr = sqr(sndx);
      const double sn2dx2 = sqrt2 * si / sn;
      const double twominus1 = 2. * sndx - sn2dx2;
      const double hs =
          32. * ix1 + 16. * twominus1 * ix2
        + (-16. + 32. * sndxsqr) * ix3 - 24. * twominus1 * ix4
        + (24. - 96. * sndxsqr) * ix5 + 60. * twominus1 * ix6
        + (-60. + 360. * sndxsqr) * ix7 - 210. * twominus1 * ix8
        + (210. - 1680. * sndxsqr) * ix9 + 945. * twominus1 * ix10
        + (9450. * sndxsqr) * ix11;
      // g = expterm gs / (sqrt(pi) sn^2), h = expterm hs / (32 sqrt(pi) sn^3)
      // and p = expterm erfcx / (2 sn)
      g_over_p = 2. * gs / (sqrt_pi * sn * erfcx);
      h_over_p = hs / (16. * sqrt_pi * sn2 * erfcx);
    }
    result r;
    r.log_likelihood = ln_p;
    r.d_ll_d_ln_sigma = sn * g_over_p;
    r.d2_ll_d_ln_sigma2 = sn * g_over_p + sn2 * (h_over_p - sqr(g_over_p));
    return r;
  }

  //! Centric ln p(I) and derivatives in ln Sigma, via llgi_exact::integral
  //! on the E^2 scale (e = I/Sigma, s = sigma_i/Sigma).
  inline
  result
  evaluate_centric(double iobs, double sigma_i, double sigma_n)
  {
    CCTBX_ASSERT(sigma_i > 0);
    CCTBX_ASSERT(sigma_n > 0);
    const double e = iobs / sigma_n, s = sigma_i / sigma_n;
    llgi_exact::integral m(e, s, 0.5, 0., true);
    const double var = m.e_j2 - m.e_j * m.e_j;
    result r;
    r.log_likelihood = m.log_z - std::log(sigma_n)
      - 0.5 * std::log(scitbx::constants::two_pi)
      - std::log(detail::sqrt_2pi * s);
    r.d_ll_d_ln_sigma = 0.5 * (m.e_j - 1.);
    r.d2_ll_d_ln_sigma2 = 0.25 * var - 0.5 * m.e_j;
    return r;
  }

  //! Arrays in, arrays out. sig_i_obs > 0 is asserted.
  class evaluate_many
  {
    public:
      af::shared<double> log_likelihood, d_ll_d_ln_sigma, d2_ll_d_ln_sigma2;

      evaluate_many(
        af::const_ref<double> const& i_obs,
        af::const_ref<double> const& sig_i_obs,
        af::const_ref<double> const& sigma,
        af::const_ref<bool> const& centric_flags)
      {
        std::size_t n = i_obs.size();
        CCTBX_ASSERT(sig_i_obs.size() == n && sigma.size() == n);
        CCTBX_ASSERT(centric_flags.size() == n);
        log_likelihood.reserve(n);
        d_ll_d_ln_sigma.reserve(n);
        d2_ll_d_ln_sigma2.reserve(n);
        for (std::size_t i = 0; i < n; i++) {
          result r = centric_flags[i]
            ? evaluate_centric(i_obs[i], sig_i_obs[i], sigma[i])
            : evaluate_acentric(i_obs[i], sig_i_obs[i], sigma[i]);
          log_likelihood.push_back(r.log_likelihood);
          d_ll_d_ln_sigma.push_back(r.d_ll_d_ln_sigma);
          d2_ll_d_ln_sigma2.push_back(r.d2_ll_d_ln_sigma2);
        }
      }
  };

  // ---- tail probabilities (E^2 scale, ports of phasertng pIntensity) ----

  //! P(E^2 observed < eosq), eosq <= 0 (sigesq > 0).
  inline
  double
  p_smaller_eosq_acen(double eosq, double sigesq)
  {
    using detail::sqr;
    if (eosq > 0.) return 1.;
    CCTBX_ASSERT(sigesq > 0.);
    const double crossover1 = 6.;
    const double crossover2 = detail::log_half_max();
    double prob = 0.;
    double erfcarg1 = (-eosq + sqr(sigesq)) / (detail::sqrt2 * sigesq);
    double erfcarg2 = -eosq / (detail::sqrt2 * sigesq);
    if (erfcarg1 < crossover1) {
      double exparg = -eosq + sqr(sigesq) / 2.;
      if (exparg > -crossover2)
        prob -= std::exp(exparg) * scitbx::math::erfc(erfcarg1) / 2.;
    }
    else {
      double exparg = -sqr(eosq / sigesq) / 2.;
      if ((erfcarg1 > -crossover1) && (exparg > -crossover2))
        prob -= std::exp(exparg) * scitbx::math::erfcx(erfcarg1) / 2.;
    }
    if (erfcarg2 < crossover1)
      prob += scitbx::math::erfc(erfcarg2) / 2.;
    if (prob < 0.) prob = 0.;
    return prob;
  }

  //! P(E^2 observed > eosq), eosq >= 0 (sigesq >= 0).
  inline
  double
  p_larger_eosq_acen(double eosq, double sigesq)
  {
    using detail::sqr;
    CCTBX_ASSERT(sigesq >= 0.);
    const double crossover1 = 6.;
    const double crossover2 = detail::log_half_max();
    double prob = 0.;
    if (eosq < 0.) return 1.;
    if (sigesq == 0.) {
      if (eosq < crossover2) prob = std::exp(-eosq);
    }
    else {
      double erfcarg1 = (-eosq + sqr(sigesq)) / (detail::sqrt2 * sigesq);
      double erfcarg2 = eosq / (detail::sqrt2 * sigesq);
      if (erfcarg1 < crossover1) {
        double exparg = -eosq + sqr(sigesq) / 2.;
        if (exparg > -crossover2)
          prob += std::exp(exparg) * scitbx::math::erfc(erfcarg1) / 2.;
      }
      else {
        double exparg = -sqr(eosq / sigesq) / 2.;
        if ((erfcarg1 > -crossover1) && (exparg > -crossover2))
          prob += std::exp(exparg) * scitbx::math::erfcx(erfcarg1) / 2.;
      }
      if (erfcarg2 < crossover1)
        prob += scitbx::math::erfc(erfcarg2) / 2.;
    }
    return prob;
  }

  //! Centric density of the observed E^2 (Wilson convolved with a Gaussian
  //! of standard deviation sigesq), by the Gauss-Legendre integral.
  inline
  double
  p_eosq_cen(double eosq, double sigesq)
  {
    CCTBX_ASSERT(sigesq >= 0.);
    if (sigesq > 0.) {
      llgi_exact::integral m(eosq, sigesq, 0.5, 0., true);
      return std::exp(m.log_z - 0.5 * std::log(scitbx::constants::two_pi)
                      - std::log(detail::sqrt_2pi * sigesq));
    }
    // Zero sigma: Wilson density (as phasertng, with the same guard)
    const double cuberteps =
      std::pow(std::numeric_limits<double>::epsilon(), 1. / 3.);
    double eosq_pos = std::max(eosq, cuberteps);
    double eo = std::sqrt(eosq_pos);
    if (eo > 1.e-6)
      return std::exp(-eosq_pos / 2.) / (eo * detail::sqrt_2pi);
    return (1. / eo - eo / 2.) / detail::sqrt_2pi;
  }

  //! P(E^2 observed < eosq), eosq <= 0 (sigesq > 0): 50-point Simpson rule.
  inline
  double
  p_smaller_eosq_cen(double eosq, double sigesq)
  {
    if (eosq > 0.) return 1.;
    CCTBX_ASSERT(sigesq > 0.);
    const int npts = 51; // odd (phasertng: 50), as Simpson needs
    double lower = std::min(eosq - 3 * sigesq, -6 * sigesq);
    double de2 = (eosq - lower) / (npts - 1);
    double prob = p_eosq_cen(lower, sigesq) + p_eosq_cen(eosq, sigesq);
    for (int i = 1; i < npts - 1; i++) {
      double wt = (i % 2) ? 4 : 2; // Simpson: odd 4, even 2
      prob += wt * p_eosq_cen(lower + i * de2, sigesq);
    }
    prob *= de2 / 3;
    return prob;
  }

  //! P(E^2 observed > eosq), eosq >= 0 (sigesq >= 0).
  inline
  double
  p_larger_eosq_cen(double eosq, double sigesq)
  {
    using detail::sqr;
    CCTBX_ASSERT(sigesq >= 0.);
    if (eosq < 0.) return 1.;
    double prob = 0., erfcarg;
    const double maxerfcarg = 9.;
    if (sigesq > 22.) { // Large measurement error Gaussian approximation
      double z = (eosq - 1.) / std::sqrt(2. + sqr(sigesq));
      erfcarg = z / std::sqrt(2.);
      if (erfcarg < maxerfcarg) prob = scitbx::math::erfc(erfcarg) / 2.;
      return prob;
    }
    erfcarg = std::sqrt(eosq / 2.);
    double upper = 4. + 8. * sigesq;
    if (erfcarg > maxerfcarg) // Danger of underflow
      prob = 0.;
    else if (sigesq < 0.00001)
      prob = scitbx::math::erfc(erfcarg);
    else {
      if (eosq > upper) // Good approximation to large outlier integral
        prob = scitbx::math::erfc(erfcarg)
             * p_eosq_cen(eosq, sigesq) / p_eosq_cen(eosq, 0.);
      else {
        const int npts = 201; // odd (phasertng: 200), as Simpson needs
        double de2 = (upper - eosq) / (npts - 1);
        prob = p_eosq_cen(eosq, sigesq) + p_eosq_cen(upper, sigesq);
        for (int i = 1; i < npts - 1; i++) {
          double wt = (i % 2) ? 4 : 2; // Simpson: odd 4, even 2
          prob += wt * p_eosq_cen(eosq + i * de2, sigesq);
        }
        prob *= de2 / 3;
        if (std::sqrt(upper / 2.) <= maxerfcarg)
          prob += scitbx::math::erfc(std::sqrt(upper / 2.))
                * p_eosq_cen(upper, sigesq) / p_eosq_cen(upper, 0.);
      }
    }
    return prob;
  }

  //! Tail probability of each observation under Wilson: the smaller-tail
  //! for I < 0, the larger-tail for I >= 0 (sig_i_obs > 0 asserted).
  class tail_probability_many
  {
    public:
      af::shared<double> probability;

      tail_probability_many(
        af::const_ref<double> const& i_obs,
        af::const_ref<double> const& sig_i_obs,
        af::const_ref<double> const& sigma,
        af::const_ref<bool> const& centric_flags)
      {
        std::size_t n = i_obs.size();
        CCTBX_ASSERT(sig_i_obs.size() == n && sigma.size() == n);
        CCTBX_ASSERT(centric_flags.size() == n);
        probability.reserve(n);
        for (std::size_t i = 0; i < n; i++) {
          CCTBX_ASSERT(sig_i_obs[i] > 0 && sigma[i] > 0);
          double eosq = i_obs[i] / sigma[i];
          double sigesq = sig_i_obs[i] / sigma[i];
          double p;
          if (i_obs[i] < 0.)
            p = centric_flags[i] ? p_smaller_eosq_cen(eosq, sigesq)
                                 : p_smaller_eosq_acen(eosq, sigesq);
          else
            p = centric_flags[i] ? p_larger_eosq_cen(eosq, sigesq)
                                 : p_larger_eosq_acen(eosq, sigesq);
          probability.push_back(p);
        }
      }
  };

}}}} // namespace cctbx::xray::targets::wilson_intensity

#endif // CCTBX_XRAY_TARGETS_WILSON_INTENSITY_H
