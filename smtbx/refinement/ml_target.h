#ifndef SMTBX_REFINEMENT_ML_TARGET_H
#define SMTBX_REFINEMENT_ML_TARGET_H

/// Maximum-likelihood refinement expressed as weighted least squares.
/** The build loop passes its accumulator a computed observable, its gradient,
    an observation and a weight, per reflection. That interface carries no
    reflection index, so per-reflection distribution parameters cannot be
    looked up downstream. Instead an *effective* observation and weight are
    substituted for the measured ones, chosen so that the normal equations
    assembled are those of the likelihood.

    For the amplitude target (Rice; cctbx/xray/targets/mlf.h) the derivative
    separates as
    \f$ \partial T/\partial F_c
        = (\partial T/\partial |F_c|)\, \overline{F_c}/|F_c| \f$,
    the real scalar \f$ \partial T/\partial |F_c| \f$ being the sum of the two
    terms mlf.h names d1 and d2. With \f$ a = \alpha k \f$,
    \f$ b = \beta k^2 \f$ and \f$ X = I_1/I_0 \f$,

    \f[ w_{\mathrm{eff}} = \frac{2a^2}{\epsilon b} \quad\text{(acentric)},
        \qquad \frac{a^2}{\epsilon b} \quad\text{(centric)} \f]
    \f[ y_{\mathrm{eff}}
        = \frac{f_o}{a} X\!\left(\frac{2 a f_o |F_c|}{\epsilon b}\right)
          \quad\text{(acentric)},
        \qquad \frac{f_o}{a}\tanh\!\left(\frac{a f_o |F_c|}{\epsilon b}\right)
          \quad\text{(centric)} \f]

    satisfy \f$ w_{\mathrm{eff}}(y_{\mathrm{eff}} - |F_c|)
                = -\partial T/\partial |F_c| \f$ exactly, so the ordinary
    Gauss-Newton assembly yields the likelihood gradient. This is the reduction
    used by Refmac: the likelihood becomes least squares against an expected
    amplitude re-derived each cycle.

    \f$ w_{\mathrm{eff}} \f$ is the positive part of the curvature only. That
    keeps the normal matrix positive definite, as required by the Cholesky
    decomposition in linear_ls::solve and by the non-negative weight that
    matrix::rank_n_update and store_design_row both assume. The omitted term is
    \f$ (2 a f_o/\epsilon b)^2 X'(u) \f$, of relative size
    \f$ O(1/|F_c|^2) \f$ against the retained one.

    The observable remains \f$ u = |F_c|^2 \f$. Since
    \f$ \partial u/\partial p = 2|F_c|\, \partial |F_c|/\partial p \f$, taking

    \f[ w_u = \frac{w_{\mathrm{eff}}}{4u}, \qquad
        y_u = 2|F_c| y_{\mathrm{eff}} - u \f]

    reproduces the \f$ |F_c| \f$-space normal equations identically,

    \f[ \sum w_u \frac{\partial u}{\partial p}
                 \frac{\partial u}{\partial p}^T
        = \sum w_{\mathrm{eff}} \frac{\partial |F_c|}{\partial p}
                                \frac{\partial |F_c|}{\partial p}^T \f]
    \f[ \sum w_u (y_u - u) \frac{\partial u}{\partial p}
        = \sum w_{\mathrm{eff}} (y_{\mathrm{eff}} - |F_c|)
               \frac{\partial |F_c|}{\partial p} \f]

    leaving the observable, the gradients, the meaning of the scale factor and
    the existing bindings unchanged.

    Twinning, a non-trivial Fc correction and non-unit measured scales are
    refused by the builder; see check_maximum_likelihood() in least_squares.h.

    References:
     - Murshudov, Vagin & Dodson (1997) Acta Cryst. D53, 240-255.
     - Pannu & Read (1996) Acta Cryst. A52, 659-668.
     - Read (1990) Acta Cryst. A46, 900-912.
*/

/* mlf.h is written to be included from cctbx/xray/boost_python/targets.cpp and
   does not include its own prerequisites, hence the two headers before it.
 */
#include <cctbx/xray/targets.h>
#include <cctbx/xray/targets/common_results.h>
#include <cctbx/xray/targets/mlf.h>

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/math/bessel.h>
#include <cmath>
#include <limits>
#include <vector>

namespace smtbx { namespace refinement { namespace least_squares {

  /// Per-reflection distribution parameters for the likelihood targets.
  /** A default-constructed instance is inactive, denoting ordinary least
      squares and leaving that path untouched. Every array is indexed by the
      position of the reflection in the observation list, so a worker needs
      only the index it already holds; all are read-only during a build and
      therefore need no per-thread copy.
   */
  template <typename FloatType>
  struct ml_data {
    /* Owning copies rather than references: instances are constructed from
       Python, where the arrays passed may be temporaries whose lifetime ends
       before the build does.
     */
    scitbx::af::shared<FloatType> alpha, beta, epsilon;
    scitbx::af::shared<bool> centric;
    /// Reflections excluded from the target sum; empty excludes none.
    /** alpha and beta are estimated on this set, so including it in the sum
        as well would make that estimate self-referential, biasing alpha up
        and beta down and hence every effective weight.
     */
    scitbx::af::shared<bool> free_flags;
    /* Activity is a flag rather than a test on alpha.size(): the arrays arrive
       as scitbx::af::const_ref, whose default constructor initialises nothing,
       so size() on a default-constructed one is indeterminate.
     */
    bool active_;
    /// Selects the intensity target over the amplitude one.
    bool intensity_;
    /* The scale the *likelihood* works at, which is not the crystallographic
       one. mlf.h forms a = alpha k, so alpha and k have to come from the same
       estimation; the builder's scale_factor is the Fo/Fc scale that R1 and
       the difference map want, and feeding that here would apply it a second
       time on top of the one alpha already carries.
     */
    FloatType scale_factor_;

    ml_data() : active_(false), intensity_(false), scale_factor_(1) {}

    ml_data(scitbx::af::const_ref<FloatType> const &alpha_,
            scitbx::af::const_ref<FloatType> const &beta_,
            scitbx::af::const_ref<FloatType> const &epsilon_,
            scitbx::af::const_ref<bool> const &centric_,
            bool intensity = false,
            scitbx::af::const_ref<bool> const &free_flags_
              = scitbx::af::const_ref<bool>(0, 0),
            FloatType scale_factor = 1)
      : alpha(alpha_.begin(), alpha_.end()),
        beta(beta_.begin(), beta_.end()),
        epsilon(epsilon_.begin(), epsilon_.end()),
        centric(centric_.begin(), centric_.end()),
        free_flags(free_flags_.begin(), free_flags_.end()),
        active_(true),
        intensity_(intensity),
        scale_factor_(scale_factor > 0 ? scale_factor : FloatType(1))
    {}

    /// Whether this is a maximum-likelihood refinement at all.
    bool active() const { return active_; }
    /// The scale alpha and beta were estimated against.
    FloatType scale_factor() const { return scale_factor_; }
    /// True for the intensity target, which convolves with the measured sigma.
    bool is_intensity() const { return intensity_; }
    /// Whether reflection i_h is held out of the target sum.
    bool is_free(std::size_t i_h) const {
      return free_flags.size() != 0 && free_flags[i_h];
    }

    /// Whether the arrays describe exactly n_reflections reflections.
    /** free_flags may also be empty, denoting that none are held out.
     */
    bool is_consistent_with(std::size_t n_reflections) const {
      if (!active()) return true;
      return alpha.size() == n_reflections
          && beta.size() == n_reflections
          && epsilon.size() == n_reflections
          && centric.size() == n_reflections
          && (free_flags.size() == 0 || free_flags.size() == n_reflections);
    }
  };

  /// The effective observation and weight for one reflection, in F^2 space.
  /** \param fo_sq the measured intensity
      \param observable the computed intensity \f$ |F_c|^2 \f$
      \param yo_eff set to the effective observation \f$ y_u \f$
      \param w_eff set to the effective weight \f$ w_u \f$
      \return false if the reflection cannot contribute, in which case the
              caller gives it zero weight. The guards match those in mlf.h, so
              a reflection refused here is one for which the target would
              return zero anyway.
   */
  /// Whether the F^2-space conversion is safe for this reflection.
  /** Both effective observations divide by the *calculated* intensity, and a
      reflection the model puts at zero makes that enormous. Shared so the two
      targets cannot disagree about which reflections they drop - the amplitude
      one was guarded first and the intensity one was not, which left MLI
      unchanged and looked like the guard having no effect at all.
   */
  template <typename FloatType>
  inline bool ml_intensity_is_usable(FloatType fo_sq, FloatType sigma_fo_sq,
                                     FloatType observable)
  {
    /* Against the larger of the measurement and its uncertainty. A purely
       relative test on fo_sq alone does not bound the weight: where Fo^2 is
       itself tiny the threshold shrinks with it and a |Fc|^2 of 1e-25 still
       clears a bar of 1e-28. That left the worst weights at 1e29 after the
       first attempt, with the 99th percentile already down at 16.7.

       sigma is the scale that does not vanish - a measured reflection has a
       real uncertainty whatever its intensity came out at - so it is what
       makes this an absolute bound in the units of the data.
     */
    FloatType const scale = std::max(fo_sq, sigma_fo_sq);
    return observable > 1e-6*scale;
  }

  template <typename FloatType>
  bool ml_effective_observation(
    FloatType fo_sq,
    FloatType sigma_fo_sq,
    FloatType observable,
    FloatType alpha,
    FloatType beta,
    FloatType scale_factor,
    FloatType epsilon,
    bool centric,
    FloatType &yo_eff,
    FloatType &w_eff)
  {
    yo_eff = observable;
    w_eff = 0;
    if (scale_factor <= 0) scale_factor = 1;
    // mlf.h's own guards, so that the two agree about which reflections count
    if (alpha <= 0 || beta <= 1e-3 || fo_sq <= 0 || observable <= 0) {
      return false;
    }
    if (epsilon <= 0) {
      return false;
    }
    FloatType const a = alpha*scale_factor;
    FloatType const b = beta*scale_factor*scale_factor;
    FloatType const eb = epsilon*b;
    FloatType const fo = std::sqrt(fo_sq);
    FloatType const modulus_fc = std::sqrt(observable);
    FloatType w_amplitude, yo_amplitude;
    if (centric) {
      w_amplitude = a*a/eb;
      yo_amplitude = (fo/a)*std::tanh(a*fo*modulus_fc/eb);
    }
    else {
      w_amplitude = 2*a*a/eb;
      yo_amplitude = (fo/a)*scitbx::math::bessel::i1_over_i0(
        2*a*fo*modulus_fc/eb);
    }
    /* Into F^2 space; see the header comment for why this is exact.

       The division is by the *calculated* intensity, and a reflection the
       model puts at zero makes it enormous: measured on jaca_A14, 193 of 6819
       reflections have |Fc|^2 below 1e-6 and the smallest is 6e-31, which
       turned the largest weight into 6e29 and the objective into 2.5e22
       against least squares' 4e-3.

       Zero is the right weight there, not infinity. The conversion exists
       because du/dp = 2|Fc| d|Fc|/dp, so as |Fc| goes to zero the row of the
       design matrix vanishes too and the reflection stops saying anything
       about the parameters through this route. The product w (x) row (x) row
       stays finite - it is w_amplitude g (x) g - but w and the row are kept
       apart, and sum w yo^2, which normalises the objective, is not protected
       by that cancellation.

       The threshold is relative to what was measured, so it means "the model
       says this reflection is absent" rather than a fixed number of electrons.
       It wants review against a structure with genuine systematic absences.
     */
    if (!ml_intensity_is_usable(fo_sq, sigma_fo_sq, observable)) {
      w_eff = 0;
      yo_eff = observable;
      return false;
    }
    w_eff = w_amplitude/(4*observable);
    yo_eff = 2*modulus_fc*yo_amplitude - observable;
    return true;
  }

  /// Intensity-based maximum likelihood: the Rice model convolved with the
  /// experimental error.
  /** This is not a change of variable. Expressing the Rice distribution in
      \f$ I_o = F_o^2 \f$ rather than \f$ F_o \f$ introduces the Jacobian
      \f$ |dF_o/dI_o| = 1/(2F_o) \f$, which is independent of \f$ F_c \f$ and
      so alters no derivative with respect to the model; a target formed that
      way refines identically to the amplitude one.

      The intensity likelihood proper convolves the model distribution with the
      experimental error, which allows a negative measured intensity to be used
      rather than rejected and requires no French-Wilson conversion:

      \f[ P(I_o) = \int_0^{\infty}
                   P_{\mathrm{Rice}}(I_t \mid F_c)\,
                   N(I_o; I_t, \sigma^2)\, dI_t \f]

      with, for acentric reflections,

      \f[ P_{\mathrm{Rice}}(I_t \mid F_c)
          = \frac{1}{\epsilon b}
            \exp\!\left(-\frac{I_t + (a F_c)^2}{\epsilon b}\right)
            I_0\!\left(\frac{2 a \sqrt{I_t} F_c}{\epsilon b}\right) \f]

      It departs from the amplitude target by an amount which depends on
      \f$ F_c \f$, grows with \f$ \sigma^2 \f$ and vanishes as
      \f$ \sigma \to 0 \f$.

      The integral is evaluated by Gauss-Hermite quadrature after the
      substitution \f$ I_t = I_o + \sqrt{2}\,\sigma x \f$, which casts it in the
      weight \f$ e^{-x^2} \f$ the rule assumes. Nodes at \f$ I_t < 0 \f$
      contribute nothing.

      Reference: Pannu & Read (1996) Acta Cryst. A52, 659-668.
   */
  namespace mli_quadrature {
    /* Twenty nodes: the integrand is smooth and unimodal in one variable, so
       this resolves it well beyond the precision carried elsewhere, at a cost
       of twenty Bessel evaluations per reflection.
     */
    inline int n_nodes() { return 20; }

    inline double const* nodes() {
      static const double x[20] = {
        -5.387480890011233, -4.603682449550744, -3.944764040115625,
        -3.347854567383216, -2.788806058428131, -2.254974002089276,
        -1.738537712116586, -1.234076215395323, -0.737473728545394,
        -0.245340708300901,  0.245340708300901,  0.737473728545394,
         1.234076215395323,  1.738537712116586,  2.254974002089276,
         2.788806058428131,  3.347854567383216,  3.944764040115625,
         4.603682449550744,  5.387480890011233 };
      return x;
    }

    /// Weights of the rule, already divided by sqrt(pi) so they sum to one.
    inline double const* weights() {
      static const double w[20] = {
        2.229393645534151e-13, 4.399340992273181e-10, 1.086069370769282e-07,
        7.802556478532064e-06, 2.283386360163540e-04, 3.243773342237863e-03,
        2.481052088746362e-02, 1.090172060200233e-01, 2.866755053628341e-01,
        4.622436696006101e-01, 4.622436696006101e-01, 2.866755053628341e-01,
        1.090172060200233e-01, 2.481052088746362e-02, 3.243773342237863e-03,
        2.283386360163540e-04, 7.802556478532064e-06, 1.086069370769282e-07,
        4.399340992273181e-10, 2.229393645534151e-13 };
      return w;
    }
  }

  /// dT/d|Fc| for the intensity target, and whether it could be computed.
  /** The gradient is exact for the convolved target; the curvature used by the
      caller is the amplitude one, which is a positive approximation. That is
      ordinary Gauss-Newton practice: an approximate Hessian changes how fast a
      refinement converges, never what it converges to, because the gradient is
      what vanishes at the answer.
   */
  template <typename FloatType>
  bool mli_d_target_d_modulus_fc(
    FloatType fo_sq,
    FloatType sigma_fo_sq,
    FloatType observable,
    FloatType alpha,
    FloatType beta,
    FloatType scale_factor,
    FloatType epsilon,
    bool centric,
    FloatType &d_target)
  {
    d_target = 0;
    if (scale_factor <= 0) scale_factor = 1;
    if (alpha <= 0 || beta <= 1e-3 || observable <= 0 || sigma_fo_sq <= 0) {
      return false;
    }
    if (epsilon <= 0) {
      return false;
    }
    FloatType const a = alpha*scale_factor;
    FloatType const b = beta*scale_factor*scale_factor;
    FloatType const eb = epsilon*b;
    FloatType const modulus_fc = std::sqrt(observable);
    FloatType const root_two_sigma =
      std::sqrt(FloatType(2))*sigma_fo_sq;

    int const n = mli_quadrature::n_nodes();
    double const *x = mli_quadrature::nodes();
    double const *w = mli_quadrature::weights();

    /* Accumulated with the largest log-density factored out: the Rice density
       underflows to zero for a reflection the model fits badly, and the ratio
       below is then 0/0 rather than the finite number it should be.
     */
    FloatType ln_max = -std::numeric_limits<FloatType>::max();
    std::vector<FloatType> ln_p(n), dln_p(n);
    std::vector<bool> usable(n, false);
    for (int i = 0; i < n; i++) {
      FloatType const it = fo_sq + root_two_sigma*static_cast<FloatType>(x[i]);
      if (it <= 0) {
        continue;             // no density below zero true intensity
      }
      FloatType const root_it = std::sqrt(it);
      FloatType const arg = 2*a*root_it*modulus_fc/eb;
      ln_p[i] = -std::log(eb) - (it + a*a*observable)/eb
                + scitbx::math::bessel::ln_of_i0(arg);
      // d ln P_Rice / d|Fc|
      dln_p[i] = -2*a*a*modulus_fc/eb
                 + (2*a*root_it/eb)*scitbx::math::bessel::i1_over_i0(arg);
      usable[i] = true;
      if (ln_p[i] > ln_max) {
        ln_max = ln_p[i];
      }
    }
    if (ln_max == -std::numeric_limits<FloatType>::max()) {
      return false;
    }
    FloatType sum = 0, sum_d = 0;
    for (int i = 0; i < n; i++) {
      if (!usable[i]) {
        continue;
      }
      FloatType const p = static_cast<FloatType>(w[i])*std::exp(ln_p[i] - ln_max);
      sum += p;
      sum_d += p*dln_p[i];
    }
    if (sum <= 0) {
      return false;
    }
    // T = -ln sum, so dT/d|Fc| = -(sum of P dlnP)/(sum of P)
    d_target = -sum_d/sum;
    return true;
  }

  /// The effective observation and weight for the intensity target, in F^2.
  /** Same two-step shape as the amplitude case: the exact gradient decides the
      observation, the amplitude curvature supplies a positive weight.
   */
  template <typename FloatType>
  bool mli_effective_observation(
    FloatType fo_sq,
    FloatType sigma_fo_sq,
    FloatType observable,
    FloatType alpha,
    FloatType beta,
    FloatType scale_factor,
    FloatType epsilon,
    bool centric,
    FloatType &yo_eff,
    FloatType &w_eff)
  {
    yo_eff = observable;
    w_eff = 0;
    FloatType d_target;
    if (!mli_d_target_d_modulus_fc(fo_sq, sigma_fo_sq, observable, alpha, beta,
                                   scale_factor, epsilon, centric, d_target))
    {
      return false;
    }
    FloatType const k = (scale_factor <= 0) ? FloatType(1) : scale_factor;
    FloatType const a = alpha*k;
    FloatType const b = beta*k*k;
    FloatType const eb = epsilon*b;
    FloatType const w_amplitude = centric ? a*a/eb : 2*a*a/eb;
    FloatType const modulus_fc = std::sqrt(observable);
    // |Fc| - (dT/d|Fc|)/w, then into F^2 space exactly as above
    FloatType const yo_amplitude = modulus_fc - d_target/w_amplitude;
    // and the same near-zero guard as the amplitude case, for the same reason
    if (!ml_intensity_is_usable(fo_sq, sigma_fo_sq, observable)) {
      w_eff = 0;
      yo_eff = observable;
      return false;
    }
    w_eff = w_amplitude/(4*observable);
    yo_eff = 2*modulus_fc*yo_amplitude - observable;
    return true;
  }

}}}

#endif // GUARD
