#ifndef CCTBX_XRAY_DISPERSION_RADIAL_H
#define CCTBX_XRAY_DISPERSION_RADIAL_H

#include <cctbx/error.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <boost/shared_ptr.hpp>
#include <algorithm>
#include <complex>
#include <cmath>

namespace cctbx { namespace xray {

  /* A radial falloff shared by f' and f''.

  f' and f'' are stored on a scatterer as constants, whereas f0 decays with
  scattering angle. The resonant scattering comes from core electrons of finite
  extent, so it decays too, and with the same shape for both parts since they
  come from the same oscillator strengths. This gives them one:

    f'_eff  = fp  R_g
    f''_eff = fdp R_g
    R_g     = 1 + sum_{k=1..n} c_{g,k} u^k

  where u is either 1 - cos(2 theta) or sin(theta)/lambda; c.f. variable().
  The first is what the physics asks for and is the default. Since
  cos(2 theta) = 1 - 2 lambda^2 s^2, the two are not merely different units: a
  polynomial in cos(2 theta) is exactly a polynomial of the same degree in s^2,
  so it spans the even powers of s and only those, whereas powers of s span the
  odd ones too, which are not polynomials in cos(2 theta) at all. Choosing
  cos(2 theta) is therefore a restriction of the model, not a reparametrisation
  of it.

  There is deliberately no k=0 term: it would be exactly degenerate with fp and
  fdp themselves, scaling f' by a constant being the same thing as changing f'.
  Both choices of u vanish at theta = 0, so R(0) = 1 and the tabulated values
  are left alone at forward scattering, which is where they are by definition
  right. With every c zero R is 1 everywhere and the whole thing is a no-op,
  which is where a refinement starts from.

  g is a group index, so that the coefficients can be shared however the caller
  wants: one group per element, or one per scatterer, or one for everything. A
  scatterer whose group is negative is left alone. Which grouping was chosen is
  no business of this class, and neither is the check that a group has anomalous
  scattering to speak of -- a group whose scatterers all have fp = fdp = 0
  contributes no gradient at all and would make the normal matrix singular. The
  caller is expected to have excluded it, where it can say why.

  The coefficients refine as ordinary independent parameters, sitting at the
  tail of the gradient vector next to BASF, EXTI and the rest; c.f.
  smtbx::refinement::constraints::dispersion_radial_parameter.

  The per-reflection members are the working storage of one linearisation. This
  object is therefore not shareable between threads: it belongs to a single
  structure factor functor, which clones it in raw_fork(). Storing it beside
  one, the way fc_correction is stored, would leave the clone that accumulates
  the gradients and the clone that is read for them being different objects.
  */
  template <typename FloatType>
  struct dispersion_radial_correction {
    typedef std::complex<FloatType> complex_type;

    dispersion_radial_correction(af::shared<int> const &group_of_scatterer,
                                 int n_groups,
                                 int n_terms,
                                 bool grad=true,
                                 FloatType s_max=0,
                                 FloatType r_min=0.1,
                                 bool in_cos_two_theta=true,
                                 FloatType wavelength=0)
      : n_terms(n_terms),
        n_groups(n_groups),
        grad_index(-1),
        grad(grad),
        s_max(s_max),
        r_min(r_min),
        in_cos_two_theta(in_cos_two_theta),
        wavelength(wavelength),
        group_of_scatterer(group_of_scatterer),
        coefficients(n_groups*n_terms, 0.),
        u_pow(n_terms, 0.),
        grad_accum(n_groups*n_terms, complex_type(0)),
        grad_obs(n_groups*n_terms, 0.)
    {
      CCTBX_ASSERT(n_terms > 0);
      CCTBX_ASSERT(n_groups > 0);
      CCTBX_ASSERT(!in_cos_two_theta || wavelength > 0);
    }

    std::size_t n_param() const {
      return static_cast<std::size_t>(n_groups)*n_terms;
    }

    /** @brief The variable the polynomial is in.

    Either 1 - cos(2 theta) or sin(theta)/lambda, and the choice is not just a
    change of units: since cos(2 theta) = 1 - 2 lambda^2 s^2, a polynomial in
    one is exactly a polynomial of the same degree in s^2, i.e. it spans the
    *even* powers of s and only those. The powers of s span more than that, the
    odd ones not being polynomials in cos(2 theta) at all.

    Both vanish at theta = 0, which is what keeps R(0) = 1: the tabulated f' and
    f'' are the forward scattering limits, so there is nothing there to correct.
    A constant term would in any case be exactly degenerate with f' and f''
    themselves -- scaling f' by a constant is changing f'.
    */
    FloatType variable(FloatType d_star_sq) const {
      if (in_cos_two_theta) {
        // 1 - cos(2 theta) = 2 sin^2(theta) = lambda^2 d*^2 / 2
        return wavelength*wavelength*d_star_sq/2;
      }
      return std::sqrt(d_star_sq)/2;
    }

    /// Position on a reflection: cache u^1..u^n and start a fresh accumulation
    void set_d_star_sq(FloatType d_star_sq) const {
      FloatType const u = variable(d_star_sq);
      FloatType p = 1;
      for (int k=0; k < n_terms; ++k) {
        p *= u;
        u_pow[k] = p;
      }
      std::fill(grad_accum.begin(), grad_accum.end(), complex_type(0));
    }

    /// R for a scatterer at the reflection set_d_star_sq was last called with
    FloatType R(std::size_t i_scatterer) const {
      int const g = group_of_scatterer[i_scatterer];
      if (g < 0) return 1;
      FloatType const *c = &coefficients[g*n_terms];
      FloatType r = 1;
      for (int k=0; k < n_terms; ++k) {
        r += c[k]*u_pow[k];
      }
      return r;
    }

    /** @brief Accumulate dFc/dc for one scatterer.

    dF_dR is dFc/dR for this scatterer, i.e. grad_fp fp + grad_fdp fdp taken
    with the gradients of the *effective* f' and f'', before they are scaled by
    R. Since dR/dc_{g,k} is u^k, the rest is a multiplication.
    */
    void accumulate(std::size_t i_scatterer, complex_type const &dF_dR) const {
      int const g = group_of_scatterer[i_scatterer];
      if (g < 0) return;
      complex_type *a = &grad_accum[g*n_terms];
      for (int k=0; k < n_terms; ++k) {
        a[k] += dF_dR*u_pow[k];
      }
    }

    /// The gradients of the observable, once base_obs has converted them
    af::const_ref<FloatType> get_gradients() const {
      return grad_obs.const_ref();
    }

    /// R at an arbitrary d*^2, for reporting and for tests
    FloatType R_at(FloatType d_star_sq, int group) const {
      if (group < 0) return 1;
      FloatType const u = variable(d_star_sq);
      FloatType const *c = &coefficients[group*n_terms];
      FloatType r = 1, p = 1;
      for (int k=0; k < n_terms; ++k) {
        p *= u;
        r += c[k]*p;
      }
      return r;
    }

    /** @brief Pull a group back if the shifts drove R to nothing.

    Nothing bounds the coefficients, and one of them against several thousand
    reflections has enough leverage to absorb whatever the model of a heavy
    atom gets wrong. Left alone it will: on a real structure a single term ran
    to R = -11 at the resolution limit, which is f' and f'' reversed in sign
    several times over. That is not dispersion falling off, it is a free
    parameter mopping up.

    R is scaled back toward 1 -- the whole group by one factor, so the shape it
    found is kept and only its size is cut -- by the least amount that puts
    R >= r_min over the whole of [0, s_max]. Scaling toward 1 rather than
    clamping each coefficient keeps the correction on the line between where
    the minimiser wanted to go and no correction at all, which is what damping
    a shift does.

    Inert while s_max is 0: a caller which does not say what its data cover is
    not asking for this.

    Returns whether anything was changed.
    */
    bool validate(int n_samples=64) {
      if (s_max <= 0) return false;
      bool changed = false;
      // s_max is the resolution limit whichever variable the polynomial is in,
      // so the sampling is done in d*^2 and R_at converts
      FloatType const d_star_sq_max = 4*s_max*s_max;
      for (int g=0; g < n_groups; ++g) {
        FloatType t = 1;
        for (int i=0; i <= n_samples; ++i) {
          FloatType const r = R_at(d_star_sq_max*i/n_samples, g);
          if (r < r_min) {
            // (1 - r_min) + t (R - 1) >= 0, with R < 1 here
            FloatType const t_i = (1 - r_min)/(1 - r);
            if (t_i < t) t = t_i;
          }
        }
        if (t < 1) {
          for (int k=0; k < n_terms; ++k) {
            coefficients[g*n_terms + k] *= t;
          }
          changed = true;
        }
      }
      return changed;
    }

    /// A copy with the same coefficients but its own working storage
    boost::shared_ptr<dispersion_radial_correction> fork() const {
      boost::shared_ptr<dispersion_radial_correction> cr(
        new dispersion_radial_correction(group_of_scatterer, n_groups, n_terms,
                                         grad, s_max, r_min, in_cos_two_theta,
                                         wavelength));
      std::copy(coefficients.begin(), coefficients.end(),
                cr->coefficients.begin());
      cr->grad_index = grad_index;
      return cr;
    }

    int n_terms, n_groups, grad_index;
    bool grad;
    /// the data's outermost sin(theta)/lambda, or 0 to leave validate() inert
    FloatType s_max;
    /// how far R is allowed to fall before validate() pulls it back
    FloatType r_min;
    /// powers of 1 - cos(2 theta) rather than of sin(theta)/lambda
    bool in_cos_two_theta;
    /// needed only by the former, which is a function of theta and not of d*
    FloatType wavelength;
    /// which group each scatterer belongs to, negative for none
    af::shared<int> group_of_scatterer;
    /// n_groups*n_terms, group-major: c_{g,k} is coefficients[g*n_terms + k]
    af::shared<FloatType> coefficients;

    // working storage of one reflection
    mutable af::shared<FloatType> u_pow;
    mutable af::shared<complex_type> grad_accum;
    mutable af::shared<FloatType> grad_obs;
  };

}} // namespace cctbx::xray

#endif // GUARD
