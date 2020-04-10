#ifndef SMTBX_STRUCTURE_FACTORS_DIRECT_STANDARD_XRAY_H
#define SMTBX_STRUCTURE_FACTORS_DIRECT_STANDARD_XRAY_H

#include <scitbx/math/imaginary.h>
#include <scitbx/math/copysign.h>
#include <scitbx/array_family/owning_ref.h>
#include <cctbx/math/cos_sin_table.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <cctbx/xray/hr_ht_cache.h>
#include <smtbx/import_cctbx.h>
#include <smtbx/import_scitbx_af.h>
#include <smtbx/error.h>

#include <boost/shared_ptr.hpp>
#include <iostream>

#define SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS                            \
  typedef FloatType float_type;                                            \
  typedef std::complex<float_type> complex_type;                           \
  typedef scitbx::sym_mat3<complex_type> grad_u_star_type;                 \
  typedef af::tiny<complex_type, 3> grad_site_type;                        \
  typedef xray::scatterer<float_type> scatterer_type;                      \
  typedef scitbx::vec3<float_type> pol_vec_type;


#define SMTBX_STRUCTURE_FACTORS_DIRECT_ONE_SCATTERER_ONE_H_USING           \
  using core<float_type>::structure_factor;                                \
  using core<float_type>::grad_site;                                       \
  using core<float_type>::grad_u_star;                                     \
  using core<float_type>::grad_u_iso;                                      \
  using core<float_type>::grad_occ;                                        \
  using core<float_type>::grad_fp_star;                                    \
  using core<float_type>::grad_fdp_star;                                   \
  using base_t::hr_ht;                                                     \
  using base_t::d_star_sq;

namespace smtbx { namespace structure_factors { namespace direct {

  using cctbx::xray::structure_factors::hr_ht_group;
  using cctbx::xray::structure_factors::hr_ht_cache;


  namespace one_scatterer_one_h {

    /** Bare bundle of the structure factor of one scatterer
        for a given Miller index and its gradient wrt to the crystallagraphic
        parameters of that scatterer.
     */
    template <typename FloatType>
    struct core
    {
      SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS;

      complex_type structure_factor;
      grad_site_type grad_site;
      complex_type grad_fp;
      complex_type grad_fdp;
      grad_u_star_type grad_u_star;
      complex_type grad_u_iso, grad_occ;
      grad_u_star_type grad_fp_star, grad_fdp_star;
    };

    /** Base class for the linearisation or the evaluation of the structure
        factor for one scatterer for a given Miller index.

        This uses the CRTP, delegating to Heir the key step of the actual
        computation.
     */
    template <typename FloatType, class ExpI2PiFunctor, class Heir>
    struct base : core<FloatType>
    {
      SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS;

      hr_ht_cache<float_type> hr_ht;
      float_type d_star_sq;
      ExpI2PiFunctor const &exp_i_2pi;

      //vectors for the incident and scattered polarization, for refinement of
      //anisotropic fp, fdp tensors
      hr_ht_cache<float_type, pol_vec_type> e1r_e1t;
      hr_ht_cache<float_type, pol_vec_type> e2r_e2t;
      float_type pol_factor;
      bool has_polarization;
      float_type m_f0;

      /** Construct the linearisation or the evaluation for the given h
          in the given space-group.

       The functor exp_i_2pi_functor will be used to compute exp( i 2pi h.x ).
       */
      base(sgtbx::space_group const &space_group,
           miller::index<> const &h,
           float_type d_star_sq,
           ExpI2PiFunctor const &exp_i_2pi_functor)

        : hr_ht(exp_i_2pi_functor, space_group, h),
          d_star_sq(d_star_sq),
          exp_i_2pi(exp_i_2pi_functor),
          has_polarization(false)

      {}

      base(sgtbx::space_group const &space_group,
           miller::index<> const &h,
           float_type d_star_sq,
           ExpI2PiFunctor const &exp_i_2pi_functor,
           pol_vec_type const &pol_inc,
           pol_vec_type const &pol_scat,
           float_type const &pol_factor)

        : hr_ht(exp_i_2pi_functor, space_group, h),
          d_star_sq(d_star_sq),
          exp_i_2pi(exp_i_2pi_functor),
          e1r_e1t(exp_i_2pi_functor, space_group, pol_inc),
          e2r_e2t(exp_i_2pi_functor, space_group, pol_scat),
          has_polarization(true),
          pol_factor(pol_factor)

      {}

      /** Compute the structure factor of the given scatterer
       as well as its gradients wrt to the crystallographic parameters of that
       scatterer if requested so.

       The argument f0 is the elastic form factor
       for that type of chemical element at the miller index at hand.
       */
      void compute(scatterer_type const &scatterer,
                   float_type f0,
                   bool compute_grad)
      {
        Heir &heir = static_cast<Heir &> (*this);

        this->structure_factor = 0;
        if (compute_grad) {
          this->grad_site = grad_site_type(0, 0, 0);
          this->grad_u_star = grad_u_star_type(0, 0, 0, 0, 0, 0);
          this->grad_fp = 0;
          this->grad_fdp = 0;
          this->grad_fp_star = grad_u_star_type(0, 0, 0, 0, 0, 0);
          this->grad_fdp_star = grad_u_star_type(0, 0, 0, 0, 0, 0);
        }

        if (scatterer.flags.use_fp_fdp_aniso()) {
          SMTBX_ASSERT(has_polarization);
          m_f0 = f0;
          heir.compute_anisotropic_part(scatterer, compute_grad);
          FloatType dummy_ff = 1; //form factor was already used in aniso part
          heir.multiply_by_isotropic_part(scatterer, dummy_ff, compute_grad);
        }
        else {
          heir.compute_anisotropic_part(scatterer, compute_grad);
          if (scatterer.flags.use_fp_fdp()) {
            complex_type ff(f0 + scatterer.fp, scatterer.fdp);
            heir.multiply_by_isotropic_part(scatterer, ff, compute_grad);
          }
          else {
            heir.multiply_by_isotropic_part(scatterer, f0, compute_grad);
          }
        }
      }
    };



    /** Key steps of the evaluation or linearisation of the structure factor
        of one scatterer for a given Miller index in any space group.
     */
    template <typename FloatType, class ExpI2PiFunctor>
    struct in_generic_space_group : base<FloatType,
                                        ExpI2PiFunctor,
                                        in_generic_space_group<FloatType,
                                                               ExpI2PiFunctor> >

    {
      typedef base<FloatType,
                   ExpI2PiFunctor,
                   in_generic_space_group<FloatType,
                                          ExpI2PiFunctor> >
              base_t;

      SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS;
      SMTBX_STRUCTURE_FACTORS_DIRECT_ONE_SCATTERER_ONE_H_USING;

      /** \copydoc one_scatterer_one_h::base */
      in_generic_space_group(sgtbx::space_group const &space_group,
                             miller::index<> const &h,
                             float_type d_star_sq,
                             ExpI2PiFunctor const &exp_i_2pi_functor)

        : base_t(space_group, h, d_star_sq, exp_i_2pi_functor)
      {}

      in_generic_space_group(sgtbx::space_group const &space_group,
                             miller::index<> const &h,
                             float_type d_star_sq,
                             ExpI2PiFunctor const &exp_i_2pi_functor,
                             pol_vec_type const &pol_inc,
                             pol_vec_type const &pol_scat,
                             float_type pol_factor)
        : base_t(space_group, h, d_star_sq, exp_i_2pi_functor, pol_inc,
            pol_scat, pol_factor)
      {}

      /** Sum of Debye-Waller * Fourier basis function,
         over the Miller indices equivalent to h, and its gradients.

         A centric space group G is decomposed as G = G' U inversion*G'
         Denoting the respective sums as S' and S'', we seek S = S' + S''
         and its gradients.

        In general, S'' = conj(S') exp(i 2pi h.t_inv)
        where t_inv is the translation part of the inversion.
        It is cached in hr_ht.f_h_inv_t

        For a non-centric space group, G = G' and S = S'

        Centring translations are not included in those sums: they just
        result in a multiplicative factor later.
       */
      void compute_anisotropic_part(scatterer_type const &scatterer,
                                    bool compute_grad)
      {
        using namespace adptbx;
        using namespace scitbx::constants;
        complex_type const i(0,1);

        // Compute S'
        for (int k=0; k < hr_ht.groups.size(); ++k) {

          /* If f' and/or f" are polarization-anisotropic, then the form factor
           * depends on hR and further we need it before calculating grad_site,
           * etc. Thus we calculate it right away.
           * */
          complex_type formfactor;
          hr_ht_group<float_type, pol_vec_type> e1_g;
          hr_ht_group<float_type, pol_vec_type> e2_g;
          if (scatterer.flags.use_fp_fdp_aniso()) {
            e1_g = base_t::e1r_e1t.groups[k];
            e2_g = base_t::e2r_e2t.groups[k];
            float_type fp =
              e1_g.hr * scatterer.fp_aniso * e2_g.hr / base_t::pol_factor;
            float_type fdp =
              e1_g.hr * scatterer.fdp_aniso * e2_g.hr / base_t::pol_factor;
            formfactor =  complex_type(base_t::m_f0 + fp, fdp);
          }


          hr_ht_group<float_type> const &g = hr_ht.groups[k];
          float_type hrx = g.hr * scatterer.site;
          complex_type f = this->exp_i_2pi(hrx + g.ht);
          if (scatterer.flags.use_u_aniso()) {
            float_type dw = debye_waller_factor_u_star(g.hr, scatterer.u_star);
            f *= dw;
            if (compute_grad) {
              if (scatterer.flags.grad_u_aniso()) {
                scitbx::sym_mat3<float_type> log_grad_u_star
                  = debye_waller_factor_u_star_gradient_coefficients<
                      float_type>(g.hr);
                complex_type grad_u_star_factor = -two_pi_sq * f;
                if (scatterer.flags.use_fp_fdp_aniso()) {
                  grad_u_star_factor *= formfactor;
                }
                for (int j=0; j<6; ++j) {
                  grad_u_star[j] += grad_u_star_factor * log_grad_u_star[j];
                }
              }
            }
          }

          if (compute_grad && scatterer.flags.grad_site()) {
            complex_type grad_site_factor = i * two_pi * f;
            if (scatterer.flags.use_fp_fdp_aniso()) {
              grad_site_factor *= formfactor;
            }
            for (int j=0; j<3; ++j) {
              float_type h_j = g.hr[j];
              grad_site[j] += grad_site_factor * h_j;
            }
          }

          if (scatterer.flags.use_fp_fdp_aniso()) {
            if(compute_grad && scatterer.flags.grad_fp_aniso()) {
              grad_fp_star[0] += 
                  e1_g.hr[0] * e2_g.hr[0] * f / base_t::pol_factor;
              grad_fp_star[1] += 
                  e1_g.hr[1] * e2_g.hr[1] * f / base_t::pol_factor;
              grad_fp_star[2] += 
                  e1_g.hr[2] * e2_g.hr[2] * f / base_t::pol_factor;
              grad_fp_star[3] +=
                  (e1_g.hr[0]*e2_g.hr[1] + e1_g.hr[1]*e2_g.hr[0])
                  * f / base_t::pol_factor;
              grad_fp_star[4] +=
                  (e1_g.hr[0]*e2_g.hr[2] + e1_g.hr[2]*e2_g.hr[0])
                  * f / base_t::pol_factor;
              grad_fp_star[5] +=
                  (e1_g.hr[1]*e2_g.hr[2] + e1_g.hr[2]*e2_g.hr[1])
                  * f / base_t::pol_factor;
            }
            if(compute_grad && scatterer.flags.grad_fdp_aniso()) {
              grad_fdp_star[0] +=
                  i * e1_g.hr[0] * e2_g.hr[0] * f / base_t::pol_factor;
              grad_fdp_star[1] +=
                  i * e1_g.hr[1] * e2_g.hr[1] * f / base_t::pol_factor;
              grad_fdp_star[2] +=
                  i * e1_g.hr[2] * e2_g.hr[2] * f / base_t::pol_factor;
              grad_fdp_star[3] +=
                  i * (e1_g.hr[0]*e2_g.hr[1] + e1_g.hr[1]*e2_g.hr[0])
                  * f / base_t::pol_factor;
              grad_fdp_star[4] +=
                  i * (e1_g.hr[0]*e2_g.hr[2] + e1_g.hr[2]*e2_g.hr[0])
                  * f / base_t::pol_factor;
              grad_fdp_star[5] +=
                  i * (e1_g.hr[1]*e2_g.hr[2] + e1_g.hr[2]*e2_g.hr[1])
                  * f / base_t::pol_factor;
            }
            f *= formfactor;
          }

          structure_factor += f;
        }

        if (hr_ht.is_centric) {
          /*
           * This is a major problem! Might need to add the conjugate after each
           * hr_ht_group is processed. The basic problem is that we're factoring
           * in the form factor while looping over each operator, then lumping 
           * everything together and forgetting the contribution of the form
           * factor. Then you take the conj of the whole thing. Instead we need
           * to apply form factor to pairs (f + conj(f)).
           * */
          // Compute S = S' + conj(S') exp(i 2pi h.t_inv) and its gradients
          structure_factor += std::conj(structure_factor) * hr_ht.f_h_inv_t;
          if (compute_grad) {
            if (scatterer.flags.grad_site()) {
              for (int j=0; j<3; ++j) {
                grad_site[j] += std::conj(grad_site[j]) * hr_ht.f_h_inv_t;
              }
            }
            if (scatterer.flags.use_u_aniso() && scatterer.flags.grad_u_aniso())
            {
              for (int j=0; j<6; ++j) {
                grad_u_star[j] += std::conj(grad_u_star[j]) * hr_ht.f_h_inv_t;
              }
            }
            if (scatterer.flags.use_fp_fdp_aniso()
                && scatterer.flags.grad_fp_aniso())
            {
              for (int j=0; j<6; ++j)
              {
                grad_fp_star[j] +=
                  std::conj(grad_fp_star[j]) * hr_ht.f_h_inv_t;
              }
            }
            if (scatterer.flags.use_fp_fdp_aniso()
                && scatterer.flags.grad_fdp_aniso())
            {
              for (int j=0; j<6; ++j)
              {
                grad_fdp_star[j] +=
                  std::conj(grad_fdp_star[j]) * hr_ht.f_h_inv_t;
              }
            }
          }
        }

        // This balances the extra (conditional) curly brace that follows
        #if (0)
        {
        #endif

        #if (   defined(__linux__) && defined(__GNUC__) \
             && __GNUC__ == 4 && __GNUC_MINOR__ == 0 && __GNUC_PATCHLEVEL__ == 0)
        /** Careful analysis with valgrind showed that the compiler seems to
            generate undefined bits in the array grad_site in this member
            function.
            We hypothesised that an overzealous optimiser was at fault, hence
            the idea of the volatile member. That it solves the problem seems
            to vindicate our intuition.
        */
        foo = grad_site.begin();
      }
      complex_type volatile *foo;
        #else
      }
        #endif

      /** The isotropic factor ff * occupancy * isotropic debye-waller
          is factored out of the sum over equivalent reflections
          and multiplied with the latter.

          ff may be either the elastic form factor f0 (FormFactorType is a real
          number type) or the sum f0 + f' + if'' (FormFactorType is a complex
          number type) including the effect of inelastic scattering.
          The former case can be done more efficiently and the implementation
          thrives to take advantage of that.
       */
      template <typename FormFactorType>
      void multiply_by_isotropic_part(scatterer_type const &scatterer,
                                      FormFactorType const &ff,
                                      bool compute_grad)
      {
        using namespace adptbx;
        using namespace scitbx::constants;

        // factor from centring translation group * factor from special position
        float_type f_iso = hr_ht.ltr_factor * scatterer.weight_without_occupancy();

        // isotropic Debye-Waller
        if (scatterer.flags.use_u_iso()) {
          float_type dw_iso = debye_waller_factor_u_iso(d_star_sq/4,
                                                        scatterer.u_iso);
          f_iso *= dw_iso;
        }

        FormFactorType ff_iso_p = ff * f_iso;

        // occupancy
        if (compute_grad) {
          if (scatterer.flags.grad_occupancy()) {
            grad_occ = ff_iso_p * structure_factor;
          }
          if (scatterer.flags.grad_fp() || scatterer.flags.grad_fdp()) {
            complex_type p = f_iso * structure_factor * scatterer.occupancy;
            if (scatterer.flags.grad_fp()) {
              base_t::grad_fp = p;
            }
            if (scatterer.flags.grad_fdp()) {
              base_t::grad_fdp = complex_type(-p.imag(), p.real());
            }
          }
        }


        // Scattering factor
        FormFactorType ff_iso = ff_iso_p * scatterer.occupancy;

        // Finish
        structure_factor *= ff_iso;

        if (!compute_grad) return;

        if (scatterer.flags.use_u_iso() && scatterer.flags.grad_u_iso()) {
          grad_u_iso = -two_pi_sq * d_star_sq * structure_factor;
        }
        if (ff_iso != complex_type(1,0)) {
          if (scatterer.flags.grad_site()) {
            for (int j=0; j<3; ++j) grad_site[j] *= ff_iso;
          }
          if (scatterer.flags.grad_u_aniso()) {
            for (int j=0; j<6; ++j) grad_u_star[j] *= ff_iso;
          }
          if (scatterer.flags.grad_fp_aniso()) {
            for (int j=0; j<6; ++j) grad_fp_star[j] *= ff_iso;
          }
          if (scatterer.flags.grad_fdp_aniso()) {
            for (int j=0; j<6; ++j) grad_fdp_star[j] *= ff_iso;
          }
        }
      }
    };


    /** Key steps of the evaluation or linearisation of the structure factor
        of one scatterer for a given Miller index
        in an origin centric space group.
     */
    template <typename FloatType, class ExpI2PiFunctor>
    struct in_origin_centric_space_group
      : base<FloatType,
             ExpI2PiFunctor,
             in_origin_centric_space_group<FloatType,
                                           ExpI2PiFunctor> >
    {
      typedef
        base<FloatType,
             ExpI2PiFunctor,
             in_origin_centric_space_group<FloatType,
                                           ExpI2PiFunctor> >
        base_t;

      SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS;
      SMTBX_STRUCTURE_FACTORS_DIRECT_ONE_SCATTERER_ONE_H_USING;

      /** \copydoc one_scatterer_one_h::base */
      in_origin_centric_space_group(sgtbx::space_group const &space_group,
                                    miller::index<> const &h,
                                    float_type d_star_sq,
                                    ExpI2PiFunctor const &exp_i_2pi_functor)

        : base_t(space_group, h, d_star_sq, exp_i_2pi_functor)
      {}

      in_origin_centric_space_group(sgtbx::space_group const &space_group,
                             miller::index<> const &h,
                             float_type d_star_sq,
                             ExpI2PiFunctor const &exp_i_2pi_functor,
                             pol_vec_type const &pol_inc,
                             pol_vec_type const &pol_scat,
                             float_type pol_factor)
        : base_t(space_group, h, d_star_sq, exp_i_2pi_functor, pol_inc,
            pol_scat, pol_factor)
      {}

      /** Sum of Debye-Waller * Fourier basis function,
       over the Miller indices equivalent to h, and its gradients.

       The origin centric spacegroup G is decomposed as G = G' U inversion*G'
       Denoting the respective sums as S' and S'', we seek S = S' + S''
       and its gradients. Then S'' = conj(S') and thererefore S is real:
       S = 2 real( S' ). As a result, one only needs to compute the real
       part of S' and its derivatives.

       Centring translations are not included in those sums: they just
       result in a multiplicative factor later, which is combined with the
       factor 2 above.
       */
      void compute_anisotropic_part(scatterer_type const &scatterer,
                                    bool compute_grad)
      {
        using namespace adptbx;
        using namespace scitbx::constants;
        complex_type const i(0,1);

        /* Compute the real part of S' and its gradients
         We only ever touch the real part of structure_factors and grad_xxx
         members. So this code should be as efficient as working with
         real numbers.
         */



        for (int k=0; k < hr_ht.groups.size(); ++k) {

          /* If f' and/or f" are polarization-anisotropic, then the form factor
           * depends on hR and further we need it before calculating grad_site,
           * etc. Thus we calculate it right away.
           * */
          complex_type formfactor;
          hr_ht_group<float_type, pol_vec_type> e1_g;
          hr_ht_group<float_type, pol_vec_type> e2_g;
          if (scatterer.flags.use_fp_fdp_aniso()) {
            e1_g = base_t::e1r_e1t.groups[k];
            e2_g = base_t::e2r_e2t.groups[k];
            float_type fp =
              e1_g.hr * scatterer.fp_aniso * e2_g.hr / base_t::pol_factor;
            float_type fdp =
              e1_g.hr * scatterer.fdp_aniso * e2_g.hr / base_t::pol_factor;
            formfactor =  complex_type(base_t::m_f0 + fp, fdp);
          }

          hr_ht_group<float_type> const &g = hr_ht.groups[k];
          float_type hrx = g.hr * scatterer.site;
          complex_type f = this->exp_i_2pi(hrx + g.ht);
          if (scatterer.flags.use_u_aniso()) {
            float_type dw = debye_waller_factor_u_star(g.hr, scatterer.u_star);
            f *= dw;
            if (compute_grad && scatterer.flags.grad_u_aniso()) {
              scitbx::sym_mat3<float_type> log_grad_u_star
                = debye_waller_factor_u_star_gradient_coefficients<
                    float_type>(g.hr);
              complex_type grad_u_star_factor = -two_pi_sq * f.real();
              if (scatterer.flags.use_fp_fdp_aniso()) {
                grad_u_star_factor *= formfactor;
              }
              for (int j=0; j<6; ++j) {
                grad_u_star[j] += grad_u_star_factor * log_grad_u_star[j];
              }
            }
          }

          if (compute_grad && scatterer.flags.grad_site()) {
            complex_type grad_site_factor = -two_pi * f.imag();
            if (scatterer.flags.use_fp_fdp_aniso()) {
              grad_site_factor *= formfactor;
            }
            for (int j=0; j<3; ++j) {
              float_type h_j = g.hr[j];
              grad_site[j] += grad_site_factor * h_j;
            }
          }

          complex_type f_r(f.real(), 0.);



          if (scatterer.flags.use_fp_fdp_aniso()) {

            if(compute_grad && scatterer.flags.grad_fp_aniso()) {
              grad_fp_star[0] += 
                  e1_g.hr[0] * e2_g.hr[0] * f_r / base_t::pol_factor;
              grad_fp_star[1] += 
                  e1_g.hr[1] * e2_g.hr[1] * f_r / base_t::pol_factor;
              grad_fp_star[2] += 
                  e1_g.hr[2] * e2_g.hr[2] * f_r / base_t::pol_factor;
              grad_fp_star[3] +=
                  (e1_g.hr[0]*e2_g.hr[1] + e1_g.hr[1]*e2_g.hr[0])
                  * f_r / base_t::pol_factor;
              grad_fp_star[4] +=
                  (e1_g.hr[0]*e2_g.hr[2] + e1_g.hr[2]*e2_g.hr[0])
                  * f_r / base_t::pol_factor;
              grad_fp_star[5] +=
                  (e1_g.hr[1]*e2_g.hr[2] + e1_g.hr[2]*e2_g.hr[1])
                  * f_r / base_t::pol_factor;
            }
            if(compute_grad && scatterer.flags.grad_fdp_aniso()) {
              grad_fdp_star[0] +=
                  i * e1_g.hr[0] * e2_g.hr[0] * f_r / base_t::pol_factor;
              grad_fdp_star[1] +=
                  i * e1_g.hr[1] * e2_g.hr[1] * f_r / base_t::pol_factor;
              grad_fdp_star[2] +=
                  i * e1_g.hr[2] * e2_g.hr[2] * f_r / base_t::pol_factor;
              grad_fdp_star[3] +=
                  i * (e1_g.hr[0]*e2_g.hr[1] + e1_g.hr[1]*e2_g.hr[0])
                  * f_r / base_t::pol_factor;
              grad_fdp_star[4] +=
                  i * (e1_g.hr[0]*e2_g.hr[2] + e1_g.hr[2]*e2_g.hr[0])
                  * f_r / base_t::pol_factor;
              grad_fdp_star[5] +=
                  i * (e1_g.hr[1]*e2_g.hr[2] + e1_g.hr[2]*e2_g.hr[1])
                  * f_r / base_t::pol_factor;
            }

            f_r *= formfactor;
          }


          structure_factor += f_r;
        }

        /* We should now multiply S' and its gradients by 2
         to get S and its gradients. Postponed to next stage for efficiency.
         */
      }

      /** \copydoc in_generic_space_group::multiply_by_isotropic_part */
      template <typename FormFactorType>
      void multiply_by_isotropic_part(scatterer_type const &scatterer,
                                      FormFactorType const &ff,
                                      bool compute_grad)
      {
        using namespace adptbx;
        using namespace scitbx::constants;

        // factor from inversion
        float_type f_iso = 2.;

        // factor from centring translation group * factor from special position
        f_iso *= hr_ht.ltr_factor * scatterer.weight_without_occupancy();

        // isotropic Debye-Waller
        if (scatterer.flags.use_u_iso()) {
          float_type dw_iso = debye_waller_factor_u_iso(d_star_sq/4,
                                                        scatterer.u_iso);
          f_iso *= dw_iso;
        }

        FormFactorType ff_iso_p = ff * f_iso;

        // occupancy
        if (compute_grad) {
          if (scatterer.flags.grad_occupancy()) {
            grad_occ = ff_iso_p * structure_factor;
          }
          if (scatterer.flags.grad_fp() || scatterer.flags.grad_fdp()) {
            complex_type p = f_iso * structure_factor * scatterer.occupancy;
            if (scatterer.flags.grad_fp()) {
              base_t::grad_fp = p;
            }
            if (scatterer.flags.grad_fdp()) {
              base_t::grad_fdp = complex_type(0, 1) * p;
            }
          }
        }

        // Scattering factor
        FormFactorType ff_iso = ff_iso_p * scatterer.occupancy;

        // Finish
        structure_factor = ff_iso * structure_factor;

        if (!compute_grad) return;

        if (scatterer.flags.use_u_iso() && scatterer.flags.grad_u_iso()) {
          grad_u_iso = -two_pi_sq * d_star_sq * structure_factor;
        }
        if (scatterer.flags.grad_site()) {
          for (int j=0; j<3; ++j) {
            grad_site[j] = ff_iso * grad_site[j];
          }
        }
        if (scatterer.flags.grad_u_aniso()) {
          for (int j=0; j<6; ++j) {
            grad_u_star[j] = ff_iso * grad_u_star[j];
          }
        }
        if (scatterer.flags.grad_fp_aniso()) {
          for (int j=0; j<6; ++j) {
            grad_fp_star[j] *= ff_iso;
          }
        }
        if (scatterer.flags.grad_fdp_aniso()) {
          for (int j=0; j<6; ++j) {
            grad_fdp_star[j] *= ff_iso;
          }
        }
      }
    };

  } // namespace one_scatterer_one_h


  namespace one_h {

    /** @brief Evaluation or linearisation of \f$F_c(h)\f$
        and of a derived observable, for any miller index h,
        as functions of crystallographic parameters.

        The observable is modelled by the type ObservableType, two examples
        of which are the class modulus and modulus_squared in this namespace.

        The ordering of the array \f$\nabla F_c(h)\f$ is a very important
        property since it has to be known by some key client code, e.g.
        the constraint framework in smtbx::refinement::contraints. This class
        guarantees that the partial derivatives are stored in the following
        order for each scatterer:
          x, y, z, u_iso, u_11, u_22, u_33, u_12, u_13, u_23, occupancy
        If a scatterer's site is not refined, then x,y,z is taken out, etc.
        Then that scheme is applied for all scatterers in the order of the
        array of xray::scatterer passed to the constructor.

        This class is a CRTP: class Heir should be a heir of this class
        and provide

        - a member exp_i_2pi of type ExpI2PiFunctor to compute
          \f$\exp i2\pi x\f$;
        - a member function raw_fork to implement member function fork.
     */
    template <typename FloatType,
              template<typename> class ObservableType,
              template<typename> class ExpI2PiFunctor,
              class Heir>
    class base
    {
    public:
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;
      typedef ObservableType<float_type> observable_type;
      typedef ExpI2PiFunctor<float_type> exp_i_2pi_functor;
      typedef boost::shared_ptr<Heir> pointer_type;
      typedef scitbx::vec3<float_type> pol_vec_type;

    protected:
      cctbx::xray::scatterer_grad_flags_counts grad_flags_counts;
      uctbx::unit_cell const &unit_cell;
      sgtbx::space_group const &space_group;
      bool origin_centric_case;
      af::ref_owning_shared< xray::scatterer<float_type> > scatterers;
      xray::scattering_type_registry const &scattering_type_registry;
      af::shared<std::size_t> scattering_type_indices;

      complex_type *grad_f_calc_cursor;
      bool has_computed_grad;

    public:
      complex_type f_calc;
      af::ref_owning_shared<complex_type> grad_f_calc;
      float_type observable;
      af::ref_owning_shared<float_type> grad_observable;

      // A second Fc at this h; e.g. for the second polarization direction when
      // polarization anisotropy is present
      complex_type f_calc_prime;
      af::ref_owning_shared<complex_type> grad_f_calc_prime;
      complex_type *grad_f_calc_prime_cursor;

    public:
      /** @brief The evaluation or linearisation of \f$F_c\f$
          for the given structure is to be computed.

          After construction, only the value of the refined parameters may
          change. That is, any of the following operations are forbidden:

            - adding or removing scatterers
            - changing the scattering type of a scatterer
            - disabling the refinement of some parameters
            - enabling the refinement of some parameters

          Not only the computation will be incorrect but illegal memory
          accesses will be prone to occur.
       */
      base(uctbx::unit_cell const &unit_cell,
           sgtbx::space_group const &space_group,
           af::shared< xray::scatterer<float_type> > const &scatterers,
           xray::scattering_type_registry const &scattering_type_registry)

        : grad_flags_counts(scatterers.const_ref()),
          unit_cell(unit_cell),
          space_group(space_group),
          origin_centric_case( space_group.is_origin_centric() ),
          scatterers(scatterers),
          scattering_type_registry(scattering_type_registry),
          scattering_type_indices(
            scattering_type_registry.unique_indices(this->scatterers) ),
          grad_f_calc(grad_flags_counts.n_parameters(),
                      af::init_functor_null<complex_type>()),
          grad_f_calc_prime(grad_flags_counts.n_parameters(),
                      af::init_functor_null<complex_type>()),
          grad_observable(grad_flags_counts.n_parameters(),
                          af::init_functor_null<float_type>()),
          has_computed_grad(false)
      {}

      /// A functor able to independently perform the same crystallographic
      /// computation as this
      /** That is to say that it shares with this the same unit cell,
          space group, scatterers and scattering factors but it has its own
          storage for f_calc and its gradients, as well as the observable and
          its gradients.
       */
      pointer_type fork() {
        Heir &heir = static_cast<Heir &>(*this);
        return pointer_type(heir.raw_fork());
      }
      /// Evaluate the structure factors
      void evaluate(miller::index<> const &h,
                    boost::optional<complex_type> const &f_mask=boost::none)
      {
        compute(h, f_mask, false);
      }

      /// Linearise the structure factors
      void linearise(miller::index<> const &h,
                     boost::optional<complex_type> const &f_mask=boost::none)
      {
        compute(h, f_mask, true);
      }

      /// Whether this is a linearisation of Fc.
      /** As opposed to a mere evaluation: for the latter, only the value
       of Fc was computed whereas for the former, its gradient was computed
       as well.
       */
      bool is_linearisation() const { return has_computed_grad; }

      /// Compute the evaluation or the linearisation
      void compute(miller::index<> const &h,
                   boost::optional<complex_type> const &f_mask=boost::none,
                   bool compute_grad=true)
      {
        float_type d_star_sq = unit_cell.d_star_sq(h);
        af::shared<float_type> elastic_form_factors
          = scattering_type_registry.unique_form_factors_at_d_star_sq(d_star_sq);

        Heir &heir = static_cast<Heir &>(*this);

        typedef one_scatterer_one_h::in_generic_space_group<
                  float_type, exp_i_2pi_functor>
                generic_linearisation_t;

        typedef one_scatterer_one_h::in_origin_centric_space_group<
                  float_type, exp_i_2pi_functor>
                origin_centric_linearisation_t;

        if (!origin_centric_case) {
           generic_linearisation_t lin_for_h(space_group, h, d_star_sq,
                                             heir.exp_i_2pi);
          compute(elastic_form_factors.ref(), lin_for_h, f_mask, compute_grad);
        }
        else {
          origin_centric_linearisation_t lin_for_h(space_group, h, d_star_sq,
                                                   heir.exp_i_2pi);
          compute(elastic_form_factors.ref(), lin_for_h, f_mask, compute_grad);
        }
      }

      void compute(miller::index<> const h,
                   pol_vec_type u_inc, // u,v convention as in Schiltz 2008
                   pol_vec_type u_scat,
                   pol_vec_type v_scat,
                   float_type pol_factor,
                   boost::optional<complex_type> const &f_mask=boost::none,
                   bool compute_grad=true)
      {
        float_type d_star_sq = unit_cell.d_star_sq(h);
        af::shared<float_type> elastic_form_factors
          = scattering_type_registry.unique_form_factors_at_d_star_sq(d_star_sq);

        Heir &heir = static_cast<Heir &> (*this);

        typedef one_scatterer_one_h::in_generic_space_group<
                  float_type, exp_i_2pi_functor>
                generic_linearisation_t;

        typedef one_scatterer_one_h::in_origin_centric_space_group<
                  float_type, exp_i_2pi_functor>
                origin_centric_linearisation_t;

        if (!origin_centric_case) {
          generic_linearisation_t lin_for_h_parallel(space_group, h,
                                                     d_star_sq, heir.exp_i_2pi,
                                                     u_inc, u_scat, pol_factor);
          generic_linearisation_t lin_for_h_perpendicular(space_group, h,
                                                     d_star_sq, heir.exp_i_2pi,
                                                     u_inc, v_scat, pol_factor);
          compute(elastic_form_factors.ref(), lin_for_h_parallel,
                  lin_for_h_perpendicular, f_mask, compute_grad);
        }
        else {
          origin_centric_linearisation_t lin_for_h_parallel(space_group, h,
                                                     d_star_sq, heir.exp_i_2pi,
                                                     u_inc, u_scat, pol_factor);
          origin_centric_linearisation_t lin_for_h_perpendicular(space_group, h,
                                                     d_star_sq, heir.exp_i_2pi,
                                                     u_inc, v_scat, pol_factor);
          compute(elastic_form_factors.ref(), lin_for_h_parallel,
                  lin_for_h_perpendicular, f_mask, compute_grad);
        }

      }


    private:
      template <class LinearisationForMillerIndex>
      void compute(af::const_ref<float_type> const &elastic_form_factors,
                   LinearisationForMillerIndex &l,
                   boost::optional<complex_type> const &f_mask,
                   bool compute_grad)
      {
        f_calc = 0;
        grad_f_calc_cursor = grad_f_calc.begin();

        for (int j=0; j < scatterers.size(); ++j) {
          xray::scatterer<> const &sc = scatterers[j];
          float_type f0 = elastic_form_factors[ scattering_type_indices[j] ];
          l.compute(sc, f0, compute_grad);

          f_calc += l.structure_factor;

          if (!compute_grad) continue;

          if (sc.flags.grad_site()) {
            for (int j=0; j<3; ++j) {
              *grad_f_calc_cursor++ = l.grad_site[j];
            }
          }
          if (sc.flags.use_u_iso() && sc.flags.grad_u_iso()) {
            *grad_f_calc_cursor++ = l.grad_u_iso;
          }
          if (sc.flags.use_u_aniso() && sc.flags.grad_u_aniso()) {
            for (int j=0; j<6; ++j) {
              *grad_f_calc_cursor++ = l.grad_u_star[j];
            }
          }
          if (sc.flags.grad_occupancy()) {
            *grad_f_calc_cursor++ = l.grad_occ;
          }
          if (sc.flags.grad_fp()) {
            *grad_f_calc_cursor++ = l.grad_fp;
          }
          if (sc.flags.grad_fdp()) {
            *grad_f_calc_cursor++ = l.grad_fdp;
          }
        }
        if (f_mask) {
          f_calc += *f_mask;
        }
        observable_type::compute(origin_centric_case,
                                 f_calc, grad_f_calc,
                                 observable, grad_observable,
                                 compute_grad);
        has_computed_grad = compute_grad;
      }

      template <class LinearisationForMillerIndex>
      void compute(af::const_ref<float_type> const &elastic_form_factors,
                   LinearisationForMillerIndex &l,
                   LinearisationForMillerIndex &l_prime,
                   boost::optional<complex_type> const &f_mask,
                   bool compute_grad)
      {
        f_calc = 0;
        f_calc_prime = 0;
        grad_f_calc_cursor = grad_f_calc.begin();
        grad_f_calc_prime_cursor = grad_f_calc_prime.begin();

        for (int j=0; j < scatterers.size(); ++j) {
          xray::scatterer<> const &sc = scatterers[j];
          bool ff_aniso = sc.flags.use_fp_fdp_aniso();
          float_type f0 = elastic_form_factors[ scattering_type_indices[j] ];
          l.compute(sc, f0, compute_grad);
          if (ff_aniso) {
            l_prime.compute(sc, 0., compute_grad);
          }

          f_calc += l.structure_factor;
          f_calc_prime += (ff_aniso ? l_prime.structure_factor : 0);

          if (!compute_grad) continue;

          if (sc.flags.grad_site()) {
            for (int j=0; j<3; ++j) {
              *grad_f_calc_cursor++ = l.grad_site[j];
              *grad_f_calc_prime_cursor++ =
                  (ff_aniso ? l_prime.grad_site[j] : 0);
            }
          }
          if (sc.flags.use_u_iso() && sc.flags.grad_u_iso()) {
            *grad_f_calc_cursor++ = l.grad_u_iso;
            *grad_f_calc_prime_cursor++ =
                (ff_aniso ? l_prime.grad_u_iso : 0);
          }
          if (sc.flags.use_u_aniso() && sc.flags.grad_u_aniso()) {
            for (int j=0; j<6; ++j) {
              *grad_f_calc_cursor++ = l.grad_u_star[j];
              *grad_f_calc_prime_cursor++ =
                  (ff_aniso ? l_prime.grad_u_star[j] : 0);
            }
          }
          if (sc.flags.grad_occupancy()) {
            *grad_f_calc_cursor++ = l.grad_occ;
            *grad_f_calc_prime_cursor++ =
                (ff_aniso ? l_prime.grad_occ : 0);
          }
          if (sc.flags.grad_fp()) {
            *grad_f_calc_cursor++ = l.grad_fp;
            *grad_f_calc_prime_cursor++ =
                (ff_aniso ? l_prime.grad_fp : 0);
          }
          if (sc.flags.grad_fdp()) {
            *grad_f_calc_cursor++ = l.grad_fdp;
            *grad_f_calc_prime_cursor++ =
                (ff_aniso ? l_prime.grad_fdp : 0);
          }
        }

        observable_type::compute(origin_centric_case,
                                 f_calc, grad_f_calc,
                                 f_calc_prime, grad_f_calc_prime,
                                 observable, grad_observable,
                                 compute_grad);
        has_computed_grad = compute_grad;
      }

    };




    /// Specialisation of class base computing trigonometric functions
    /// with a custom functor of type ExpI2PiFunctor
    template <typename FloatType,
              template<typename> class ObservableType,
              template<typename> class ExpI2PiFunctor>
    class custom_trigonometry
      : public base<FloatType, ObservableType, ExpI2PiFunctor,
                    custom_trigonometry<
                      FloatType, ObservableType, ExpI2PiFunctor> >
    {
    public:
      typedef base<FloatType, ObservableType, ExpI2PiFunctor,
                   custom_trigonometry<
                     FloatType, ObservableType, ExpI2PiFunctor> >
              base_t;

      typedef FloatType float_type;
      ExpI2PiFunctor<float_type> const &exp_i_2pi;

      custom_trigonometry(
        uctbx::unit_cell const &unit_cell,
        sgtbx::space_group const &space_group,
        af::shared< xray::scatterer<float_type> > const &scatterers,
        xray::scattering_type_registry const &scattering_type_registry,
        ExpI2PiFunctor<float_type> const &exp_i_2pi)

        : base_t(unit_cell, space_group, scatterers, scattering_type_registry),
          exp_i_2pi(exp_i_2pi)
      {}

      custom_trigonometry *raw_fork() {
        return new custom_trigonometry(this->unit_cell,
                                       this->space_group,
                                       this->scatterers.array(),
                                       this->scattering_type_registry,
                                       exp_i_2pi);
      }
    };


    /// Specialisation of class base computing trigonometric functions
    /// with the C++ standard library
    template <typename FloatType,
              template<typename> class ObservableType>
    class std_trigonometry
      : public base<FloatType, ObservableType,
                    cctbx::math::cos_sin_exact,
                    std_trigonometry<FloatType, ObservableType> >
    {
    public:
      typedef base<FloatType, ObservableType,
                   cctbx::math::cos_sin_exact,
                   std_trigonometry<FloatType, ObservableType> >
              base_t;

      typedef FloatType float_type;
      cctbx::math::cos_sin_exact<float_type> exp_i_2pi;

      std_trigonometry(
        uctbx::unit_cell const &unit_cell,
        sgtbx::space_group const &space_group,
        af::shared< xray::scatterer<float_type> > const &scatterers,
        xray::scattering_type_registry const &scattering_type_registry)

        : base_t(unit_cell, space_group, scatterers, scattering_type_registry)
      {}

      std_trigonometry *raw_fork() {
        return new std_trigonometry(this->unit_cell,
                                    this->space_group,
                                    this->scatterers.array(),
                                    this->scattering_type_registry);
      }
    };

    template <typename FloatType>
    struct modulus_squared
    {
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;

      static
      void compute(bool origin_centric_case,
                   complex_type f_calc,
                   af::const_ref<complex_type> const &grad_f_calc,
                   float_type &observable,
                   af::ref<float_type> const &grad_observable,
                   bool compute_grad)
      {
        if (origin_centric_case && f_calc.imag() == 0) {
          observable = f_calc.real() * f_calc.real();
        }
        else {
          observable = std::norm(f_calc);
        }

        if (!compute_grad) return;

        if (!origin_centric_case) {
          for (int j=0; j < grad_f_calc.size(); ++j) {
            grad_observable[j] = 2 * (  f_calc.real() * grad_f_calc[j].real()
                                      + f_calc.imag() * grad_f_calc[j].imag() );
          }
        }
        else {
          if (f_calc.imag() == 0) {
            for (int j=0; j < grad_f_calc.size(); ++j) {
                grad_observable[j] =  2 * f_calc.real() * grad_f_calc[j].real();
            }
          }
          else {
            for (int j=0; j < grad_f_calc.size(); ++j) {
              float_type &grad = grad_observable[j];
              grad = f_calc.real() * grad_f_calc[j].real();
              if (grad_f_calc[j].imag()) {
                grad += f_calc.imag() * grad_f_calc[j].imag();
              }
              grad *= 2;
            }
          }
        }
      }

      static
      void compute(bool origin_centric_case,
                   complex_type f_calc,
                   af::const_ref<complex_type> const &grad_f_calc,
                   complex_type f_calc_prime,
                   af::const_ref<complex_type> const &grad_f_calc_prime,
                   float_type &observable,
                   af::ref<float_type> const &grad_observable,
                   bool compute_grad)
      {
        if (origin_centric_case && f_calc.imag()==0 && f_calc_prime.imag()==0) {
          observable = (  f_calc.real() * f_calc.real()
                        + f_calc_prime.real() * f_calc_prime.real() );
        }
        else {
          observable = std::norm(f_calc) + std::norm(f_calc_prime);
        }

        if (!compute_grad) return;

        if (!origin_centric_case) {
          for (int j=0; j < grad_f_calc.size(); ++j) {
            grad_observable[j] = 2 * (  f_calc.real() * grad_f_calc[j].real()
                                      + f_calc.imag() * grad_f_calc[j].imag()
                                      + f_calc_prime.real()
                                          * grad_f_calc_prime[j].real()
                                      + f_calc_prime.imag()
                                          * grad_f_calc_prime[j].imag() );
          }
        }
        else {
          if (f_calc.imag()==0 && f_calc_prime.imag()==0) {
            for (int j=0; j < grad_f_calc.size(); ++j) {
              grad_observable[j] = 2 * (  f_calc.real() * grad_f_calc[j].real()
                                        + f_calc_prime.real()
                                            * grad_f_calc_prime[j].real() );
            }
          }
          else {
            for (int j=0; j < grad_f_calc.size(); ++j) {
              float_type &grad = grad_observable[j];
              grad = f_calc.real() * grad_f_calc[j].real();
              grad += f_calc_prime.real() * grad_f_calc_prime[j].real();
              if (f_calc.imag() && grad_f_calc[j].imag()) {
                grad += f_calc.imag() * grad_f_calc[j].imag();
              }
              if (f_calc_prime.imag() && grad_f_calc_prime.imag()) {
                grad += f_calc_prime.imag() * grad_f_calc_prime[j].imag();
              }
              grad *= 2;
            }
          }
        }
      }
    };


    template <typename FloatType>
    struct modulus
    {
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;

      static
      void compute(bool origin_centric_case,
                   complex_type f_calc,
                   af::const_ref<complex_type> const &grad_f_calc,
                   float_type &observable,
                   af::ref<float_type> const &grad_observable,
                   bool compute_grad)
      {
        if (origin_centric_case && f_calc.imag() == 0) {
          observable = std::abs(f_calc.real());
        }
        else {
          observable = std::abs(f_calc);
        }

        if (!compute_grad) return;

        if (!origin_centric_case) {
          float_type observable_inverse = 1./observable;
          for (int j=0; j < grad_f_calc.size(); ++j) {
            grad_observable[j] =   f_calc.real() * grad_f_calc[j].real()
                                 + f_calc.imag() * grad_f_calc[j].imag();
            grad_observable[j] *= observable_inverse;
          }
        }
        else {
          if (f_calc.imag() == 0) {
            for (int j=0; j < grad_f_calc.size(); ++j) {
              grad_observable[j] = f_calc.real() > 0 ?  grad_f_calc[j].real()
                                                     : -grad_f_calc[j].real();
            }
          }
          else {
            float_type observable_inverse = 1./observable;
            for (int j=0; j < grad_f_calc.size(); ++j) {
              float_type &grad = grad_observable[j];
              grad = f_calc.real() * grad_f_calc[j].real();
              if (grad_f_calc[j].imag()) {
                grad += f_calc.imag() * grad_f_calc[j].imag();
              }
              grad *= observable_inverse;
            }
          }
        }
      }
    };
  } // namespace one_h


}}} // smtbx::structure_factors::direct

#endif // GUARD
