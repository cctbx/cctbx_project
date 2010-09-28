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

#include <iostream>

namespace smtbx { namespace structure_factors { namespace direct {

  using cctbx::xray::structure_factors::hr_ht_group;
  using cctbx::xray::structure_factors::hr_ht_cache;


  namespace one_scatterer_one_h_linearisation {

    #define SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS                            \
      typedef FloatType float_type;                                            \
      typedef std::complex<float_type> complex_type;                           \
      typedef scitbx::sym_mat3<complex_type> grad_u_star_type;                 \
      typedef af::tiny<complex_type, 3> grad_site_type;                        \
      typedef xray::scatterer<float_type> scatterer_type;

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
      grad_u_star_type grad_u_star;
      complex_type grad_u_iso, grad_occ;
    };

    /** Base class for the linearisation of the structure factor
        for one scatterer for a given Miller index.
     */
    template <typename FloatType, class ExpI2PiFunctor, class Heir>
    struct base : core<FloatType>
    {
      SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS;

      hr_ht_cache<float_type> hr_ht;
      float_type d_star_sq;
      ExpI2PiFunctor const &exp_i_2pi;

      /** Construct the linearisation for the given h in the given space-group.

       The functor exp_i_2pi_functor will be used to compute exp( i 2pi h.x ).
       */
      base(sgtbx::space_group const &space_group,
           miller::index<> const &h,
           float_type d_star_sq,
           ExpI2PiFunctor const &exp_i_2pi_functor)

        : hr_ht(exp_i_2pi_functor, space_group, h),
          d_star_sq(d_star_sq),
          exp_i_2pi(exp_i_2pi_functor)
      {}

      /** Compute the structure factor of the given scatterer
       and its gradients wrt to the crystallographic parameters of that
       scatterer.

       The argument f0 is the elastic form factor
       for that type of chemical element at the miller index at hand.
       */
      void compute(scatterer_type const &scatterer, float_type f0)
      {
        Heir &heir = static_cast<Heir &> (*this);

        this->structure_factor = 0;
        this->grad_site = grad_site_type(0, 0, 0);
        this->grad_u_star = grad_u_star_type(0, 0, 0, 0, 0, 0);

        heir.compute_anisotropic_part(scatterer);
        if (scatterer.flags.use_fp_fdp()) {
          complex_type ff(f0 + scatterer.fp, scatterer.fdp);
          heir.multiply_by_isotropic_part(scatterer, ff);
        }
        else {
          heir.multiply_by_isotropic_part(scatterer, f0);
        }
      }
    };


    #define SMTBX_STRUCTURE_FACTORS_DIRECT_ONE_SCATTERER_ONE_H_USING           \
      using core<float_type>::structure_factor;                                \
      using core<float_type>::grad_site;                                       \
      using core<float_type>::grad_u_star;                                     \
      using core<float_type>::grad_u_iso;                                      \
      using core<float_type>::grad_occ;                                        \
      using base_t::hr_ht;                                                     \
      using base_t::d_star_sq;


    /** Key steps of the linearisation of the structure factor
        of one scatterer for a given Miller index in any space group.
     */
    template <typename FloatType, class ExpI2PiFunctor, bool compute_grad>
    struct in_generic_space_group : base<FloatType,
                                        ExpI2PiFunctor,
                                        in_generic_space_group<FloatType,
                                                               ExpI2PiFunctor,
                                                               compute_grad> >

    {
      typedef base<FloatType,
                   ExpI2PiFunctor,
                   in_generic_space_group<FloatType,
                                          ExpI2PiFunctor,
                                          compute_grad> >
              base_t;

      SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS;
      SMTBX_STRUCTURE_FACTORS_DIRECT_ONE_SCATTERER_ONE_H_USING;

      /** \copydoc one_scatterer_one_h_linearisation::base */
      in_generic_space_group(sgtbx::space_group const &space_group,
                             miller::index<> const &h,
                             float_type d_star_sq,
                             ExpI2PiFunctor const &exp_i_2pi_functor)

        : base_t(space_group, h, d_star_sq, exp_i_2pi_functor)
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
      void compute_anisotropic_part(scatterer_type const &scatterer)
      {
        using namespace adptbx;
        using namespace scitbx::constants;
        scitbx::math::imaginary_unit_t i;

        // Compute S'
        for (int k=0; k < hr_ht.groups.size(); ++k) {
          hr_ht_group<float_type> const &g = hr_ht.groups[k];
          float_type hrx = g.hr * scatterer.site;
          complex_type f = this->exp_i_2pi(hrx + g.ht);
          if (scatterer.flags.use_u_aniso()) {
            float_type dw = debye_waller_factor_u_star(g.hr, scatterer.u_star);
            f *= dw;
            if (!compute_grad) continue;
            if (scatterer.flags.grad_u_aniso()) {
              scitbx::sym_mat3<float_type> log_grad_u_star
                = debye_waller_factor_u_star_gradient_coefficients<float_type>(
                                                                          g.hr);
              complex_type grad_u_star_factor = -two_pi_sq * f;
              for (int j=0; j<6; ++j) {
                grad_u_star[j] += grad_u_star_factor * log_grad_u_star[j];
              }
            }
          }
          structure_factor += f;
          if (compute_grad && scatterer.flags.grad_site()) {
            complex_type grad_site_factor = i * two_pi * f;
            for (int j=0; j<3; ++j) {
              float_type h_j = g.hr[j];
              grad_site[j] += grad_site_factor * h_j;
            }
          }
        }

        if (hr_ht.is_centric) {
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
          }
        }

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
                                      FormFactorType const &ff)
      {
        using namespace adptbx;
        using namespace scitbx::constants;

        // factor from centring translation group
        float_type f_iso = hr_ht.ltr_factor;

        // factor from special position
        float_type w = scatterer.weight_without_occupancy();
        if (w != 1) f_iso *= w;

        // isotropic Debye-Waller
        if (scatterer.flags.use_u_iso()) {
          float_type dw_iso = debye_waller_factor_u_iso(d_star_sq/4,
                                                        scatterer.u_iso);
          f_iso *= dw_iso;
        }

        // occupancy
        if (compute_grad && scatterer.flags.grad_occupancy()) {
          grad_occ = ff * f_iso * structure_factor;
        }
        f_iso *= scatterer.occupancy;

        // Scattering factor
        FormFactorType ff_iso = ff * f_iso;

        // Finish
        structure_factor *= ff_iso;

        if (!compute_grad) return;

        if (scatterer.flags.use_u_iso() && scatterer.flags.grad_u_iso()) {
          grad_u_iso = -two_pi_sq * d_star_sq * structure_factor;
        }
        if (scatterer.flags.grad_site()) {
          for (int j=0; j<3; ++j) grad_site[j] *= ff_iso;
        }
        if (scatterer.flags.grad_u_aniso()) {
          for (int j=0; j<6; ++j) grad_u_star[j] *= ff_iso;
        }
      }
    };


    /** Key steps of the linearisation of the structure factor
        of one scatterer for a given Miller index
        in an origin centric space group.
     */
    template <typename FloatType, class ExpI2PiFunctor, bool compute_grad>
    struct in_origin_centric_space_group
      : base<FloatType,
             ExpI2PiFunctor,
             in_origin_centric_space_group<FloatType,
                                           ExpI2PiFunctor,
                                           compute_grad> >
    {
      typedef
        base<FloatType,
             ExpI2PiFunctor,
             in_origin_centric_space_group<FloatType,
                                           ExpI2PiFunctor,
                                           compute_grad> >
        base_t;

      SMTBX_STRUCTURE_FACTORS_DIRECT_TYPEDEFS;
      SMTBX_STRUCTURE_FACTORS_DIRECT_ONE_SCATTERER_ONE_H_USING;

      /** \copydoc one_scatterer_one_h_linearisation::base */
      in_origin_centric_space_group(sgtbx::space_group const &space_group,
                                    miller::index<> const &h,
                                    float_type d_star_sq,
                                    ExpI2PiFunctor const &exp_i_2pi_functor)

        : base_t(space_group, h, d_star_sq, exp_i_2pi_functor)
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
      void compute_anisotropic_part(scatterer_type const &scatterer)
      {
        using namespace adptbx;
        using namespace scitbx::constants;

        /* Compute the real part of S' and its gradients
         We only ever touch the real part of structure_factors and grad_xxx
         members. So this code should be as efficient as working with
         real numbers.
         */
        for (int k=0; k < hr_ht.groups.size(); ++k) {
          hr_ht_group<float_type> const &g = hr_ht.groups[k];
          float_type hrx = g.hr * scatterer.site;
          complex_type f = this->exp_i_2pi(hrx + g.ht);
          float_type fa = f.real(), fb = f.imag();
          if (scatterer.flags.use_u_aniso()) {
            float_type dw = debye_waller_factor_u_star(g.hr, scatterer.u_star);
            fa *= dw;
            if (!compute_grad) continue;
            if (scatterer.flags.grad_site()) fb *= dw;
            if (scatterer.flags.grad_u_aniso()) {
              scitbx::sym_mat3<float_type> log_grad_u_star
              = debye_waller_factor_u_star_gradient_coefficients<float_type>(
                                                                          g.hr);
              float_type grad_u_star_factor = -two_pi_sq * fa;
              for (int j=0; j<6; ++j) {
                grad_u_star[j] += grad_u_star_factor * log_grad_u_star[j];
              }
            }
          }
          structure_factor += fa;
          float_type grad_site_factor = -two_pi * fb;
          if (compute_grad && scatterer.flags.grad_site()) {
            for (int j=0; j<3; ++j) grad_site[j] += grad_site_factor * g.hr[j];
          }
        }

        /* We should now multiply S' and its gradients by 2
         to get S and its gradients. Postponed to next stage for efficiency.
         */
      }

      /** \copydoc in_generic_space_group::multiply_by_isotropic_part */
      template <typename FormFactorType>
      void multiply_by_isotropic_part(scatterer_type const &scatterer,
                                      FormFactorType const &ff)
      {
        using namespace adptbx;
        using namespace scitbx::constants;

        // factor from inversion
        float_type f_iso = 2.;

        // factor from centring translation group
        f_iso *= hr_ht.ltr_factor;

        // factor from special position
        float_type w = scatterer.weight_without_occupancy();
        if (w != 1) f_iso *= w;

        // isotropic Debye-Waller
        if (scatterer.flags.use_u_iso()) {
          float_type dw_iso = debye_waller_factor_u_iso(d_star_sq/4,
                                                        scatterer.u_iso);
          f_iso *= dw_iso;
        }

        // occupancy
        if (compute_grad && scatterer.flags.grad_occupancy()) {
          grad_occ = ff * f_iso * structure_factor.real();
        }
        f_iso *= scatterer.occupancy;

        // Scattering factor
        FormFactorType ff_iso = ff * f_iso;

        // Finish
        structure_factor = ff_iso * structure_factor.real();

        if (!compute_grad) return;

        if (scatterer.flags.use_u_iso() && scatterer.flags.grad_u_iso()) {
          grad_u_iso = -two_pi_sq * d_star_sq * structure_factor;
        }
        if (scatterer.flags.grad_site()) {
          for (int j=0; j<3; ++j) {
            grad_site[j] = ff_iso * grad_site[j].real();
          }
        }
        if (scatterer.flags.grad_u_aniso()) {
          for (int j=0; j<6; ++j) {
            grad_u_star[j] = ff_iso * grad_u_star[j].real();
          }
        }
      }
    };

  } // namespace one_scatterer_one_h_linearisation


  namespace one_h_linearisation {

    /** @brief Linearisation of \f$F_c(h)\f$ and of a derived observable,
               for any miller index h,
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
        and provide a member exp_i_2pi of type ExpI2PiFunctor to compute
        \f$\exp i2\pi x\f$.
     */
    template <typename FloatType, bool compute_grad,
              template<typename, bool> class ObservableType,
              template<typename> class ExpI2PiFunctor,
              class Heir>
    class base
    {
    public:
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;
      typedef ObservableType<float_type, compute_grad> observable_type;
      typedef ExpI2PiFunctor<float_type> exp_i_2pi_functor;

    protected:
      cctbx::xray::scatterer_grad_flags_counts grad_flags_counts;
      uctbx::unit_cell const &unit_cell;
      sgtbx::space_group const &space_group;
      bool origin_centric_case;
      af::ref_owning_shared< xray::scatterer<float_type> > scatterers;
      xray::scattering_type_registry const &scattering_type_registry;
      af::shared<std::size_t> scattering_type_indices;

      complex_type *grad_f_calc_cursor;

    public:
      complex_type f_calc;
      af::ref_owning_shared<complex_type> grad_f_calc;
      float_type observable;
      af::ref_owning_shared<float_type> grad_observable;

    public:
      /** @brief The linearisation of \f$F_c\f$ for the given structure
                 is to be computed.

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
          grad_observable(grad_flags_counts.n_parameters(),
                          af::init_functor_null<float_type>())
      {}

      /// Compute the linearisation
      void compute(
        miller::index<> const &h,
        boost::optional<complex_type> const &f_mask=boost::optional<complex_type>())
      {
        float_type d_star_sq = unit_cell.d_star_sq(h);
        af::shared<float_type> elastic_form_factors
          = scattering_type_registry.unique_form_factors_at_d_star_sq(d_star_sq);

        Heir &heir = static_cast<Heir &>(*this);

        typedef one_scatterer_one_h_linearisation::in_generic_space_group<
                  float_type, exp_i_2pi_functor, compute_grad>
                generic_linearisation_t;

        typedef one_scatterer_one_h_linearisation::in_origin_centric_space_group<
                  float_type, exp_i_2pi_functor, compute_grad>
                origin_centric_linearisation_t;

        if (!origin_centric_case) {
           generic_linearisation_t lin_for_h(space_group, h, d_star_sq,
                                             heir.exp_i_2pi);
          compute(elastic_form_factors.ref(), lin_for_h, f_mask);
        }
        else {
          origin_centric_linearisation_t lin_for_h(space_group, h, d_star_sq,
                                                   heir.exp_i_2pi);
          compute(elastic_form_factors.ref(), lin_for_h, f_mask);
        }
      }

    private:
      template <class LinearisationForMillerIndex>
      void compute(af::const_ref<float_type> const &elastic_form_factors,
                   LinearisationForMillerIndex &l,
                   boost::optional<complex_type> const &f_mask)
      {
        f_calc = 0;
        grad_f_calc_cursor = grad_f_calc.begin();

        for (int j=0; j < scatterers.size(); ++j) {
          xray::scatterer<> const &sc = scatterers[j];
          float_type f0 = elastic_form_factors[ scattering_type_indices[j] ];
          l.compute(sc, f0);

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
        }
        if (f_mask) f_calc += *f_mask;
        observable_type::compute(origin_centric_case,
                                 f_calc, grad_f_calc,
                                 observable, grad_observable);
      }
    };


    template <typename FloatType, bool compute_grad,
              template<typename, bool> class ObservableType,
              template<typename> class ExpI2PiFunctor>
    class custom_trigonometry
      : public base<FloatType, compute_grad, ObservableType, ExpI2PiFunctor,
                    custom_trigonometry<
                      FloatType, compute_grad, ObservableType, ExpI2PiFunctor> >
    {
    public:
      typedef base<FloatType, compute_grad, ObservableType, ExpI2PiFunctor,
                   custom_trigonometry<
                     FloatType, compute_grad, ObservableType, ExpI2PiFunctor> >
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
    };


    template <typename FloatType, bool compute_grad,
              template<typename, bool> class ObservableType>
    class std_trigonometry
      : public base<FloatType, compute_grad, ObservableType,
                    cctbx::math::cos_sin_exact,
                    std_trigonometry<
                      FloatType, compute_grad, ObservableType>
                    >
    {
    public:
      typedef base<FloatType, compute_grad, ObservableType,
                   cctbx::math::cos_sin_exact,
                   std_trigonometry<
                     FloatType, compute_grad, ObservableType>
                   >
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
    };


    template <typename FloatType, bool compute_grad>
    struct modulus_squared
    {
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;

      static
      void compute(bool origin_centric_case,
                   complex_type f_calc,
                   af::const_ref<complex_type> const &grad_f_calc,
                   float_type &observable,
                   af::ref<float_type> const &grad_observable)
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
    };


    template <typename FloatType, bool compute_grad>
    struct modulus
    {
      typedef FloatType float_type;
      typedef std::complex<float_type> complex_type;

      static
      void compute(bool origin_centric_case,
                   complex_type f_calc,
                   af::const_ref<complex_type> const &grad_f_calc,
                   float_type &observable,
                   af::ref<float_type> const &grad_observable)
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
  } // namespace one_h_linearisation


}}} // smtbx::structure_factors::direct

#endif // GUARD
