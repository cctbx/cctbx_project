#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_H

#include <cctbx/xray/twin_component.h>
#include <cctbx/xray/dispersion_radial.h>
#include <scitbx/array_family/ref_reductions.h>
#include <boost/shared_ptr.hpp>

namespace smtbx {
  namespace refinement {
    namespace least_squares {
      using namespace cctbx;
      using namespace cctbx::xray;
      /* Need inheritance to achive more flexibility */
      template <typename FloatType>
      struct f_calc_function_base {
        typedef std::complex<FloatType> complex_t;
        virtual ~f_calc_function_base() {}

        virtual void compute(
          miller::index<> const& h,
          boost::optional<complex_t > const& f_mask = boost::none,
          twin_fraction<FloatType> const* fraction = 0,
          bool compute_grad = true) = 0;

        void compute(
          miller::index<> const& h,
          twin_fraction<FloatType> const& fraction,
          complex_t const& f_mask,
          bool compute_grad = true)
        {
          compute(h, f_mask, &fraction, compute_grad);
        }

        /// Evaluate the structure factors
        void evaluate(miller::index<> const& h)
        {
          compute(h, boost::none, 0, false);
        }

        /// Evaluate the structure factors
        void evaluate(miller::index<> const& h,
          complex_t const& f_mask)
        {
          compute(h, f_mask, 0, false);
        }

        /// Linearise the structure factors
        void linearise(miller::index<> const& h)
        {
          compute(h, boost::none, 0, true);
        }

        /// Linearise the structure factors
        void linearise(miller::index<> const& h,
          complex_t const& f_mask)
        {
          compute(h, f_mask, 0, true);
        }

        virtual boost::shared_ptr<f_calc_function_base> fork() const = 0;

        virtual FloatType get_observable() const = 0;
        virtual complex_t get_f_calc() const = 0;
        virtual af::const_ref<complex_t > get_grad_f_calc() const = 0;
        virtual af::const_ref<FloatType> get_grad_observable() const = 0;
        /* returns true if grads are for all and not independent only params */
        virtual bool raw_gradients() const { return true; }
        /* The radial correction of f' and f'', or null if there is none.

        It is owned by the structure factor functor rather than sitting beside
        the normal equations the way fc_correction does, because its gradients
        are accumulated during the computation of Fc. Whoever reads them has to
        read them from the very functor that computed them, which under
        threading is a fork of the original -- hence going through here rather
        than holding a pointer of one's own.
        */
        virtual cctbx::xray::dispersion_radial_correction<FloatType> const*
        get_dispersion_correction() const
        {
          return 0;
        }
      };

      /** @brief Put the radial f'/f'' correction's gradients in their slots.

      They sit at the tail of the gradient vector next to BASF, EXTI and the
      rest. The sparse Jacobian product which fills the rest of the vector
      leaves them at zero -- they belong to no scatterer -- so this assigns
      rather than adds. The correction has to be the one belonging to this very
      f_calc_function, which under threading is a fork of the original.
      */
      template <typename FloatType>
      void write_dispersion_gradients(
        f_calc_function_base<FloatType> const &f_calc_function,
        af::shared<FloatType> &gradients)
      {
        cctbx::xray::dispersion_radial_correction<FloatType> const *dc =
          f_calc_function.get_dispersion_correction();
        if (dc == 0 || !dc->grad) {
          return;
        }
        af::const_ref<FloatType> dg = dc->get_gradients();
        SMTBX_ASSERT(dc->grad_index >= 0
          && dc->grad_index + dg.size() <= gradients.size());
        for (std::size_t gi = 0; gi < dg.size(); gi++) {
          gradients[dc->grad_index + gi] = dg[gi];
        }
      }

      /* A thin wrapper around the concrete implementation */
      template <typename FloatType,
        class OneMillerIndexFcalc>
      struct f_calc_function_default : public f_calc_function_base<FloatType> {
        typedef std::complex<FloatType> complex_t;
        typedef f_calc_function_base<FloatType> f_calc_function_base_t;

        f_calc_function_default(boost::shared_ptr<OneMillerIndexFcalc> f_calc_function)
          : f_calc_function(f_calc_function)
        {}

        virtual void compute(
          miller::index<> const& h,
          boost::optional<complex_t > const& f_mask = boost::none,
          twin_fraction<FloatType> const* fraction = 0,
          bool compute_grad = true)
        {
          f_calc_function->compute(h, f_mask, compute_grad);
        }

        virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
          return boost::shared_ptr<f_calc_function_base_t>(
            new f_calc_function_default(f_calc_function->fork()));
        }

        virtual FloatType get_observable() const {
          return f_calc_function->get_observable();
        }
        virtual complex_t get_f_calc() const {
          return f_calc_function->f_calc;
        }
        virtual af::const_ref<complex_t> get_grad_f_calc() const {
          return f_calc_function->grad_f_calc.const_ref();
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          return f_calc_function->get_grad_observable().const_ref();
        }
        virtual cctbx::xray::dispersion_radial_correction<FloatType> const*
        get_dispersion_correction() const
        {
          return f_calc_function->disp_cr.get();
        }

        boost::shared_ptr<OneMillerIndexFcalc> f_calc_function;
      };


      /*  A thin wrapper around concrete implementation to enable caching of
      the results for symmetry related indices.
       */
      template <typename FloatType>
      struct f_calc_function_with_cache : public f_calc_function_base<FloatType>
      {
        typedef std::complex<FloatType> complex_t;
        typedef f_calc_function_base<FloatType> f_calc_function_base_t;
        struct f_calc_function_result {
          f_calc_function_result(
            FloatType const& observable,
            complex_t const& f_calc,
            af::const_ref<complex_t> const& grad_f_calc,
            af::const_ref<FloatType> const& grad_observable)
            :
            observable(observable),
            f_calc(f_calc),
            grad_f_calc(grad_f_calc.begin(), grad_f_calc.end()),
            grad_observable(grad_observable.begin(), grad_observable.end())
          {}

          f_calc_function_result(
            FloatType const& observable,
            complex_t const& f_calc)
            :
            observable(observable),
            f_calc(f_calc),
            grad_observable()
          {}

          FloatType const observable;
          complex_t const f_calc;
          af::shared<complex_t> grad_f_calc;
          af::shared<FloatType> grad_observable;
        };

        f_calc_function_with_cache(
          boost::shared_ptr<f_calc_function_base_t> f_calc_function,
            bool use_cache = false)
          : f_calc_function(f_calc_function),
          use_cache(use_cache),
          length_sq(0)
        {
          /* A cache hit skips compute() altogether, which would leave the
             radial correction's gradient accumulator holding whatever the last
             reflection that missed put there, while the observable came from
             the cache. Caching the correction's gradients alongside the rest
             would fix it; refusing the combination is what is warranted until
             something asks for it, this wrapper being unreachable from Python.
           */
          SMTBX_ASSERT(!use_cache
            || f_calc_function->get_dispersion_correction() == 0);
        }

        virtual void compute(
          miller::index<> const& h,
          boost::optional<complex_t > const& f_mask = boost::none,
          twin_fraction<FloatType> const* fraction = 0,
          bool compute_grad = true)
        {
          if (!use_cache) {
            f_calc_function->compute(h, f_mask, fraction, compute_grad);
            observable = f_calc_function->get_observable();
            grad_f_calc = f_calc_function->get_grad_f_calc();
            grad_observable = f_calc_function->get_grad_observable();
            f_calc = f_calc_function->get_f_calc();
          }
          else {
            FloatType h_length_sq = h.length_sq();
            if (h_length_sq != length_sq) {
              cache.clear();
              length_sq = h_length_sq;
            }
            typename cache_t::iterator iter = cache.find(h);
            if (iter == cache.end()) {
              f_calc_function->compute(h, f_mask, fraction, compute_grad);
              observable = f_calc_function->get_observable();
              grad_f_calc = f_calc_function->get_grad_f_calc();
              grad_observable = f_calc_function->get_grad_observable();
              f_calc = f_calc_function->get_f_calc();
              cache.insert(
                std::pair<miller::index<>, f_calc_function_result>(
                  h, f_calc_function_result(
                    observable,
                    f_calc,
                    grad_f_calc,
                    grad_observable)));
            }
            else {
              observable = iter->second.observable;
              f_calc = iter->second.f_calc;
              grad_f_calc = iter->second.grad_f_calc.const_ref();
              grad_observable = iter->second.grad_observable.const_ref();
            }
          }
        }

        void compute(miller::index<> const& h,
          bool compute_grad = true)
        {
          compute(h, /*f_mask=*/ boost::none, compute_grad);
        }

        virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
          return boost::shared_ptr<f_calc_function_base_t>(
            new f_calc_function_with_cache(f_calc_function->fork(),
              use_cache));
        }

        virtual FloatType get_observable() const {
          return observable;
        }
        virtual complex_t get_f_calc() const {
          return f_calc;
        }
        virtual af::const_ref<complex_t> get_grad_f_calc() const {
          return grad_f_calc;
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          return grad_observable;
        }
        virtual cctbx::xray::dispersion_radial_correction<FloatType> const*
        get_dispersion_correction() const
        {
          // safe only because the constructor refuses use_cache with one
          return f_calc_function->get_dispersion_correction();
        }

        typedef std::map<miller::index<>, f_calc_function_result> cache_t;

        boost::shared_ptr<f_calc_function_base_t> f_calc_function;
        FloatType observable;
        af::const_ref<complex_t> grad_f_calc;
        af::const_ref<FloatType> grad_observable;
        complex_t f_calc;
        bool use_cache;
        FloatType length_sq;
        cache_t cache;
      };

    }
  }
}


#endif // GUARD
