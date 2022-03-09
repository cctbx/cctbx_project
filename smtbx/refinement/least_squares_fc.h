#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_H

#include <cctbx/xray/twin_component.h>
#include <scitbx/array_family/ref_reductions.h>
#include <boost/shared_ptr.hpp>

namespace smtbx {
  namespace refinement {
    namespace least_squares {
      using namespace cctbx::xray;
      /* Need inheritance to achive more flexibility */
      template <typename FloatType>
      struct f_calc_function_base {

        virtual ~f_calc_function_base() {}

        virtual void compute(
          miller::index<> const& h,
          boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
          twin_fraction<FloatType> const* fraction = 0,
          bool compute_grad = true) = 0;

        void compute(
          miller::index<> const& h,
          twin_fraction<FloatType> const& fraction,
          std::complex<FloatType> const& f_mask,
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
          std::complex<FloatType> const& f_mask)
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
          std::complex<FloatType> const& f_mask)
        {
          compute(h, f_mask, 0, true);
        }

        virtual boost::shared_ptr<f_calc_function_base> fork() const = 0;

        virtual FloatType get_observable() const = 0;
        virtual std::complex<FloatType> get_f_calc() const = 0;
        virtual af::const_ref<FloatType> get_grad_observable() const = 0;
        /* returns true if grads are for all and not independent only params */
        virtual bool raw_gradients() const { return true; }
      };

      /* A thin wrapper around the concrete implementation */
      template <typename FloatType,
        class OneMillerIndexFcalc>
      struct f_calc_function_default : public f_calc_function_base<FloatType> {
        typedef f_calc_function_base<FloatType> f_calc_function_base_t;

        f_calc_function_default(boost::shared_ptr<OneMillerIndexFcalc> f_calc_function)
          : f_calc_function(f_calc_function)
        {}

        virtual void compute(
          miller::index<> const& h,
          boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
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
          return f_calc_function->observable;
        }
        virtual std::complex<FloatType> get_f_calc() const {
          return f_calc_function->f_calc;
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          return f_calc_function->grad_observable.const_ref();
        }

        boost::shared_ptr<OneMillerIndexFcalc> f_calc_function;
      };


      /*  A thin wrapper around concrete implementation to enable caching of
      the results for symmetry related indices.
       */
      template <typename FloatType>
      struct f_calc_function_with_cache : public f_calc_function_base<FloatType>
      {
        typedef f_calc_function_base<FloatType> f_calc_function_base_t;
        struct f_calc_function_result {
          f_calc_function_result(
            FloatType const& observable,
            std::complex<FloatType> const& f_calc,
            af::const_ref<FloatType> const& grad_observable)
            :
            observable(observable),
            f_calc(f_calc),
            grad_observable(grad_observable.begin(), grad_observable.end())
          {}

          f_calc_function_result(
            FloatType const& observable,
            std::complex<FloatType> const& f_calc)
            :
            observable(observable),
            f_calc(f_calc),
            grad_observable()
          {}

          FloatType const observable;
          std::complex<FloatType> const f_calc;
          af::shared<FloatType> const grad_observable;
        };

        f_calc_function_with_cache(
          boost::shared_ptr<f_calc_function_base_t> f_calc_function, bool use_cache = false)
          : f_calc_function(f_calc_function),
          use_cache(use_cache),
          length_sq(0)
        {}

        virtual void compute(
          miller::index<> const& h,
          boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
          twin_fraction<FloatType> const* fraction = 0,
          bool compute_grad = true)
        {
          if (!use_cache) {
            f_calc_function->compute(h, f_mask, fraction, compute_grad);
            observable = f_calc_function->get_observable();
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
              grad_observable = f_calc_function->get_grad_observable();
              f_calc = f_calc_function->get_f_calc();
              cache.insert(
                std::pair<miller::index<>, f_calc_function_result>(
                  h, f_calc_function_result(
                    observable,
                    f_calc,
                    grad_observable)));
            }
            else {
              observable = iter->second.observable;
              f_calc = iter->second.f_calc;
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
        virtual std::complex<FloatType> get_f_calc() const {
          return f_calc;
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          return grad_observable;
        }

        typedef std::map<miller::index<>, f_calc_function_result> cache_t;

        boost::shared_ptr<f_calc_function_base_t> f_calc_function;
        FloatType observable;
        af::const_ref<FloatType> grad_observable;
        std::complex<FloatType> f_calc;
        bool use_cache;
        FloatType length_sq;
        cache_t cache;
      };

    }
  }
}


#endif // GUARD
