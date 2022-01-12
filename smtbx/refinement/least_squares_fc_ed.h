#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_H


#include <smtbx/refinement/least_squares_fc.h>

namespace smtbx {
  namespace refinement {
    namespace least_squares {

      template <typename FloatType>
      struct f_calc_function_ed : public f_calc_function_base<FloatType>
      {
        typedef f_calc_function_base<FloatType> f_calc_function_base_t;

        f_calc_function_ed(boost::shared_ptr<f_calc_function_base<FloatType> > f_calc_function)
          : f_calc_function(f_calc_function)
        {}

        virtual void compute(
          miller::index<> const& h,
          boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
          bool compute_grad = true)
        {
          f_calc_function->compute(h, f_mask, compute_grad);
        }

        virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
          return boost::shared_ptr<f_calc_function_base_t>(
            new f_calc_function_ed(f_calc_function->fork()));
        }

        virtual FloatType get_observable() const {
          return f_calc_function->get_observable();
        }
        virtual std::complex<FloatType> get_f_calc() const {
          return f_calc_function->get_f_calc();
        }
        virtual af::const_ref<FloatType> get_grad_observable() const {
          return f_calc_function->get_grad_observable();
        }

        boost::shared_ptr<f_calc_function_base_t> f_calc_function;
      };


    }
  }
}


#endif // GUARD
