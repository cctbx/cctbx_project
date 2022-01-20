#pragma once
#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_TWINNING_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_TWINNING_H

/// Crystallographic least-squares

#include <scitbx/sparse/matrix.h>

#include <cctbx/xray/observations.h>
#include <cctbx/miller/lookup_utils.h>

#include <smtbx/error.h>
#include <smtbx/refinement/least_squares_fc.h>


namespace smtbx {
  namespace refinement {
    namespace least_squares {
      template <typename FloatType>
      class twinning_processor {
      public:
        typedef typename cctbx::xray::observations<FloatType>::iterator itr_t;
        typedef typename cctbx::xray::twin_fraction<FloatType> twf_t;
        typedef typename cctbx::xray::observations<FloatType>::index_twin_component twc_t;
        typedef std::complex<FloatType> complex_type;
        typedef std::map<cctbx::miller::index<>, complex_type,
          cctbx::miller::fast_less_than<> > lookup_t;

        twinning_processor(
          cctbx::xray::observations<FloatType> const& reflections,
          af::const_ref<complex_type> const& f_mask,
          bool compute_grad,
          scitbx::sparse::matrix<FloatType> const&
            jacobian_transpose_matching_grad_fc)
          : reflections(reflections),
          compute_grad(compute_grad),
          jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc)
        {
          if (f_mask.size() > 0 && reflections.has_twin_components()) {
            scitbx::af::shared<miller::index<> > const& indices = reflections.indices();
            for (int i_h = 0, i=0; i_h < reflections.size(); i_h++, i++) {
              mi_lookup[indices[i_h]] = f_mask[i];
              itr_t itr = reflections.iterate(i_h);
              while (itr.has_next()) {
                twc_t twc = itr.next();
                i++;
                SMTBX_ASSERT(i < f_mask.size())(i)(f_mask.size());
                mi_lookup[twc.h] = f_mask[i];
              }
            }
          }
        }

        FloatType process(int i_h,
          f_calc_function_base<FloatType>& f_calc_function,
          af::shared<FloatType>& gradients) const
        {
          FloatType obs = f_calc_function.get_observable();
          if (reflections.has_twin_components()) {
            itr_t itr = reflections.iterate(i_h);
            FloatType measured_part = obs,
              identity_part = 0,
              obs_scale = reflections.scale(i_h);
            obs *= obs_scale;
            const twf_t* measured_fraction = reflections.fraction(i_h);
            if (compute_grad) {
              gradients *= obs_scale;
              if (measured_fraction == 0) {
                identity_part = measured_part;
              }
            }
            std::size_t twc_cnt = 0;
            while (itr.has_next()) {
              twc_t twc = itr.next();
              boost::optional<complex_type> f_mask = boost::none;
              if (mi_lookup.size()) {
                typename lookup_t::const_iterator l = mi_lookup.find(twc.h);
                SMTBX_ASSERT(l != mi_lookup.end())(twc.h.as_string());
                f_mask = l->second;
              }
              f_calc_function.compute(twc.h, f_mask, compute_grad);
              obs += twc.scale() * f_calc_function.get_observable();
              if (compute_grad) {
                if (f_calc_function.raw_gradients()) {
                  af::shared<FloatType> tmp_gradients =
                    jacobian_transpose_matching_grad_fc * f_calc_function.get_grad_observable();
                  gradients += twc.scale() * tmp_gradients;
                }
                else {
                  af::shared<FloatType> tmp_gradients(
                    f_calc_function.get_grad_observable().begin(),
                    f_calc_function.get_grad_observable().end());
                  gradients += twc.scale() * tmp_gradients;
                }
                if (twc.fraction != 0) {
                  if (twc.fraction->grad) {
                    SMTBX_ASSERT(!(twc.fraction->grad_index < 0 ||
                      twc.fraction->grad_index >= gradients.size()));
                    gradients[twc.fraction->grad_index] += f_calc_function.get_observable();
                  }
                }
                else {
                  identity_part += f_calc_function.get_observable();
                }
                twc_cnt++;
              }
            }
            if (compute_grad) {
              // consider multiple reflections with the 'prime' scale
              itr.reset();
              while (itr.has_next()) {
                twc_t twc = itr.next();
                if (twc.fraction != 0 && twc.fraction->grad) {
                  gradients[twc.fraction->grad_index] -= identity_part;
                }
              }
              if (twc_cnt != 0 && measured_fraction != 0 && measured_fraction->grad) {
                SMTBX_ASSERT(!(measured_fraction->grad_index < 0 ||
                  measured_fraction->grad_index >= gradients.size()));
                gradients[measured_fraction->grad_index] +=
                  measured_part - identity_part;
              }
            }
          }
          return obs;
        }
      private:
        cctbx::xray::observations<FloatType> const& reflections;
        bool compute_grad;
        scitbx::sparse::matrix<FloatType> const&
          jacobian_transpose_matching_grad_fc;
        lookup_t mi_lookup;
      };
    }
  }
}


#endif // GUARD
