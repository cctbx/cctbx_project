#pragma once
#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_TWINNING_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_TWINNING_H

/// Crystallographic least-squares

#include <scitbx/sparse/matrix.h>

#include <cctbx/xray/observations.h>

#include <smtbx/error.h>
#include <smtbx/refinement/least_squares_fc.h>


namespace smtbx {
  namespace refinement {
    namespace least_squares {
      template <typename FloatType>
      struct MaskData {
        typedef std::complex<FloatType> complex_type;

        MaskData(af::const_ref<complex_type> const& f_mask)
          : f_mask(f_mask)
        {}

        MaskData(cctbx::xray::observations<FloatType> const& reflections,
          sgtbx::space_group const& space_group,
          bool anomalous_flag,
          af::const_ref<miller::index<> > const& indices,
          af::const_ref<complex_type> const& f_mask)
          : f_mask(f_mask)
        {
          mi_lookup = miller::lookup_utils::lookup_tensor<FloatType>(
            indices, space_group, anomalous_flag);
        }

        complex_type find(miller::index<> const& h) const {
          long index = mi_lookup.find_hkl(h);
          SMTBX_ASSERT(index >= 0 && index < f_mask.size())(h.as_string());
          return f_mask[index];
        }

        int size() const {
          return f_mask.size();
        }

        const complex_type& get(int i) const {
          return f_mask[i];
        }

      private:
        af::const_ref<complex_type> f_mask;
        miller::lookup_utils::lookup_tensor<FloatType> mi_lookup;
      };

      template <typename FloatType>
      class twinning_processor {
      public:
        typedef typename cctbx::xray::observations<FloatType>::iterator itr_t;
        typedef typename cctbx::xray::twin_fraction<FloatType> twf_t;
        typedef typename cctbx::xray::observations<FloatType>::index_twin_component twc_t;
        typedef std::complex<FloatType> complex_type;
        
        twinning_processor(
          cctbx::xray::observations<FloatType> const& reflections,
          MaskData<FloatType> const& f_mask_data,
          bool compute_grad,
          scitbx::sparse::matrix<FloatType> const&
            jacobian_transpose_matching_grad_fc)
          : reflections(reflections),
          f_mask_data(f_mask_data),
          compute_grad(compute_grad),
          jacobian_transpose_matching_grad_fc(jacobian_transpose_matching_grad_fc)
        {}

        FloatType process(int i_h,
          f_calc_function_base<FloatType>& f_calc_function,
          af::shared<FloatType>& gradients) const
        {
          FloatType obs = f_calc_function.get_observable();
          // this is the same for all components - not much useful with real twinning
          const twin_fraction<FloatType>* fraction = reflections.fraction(i_h);
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
              if (f_mask_data.size() > 0) {
                f_mask = f_mask_data.find(twc.h);
              }
              f_calc_function.compute(twc.h, f_mask, fraction, compute_grad);
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
        MaskData<FloatType> const& f_mask_data;
        bool compute_grad;
        scitbx::sparse::matrix<FloatType> const&
          jacobian_transpose_matching_grad_fc;
      };
    }
  }
}


#endif // GUARD
