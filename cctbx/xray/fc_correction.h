#pragma once

#include <cctbx/miller.h>
#include <scitbx/array_family/ref.h>
#include <boost/shared_ptr.hpp>

namespace cctbx {
  namespace xray {

    // Fc correction base, EXTI/SWAT
    template <typename FloatType>
    struct fc_correction {
      typedef fc_correction<FloatType> fc_correction_t;

      fc_correction(bool grad=false)
        : grad(grad)
      {}
      virtual ~fc_correction() {}
      // return multiplier for Fc_sq
      virtual FloatType compute(miller::index<> const& h,
        FloatType fc_sq,
        bool compute_gradient) const = 0;
      virtual FloatType get_grad_Fc_multiplier() const = 0;
      // starting grad index
      virtual int get_grad_index() const = 0;
      virtual af::const_ref<FloatType> get_gradients() const = 0;
      virtual boost::shared_ptr<fc_correction_t> fork() const = 0;
      bool grad;
    };

    // dummy fc correcction implementation
    template <typename FloatType>
    struct dummy_fc_correction : public fc_correction<FloatType> {
      typedef fc_correction<FloatType> fc_correction_t;

      dummy_fc_correction() {}
      virtual FloatType compute(miller::index<> const& h,
        FloatType fc_sq,
        bool compute_gradient) const
      {
        return 1;
      }

      virtual FloatType get_grad_Fc_multiplier() const { return 1; }
      virtual int get_grad_index() const { return -1; }
      virtual af::const_ref<FloatType> get_gradients() const {
        CCTBX_NOT_IMPLEMENTED();
        af::shared<FloatType> rv;
        return rv.const_ref();
      }
      virtual boost::shared_ptr<fc_correction_t> fork() const {
        return boost::shared_ptr<fc_correction_t>(
          new dummy_fc_correction());
      }
    };
  }
} // namespace cctbx::xray
