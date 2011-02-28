#ifndef CCTBX_XRAY_TARGET_COMMON_RESULTS_H
#define CCTBX_XRAY_TARGET_COMMON_RESULTS_H

#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <tbxx/error_utils.hpp>
#include <boost/optional.hpp>
#include <complex>

namespace cctbx { namespace xray {

//! Structure-factor target functions.
namespace targets {

  struct common_results
  {
    protected:
      af::shared<double> target_per_reflection_;
      boost::optional<double> target_work_;
      boost::optional<double> target_test_;
      af::shared<std::complex<double> > gradients_work_;
      af::shared<scitbx::vec3<double> > hessians_work_;
      public:

    common_results() {}

    common_results(
      std::size_t n_refl)
    :
      target_per_reflection_(n_refl, af::init_functor_null<double>())
    {}

    common_results(
      af::shared<double> const& target_per_reflection,
      boost::optional<double> const& target_work,
      boost::optional<double> const& target_test,
      af::shared<std::complex<double> > const& gradients_work)
    :
      target_per_reflection_(target_per_reflection),
      target_work_(target_work),
      target_test_(target_test),
      gradients_work_(gradients_work)
    {
      if (   target_per_reflection.size() != 0
          && gradients_work.size() != 0) {
        TBXX_ASSERT(gradients_work.size() <= target_per_reflection.size());
      }
    }

    common_results(
      af::shared<double> const& target_per_reflection,
      boost::optional<double> const& target_work,
      boost::optional<double> const& target_test,
      af::shared<std::complex<double> > const& gradients_work,
      af::shared<scitbx::vec3<double> > const& hessians_work)
    :
      target_per_reflection_(target_per_reflection),
      target_work_(target_work),
      target_test_(target_test),
      gradients_work_(gradients_work),
      hessians_work_(hessians_work)
    {
      if (target_per_reflection.size() != 0) {
        if (gradients_work.size() != 0) {
          TBXX_ASSERT(
            gradients_work.size() <= target_per_reflection.size());
        }
        if (hessians_work.size() != 0) {
          TBXX_ASSERT(
            hessians_work.size() <= target_per_reflection.size());
        }
      }
      if (   gradients_work.size() != 0
          && hessians_work.size() != 0) {
        TBXX_ASSERT(hessians_work.size() == gradients_work.size());
      }
    }

    af::shared<double> const&
    target_per_reflection() const { return target_per_reflection_; }

    boost::optional<double>
    target_work() const { return target_work_; }

    boost::optional<double>
    target_test() const { return target_test_; }

    //! da, db
    af::shared<std::complex<double> > const&
    gradients_work() { return gradients_work_; }

    //! daa, dbb, dab
    af::shared<scitbx::vec3<double> > const&
    hessians_work() { return hessians_work_; }
  };

}}} // namespace cctbx::xray::targets

#endif // GUARD
