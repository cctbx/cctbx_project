#ifndef SCITBX_MATH_PHASE_ERROR_H
#define SCITBX_MATH_PHASE_ERROR_H

#include <scitbx/constants.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/misc_functions.h>

namespace scitbx { namespace math {

  template <typename FloatType>
  FloatType
  signed_phase_error(
    FloatType const& phi1,
    FloatType const& phi2,
    bool deg=false)
  {
    FloatType pi_sc;
    if (deg) pi_sc = 180;
    else     pi_sc = constants::pi;
    FloatType e = std::fmod(phi2-phi1, 2*pi_sc);
    if      (e < -pi_sc) e += 2*pi_sc;
    else if (e >  pi_sc) e -= 2*pi_sc;
    return e;
  }

  template <typename FloatType>
  af::shared<FloatType>
  signed_phase_error(
    af::const_ref<FloatType> const& phi1,
    af::const_ref<FloatType> const& phi2,
    bool deg=false)
  {
    SCITBX_ASSERT(phi1.size() == phi2.size());
    af::shared<FloatType> result((af::reserve(phi1.size())));
    for(std::size_t i=0;i<phi1.size();i++) {
      result.push_back(signed_phase_error(phi1[i], phi2[i], deg));
    }
    return result;
  }

  template <typename FloatType>
  FloatType
  phase_error(
    FloatType const& phi1,
    FloatType const& phi2,
    bool deg=false)
  {
    return fn::absolute(signed_phase_error(phi1, phi2, deg));
  }

  template <typename FloatType>
  af::shared<FloatType>
  phase_error(
    af::const_ref<FloatType> const& phi1,
    af::const_ref<FloatType> const& phi2,
    bool deg=false)
  {
    SCITBX_ASSERT(phi1.size() == phi2.size());
    af::shared<FloatType> result((af::reserve(phi1.size())));
    for(std::size_t i=0;i<phi1.size();i++) {
      result.push_back(phase_error(phi1[i], phi2[i], deg));
    }
    return result;
  }

  template <typename FloatType>
  FloatType
  nearest_phase(
    FloatType const& reference,
    FloatType const& other,
    bool deg=false)
  {
    return reference + signed_phase_error(reference, other, deg);
  }

  template <typename FloatType>
  af::shared<FloatType>
  nearest_phase(
    af::const_ref<FloatType> const& reference,
    af::const_ref<FloatType> const& other,
    bool deg=false)
  {
    SCITBX_ASSERT(other.size() == reference.size());
    af::shared<FloatType> result((af::reserve(reference.size())));
    for(std::size_t i=0;i<reference.size();i++) {
      result.push_back(nearest_phase(reference[i], other[i], deg));
    }
    return result;
  }

}} // namespace scitbx::math

#endif // SCITBX_MATH_PHASE_ERROR_H
