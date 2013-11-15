#ifndef CCTBX_MILLER_PHASE_TRANSFER_H
#define CCTBX_MILLER_PHASE_TRANSFER_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/math/utils.h>

namespace cctbx { namespace miller {

  /// Transfer the phase of phase_source[i] onto amplitude_source[i]
  /**
    The sign of phase_source[i] is kept.

    The phase is actually the nearest phase compatible with space-group
    restriction for centric reflections.
  */
  template <typename FloatType>
  af::shared<std::complex<FloatType> >
  phase_transfer(
    sgtbx::space_group const& space_group,
    af::const_ref<index<> > const& miller_indices,
    af::const_ref<FloatType> const& amplitude_source,
    af::const_ref<std::complex<FloatType> > const& phase_source,
    FloatType const& epsilon=1.e-10)
  {
    CCTBX_ASSERT(amplitude_source.size() == miller_indices.size());
    CCTBX_ASSERT(phase_source.size() == miller_indices.size());
    af::shared<std::complex<FloatType> >
      result(af::reserve(miller_indices.size()));
    for(std::size_t i=0;i<miller_indices.size();i++) {
      std::complex<FloatType> p = phase_source[i];
      if (   scitbx::fn::absolute(p.real()) < epsilon
          && scitbx::fn::absolute(p.imag()) < epsilon) {
        result.push_back(0); // no phase information
      }
      else {
        index<>
#if !defined(CCTBX_SGTBX_PHASE_INFO_APPLE_LLVM2335_WORKAROUND)
        const&
#endif
        h = miller_indices[i];
        result.push_back(amplitude_source[i] * scitbx::math::unit_complex(
          space_group.phase_restriction(h).nearest_valid_phase(std::arg(p))));
      }
    }
    return result;
  }

  /// Another overload of phase_transfer passing the source phases directly
  template <typename FloatType>
  af::shared<std::complex<FloatType> >
  phase_transfer(
    sgtbx::space_group const& space_group,
    af::const_ref<index<> > const& miller_indices,
    af::const_ref<FloatType> const& amplitude_source,
    af::const_ref<FloatType> const& phase_source,
    bool deg=false)
  {
    CCTBX_ASSERT(amplitude_source.size() == miller_indices.size());
    CCTBX_ASSERT(phase_source.size() == miller_indices.size());
    af::shared<std::complex<FloatType> >
      result(af::reserve(miller_indices.size()));
    for(std::size_t i=0;i<miller_indices.size();i++) {
      FloatType p = phase_source[i];
      if (deg) p *= scitbx::constants::pi_180;
      result.push_back(amplitude_source[i] * scitbx::math::unit_complex(
        space_group.phase_restriction(
          miller_indices[i]).nearest_valid_phase(p)));
    }
    return result;
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_PHASE_TRANSFER_H
