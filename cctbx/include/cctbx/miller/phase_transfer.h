/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_PHASE_TRANSFER_H
#define CCTBX_MILLER_PHASE_TRANSFER_H

#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace miller {

  template <typename AmplitudeType,
            typename FloatType>
  af::shared<std::complex<FloatType> >
  phase_transfer(
    sgtbx::space_group const& space_group,
    af::const_ref<index<> > const& miller_indices,
    af::const_ref<AmplitudeType> const& amplitude_source,
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
        result.push_back(std::polar(
          std::abs(amplitude_source[i]),
          space_group.phase_restriction(miller_indices[i])
            .nearest_valid_phase(std::arg(p))));
      }
    }
    return result;
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_PHASE_TRANSFER_H
