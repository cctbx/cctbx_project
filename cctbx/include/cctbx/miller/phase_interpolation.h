/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Dec: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_PHASE_INTERPOLATION_H
#define CCTBX_MILLER_PHASE_INTERPOLATION_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/constants.h>
#include <cctbx/error.h>
#include <complex>

namespace cctbx { namespace miller {

  inline
  double
  phase_error(double phase_difference, bool deg=false)
  {
    double pi_sc;
    if (deg) pi_sc = 180;
    else     pi_sc = scitbx::constants::pi;
    double result = std::fmod(phase_difference, 2 * pi_sc);
    if      (result < -pi_sc) result += 2 * pi_sc;
    else if (result >  pi_sc) result -= 2 * pi_sc;
    return result;
  }

  inline
  double
  phase_interpolation(
    bool centric_flag,
    double reference_amplitude,
    double amplitude_1,
    double phase_1,
    double amplitude_2,
    double phase_2,
    bool deg=false,
    double epsilon=1.e-10)
  {
    double delta_1 = scitbx::fn::absolute(amplitude_1 - reference_amplitude);
    double delta_2 = scitbx::fn::absolute(amplitude_2 - reference_amplitude);
    if (centric_flag) {
      if (delta_2 < delta_1) return phase_2;
      return phase_1;
    }
    double sum_delta = delta_1 + delta_2;
    phase_2 = phase_1 + phase_error(phase_2 - phase_1, deg);
    if (sum_delta < epsilon) return (phase_1 + phase_2) / 2;
    return phase_1 + (phase_2 - phase_1) * delta_1 / sum_delta;
  }

  template <typename FloatType>
  af::shared<FloatType>
  phase_interpolation(
    af::const_ref<bool> const& centric_flags,
    af::const_ref<FloatType> const& reference_amplitudes,
    af::const_ref<std::complex<FloatType> > const& f_1,
    af::const_ref<std::complex<FloatType> > const& f_2,
    bool deg=false,
    double epsilon=1.e-10)
  {
    CCTBX_ASSERT(reference_amplitudes.size() == centric_flags.size());
    CCTBX_ASSERT(f_1.size() == centric_flags.size());
    CCTBX_ASSERT(f_2.size() == centric_flags.size());
    af::shared<FloatType> result(af::reserve(centric_flags.size()));
    for(std::size_t i=0;i<centric_flags.size();i++) {
      double phase = phase_interpolation(
        centric_flags[i],
        reference_amplitudes[i],
        std::abs(f_1[i]),
        std::arg(f_1[i]),
        std::abs(f_2[i]),
        std::arg(f_2[i]),
        false,
        epsilon);
      if (deg) phase /= scitbx::constants::pi_180;
      result.push_back(phase);
    }
    return result;
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_PHASE_INTERPOLATION_H
