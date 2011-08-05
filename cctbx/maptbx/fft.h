#ifndef CCTBX_MAPTBX_FFT_H
#define CCTBX_MAPTBX_FFT_H

#include <cctbx/maptbx/structure_factors.h>
#include <cctbx/maptbx/copy.h>
#include <scitbx/fftpack/real_to_complex_3d.h>
#include <cctbx/error.h>
#include <scitbx/vec3.h>
#include <scitbx/constants.h>

#include <complex>

namespace cctbx { namespace maptbx {

  //! Computes FFT map from Fourier coefficients.
  /*! See also:
        scitbx::fftpack::real_to_complex_3d
        structure_factors::to_map
   */
  template <typename FloatType>
  af::versa<FloatType, af::c_grid_padded<3> >
  fft_to_real_map(
    sgtbx::space_group const& space_group,
    af::tiny<int, 3> const& n_real,
    af::const_ref<miller::index<> > const& miller_indices,
    af::const_ref<std::complex<FloatType> > const& data)
  {
    // determine padding
    scitbx::fftpack::real_to_complex_3d<FloatType> rfft(n_real);
    // allocate map and copy structure factors
    structure_factors::to_map<FloatType> to_map(
      space_group,
      /* anomalous_flag */ false,
      miller_indices,
      data,
      rfft.n_real(),
      af::c_grid_padded<3>(rfft.n_complex(), rfft.n_complex()),
      /* conjugate_flag */ true,
      /* treat_restricted */ true);
    // tedious but trivial conversion to suitable array reference
    af::ref<std::complex<FloatType>, af::c_grid<3> > complex_map(
      to_map.complex_map().begin(), af::c_grid<3>(rfft.n_complex()));
    // actual FFT
    rfft.backward(complex_map);
    // cast from complex to real
    return af::versa<FloatType, af::c_grid_padded<3> >(
      to_map.complex_map().handle(),
      af::c_grid_padded<3>(rfft.m_real(), rfft.n_real()));
  }

  //! Computes unpadded FFT map from Fourier coefficients.
  template <typename FloatType>
  af::versa<FloatType, af::c_grid<3> >
  fft_to_real_map_unpadded(
    sgtbx::space_group const& space_group,
    af::tiny<int, 3> const& n_real,
    af::const_ref<miller::index<> > const& miller_indices,
    af::const_ref<std::complex<FloatType> > const& data)
  {
    af::versa<FloatType, af::c_grid_padded<3> > map = fft_to_real_map(
      space_group, n_real, miller_indices, data);
    unpad_in_place(map.begin(), map.accessor().all(), map.accessor().focus());
    return af::versa<FloatType, af::c_grid<3> >(
      map, af::c_grid<3>(map.accessor().focus()));
  }

  template <typename FloatType>
  std::complex<FloatType>
  direct_summation_at_point(
    af::const_ref<miller::index<> > const& miller_indices,
    af::const_ref<std::complex<FloatType> > const& data,
    scitbx::vec3<FloatType> site_frac)
  {
    CCTBX_ASSERT(data.size() == miller_indices.size());
    std::complex<FloatType> sum(0.,0.);
    using scitbx::constants::pi;
    for (unsigned i = 0; i < data.size(); i++) {
      miller::index<> hkl = miller_indices[i];
      std::complex<FloatType> F_hkl = data[i];
      std::complex<FloatType> exponent(0, -2 * pi * (hkl*site_frac));
      sum += F_hkl * std::exp(exponent);
    }
    return sum;
  }

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_FFT_H
