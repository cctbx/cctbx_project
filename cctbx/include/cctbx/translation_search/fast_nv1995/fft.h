/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     Oct 2002: Modified copy of phenix/fast_translation/fft.h (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_TRANSLATION_SEARCH_FAST_NV1995_FFT_H
#define CCTBX_TRANSLATION_SEARCH_FAST_NV1995_FFT_H

#include <scitbx/fftpack/complex_to_complex.h>
#include <scitbx/fftpack/real_to_complex.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/tiny_reductions.h>
#include <cctbx/import_scitbx_af.h>
#include <vector>

namespace cctbx { namespace translation_search { namespace fast_nv1995_detail {

  namespace fftpack = scitbx::fftpack;

  template <typename InRealType, typename OutRealType>
  struct shrinking_complex_to_real_fft
  {
    typedef InRealType in_real_type;
    typedef std::complex<in_real_type> in_complex_type;
    typedef OutRealType out_real_type;
    typedef std::complex<out_real_type> out_complex_type;

    template <typename AccessorType>
    static
    af::ref<out_real_type, AccessorType>
    run(af::ref<in_complex_type, AccessorType> ccc,
        af::int3 const& n_real_in,
        af::int3 const& n_real_out)
    {
      using af::int3;

      std::size_t fred = sizeof(in_real_type) / sizeof(out_real_type);
      int3 sampling = n_real_in / n_real_out;
      AccessorType const& dim_in = ccc.accessor();

      fftpack::complex_to_complex<out_real_type> fft_x(n_real_in[0]);
      fftpack::complex_to_complex<out_real_type> fft_y(n_real_in[1]);
      fftpack::real_to_complex<out_real_type>    fft_z(n_real_in[2]);

      std::vector<out_complex_type>
        c_seq(af::max(int3(fft_x.n(), fft_y.n(), fft_z.n_complex())));

      af::ref<out_complex_type, AccessorType>
      rcc(reinterpret_cast<out_complex_type*>(ccc.begin()),
        AccessorType(n_real_out[0], dim_in[1], dim_in[2] * fred));

      af::ref<out_complex_type, AccessorType>
      rrc(rcc.begin(),
        AccessorType(n_real_out[0], n_real_out[1], dim_in[2] * fred));

      af::ref<out_real_type, AccessorType>
      rrr(reinterpret_cast<out_real_type*>(ccc.begin()),
        AccessorType(n_real_out));

      for (std::size_t iz = 0; iz < fft_z.n_complex(); iz++) {
        for (std::size_t iy = 0; iy < n_real_in[1]; iy++) {
          for (std::size_t ix = 0; ix < n_real_in[0]; ix++) {
            c_seq[ix] = out_complex_type(ccc(ix, iy, iz));
          }
          fft_x.backward(c_seq.begin()); // Transform along x (slow direction)
          for (std::size_t ix = 0; ix < n_real_out[0]; ix++) {
            rcc(ix, iy, iz) = c_seq[ix * sampling[0]];
          }
        }
        for (std::size_t ix = 0; ix < n_real_out[0]; ix++) {
          for (std::size_t iy = 0; iy < n_real_in[1]; iy++) {
            c_seq[iy] = rcc(ix, iy, iz);
          }
          fft_y.backward(c_seq.begin());// Transform along y (medium direction)
          for (std::size_t iy = 0; iy < n_real_out[1]; iy++) {
            rrc(ix, iy, iz) = c_seq[iy * sampling[1]];
          }
        }
      }
      for (std::size_t ix = 0; ix < n_real_out[0]; ix++) {
        for (std::size_t iy = 0; iy < n_real_out[1]; iy++) {
          out_real_type*
            r_seq = reinterpret_cast<out_real_type*>(&rrc(ix, iy, 0));
          fft_z.backward(r_seq); // Transform along z (fast direction)
          for (std::size_t iz = 0; iz < n_real_out[2]; iz++) {
            rrr(ix, iy, iz) = r_seq[iz * sampling[2]];
          }
        }
      }
      return rrr;
    }
  };

}}} // namespace cctbx::translation_search::fast_nv1995_detail

#endif // CCTBX_TRANSLATION_SEARCH_FAST_NV1995_FFT_H
