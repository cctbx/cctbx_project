#include <scitbx/fftpack/real_to_complex.h>

#include <rstbx/dps_core/directional_fft.h>

rstbx::Directional_FFT::Directional_FFT (
  const Direction& angle, const veclist_t& xyzdata,
  const double& granularity, const double& amax,
  const sztype& F0_specific_cutoff):
    xy(xyzdata.const_ref()),
    has_power_spectrum(false),
    has_kval(false),
    F0_specific_cutoff(F0_specific_cutoff){

      shareddouble_t   projections(xy.size(), af::init_functor_null<double>() );
      double* ps = projections.begin();

      for (sztype k=0; k < projections.size(); ++k) {
        // unroll the dot product; 1% increase in efficiency
        ps[k] = xy[k][0]*angle.dvec[0]+xy[k][1]*angle.dvec[1]+xy[k][2]*angle.dvec[2];
      }

      typedef shareddouble_t::const_iterator For;
      For b(projections.begin());
      For e(projections.end());
             pmin = *(std::min_element<For>(b,e));
      double pmax = *(std::max_element<For>(b,e));

      m = 1 + (sztype)((pmax - pmin) * granularity * amax );
      m = primecheck(m);
      m = std::max<int>(2,m); //Guard against extremely rare divide by zero with m==1

      delta_p = (pmax - pmin) / (m-1);

      scitbx::fftpack::real_to_complex<double> FFT(m);

      flexdouble frequency_distribution( scitbx::af::flex_grid<>(FFT.m_real()));
      double* f_d_ptr = frequency_distribution.begin();

      for (sztype k=0; k < projections.size(); ++k) {
        f_d_ptr[  (int)(((ps[k]-pmin)/delta_p) + 0.5)  ]+=1.0;
      }

      FFT.forward(frequency_distribution.begin());
      fft_n_complex = FFT.n_complex();
      fft_result = flexcomplex(frequency_distribution.handle(),
                         scitbx::af::flex_grid<>((fft_n_complex)) );
}

rstbx::flexdouble
rstbx::Directional_FFT::power_spectrum() {
      if (!has_power_spectrum) {
        pspectrum = flexdouble(fft_n_complex, af::init_functor_null<double>() );
        double* psptr = pspectrum.begin();
        std::complex<double>* resptr = fft_result.begin();
        for (sztype k=0; k<fft_n_complex; ++k) {
          psptr[k]=std::abs(resptr[k]);
        }
        has_power_spectrum = true;
      }
      return pspectrum;
}

rstbx::sztype
rstbx::Directional_FFT::kmax() {
      if (!has_kval) {
        kval();
      }
      return p_kmax;
}

double
rstbx::Directional_FFT::kval() {
      if (!has_kval) {
        p_kmax = 0;
        p_kval = 0.0;
        power_spectrum();
        scitbx::af::ref<double> ps(pspectrum.begin(),pspectrum.size());
        for (sztype k = F0_specific_cutoff; k<ps.size(); ++k){
          if ( ps[k] > p_kval ) { p_kval=ps[k]; p_kmax=k;}
        }
        has_kval = true;
      }
      return p_kval;
}

double
rstbx::Directional_FFT::kval0() {
      if (!has_kval) { kval(); }
      scitbx::af::ref<double> ps(pspectrum.begin(),pspectrum.size());
      { return ps[0]; }
}

double
rstbx::Directional_FFT::kval2() {
      if (!has_kval) { kval(); }
      scitbx::af::ref<double> ps(pspectrum.begin(),pspectrum.size());
      if (p_kmax * 2 < ps.size()) { return ps[p_kmax * 2]; }
      return 0.0;
}

double
rstbx::Directional_FFT::kval3() {
      if (!has_kval) { kval(); }
      scitbx::af::ref<double> ps(pspectrum.begin(),pspectrum.size());
      if (p_kmax * 3 < ps.size()) { return ps[p_kmax * 3]; }
      return 0.0;
}

rstbx::sztype
rstbx::primecheck(const sztype& input_m){
  sztype residual(input_m);
  sztype return_m(input_m);
  while (residual%2==0){ residual/=2; }
  while (residual%3==0){ residual/=3; }
  while (residual%5==0){ residual/=5; }
  if (residual>1) { return_m = primecheck(return_m+1); }
  return return_m;
}
