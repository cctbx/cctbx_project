#include <scitbx/fftpack/real_to_complex.h>

#include <rstbx/dps_core/directional_fft.h>

rstbx::Directional_FFT::Directional_FFT (
  Direction& angle, veclist_t& xyzdata,
  const double& granularity, const double& amax):
    xy(xyzdata.ref()),
    has_power_spectrum(false),
    has_kval(false),
    amax(amax),
    granularity(granularity),
    angle(angle){

      shareddouble_t   projections(xy.size());
      shareddouble_ref ps(projections.ref());//optimization

      for (sztype k=0; k<ps.size(); ++k) {
        ps[k] = xy[k]*angle.dvec;
      }

      typedef shareddouble_t::const_iterator For;
      For b(projections.begin());
      For e(projections.end());
             pmin = *(std::min_element<For>(b,e));
      double pmax = *(std::max_element<For>(b,e));

      m = 1 + (sztype)((pmax - pmin) * granularity * amax );
      m = primecheck(m);

      delta_p = (pmax - pmin) / (m-1);

      scitbx::fftpack::real_to_complex<double> FFT(m);

      flexdouble frequency_distribution( scitbx::af::flex_grid<>(FFT.m_real()));
      double* f_d_ptr = frequency_distribution.begin();

      for (sztype k=0; k<ps.size(); ++k) {
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
        pspectrum = flexdouble(fft_n_complex);
        for (sztype k=0; k<fft_n_complex; ++k) {
          pspectrum[k]=std::abs(fft_result[k]);
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
        for (sztype k = F0_cutoff; k<ps.size(); ++k){
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

void
rstbx::Directional_FFT::extract_directional_properties(const bool PS){
  angle.kmax = kmax();
  angle.kval = kval();
  angle.kval0 = kval0();
  angle.kval2 = kval2();
  angle.kval3 = kval3();
  angle.m = m;
  angle.delta_p = delta_p;
  angle.pmin = pmin;
  angle.uc_length = kmax()/ (m*delta_p);
  if (PS) {
    angle.ff = power_spectrum();
  }
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
