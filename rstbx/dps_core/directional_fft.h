#ifndef RSTBX_DIRECTIONAL_FFT_H
#define RSTBX_DIRECTIONAL_FFT_H

#include <algorithm>
#include <scitbx/fftpack/real_to_complex.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid.h>

#include <rstbx/dps_core/direction.h>

namespace af = scitbx::af;
namespace rstbx {

typedef std::size_t                          sztype;
typedef af::shared<scitbx::vec3<double> >    veclist_t;
typedef af::const_ref<scitbx::vec3<double> > veclist_ref;
typedef af::shared<double>                   shareddouble_t;
typedef af::ref<double>                      shareddouble_ref;
typedef af::flex_double                      flexdouble;
typedef af::flex_complex_double              flexcomplex;
typedef af::flex_int                         flexint;

static const sztype F0_cutoff =  3; // smallest channel for direct-space maximum

sztype
primecheck(const sztype& input_m);

struct Directional_FFT {
  Directional_FFT (const Direction& angle, const veclist_t& xyzdata,
                   const double& granularity, const double& amax,
                   const sztype& = F0_cutoff);
  flexdouble
  power_spectrum();

  //result values to be publicly accessible after FFT:
  flexcomplex fft_result;
  sztype m;           // size of FFT domain
  double delta_p;     // reciprocal-space interval
  double pmin;        // minimum projection
  sztype kmax();      // position of the power spectrum's maximum value
  double kval();      // the power spectrum's maximum value
  double kval0();      // the power spectrum's maximum value
  double kval2();      // the power spectrum's maximum value
  double kval3();      // the power spectrum's maximum value

 private:
  veclist_ref xy; //reference to the data; invalid after constructor call
  sztype fft_n_complex;
  bool has_power_spectrum, has_kval;
  flexdouble pspectrum; //private cached copy of power spectrum
  sztype p_kmax;
  double p_kval;

  sztype F0_specific_cutoff; // smallest channel for direct-space maximum

};

} //namespace

#endif //RSTBX_DIRECTIONAL_FFT_H
